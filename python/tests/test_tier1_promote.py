"""Tier-1 promotions out of the lossless passthrough into typed, axis-bound fields:
  - color palettes  uns['<key>_colors']        -> a label field over the factor axis (category order)
  - PCA variance     uns['pca']['variance*']    -> measures over the pca coordinate axis
  - RNA-velocity     layers/obsp/var fit_*      -> measures + a relation (free via the generic paths)
All round-trip through the store (ordinary fields) and regenerate their native uns/layers on write-back.

Run: PYTHONPATH=python/src python3 python/tests/test_tier1_promote.py
"""
import os
import tempfile
import warnings

import numpy as np

import lstar

warnings.filterwarnings("ignore")


def _store():
    return os.path.join(tempfile.mkdtemp(), "t1.lstar.zarr")


def _adata():
    import anndata as ad
    import pandas as pd
    import scanpy as sc
    import scipy.sparse as sp
    rng = np.random.default_rng(0); n, g = 80, 24
    a = ad.AnnData(X=sp.csr_matrix(rng.poisson(0.5, (n, g)).astype("float32")),
                   obs=pd.DataFrame({"leiden": pd.Categorical([str(i % 3) for i in range(n)])},
                                    index=[f"c{i}" for i in range(n)]),
                   var=pd.DataFrame(index=[f"g{j}" for j in range(g)]))
    sc.pp.normalize_total(a); sc.pp.log1p(a); sc.pp.pca(a, n_comps=8)
    a.uns["leiden_colors"] = np.array(["#1f77b4", "#ff7f0e", "#2ca02c"])   # one per leiden category
    return a


def test_color_palette_promotion():
    a = _adata()
    ds = lstar.read_anndata(a)
    cf = ds.field("leiden_colors")
    assert cf.subtype == "color" and cf.span == ["leiden"]                # bound to the factor axis
    assert len(np.asarray(cf.values)) == len(ds.axis("leiden"))           # one color per category
    assert "leiden_colors" not in ds.aux.get("anndata.uns", {})           # promoted out of the tail
    assert not lstar.validate(ds)

    p = _store(); lstar.write(ds, p)
    a2 = lstar.write_anndata(lstar.read(p))
    assert list(a2.uns["leiden_colors"]) == list(a.uns["leiden_colors"])  # regenerated on write-back
    print("colors: uns['leiden_colors'] -> field over the leiden factor axis; round-trips + regenerates")


def test_pca_variance_promotion():
    a = _adata()
    ds = lstar.read_anndata(a)
    vr = ds.field("pca_variance_ratio")
    assert vr.subtype == "pca_var" and vr.span == ["pca"]                 # measure over the pca axis
    assert np.allclose(np.asarray(vr.values), a.uns["pca"]["variance_ratio"])  # input not mutated
    assert "variance_ratio" not in ds.aux.get("anndata.uns", {}).get("pca", {})  # promoted out of tail
    p = _store(); lstar.write(ds, p)
    a2 = lstar.write_anndata(lstar.read(p))
    assert np.allclose(a2.uns["pca"]["variance_ratio"], a.uns["pca"]["variance_ratio"])
    print("pca: uns['pca']['variance_ratio'] -> measure over the pca axis; round-trips + regenerates")


def test_velocity_free_path():
    import anndata as ad
    import pandas as pd
    import scipy.sparse as sp
    rng = np.random.default_rng(1); n, g = 40, 15
    a = ad.AnnData(X=sp.random(n, g, density=0.3, format="csr", random_state=1),
                   obs=pd.DataFrame(index=[f"c{i}" for i in range(n)]),
                   var=pd.DataFrame(index=[f"g{j}" for j in range(g)]))
    a.layers["spliced"] = sp.random(n, g, density=0.3, format="csr", random_state=2)
    a.layers["unspliced"] = sp.random(n, g, density=0.3, format="csr", random_state=3)
    a.var["fit_likelihood"] = rng.random(g).astype("float32")
    a.obsp["velocity_graph"] = sp.random(n, n, density=0.2, format="csr", random_state=4)

    ds = lstar.read_anndata(a)
    # velocity comes free: layers -> measures, var fit_* -> measure over genes, graph -> relation
    assert ds.field("spliced").role == "measure" and ds.field("spliced").span == ["cells", "genes"]
    assert ds.field("unspliced").role == "measure"
    assert ds.field("fit_likelihood").span == ["genes"]
    assert ds.field("velocity_graph").role == "relation" and ds.field("velocity_graph").span == ["cells", "cells"]
    assert not lstar.validate(ds)
    p = _store(); lstar.write(ds, p); a2 = lstar.write_anndata(lstar.read(p))
    assert "spliced" in a2.layers and "unspliced" in a2.layers
    assert np.allclose(a2.var["fit_likelihood"].values.astype(float), a.var["fit_likelihood"].values.astype(float))
    print("velocity: spliced/unspliced layers + var fit_* + velocity_graph captured & round-trip (free path)")


if __name__ == "__main__":
    test_color_palette_promotion()
    test_pca_variance_promotion()
    test_velocity_free_path()
