"""Lossless passthrough: a format's untyped long-tail (AnnData `uns`, Seurat `@misc`) round-trips
**verbatim** through a self-describing `aux/` subtree -- nested dicts/lists, scalars, numeric/string
arrays, and numpy structured arrays (`rank_genes_groups`) all survive Py->L*->Py exactly, instead of
being recorded name-only in `dropped`.

Run: PYTHONPATH=python/src python3 python/tests/test_aux.py
"""
import os
import tempfile

import numpy as np

import lstar
from lstar.passthrough import from_store, to_store


def _store():
    return os.path.join(tempfile.mkdtemp(), "aux.lstar.zarr")


def _rich_uns():
    return {
        "log1p": {"base": None},
        "iroot": 42,
        "a_flag": True,
        "a_list": [1, 2, 3],
        "pca": {"variance_ratio": np.array([0.5, 0.3, 0.2]),
                "params": {"n_comps": 50, "zero_center": True}},
        "leiden_colors": np.array(["#ff0000", "#00ff00", "#0000ff"]),
        "neighbors": {"connectivities_key": "connectivities",
                      "params": {"n_neighbors": 15, "method": "umap", "metric": "euclidean"}},
        "dendrogram_leiden": {"linkage": np.array([[0., 1., 0.5, 2.], [2., 3., 0.7, 3.]])},
        "rank_genes_groups": {
            "params": {"groupby": "leiden", "method": "t-test", "reference": "rest"},
            "names": np.array([("g1", "g2"), ("g3", "g4")], dtype=[("A", "U10"), ("B", "U10")]),
            "scores": np.array([(1.5, 2.5), (3.5, 4.5)], dtype=[("A", "f4"), ("B", "f4")]),
        },
    }


def _eq(a, b):
    if isinstance(a, dict):
        return isinstance(b, dict) and list(a) == list(b) and all(_eq(a[k], b[k]) for k in a)
    if isinstance(a, list):
        return isinstance(b, list) and len(a) == len(b) and all(_eq(x, y) for x, y in zip(a, b))
    if isinstance(a, np.ndarray):
        b = np.asarray(b)
        if a.dtype.names:                          # structured array: compare field by field
            return a.dtype.names == b.dtype.names and all(_eq(a[n], b[n]) for n in a.dtype.names)
        if a.dtype.kind in ("f",):
            return a.shape == b.shape and np.allclose(a, b)
        return a.shape == b.shape and (a == b).all()
    return a == b


def test_serializer_roundtrip():
    obj = _rich_uns()
    tree, arrays = to_store(obj)
    back = from_store(tree, arrays)
    assert _eq(obj, back), "serializer round-trip mismatch"
    # the tree is pure JSON (array leaves are references, not inline data)
    import json
    json.dumps(tree)
    print("aux serializer: nested dict/list/scalars + numeric/string/structured arrays round-trip exact")


def test_store_roundtrip():
    ds = lstar.Dataset(kind="sample")
    ds.add_axis("cells", [f"c{i}" for i in range(3)])
    ds.aux["anndata.uns"] = _rich_uns()
    p = _store(); lstar.write(ds, p)
    ds2 = lstar.read(p)
    assert "anndata.uns" in ds2.aux
    assert _eq(ds.aux["anndata.uns"], ds2.aux["anndata.uns"])
    # a rank_genes_groups structured array survives with its dtype + values
    rg = ds2.aux["anndata.uns"]["rank_genes_groups"]
    assert rg["names"].dtype.names == ("A", "B")
    assert list(rg["names"]["A"]) == ["g1", "g3"]
    print("aux store round-trip: uns preserved verbatim through the L* store (incl. structured arrays)")


def test_anndata_uns_passthrough():
    import anndata as ad
    import pandas as pd
    import scipy.sparse as sp
    from lstar import read_anndata, write_anndata, write, read

    a = ad.AnnData(X=sp.random(4, 3, density=0.4, format="csr", random_state=0),
                   obs=pd.DataFrame(index=[f"c{i}" for i in range(4)]))
    a.uns["pca"] = {"variance_ratio": np.array([0.6, 0.4]), "params": {"n_comps": 2}}
    a.uns["leiden_colors"] = np.array(["#aaa", "#bbb"])
    a.uns["log1p"] = {"base": None}

    ds = read_anndata(a)
    assert "anndata.uns" in ds.aux                  # captured, not just named in `dropped`
    a2 = write_anndata(read(_write(ds)))            # round-trip through the store
    assert np.allclose(a2.uns["pca"]["variance_ratio"], a.uns["pca"]["variance_ratio"])
    assert a2.uns["pca"]["params"]["n_comps"] == 2
    assert list(a2.uns["leiden_colors"]) == ["#aaa", "#bbb"]
    assert "log1p" in a2.uns
    print("anndata uns passthrough: params + colors + nested dicts reproduced on write-back")


def _write(ds):
    p = _store(); lstar.write(ds, p); return p


if __name__ == "__main__":
    test_serializer_roundtrip()
    test_store_roundtrip()
    test_anndata_uns_passthrough()
