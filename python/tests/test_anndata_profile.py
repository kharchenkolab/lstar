"""M2: AnnData profile round-trip.

  AnnData --read_anndata--> L* Dataset --write_anndata--> AnnData     (profile round trip)
  AnnData --read_anndata--> L* --zarr--> L* --write_anndata--> AnnData (full pipeline)

Both must preserve X, layers, obs, var, obsm, varm, obsp on the shared-vocabulary core.

Run: PYTHONPATH=python/src python3 python/tests/test_anndata_profile.py
"""
import os
import tempfile

import numpy as np
import scipy.sparse as sp

import anndata as ad
import pandas as pd

from lstar import read_anndata, write_anndata, write, read


def make_adata():
    n, g, pc = 80, 40, 8
    rng = np.random.default_rng(1)
    X = sp.random(n, g, density=0.2, format="csr", random_state=1)
    counts = sp.random(n, g, density=0.2, format="csr", random_state=2)
    obs = pd.DataFrame(
        {"leiden": pd.Categorical(["c%d" % (i % 4) for i in range(n)]),
         "n_umi": rng.integers(100, 1000, n).astype(float)},
        index=["cell%d" % i for i in range(n)])
    var = pd.DataFrame(
        {"dispersion": rng.standard_normal(g).astype(float),
         "highly_variable": (rng.random(g) > 0.5)},
        index=["g%d" % j for j in range(g)])
    a = ad.AnnData(X=X, obs=obs, var=var, layers={"counts": counts})
    a.obsm["X_pca"] = rng.standard_normal((n, pc)).astype("float32")
    a.obsm["X_umap"] = rng.standard_normal((n, 2)).astype("float32")
    a.varm["PCs"] = rng.standard_normal((g, pc)).astype("float32")
    a.obsp["connectivities"] = sp.random(n, n, density=0.1, format="csr", random_state=3)
    a.uns["something"] = {"note": "not in the shared vocabulary"}
    return a


def _eq(x, y):
    if sp.issparse(x) or sp.issparse(y):
        x = x.toarray() if sp.issparse(x) else np.asarray(x)
        y = y.toarray() if sp.issparse(y) else np.asarray(y)
    return np.allclose(np.asarray(x, dtype=float), np.asarray(y, dtype=float))


def check(a, a2):
    assert list(a2.obs_names) == list(a.obs_names)
    assert list(a2.var_names) == list(a.var_names)
    assert _eq(a2.X, a.X)
    assert _eq(a2.layers["counts"], a.layers["counts"])
    assert _eq(a2.obsm["X_pca"], a.obsm["X_pca"])
    assert _eq(a2.obsm["X_umap"], a.obsm["X_umap"])
    assert _eq(a2.varm["PCs"], a.varm["PCs"])
    assert _eq(a2.obsp["connectivities"], a.obsp["connectivities"])
    assert (a2.obs["leiden"].astype(str).values == a.obs["leiden"].astype(str).values).all()
    assert np.allclose(a2.obs["n_umi"].values.astype(float), a.obs["n_umi"].values.astype(float))
    assert np.allclose(a2.var["dispersion"].values.astype(float), a.var["dispersion"].values.astype(float))
    assert (a2.var["highly_variable"].astype(str).values == a.var["highly_variable"].astype(str).values).all()


def run():
    a = make_adata()
    ds = read_anndata(a)

    # shared-vocabulary signatures
    assert ds.field("X").role == "measure"
    assert ds.field("counts").role == "measure" and ds.field("counts").state == "raw"
    assert ds.field("pca").role == "embedding" and "pca" in ds.axes
    assert ds.field("pca_loadings").role == "loading" and ds.field("pca_loadings").span == ["genes", "pca"]
    assert ds.field("umap").role == "embedding"
    assert ds.field("leiden").role == "label"
    assert ds.field("connectivities").role == "relation"

    # loss is recorded, not silent: the uns extra is dropped at read time and noted
    assert "uns/something" in ds.dropped

    # (1) profile-only round trip
    check(a, write_anndata(ds))

    # (2) full pipeline through the zarr store
    p = os.path.join(tempfile.mkdtemp(), "a.lstar.zarr")
    write(ds, p)
    ds2 = read(p)
    assert "uns/something" in ds2.dropped  # the loss record survives serialization
    a3 = write_anndata(ds2)
    check(a, a3)

    # (3) variable-length round trip: AnnData -> L* -> AnnData, repeated, is a fixed point, so a
    # conversion chain of any length returns to the original native format unchanged.
    cur = a
    for _ in range(4):
        cur = write_anndata(read_anndata(cur))
        check(a, cur)

    print("anndata profile OK: round-trip via profile and via zarr; fixed point over 4 cycles; "
          "%d fields, %d axes" % (len(ds.fields), len(ds.axes)))


def test_anndata_roundtrip():
    run()


if __name__ == "__main__":
    run()
