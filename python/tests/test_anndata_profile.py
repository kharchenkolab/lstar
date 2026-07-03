"""AnnData profile round-trip -- grounded in **real** data (pbmc68k_reduced: real X + .raw counts over a
divergent gene set, real obs categoricals/numerics, real obsm pca/umap, varm PCs, obsp graphs, uns).

  AnnData --read_anndata--> L* --write_anndata--> AnnData            (profile round trip)
  AnnData --read_anndata--> L* --zarr--> L* --write_anndata--> AnnData (full pipeline)
  repeated conversion is a fixed point (stable over cycles).

Run: PYTHONPATH=python/src python3 python/tests/test_anndata_profile.py
"""
import os
import sys
import tempfile
import warnings

import numpy as np
import scipy.sparse as sp

sys.path.insert(0, os.path.dirname(__file__))
import corpus  # noqa: E402

from lstar import read, read_anndata, write, write_anndata  # noqa: E402

warnings.filterwarnings("ignore")


def _dense(x):
    return x.toarray() if sp.issparse(x) else np.asarray(x)


def _eq(x, y):
    # equal_nan: real data carries NaN (e.g. undefined PCA loadings in varm['PCs']); NaN is a valid
    # float that round-trips faithfully through the dense encoding, so NaN==NaN counts as equal here.
    return np.allclose(_dense(x).astype(float), _dense(y).astype(float), rtol=1e-5, atol=1e-6, equal_nan=True)


def check(a, a2):
    assert list(a2.obs_names) == list(a.obs_names)
    assert list(a2.var_names) == list(a.var_names)
    assert _eq(a2.X, a.X)
    assert a2.raw is not None and _eq(a2.raw.X, a.raw.X)               # real .raw counts (divergent genes)
    for k in a.obsm:                                                    # X_pca, X_umap
        assert _eq(a2.obsm[k], a.obsm[k]), k
    for k in a.varm:                                                    # PCs
        assert _eq(a2.varm[k], a.varm[k]), k
    for k in a.obsp:                                                    # distances, connectivities
        assert _eq(a2.obsp[k], a.obsp[k]), k
    for c in a.obs.columns:                                             # categoricals + numerics
        assert list(a2.obs[c].astype(str)) == list(a.obs[c].astype(str)), c


def run():
    a = corpus.pbmc68k_reduced()
    if a is None:
        print("  SKIP test_anndata_profile (corpus unavailable)"); return
    ds = read_anndata(a)

    # shared-vocabulary signatures on real data
    assert ds.field("X").role == "measure" and ds.field("X").state == "scaled"   # scaled .X, inferred from content
    assert ds.field("raw").state == "lognorm"                                    # pbmc68k .raw is log-normalized, not raw counts
    assert ds.field("pca").role == "embedding" and "pca" in ds.axes
    assert ds.field("bulk_labels").role == "label" and ds.axis("bulk_labels").role == "factor"
    assert ds.field("connectivities").role == "relation"

    # (1) profile-only round trip
    check(a, write_anndata(ds))

    # (2) full pipeline through the zarr store; the untyped uns tail survives serialization
    p = os.path.join(tempfile.mkdtemp(), "a.lstar.zarr")
    write(ds, p)
    ds2 = read(p)
    assert "anndata.uns" in ds2.aux                                    # passthrough survived the store
    check(a, write_anndata(ds2))

    # (3) repeated conversion is a fixed point
    cur = a
    for _ in range(3):
        cur = write_anndata(read_anndata(cur))
        check(a, cur)

    print("anndata profile (real pbmc68k): round-trip via profile + via zarr; fixed point over 3 cycles; "
          "%d fields, %d axes" % (len(ds.fields), len(ds.axes)))


def test_anndata_roundtrip():
    run()


def test_native_direct_backend_state_parity(tmp_path):
    """The native (anndata) and direct (h5py) `.h5ad` backends must infer the SAME measure state + names,
    so a store's `state` (and thus viewer-prep success) doesn't depend on whether `anndata` is installed
    (audit T1.2). Direct previously left X `state=None`, named `X`; now it shares `_infer_state`."""
    import anndata as ad
    from lstar.profiles.anndata_direct import read_h5ad_direct
    rng = np.random.default_rng(0)
    Xraw = sp.random(30, 10, density=0.4, format="csr", random_state=0)
    Xraw.data = (rng.poisson(3, Xraw.data.shape) + 1).astype(np.float32)          # genuinely-raw counts
    a = ad.AnnData(X=Xraw)
    a.layers["lognorm"] = a.X.copy(); a.layers["lognorm"].data = np.log1p(a.layers["lognorm"].data)
    h5 = str(tmp_path / "t.h5ad"); a.write_h5ad(h5)
    meas = lambda ds: {n: (f.role, f.state) for n, f in ds.fields.items() if f.role == "measure"}
    native, direct = meas(read_anndata(a)), meas(read_h5ad_direct(h5))
    assert native == direct, f"backend state/name drift: native={native} direct={direct}"
    assert native.get("counts") == ("measure", "raw")                            # raw X renamed to counts on both


if __name__ == "__main__":
    run()
