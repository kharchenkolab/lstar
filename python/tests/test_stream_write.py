"""Streaming write: a sparse measure is written to the L* store block-by-block from a backed/lazy
source, never materialized -- the basis for bounded-memory h5ad->L* and L*->L* conversions. The
streamed store must be byte-identical to the eager one, and the on-disk h5ad sparse format must be
recognized across versions (modern `encoding-type` and legacy `h5sparse_format`).
"""
import os
import tempfile

import numpy as np
import scipy.sparse as sp

import lstar
from lstar.profiles.anndata import read_anndata, convert_anndata, _BackedH5Sparse


def _sparse_eq(a, b):
    return a.shape == b.shape and (sp.csr_matrix(a) != sp.csr_matrix(b)).nnz == 0


def test_streaming_h5ad_equals_eager():
    import anndata as ad
    d = tempfile.mkdtemp()
    X = sp.random(120, 50, density=0.15, format="csr", random_state=0)
    X.data = np.round(X.data * 9 + 1).astype("f4")
    a = ad.AnnData(X=X)
    a.layers["counts"] = sp.csr_matrix(X)
    a.obsm["X_pca"] = np.random.RandomState(1).randn(120, 6).astype("f4")
    a.obs["leiden"] = np.array([f"c{i % 5}" for i in range(120)])
    h5 = os.path.join(d, "t.h5ad")
    a.write_h5ad(h5)

    eager = os.path.join(d, "e.lstar.zarr")
    stream = os.path.join(d, "s.lstar.zarr")
    lstar.write(read_anndata(ad.read_h5ad(h5)), eager)     # in-memory: the whole matrix is resident
    convert_anndata(h5, stream)                            # backed read + block-by-block streamed write

    de, ds = lstar.read(eager), lstar.read(stream)
    assert sorted(de.fields) == sorted(ds.fields)
    assert not lstar.validate(ds)
    for nm in ("X", "counts"):
        assert de.field(nm).encoding in ("csr", "csc")
        assert _sparse_eq(de.field(nm).values, ds.field(nm).values), nm
    assert np.allclose(np.asarray(de.field("pca").values), np.asarray(ds.field("pca").values))
    assert (np.asarray(de.field("leiden").values) == np.asarray(ds.field("leiden").values)).all()
    print("streaming h5ad->L* == eager (X/counts/pca/leiden identical)")


def test_streaming_lstar_to_lstar():
    # Read a store lazily (field values become LazyCSX streaming sources), then re-write it streamed
    # with a different chunking/compression -- a bounded-memory recompress. Result must equal source.
    import numcodecs
    d = tempfile.mkdtemp()
    X = sp.csc_matrix(sp.random(80, 40, density=0.2, format="csc", random_state=2))
    X.data = X.data * 5 + 1
    base = lstar.Dataset(kind="sample")
    base.add_axis("cells", [f"c{i}" for i in range(80)])
    base.add_axis("genes", [f"g{i}" for i in range(40)])
    base.add_field("counts", X, role="measure", span=["cells", "genes"], state="raw")
    src = os.path.join(d, "src.lstar.zarr")
    lstar.write(base, src, chunk_elems=500)

    lazy = lstar.read(src, lazy=True)                      # values are LazyCSX (a streaming source)
    out = os.path.join(d, "out.lstar.zarr")
    lstar.write(lazy, out, stream=True, compressor=numcodecs.GZip(5))
    assert _sparse_eq(lstar.read(src).field("counts").values, lstar.read(out).field("counts").values)
    print("streaming L*->L* (lazy read -> streamed recompress) == original")


def test_legacy_h5ad_sparse_recognized():
    # A pre-0.7 h5ad uses h5sparse_format/h5sparse_shape attrs; the streaming source must read it.
    import h5py
    d = tempfile.mkdtemp()
    p = os.path.join(d, "legacy.h5")
    M = sp.random(30, 20, density=0.25, format="csr", random_state=3)
    M.data = np.round(M.data * 9 + 1).astype("f4")
    with h5py.File(p, "w") as f:
        g = f.create_group("X")
        g.create_dataset("data", data=M.data)
        g.create_dataset("indices", data=M.indices)
        g.create_dataset("indptr", data=M.indptr)
        g.attrs["h5sparse_format"] = "csr"
        g.attrs["h5sparse_shape"] = np.array([30, 20])
    s = _BackedH5Sparse(p, "X")
    assert s.fmt == "csr" and s.shape == (30, 20) and s.nnz == M.nnz
    rebuilt = sp.vstack([sub for _, _, sub in s.blocks(7)]).tocsr()   # reassemble from blocks
    assert _sparse_eq(rebuilt, M)
    s.close()
    print("legacy h5ad sparse (h5sparse_format) recognized + streamed")


if __name__ == "__main__":
    test_streaming_h5ad_equals_eager()
    test_streaming_lstar_to_lstar()
    test_legacy_h5ad_sparse_recognized()
