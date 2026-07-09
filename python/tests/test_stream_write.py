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
from lstar.profiles.anndata import (read_anndata, write_anndata, convert_anndata,
                                    write_anndata_streamed, _BackedH5Sparse)


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


def test_streaming_lazy_cross_format():
    # Bounded-memory conversion ACROSS Zarr formats: read a store lazily (LazyCSX streaming source) and
    # stream-write it to the OTHER on-disk format. Guards the asymmetric lazy path (v2->v3 and v3->v2);
    # the symmetric case above only exercises same-format. The lazy source must genuinely stream (not
    # silently materialize), the output must be a genuine v2/v3 store, and values must survive.
    from lstar.lazy import is_stream_source
    d = tempfile.mkdtemp()
    X = sp.csc_matrix(sp.random(90, 45, density=0.2, format="csc", random_state=5)); X.data = X.data * 5 + 1
    base = lstar.Dataset(kind="sample")
    base.add_axis("cells", [f"c{i}" for i in range(90)])
    base.add_axis("genes", [f"g{i}" for i in range(45)])
    base.add_field("counts", X, role="measure", span=["cells", "genes"], state="raw")
    src = {f: os.path.join(d, f"src_{f}.lstar.zarr") for f in ("v2", "v3")}
    for f in ("v2", "v3"):
        lstar.write(base, src[f], chunk_elems=300, format=f)   # chunked so streaming spans >1 block

    marker = {"v2": ".zmetadata", "v3": "zarr.json"}            # a format-genuine top-level file
    for src_fmt, out_fmt in (("v2", "v3"), ("v3", "v2")):
        lazy = lstar.read(src[src_fmt], lazy=True)
        assert is_stream_source(lazy.field("counts").values)   # a streaming proxy, not a materialized matrix
        out = os.path.join(d, f"out_{src_fmt}_{out_fmt}.lstar.zarr")
        lstar.write(lazy, out, stream=True, format=out_fmt)     # bounded-memory write to the other format
        assert os.path.exists(os.path.join(out, marker[out_fmt])), f"{out_fmt} store not genuine"
        assert not lstar.validate(lstar.read(out))
        assert _sparse_eq(X, lstar.read(out).field("counts").values), f"{src_fmt}->{out_fmt}"
    print("streaming lazy cross-format (v2->v3 and v3->v2) == original")


def test_streamed_h5ad_write_equals_eager():
    # The reverse direction: write an L* dataset to an .h5ad with bounded memory (small parts via
    # anndata, big sparse measures streamed block-by-block via h5py). Must equal the eager write,
    # preserving each measure's native orientation (X csr, counts csc) and a .raw over a larger gene
    # set. Streamed from an on-disk store read lazily -- the actual bounded path.
    import anndata as ad
    import pandas as pd
    d = tempfile.mkdtemp()
    X = sp.random(140, 60, density=0.12, format="csr", random_state=0); X.data = np.round(X.data * 9 + 1).astype("f4")
    counts = sp.csc_matrix(sp.random(140, 60, density=0.2, format="csc", random_state=1)); counts.data = np.round(counts.data * 9 + 1).astype("f4")
    rawX = sp.random(140, 95, density=0.1, format="csr", random_state=2); rawX.data = np.round(rawX.data * 9 + 1).astype("f4")
    a = ad.AnnData(X=X)
    a.layers["counts"] = counts
    a.obsm["X_pca"] = np.random.RandomState(3).randn(140, 5).astype("f4")
    a.obs["leiden"] = pd.Categorical([f"c{i % 6}" for i in range(140)])
    a.raw = ad.AnnData(X=rawX, obs=pd.DataFrame(index=a.obs_names),
                       var=pd.DataFrame(index=[f"g{i}" for i in range(95)]))
    h5 = os.path.join(d, "src.h5ad"); a.write_h5ad(h5)

    # round-trip the source through an on-disk L* store, then stream that store back out to h5ad
    store = os.path.join(d, "s.lstar.zarr")
    convert_anndata(h5, store, chunk_elems=500)
    eager = write_anndata(lstar.read(store)); eager_h5 = os.path.join(d, "eager.h5ad"); eager.write_h5ad(eager_h5)
    stream_h5 = os.path.join(d, "stream.h5ad")
    lstar.convert_to_h5ad(store, stream_h5, chunk_elems=500)

    e, s = ad.read_h5ad(eager_h5), ad.read_h5ad(stream_h5)
    assert _sparse_eq(e.X, s.X), "X"
    assert type(s.X).__name__ == "csr_matrix" and type(s.layers["counts"]).__name__ == "csc_matrix"  # orientation kept
    assert _sparse_eq(e.layers["counts"], s.layers["counts"]), "counts"
    assert _sparse_eq(e.raw.X, s.raw.X) and s.raw.X.shape == (140, 95), "raw"
    assert np.allclose(e.obsm["X_pca"], s.obsm["X_pca"])
    assert (np.asarray(e.obs["leiden"]) == np.asarray(s.obs["leiden"])).all()
    print("streamed L*->h5ad == eager (X csr / counts csc / raw over larger gene set / pca / leiden)")


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
    test_streaming_lazy_cross_format()
    test_streamed_h5ad_write_equals_eager()
    test_legacy_h5ad_sparse_recognized()
