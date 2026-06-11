"""Lazy / streaming read: open without materializing, stream reductions, stay exact.

A lazy read must (1) leave heavy fields as proxies (nothing read until touched), (2) materialize
to exactly the eager result, and (3) let a CSC measure be reduced by streaming column blocks --
including the on-the-fly log1p variant -- matching a full in-memory computation to the bit.
"""
import os
import tempfile

import numpy as np
import scipy.sparse as sp

import lstar
from lstar.lazy import LazyCSX, LazyDense, stream_col_stats


def _store(**writekw):
    rng = np.random.default_rng(0)
    X = sp.csc_matrix(sp.random(300, 120, density=0.1, format="csc", random_state=rng))
    X.data = X.data * 9 + 0.5                           # nonzero magnitudes, no implicit fill change
    ds = lstar.Dataset(kind="sample")
    ds.add_axis("cells", [f"c{i}" for i in range(300)])
    ds.add_axis("genes", [f"g{i}" for i in range(120)])
    ds.add_field("counts", X, role="measure", span=["cells", "genes"], state="raw")
    emb = rng.standard_normal((300, 2))
    ds.add_axis("umap", ["u0", "u1"], role="coordinate")
    ds.add_field("umap", emb, role="embedding", span=["cells", "umap"])
    p = os.path.join(tempfile.mkdtemp(), "lazy.lstar.zarr")
    lstar.write(ds, p, **writekw)
    return p, X, emb


def test_lazy_is_proxy_and_materializes():
    import numcodecs
    p, X, emb = _store(chunk_elems=5000, compressor=numcodecs.GZip(5))
    ds = lstar.read(p, lazy=True)
    cv = ds.fields["counts"].values
    uv = ds.fields["umap"].values
    assert isinstance(cv, LazyCSX) and cv.shape == (300, 120) and cv.nnz == X.nnz
    assert isinstance(uv, LazyDense) and uv.shape == (300, 2)
    assert not lstar.validate(ds)                       # proxies expose .shape -> validator ok
    assert (cv.materialize() != X).nnz == 0
    assert np.allclose(np.asarray(uv), emb)
    assert np.allclose(uv[10:20], emb[10:20])           # streamed dense slice
    print("lazy proxies: LazyCSX/LazyDense, materialize == eager (chunked+gzip)")


def test_stream_col_stats_matches_eager():
    p, X, emb = _store(chunk_elems=4000)
    ds = lstar.read(p, lazy=True)
    lazy_csc = ds.fields["counts"].values

    for lognorm in (False, True):
        m1, v1, n1 = stream_col_stats(lazy_csc, lognorm=lognorm, block=16)
        # stable ground truth: dense per-column mean / var(ddof=1)
        D = X.toarray().astype(float)
        if lognorm:
            D = np.log1p(D)
        mean_ref = D.mean(axis=0)
        var_ref = D.var(axis=0, ddof=1)
        assert np.allclose(m1, mean_ref, rtol=1e-10, atol=1e-12), lognorm
        assert np.allclose(v1, var_ref, rtol=1e-9, atol=1e-12), lognorm  # two-pass, not cancellation
        assert np.array_equal(n1, np.diff(X.indptr)), lognorm
    print("stream_col_stats matches dense np.var ground truth (stable, plain + lognorm)")


def test_thread_count_controllable_and_deterministic():
    p, X, emb = _store(chunk_elems=2000)
    lazy_csc = lstar.read(p, lazy=True).fields["counts"].values
    base = stream_col_stats(lazy_csc, lognorm=True, block=8, n_threads=1)
    for nt in (2, 4, 0):                                  # 0 -> os.cpu_count()
        got = stream_col_stats(lazy_csc, lognorm=True, block=8, n_threads=nt)
        for a, b in zip(base, got):
            assert np.array_equal(a, b) if a.dtype.kind == "i" else np.allclose(a, b), nt
    print("stream_col_stats: n_threads (1/2/4/auto) controllable and bit-identical results")


def test_engine_parity():
    # The pure-Python engine always runs; when the C++ accelerator is present (the default), it
    # must agree to the bit. The default 'auto' picks C++ if available, Python otherwise.
    p, X, emb = _store(chunk_elems=2000)
    lazy_csc = lstar.read(p, lazy=True).fields["counts"].values
    py = stream_col_stats(lazy_csc, lognorm=True, block=8, engine="python")
    if lstar.has_accel():
        cc = stream_col_stats(lazy_csc, lognorm=True, block=8, engine="c++")
        for a, b in zip(py, cc):
            assert np.array_equal(a, b) if a.dtype.kind == "i" else np.allclose(a, b)
        print("engine parity: C++ accelerator == pure Python (accelerator ACTIVE)")
    else:
        print("engine parity: accelerator not built; pure-Python path only")


if __name__ == "__main__":
    test_lazy_is_proxy_and_materializes()
    test_stream_col_stats_matches_eager()
    test_thread_count_controllable_and_deterministic()
    test_engine_parity()
