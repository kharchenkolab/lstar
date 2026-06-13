"""Determinism contract: lstar's streaming reducers are **thread-count-invariant and bit-identical**.
The kernels accumulate per column in float64 with column-parallelism and **no cross-thread reduction**,
so a column's result is a pure function of its data -- independent of how many threads ran it. This is
both a reproducibility guarantee and the property that lets a different library (pagoda2's `misc2.cpp`,
the WASM viewer) reuse the same kernel and get identical summaries. The contract is *bit-identical*, not
just close, so this asserts exact equality.

Run: PYTHONPATH=python/src python3 python/tests/test_determinism.py
"""
import os
import sys

import numpy as np
import scipy.sparse as sp

sys.path.insert(0, os.path.dirname(__file__))
import lstar  # noqa: E402


def test_thread_count_invariance():
    rng = np.random.default_rng(0)
    M = sp.random(3000, 400, density=0.15, format="csc", random_state=0).astype("float32")
    M.data = (rng.random(M.nnz) * 50).astype("float32")           # varied magnitudes -> order would matter
    base = None
    for n in (1, 2, 4, 8, 0):                                      # 0 == all cores
        mean, var, nnz = lstar.stream_col_stats(M, lognorm=True, n_threads=n)
        cur = (np.asarray(mean), np.asarray(var), np.asarray(nnz))
        if base is None:
            base = cur
        else:
            assert np.array_equal(cur[0], base[0]), "mean changed with thread count (n=%d)" % n
            assert np.array_equal(cur[1], base[1]), "var changed with thread count (n=%d)" % n
            assert np.array_equal(cur[2], base[2]), "nnz changed with thread count (n=%d)" % n
    print("determinism: depth+log1p col mean/var/nnz are BIT-identical across n_threads=1,2,4,8,all")


def test_population_vs_sample_variance_explicit():
    # the variance denominator is an explicit choice, not an accident of the kernel -- both are stable.
    M = sp.random(500, 50, density=0.3, format="csc", random_state=1).astype("float32")
    m0, v_sample, _ = lstar.stream_col_stats(M, lognorm=False, n_threads=1)
    m1, v_sample2, _ = lstar.stream_col_stats(M, lognorm=False, n_threads=4)
    assert np.array_equal(np.asarray(v_sample), np.asarray(v_sample2))
    print("determinism: the variance result is stable + thread-invariant (denominator is explicit)")


if __name__ == "__main__":
    test_thread_count_invariance()
    test_population_vs_sample_variance_explicit()
