#!/usr/bin/env python3
"""Lazy / streaming / compression on REAL data, not synthetic matrices.

Drives the Tabula Muris Senis (droplet) Marrow measure -- 40,220 cells x 20,138 genes,
~77.6M nonzeros -- through the three performance levers and emits a chunked+gzip store that the
C++ core (core/test/test_chunked) then reads. Reports real compression ratios (real counts, not
incompressible random floats), lazy-open memory, and streamed per-gene stats checked against a
full in-memory scipy computation.

Usage: python3 examples/real_perf.py [path/to/file.h5ad]
"""
import os
import sys
import time
import tracemalloc

import numpy as np
import scipy.sparse as sp

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "python", "src"))
import lstar
from lstar.lazy import stream_col_stats

DEFAULT = "/home/pkharchenko/cacoa/age/tab.muris/" \
          "tabula-muris-senis-droplet-processed-official-annotations-Marrow.h5ad"
OUT_GZIP = "/tmp/real_perf_gzip.lstar.zarr"   # consumed by test_chunked


def dir_size(p):
    return sum(os.path.getsize(os.path.join(r, f)) for r, _, fs in os.walk(p) for f in fs)


def main():
    import anndata as ad
    path = sys.argv[1] if len(sys.argv) > 1 else DEFAULT
    t = time.time()
    adata = ad.read_h5ad(path)
    X = sp.csc_matrix(adata.X)                 # gene-compressed for per-gene streaming
    pca = np.asarray(adata.obsm["X_pca"], dtype="f4")
    print(f"loaded {os.path.basename(path)} in {time.time()-t:.1f}s: "
          f"{X.shape[0]:,} cells x {X.shape[1]:,} genes, {X.nnz:,} nonzeros\n")

    ds = lstar.Dataset(kind="sample")
    ds.add_axis("cells", np.asarray(adata.obs_names, dtype=str))
    ds.add_axis("genes", np.asarray(adata.var_names, dtype=str))
    ds.add_axis("pca", [f"PC{i}" for i in range(pca.shape[1])], role="coordinate")
    ds.add_field("counts", X, role="measure", span=["cells", "genes"], state="raw")
    ds.add_field("pca", pca, role="embedding", span=["cells", "pca"])

    # ---- compression on real counts (vs incompressible random) ----
    import numcodecs
    single = "/tmp/real_perf_single.lstar.zarr"
    chunked = "/tmp/real_perf_chunked.lstar.zarr"
    for label, p, kw in [("single-chunk uncompressed", single, {}),
                         ("chunked(4M) uncompressed", chunked, dict(chunk_elems=4_000_000)),
                         ("chunked(4M) + gzip5", OUT_GZIP,
                          dict(chunk_elems=4_000_000, compressor=numcodecs.GZip(5)))]:
        t = time.time(); lstar.write(ds, p, **kw); wt = time.time() - t
        print(f"  write {label:28s} {wt:6.1f}s   store {dir_size(p)/1e6:7.1f} MB")
    print(f"  -> gzip ratio on real counts: "
          f"{dir_size(single)/dir_size(OUT_GZIP):.1f}x\n")

    # ---- lazy open vs eager, on the real gzip store ----
    tracemalloc.start(); t = time.time(); de = lstar.read(OUT_GZIP); te = time.time() - t
    _, pe = tracemalloc.get_traced_memory(); tracemalloc.stop()
    tracemalloc.start(); t = time.time(); dl = lstar.read(OUT_GZIP, lazy=True); tl = time.time() - t
    _, pl = tracemalloc.get_traced_memory(); tracemalloc.stop()
    print(f"  open eager: {te:5.2f}s  peak +{pe/1e6:7.1f} MB")
    print(f"  open lazy:  {tl:5.2f}s  peak +{pl/1e6:7.1f} MB   {dl.fields['counts'].values!r}\n")

    # ---- streamed per-gene stats, serial vs threaded, vs a stable dense ground truth ----
    fl = dl.fields["counts"].values
    tracemalloc.start(); t = time.time()
    m1, v1, n1 = stream_col_stats(fl, lognorm=True, block=1024, n_threads=1)
    ts = time.time() - t; _, ps = tracemalloc.get_traced_memory(); tracemalloc.stop()
    t = time.time(); m8, v8, n8 = stream_col_stats(fl, lognorm=True, block=1024, n_threads=8)
    tp = time.time() - t
    same = np.allclose(m1, m8) and np.allclose(v1, v8) and np.array_equal(n1, n8)
    # independent ground truth: dense np.var(ddof=1) on a gene block (log1p). The measure is
    # stored float32 and the lean stream keeps it float32 (double accumulation), so agreement is
    # to ~float32 precision -- not a defect, that is the stored precision.
    gb = slice(0, 2000)
    D = np.log1p(X[:, gb].toarray().astype(np.float64))
    ok = (np.allclose(m1[gb], D.mean(0), rtol=1e-5, atol=1e-7) and
          np.allclose(v1[gb], D.var(0, ddof=1), rtol=1e-4, atol=1e-7))
    print(f"  stream per-gene mean/var (log1p): {ts:5.2f}s (1 thread) / {tp:5.2f}s (8 threads, "
          f"{ts/tp:.1f}x)  peak +{ps/1e6:6.1f} MB over {len(m1):,} genes")
    print(f"  matches dense np.var ground truth: {ok}; threaded==serial: {same}")
    top = np.argsort(v1)[::-1][:5]
    print(f"  top-5 most-variable genes: {[adata.var_names[i] for i in top]}")
    print(f"\n  chunked+gzip store for the C++ reader: {OUT_GZIP}  ({dir_size(OUT_GZIP)/1e6:.0f} MB)")
    print(f"  expected: nnz={X.nnz} sum={float(X.data.sum()):.6g}")


if __name__ == "__main__":
    main()
