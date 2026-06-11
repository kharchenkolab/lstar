#!/usr/bin/env python3
"""Lazy open + streaming reduction + compression, on a realistically-sized measure.

Shows the three Python performance levers:
  1. lazy open      -- read(path, lazy=True) opens the store without materializing the heavy
                       arrays (constant, tiny memory) vs eager (whole matrix resident).
  2. streaming      -- stream_col_stats reduces a CSC measure by column block, so per-gene
                       variance/HVG stats run in bounded memory, never holding the full matrix.
  3. compression    -- gzip chunks shrink the store several-fold; chunking is what lets the
                       lazy path read only the touched blocks.

Run: python3 examples/lazy_streaming_demo.py
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


def dir_size(p):
    return sum(os.path.getsize(os.path.join(r, f)) for r, _, fs in os.walk(p) for f in fs)


def main():
    import numcodecs
    cells, genes, density = 20000, 8000, 0.05
    rng = np.random.default_rng(0)
    X = sp.csc_matrix(sp.random(cells, genes, density=density, format="csc", random_state=rng))
    X.data = X.data * 20 + 0.5
    ds = lstar.Dataset(kind="sample")
    ds.add_axis("cells", [f"c{i}" for i in range(cells)])
    ds.add_axis("genes", [f"g{i}" for i in range(genes)])
    ds.add_field("counts", X, role="measure", span=["cells", "genes"], state="raw")
    print(f"synthetic measure: {cells} x {genes} CSC, {X.nnz:,} nonzeros\n")

    plain = "/tmp/lazy_demo_plain.lstar.zarr"
    gz = "/tmp/lazy_demo_gzip.lstar.zarr"
    lstar.write(ds, plain, chunk_elems=1_000_000)
    lstar.write(ds, gz, chunk_elems=1_000_000, compressor=numcodecs.GZip(5))
    print(f"store size:  uncompressed {dir_size(plain)/1e6:6.1f} MB   "
          f"gzip-5 {dir_size(gz)/1e6:6.1f} MB   ({dir_size(plain)/dir_size(gz):.1f}x smaller)\n")

    # eager vs lazy open
    tracemalloc.start(); t = time.time(); de = lstar.read(gz); te = time.time() - t
    _, peak_e = tracemalloc.get_traced_memory(); tracemalloc.stop()
    tracemalloc.start(); t = time.time(); dl = lstar.read(gz, lazy=True); tl = time.time() - t
    _, peak_l = tracemalloc.get_traced_memory(); tracemalloc.stop()
    print(f"open eager:  {te:5.2f}s  peak +{peak_e/1e6:6.1f} MB   (full matrix resident)")
    print(f"open lazy:   {tl:5.2f}s  peak +{peak_l/1e6:6.1f} MB   {dl.fields['counts'].values!r}\n")

    # streaming per-gene stats (lognorm), constant memory, vs eager reference
    fe = de.fields["counts"].values
    fl = dl.fields["counts"].values
    tracemalloc.start(); t = time.time()
    m1, v1, n1 = stream_col_stats(fl, lognorm=True, block=512)
    ts = time.time() - t; _, peak_s = tracemalloc.get_traced_memory(); tracemalloc.stop()
    m2, v2, n2 = stream_col_stats(fe, lognorm=True, block=512)
    ok = np.allclose(m1, m2) and np.allclose(v1, v2) and np.array_equal(n1, n2)
    print(f"stream per-gene mean/var (log1p): {ts:5.2f}s  peak +{peak_s/1e6:6.1f} MB  "
          f"over {len(m1):,} genes  [matches eager: {ok}]")
    print(f"  -> reduces a {X.nnz:,}-nonzero measure without ever materializing it densely")


if __name__ == "__main__":
    main()
