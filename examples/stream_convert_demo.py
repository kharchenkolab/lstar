#!/usr/bin/env python3
"""Bounded-memory format conversion: stream large matrices between an .h5ad and an L* store
without ever loading them whole -- in *both* directions.

The ordinary conversion (`read_anndata(ad.read_h5ad(path))` then `lstar.write`) holds the whole
expression matrix in RAM. For a large atlas that can be many gigabytes. The streaming path instead
reads the source matrices straight from disk and writes them block-by-block, so peak memory is
bounded by one block -- the conversion runs on a laptop. The result is byte-identical to the eager
one; it's a memory-vs-speed trade, so it's opt-in.

This demo measures peak memory (via tracemalloc) for both:
  (1) h5ad -> L* store   via `lstar.convert_anndata`   (backed read + streamed write)
  (2) L* store -> h5ad   via `lstar.convert_to_h5ad`   (lazy read + streamed write)

Usage: python3 examples/stream_convert_demo.py [path/to/file.h5ad]
  With no argument, a synthetic h5ad is generated so the demo runs anywhere. The default real
  dataset is a local path on the author's machine (Tabula Muris Marrow, 77.6M nonzeros).
"""
import os
import sys
import time
import tracemalloc

import numpy as np
import scipy.sparse as sp

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "python", "src"))
import lstar
from lstar.profiles.anndata import read_anndata, convert_anndata

DEFAULT = "/home/pkharchenko/cacoa/age/tab.muris/" \
          "tabula-muris-senis-droplet-processed-official-annotations-Marrow.h5ad"


def main():
    import anndata as ad

    # Resolve an input .h5ad: the CLI arg, else the local default, else a generated synthetic one.
    path = sys.argv[1] if len(sys.argv) > 1 else (DEFAULT if os.path.exists(DEFAULT) else None)
    if path is None:
        rng = np.random.default_rng(0)
        X = sp.random(20000, 8000, density=0.05, format="csr", random_state=rng)
        X.data = np.round(X.data * 9 + 1).astype("f4")
        a = ad.AnnData(X=sp.csr_matrix(X))
        a.obs["leiden"] = np.array([f"c{i % 8}" for i in range(20000)])
        path = "/tmp/stream_demo.h5ad"
        a.write_h5ad(path)
        print("(no input given; generated a synthetic h5ad)")
    print(f"input: {path}  ({os.path.getsize(path)/1e6:.0f} MB on disk)\n")

    eager_store = "/tmp/convert_eager.lstar.zarr"
    stream_store = "/tmp/convert_stream.lstar.zarr"

    # --- eager conversion: load the whole AnnData, then write. Peak includes the full matrix. ---
    tracemalloc.start(); t = time.time()
    lstar.write(read_anndata(ad.read_h5ad(path)), eager_store)
    te = time.time() - t
    _, peak_e = tracemalloc.get_traced_memory(); tracemalloc.stop()
    print(f"eager      (load all, then write):  {te:5.1f}s   peak +{peak_e/1e6:7.1f} MB")

    # --- streaming conversion: backed read + block-by-block write. Peak ~ one block. ---
    tracemalloc.start(); t = time.time()
    convert_anndata(path, stream_store, chunk_elems=2_000_000)   # ~2M nonzeros per streamed block
    ts = time.time() - t
    _, peak_s = tracemalloc.get_traced_memory(); tracemalloc.stop()
    print(f"streaming  (convert_anndata):        {ts:5.1f}s   peak +{peak_s/1e6:7.1f} MB"
          f"   ({peak_e/peak_s:.0f}x less memory)")

    # --- confirm the two stores hold the identical expression matrix ---
    e, s = lstar.read(eager_store), lstar.read(stream_store)
    measures = [k for k, f in e.fields.items()
                if f.role == "measure" and len(f.span) == 2 and f.encoding in ("csr", "csc")]
    mk = measures[0]
    ev, svv = e.field(mk).values, s.field(mk).values
    same = ev.shape == svv.shape and (sp.csr_matrix(ev) != sp.csr_matrix(svv)).nnz == 0
    print(f"\nmeasure '{mk}': {ev.nnz:,} nonzeros; streamed store identical to eager = {same}")
    print("  -> (1) a large h5ad converts to L* in bounded memory; the data is unchanged.\n")

    # --- (2) the reverse: L* store -> h5ad, eager vs streamed ---
    from lstar.profiles.anndata import write_anndata
    eager_h5 = "/tmp/convert_eager.h5ad"
    stream_h5 = "/tmp/convert_stream.h5ad"

    # eager: materialize the whole store, build an AnnData, write it. Peak includes the full matrix.
    tracemalloc.start(); t = time.time()
    write_anndata(lstar.read(stream_store)).write_h5ad(eager_h5)
    te2 = time.time() - t
    _, peak_e2 = tracemalloc.get_traced_memory(); tracemalloc.stop()
    print(f"eager      (load store, then write h5ad): {te2:5.1f}s   peak +{peak_e2/1e6:7.1f} MB")

    # streaming: lazy read + block-by-block write into the h5ad's sparse groups. Peak ~ one block.
    tracemalloc.start(); t = time.time()
    lstar.convert_to_h5ad(stream_store, stream_h5, chunk_elems=2_000_000)
    ts2 = time.time() - t
    _, peak_s2 = tracemalloc.get_traced_memory(); tracemalloc.stop()
    print(f"streaming  (convert_to_h5ad):             {ts2:5.1f}s   peak +{peak_s2/1e6:7.1f} MB"
          f"   ({peak_e2/peak_s2:.0f}x less memory)")

    import anndata as ad2
    eh, sh = ad2.read_h5ad(eager_h5), ad2.read_h5ad(stream_h5)
    same2 = (sp.csr_matrix(eh.X) != sp.csr_matrix(sh.X)).nnz == 0
    print(f"\nh5ad X identical (eager vs streamed) = {same2}")
    print("  -> (2) a large L* store converts to h5ad in bounded memory; the data is unchanged.")


if __name__ == "__main__":
    main()
