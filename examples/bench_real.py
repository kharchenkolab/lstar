#!/usr/bin/env python3
"""Exercise the L* pipeline on a realistically-sized AnnData store.

Toy fixtures (100x50) prove correctness but say nothing about whether the
ingest/serialize/round-trip paths hold up at scale. This drives a real h5ad
through read_anndata -> zarr write -> zarr read -> validate and reports the
timings, the on-disk store size, and that the heavy fields survive byte-faithfully.

Usage:
    python3 examples/bench_real.py [path/to/file.h5ad]

Default dataset is Tabula Muris Senis (droplet) Marrow: 40,220 cells x 20,138 genes.
"""
import os
import sys
import time

import numpy as np
import scipy.sparse as sp

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "python", "src"))
import lstar
from lstar.profiles.anndata import read_anndata

DEFAULT = "/home/pkharchenko/cacoa/age/tab.muris/" \
          "tabula-muris-senis-droplet-processed-official-annotations-Marrow.h5ad"


def dir_size(path):
    total = 0
    for root, _, files in os.walk(path):
        for f in files:
            total += os.path.getsize(os.path.join(root, f))
    return total


def human(n):
    for unit in ("B", "KB", "MB", "GB"):
        if n < 1024:
            return f"{n:.1f} {unit}"
        n /= 1024
    return f"{n:.1f} TB"


def main():
    import anndata as ad

    path = sys.argv[1] if len(sys.argv) > 1 else DEFAULT
    out = "/tmp/bench_real.lstar.zarr"
    print(f"dataset: {path}\n  on-disk h5ad: {human(os.path.getsize(path))}")

    t = time.time()
    adata = ad.read_h5ad(path)
    print(f"  read_h5ad: {time.time()-t:.1f}s  X={adata.shape} {type(adata.X).__name__}")
    X = adata.X
    nnz = X.nnz if sp.issparse(X) else X.size
    xsum = float(X.sum())
    print(f"  X: nnz={nnz:,}  sum={xsum:.6g}  obsm={list(adata.obsm)}  obs cols={len(adata.obs.columns)}")

    t = time.time()
    ds = read_anndata(adata)
    print(f"\nread_anndata -> L*: {time.time()-t:.2f}s")
    print(f"  axes:   {', '.join(f'{a.name}({a.length})' for a in ds.axes.values())}")
    print(f"  fields: {', '.join(ds.fields)}")
    if ds.dropped:
        print(f"  dropped (recorded): {ds.dropped}")

    t = time.time()
    lstar.write(ds, out)
    wt = time.time() - t
    sz = dir_size(out)
    print(f"\nwrite L* zarr: {wt:.2f}s  store={human(sz)}  ({nnz*12/sz:.2f}x raw csr bytes/store)")

    t = time.time()
    ds2 = lstar.read(out)
    print(f"read L* zarr:  {time.time()-t:.2f}s")

    print("\nvalidate:", end=" ")
    errs = lstar.validate(ds2)
    print("clean" if not errs else f"{len(errs)} issue(s): {errs}")

    # Verify the heavy fields survived byte-faithfully.
    print("\nfidelity check:")
    f = ds2.fields["X"] if "X" in ds2.fields else ds2.fields[[k for k in ds2.fields if k != "X"][0]]
    measure = next(k for k, fl in ds2.fields.items() if fl.role == "measure")
    fl = ds2.fields[measure]
    v = fl.values
    got_nnz = v.data.size if hasattr(v, "data") and hasattr(v, "indptr") else np.asarray(v).size
    got_sum = float(v.sum())
    print(f"  measure '{measure}': nnz={got_nnz:,} (orig {nnz:,})  sum={got_sum:.6g} (orig {xsum:.6g})  "
          f"-> {'OK' if got_nnz == nnz and abs(got_sum-xsum) < 1e-3*abs(xsum) else 'MISMATCH'}")
    for name, fld in ds2.fields.items():
        if fld.role in ("embedding", "loading"):
            arr = np.asarray(fld.values)
            print(f"  {fld.role} '{name}': shape={arr.shape} dtype={arr.dtype}")
    print(f"\nstore: {out}  ({human(sz)})")


if __name__ == "__main__":
    main()
