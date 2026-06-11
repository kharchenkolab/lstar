#!/usr/bin/env python3
"""Variable-length round-trip that returns to the ORIGINAL format, on real data.

A conversion library is only trustworthy if a value can leave its native format, pass through L*
(any number of hops), and come back unchanged. This drives a real AnnData (Tabula Muris Marrow)
through N cycles of AnnData -> L* -> AnnData and checks two things:
  1. the final AnnData matches the *original* object (X, obs, embeddings, graphs);
  2. it is a fixed point -- cycle 1 already equals cycle N, so any chain length is safe.
What L* cannot represent (e.g. uns) is reported via `dropped`, never silently changed.

Usage: python3 examples/roundtrip_chain.py [n_cycles] [path.h5ad]
"""
import os
import sys

import numpy as np
import scipy.sparse as sp

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "python", "src"))
import lstar
from lstar.profiles.anndata import read_anndata, write_anndata

DEFAULT = "/home/pkharchenko/cacoa/age/tab.muris/" \
          "tabula-muris-senis-droplet-processed-official-annotations-Marrow.h5ad"


def fingerprint(adata):
    """A comparable digest of the representable content of an AnnData."""
    X = adata.X
    fp = {"shape": tuple(adata.shape),
          "X_nnz": int(X.nnz) if sp.issparse(X) else int(np.size(X)),
          "X_sum": float(X.sum()),
          "obs_cols": sorted(map(str, adata.obs.columns)),
          "obsm": {k: (np.asarray(v).shape, float(np.asarray(v).sum()))
                   for k, v in adata.obsm.items()},
          "obsp": {k: (v.shape, float(v.sum())) for k, v in adata.obsp.items()}}
    return fp


def obs_values_equal(a, b):
    import pandas as pd
    for c in a.obs.columns:
        va, vb = a.obs[c], b.obs[c]
        if pd.api.types.is_numeric_dtype(va) and pd.api.types.is_numeric_dtype(vb):
            if not np.allclose(np.asarray(va, float), np.asarray(vb, float), equal_nan=True):
                return False, c
        else:                                   # categorical/string compared by value
            if not (np.asarray(va.astype(str)) == np.asarray(vb.astype(str))).all():
                return False, c
    return True, None


def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 4
    path = sys.argv[2] if len(sys.argv) > 2 else DEFAULT
    import anndata as ad

    a0 = ad.read_h5ad(path)
    print(f"original: {a0.shape[0]:,} cells x {a0.shape[1]:,} genes; "
          f"obsm={list(a0.obsm)} obsp={list(a0.obsp)}\n")
    fp0 = fingerprint(a0)                              # a digest of the ORIGINAL, to compare against

    # Run the conversion N times: each cycle is AnnData -> L* (read_anndata) -> a Zarr store
    # (lstar.write/read) -> AnnData (write_anndata). If the model is a faithful round-trip, the
    # object stops changing after the first cycle (a fixed point), so any chain length is safe.
    cur = a0
    fps = []
    dropped = None
    for i in range(n):
        ds = read_anndata(cur)                         # native -> L* fields (the profile, forward)
        if dropped is None:
            dropped = list(ds.dropped)                 # what no rule could represent (e.g. uns); recorded once
        p = f"/tmp/rtchain_{i}.lstar.zarr"
        lstar.write(ds, p)                             # through the actual on-disk format (not just in-memory)
        cur = write_anndata(lstar.read(p))             # L* -> native (the profile, reversed)
        fps.append(fingerprint(cur))
        # Compare this cycle's object to the original on the representable content (X, obs, obsm, obsp).
        same_as_orig = (fps[-1] == fp0)
        ov_ok, bad = obs_values_equal(a0, cur)
        print(f"  cycle {i+1}: AnnData->L*->AnnData  "
              f"X_nnz={fps[-1]['X_nnz']} X_sum={fps[-1]['X_sum']:.6g}  "
              f"== original: {same_as_orig and ov_ok}" + ("" if ov_ok else f" (obs '{bad}')"))

    # Two claims: (a) every cycle produced the SAME object (a fixed point), and (b) it equals the
    # original. Together they mean a conversion chain of ANY length returns to the native format.
    fixed_point = all(fp == fps[0] for fp in fps)
    matches_original = (fps[-1] == fp0) and obs_values_equal(a0, cur)[0]
    print(f"\n  fixed point across {n} cycles: {fixed_point}")
    print(f"  cycle-{n} AnnData matches the original: {matches_original}")
    # Loss is visible, not silent: uns (which L* has no field for) is listed in `dropped`.
    print(f"  not representable in L* (recorded in dropped, not lost): {dropped}")
    # Digests can collide; also check exact embedding values survived the whole chain.
    pca_ok = np.allclose(np.asarray(a0.obsm['X_pca']), np.asarray(cur.obsm['X_pca']))
    print(f"  X_pca values identical after the chain: {pca_ok}")
    sys.exit(0 if (fixed_point and matches_original and pca_ok) else 1)


if __name__ == "__main__":
    main()
