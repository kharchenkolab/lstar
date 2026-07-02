#!/usr/bin/env python3
"""Reproduce / regression-gate the macOS arm64 dual-libomp crash in extend_for_viewer.

lstar's wheel bundles its own libomp. If another installed package (scikit-learn, numba, some scipy
builds) also loads a libomp, two llvm-openmp runtimes coexist in one process; on macOS arm64 the __kmp_*
worker/suspend machinery cross-binds between the copies and extend_for_viewer's OpenMP kernels crash with
EXC_BAD_ACCESS (SIGSEGV). See build_tools/cibw_macos_repair.sh for the mechanism and the packaging fix.

This script FORCES the two-runtime condition (imports scikit-learn if present, then spins up BLAS/OpenMP)
and runs the exact kernel path that crashed. Exit 0 = no crash. A crash is an uncatchable SIGSEGV, so the
process dies with a nonzero signal and any CI job running this fails -- which is the regression gate.

Standalone (e.g. on the reporting Mac):
    pip install lstar-sc scikit-learn        # or: pip install <branch wheel> scikit-learn
    python build_tools/macos_libomp_repro.py
Expected: SIGSEGV on an unfixed wheel; "OK" on a fixed one.
"""
import sys

# 1) Load a SECOND OpenMP runtime first, the way a real scientific env does. scikit-learn's macOS wheels
#    bundle their own libomp; importing it recreates the two-runtime process the bug needs. Absent sklearn
#    this degrades to a single-runtime smoke test (still worth running -- it exercises the kernels).
_second = []
for mod in ("sklearn.linear_model", "sklearn.cluster", "scipy.linalg"):
    try:
        __import__(mod)
        _second.append(mod)
    except Exception:
        pass

import numpy as np
# actually resolve + start the BLAS/OpenMP thread pool (mapping the dylib isn't enough to trigger the bug)
np.linalg.svd(np.random.default_rng(0).standard_normal((96, 96)), full_matrices=False)

import scipy.sparse as sp
import lstar


def _synthetic(nc=400, ng=150, K=6, seed=0):
    rng = np.random.default_rng(seed)
    X = sp.random(nc, ng, density=0.2, format="csc", random_state=seed)
    X.data = (rng.poisson(3, size=X.data.shape) + 1).astype(np.int32)   # raw integer counts
    ds = lstar.Dataset(kind="sample")
    ds.add_axis("cells", [f"c{i}" for i in range(nc)])
    ds.add_axis("genes", [f"g{j}" for j in range(ng)])
    ds.add_axis("umap", ["umap1", "umap2"])
    ds.add_field("counts", X, role="measure", span=["cells", "genes"], state="raw")
    ds.add_field("umap", rng.normal(size=(nc, 2)), role="embedding", span=["cells", "umap"])
    ds.add_field("leiden", rng.integers(0, K, size=nc).astype(str), role="label", span=["cells"])
    return ds


def main():
    print("lstar %s | accel=%s | second OpenMP runtime: %s"
          % (lstar.__version__, lstar.has_accel(), _second or "(none found)"))
    if not lstar.has_accel():
        print("WARNING: compiled accelerator absent -- OpenMP path not exercised (repro is a no-op)")
    # the cross-binding crash is load/schedule-order sensitive; repeat to expose it reliably.
    for i in range(6):
        ds = _synthetic(seed=i)
        lstar.extend_for_viewer(ds)               # raw counts + grouping -> the OpenMP kernels
        assert "od_score" in ds.fields and "stats_leiden_sum" in ds.fields
    print("extend_for_viewer OK x6 -- no dual-libomp crash")
    return 0


if __name__ == "__main__":
    sys.exit(main())
