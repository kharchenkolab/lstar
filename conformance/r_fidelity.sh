#!/usr/bin/env bash
# R-path fidelity: a float32 value dtype and a graph relation's directed/weighted flags must survive a
# Python -> R (lstar_read/lstar_write) -> Python round-trip. Before, the R binding widened every non-raw
# measure f4->f8 (audit T2.2) and never surfaced directed/weighted, so they were dropped (audit T1.4).
set -e
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
export PYTHONPATH="$ROOT/python/src${PYTHONPATH:+:$PYTHONPATH}"
export LSTAR_RLIB="${LSTAR_RLIB:-$ROOT/.Rlib}"
T="$(mktemp -d)"; trap 'rm -rf "$T"' EXIT

python3 - "$T/s0.lstar.zarr" <<'PY'
import sys, numpy as np, scipy.sparse as sp, lstar
ds = lstar.Dataset(kind="sample")
ds.add_axis("cells", [f"c{i}" for i in range(30)]); ds.add_axis("genes", [f"g{j}" for j in range(10)])
X = sp.random(30, 10, density=0.4, format="csc", random_state=0); X.data = X.data.astype(np.float32)
ds.add_field("logX", X, role="measure", span=["cells", "genes"], state="lognorm")     # float32 measure
G = sp.random(30, 30, density=0.2, format="csr", random_state=1); G.data = G.data.astype(np.float32)
ds.add_field("knn", G, role="relation", span=["cells", "cells"], directed=True, weighted=True)
lstar.write(ds, sys.argv[1]); print("  [py] wrote a float32 measure + a directed/weighted graph")
PY

Rscript -e '.libPaths(c(Sys.getenv("LSTAR_RLIB"), .libPaths())); suppressMessages(library(lstar))
a <- commandArgs(TRUE); lstar_write(lstar_read(a[1]), a[2]); cat("  [R ] lstar_read -> lstar_write\n")' \
  "$T/s0.lstar.zarr" "$T/s1.lstar.zarr"

python3 - "$T/s1.lstar.zarr" <<'PY'
import sys, lstar
d = lstar.read(sys.argv[1])
dt = str(d.field("logX").values.dtype)
assert dt == "float32", "T2.2: value dtype not preserved through R (got %s)" % dt
assert d.field("knn").directed is True and d.field("knn").weighted is True, "T1.4: directed/weighted lost through R"
print("  [py] float32 dtype + directed/weighted survived the R round-trip")
PY
echo "R-path fidelity PASSED."
