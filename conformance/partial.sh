#!/usr/bin/env bash
# Partial coverage round-trips through the R path, including when `index_axis` is a **derived union**
# axis of a collection (the pagoda2/Conos case: a facet measured on a subset of the joint union cells).
# Py -> store -> R -> store -> Py preserves coverage/index/index_axis and the partial values.
set -eo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
RLIB="$ROOT/.Rlib"
S0=/tmp/partc0.lstar.zarr
S1=/tmp/partc1.lstar.zarr

PYTHONPATH="$ROOT/python/src" python3 - "$S0" <<'PY'
import sys, warnings; warnings.filterwarnings("ignore")
import numpy as np, scipy.sparse as sp, lstar
ds = lstar.Dataset(kind="collection")
ds.add_axis("cells", [f"u{i}" for i in range(10)], origin="derived", role="observation")  # DERIVED union
ds.add_axis("genes", [f"g{i}" for i in range(4)])
idx = np.array([0, 2, 4, 6, 8])                              # a facet on 5 of the 10 joint cells
ds.add_field("adt", sp.random(5, 4, density=0.5, format="csc").astype("float32"),
             role="measure", span=["cells", "genes"], state="raw", index=idx, index_axis="cells")
assert not lstar.validate(ds)
lstar.write(ds, sys.argv[1])
print("  [py] wrote a partial measure over a DERIVED union cells axis (5 of 10)")
PY

Rscript -e '.libPaths(c("'"$RLIB"'", .libPaths())); suppressMessages(library(lstar))
ds <- lstar_read("'"$S0"'"); f <- ds$fields[["adt"]]
stopifnot(identical(f$coverage, "partial"), identical(f$index_axis, "cells"),
          length(f$index) == 5, identical(ds$axes$cells$origin, "derived"), nrow(f$values) == 5)
lstar_write(ds, "'"$S1"'")
cat("  [R ] read partial coverage over the derived axis + rewrote\n")' 2>&1 | grep -E "^  \[R"

PYTHONPATH="$ROOT/python/src" python3 - "$S1" <<'PY'
import sys, warnings; warnings.filterwarnings("ignore")
import numpy as np, lstar
ds = lstar.read(sys.argv[1]); f = ds.field("adt")
assert f.coverage == "partial" and f.index_axis == "cells"
assert list(np.asarray(f.index)) == [0, 2, 4, 6, 8] and f.values.shape == (5, 4)
assert ds.axes["cells"].origin == "derived" and not lstar.validate(ds)
print("  [py] partial coverage over a derived union axis survived Py -> R -> Py")
PY
echo "partial-coverage conformance PASSED."
