#!/usr/bin/env bash
# Arity-3 fields round-trip across languages. A 3-axis tensor (a cell-cell communication
# sender×receiver×lr_pair score) must survive Py -> store -> R -> store -> Py with its 3rd axis intact
# and values exact. The R dense reconstruction used to collapse arity-3 to 2-D (drop the 3rd axis); this
# guards the n-D fix in `.lstar_assemble` + the write payload.
# Origin coverage: Py-authored ✓ | R-authored ✓ (each cross-read by the other language) — see conformance/README.md
set -eo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
RLIB="$ROOT/.Rlib"
S0=/tmp/a3_0.lstar.zarr
S1=/tmp/a3_1.lstar.zarr

PYTHONPATH="$ROOT/python/src" python3 - "$S0" <<'PY'
import sys, warnings; warnings.filterwarnings("ignore")
import numpy as np, lstar
ds = lstar.Dataset(kind="sample")
ds.add_axis("senders", ["T", "B", "Mono", "NK"]); ds.add_axis("receivers", ["T", "B", "Mono", "NK"])
ds.add_axis("lr_pairs", ["L1_R1", "L2_R2", "L3_R3"])
T = (np.arange(4 * 4 * 3, dtype="float32") + 0.5).reshape(4, 4, 3)   # distinct values -> order matters
ds.add_field("ccc", T, role="measure", span=["senders", "receivers", "lr_pairs"], subtype="communication")
assert not lstar.validate(ds)
lstar.write(ds, sys.argv[1])
print("  [py] wrote a (4×4×3) communication tensor")
PY

Rscript -e '.libPaths(c("'"$RLIB"'", .libPaths())); suppressMessages(library(lstar))
ds <- lstar_read("'"$S0"'"); f <- ds$fields[["ccc"]]
stopifnot(identical(dim(f$values), c(4L, 4L, 3L)),                       # 3rd axis NOT collapsed
          identical(as.character(f$span), c("senders","receivers","lr_pairs")))
lstar_write(ds, "'"$S1"'")
cat("  [R ] read arity-3 (4×4×3), 3rd axis intact, rewrote\n")' 2>&1 | grep -E "^  \[R"

PYTHONPATH="$ROOT/python/src" python3 - "$S1" <<'PY'
import sys, warnings; warnings.filterwarnings("ignore")
import numpy as np, lstar
T = (np.arange(4 * 4 * 3, dtype="float32") + 0.5).reshape(4, 4, 3)
f = lstar.read(sys.argv[1]).field("ccc")
v = np.asarray(f.values)
assert v.shape == (4, 4, 3) and np.allclose(v, T), (v.shape,)
assert f.span == ["senders", "receivers", "lr_pairs"]
print("  [py] arity-3 tensor survived Py -> R -> Py: shape + values exact")
PY
# R-AUTHORED arity-3 (origin coverage): R builds a 3-D tensor from scratch -> Python reads it value-equal.
# The Py-origin flow above never exercises the R *writer's* n-D path on an R-native array (only as a rewriter).
S2=/tmp/a3_2.lstar.zarr
Rscript -e '.libPaths(c("'"$RLIB"'", .libPaths())); suppressMessages(library(lstar))
T <- array(as.numeric(1:24), dim=c(2L,3L,4L))                          # col-major fill -> distinct values
ds <- structure(list(kind="sample", spec_version="0.1", profiles=character(0), dropped=character(0),
  axes=list(a=list(labels=paste0("a",1:2),origin="observed",role="feature"),
            b=list(labels=paste0("b",1:3),origin="observed",role="feature"),
            c=list(labels=paste0("c",1:4),origin="observed",role="feature")),
  fields=list(tns=list(values=T, role="measure", span=c("a","b","c"), encoding="dense"))),
  class="lstar_dataset")
lstar_write(ds, "'"$S2"'")
rb <- field_value(lstar_read("'"$S2"'"), "tns")
stopifnot(identical(dim(rb), c(2L,3L,4L)), rb[2,3,4]==24, rb[1,1,1]==1)
cat("  [R ] authored a (2×3×4) tensor from scratch, rewrote\n")' 2>&1 | grep -E "^  \[R"

PYTHONPATH="$ROOT/python/src" python3 - "$S2" <<'PY'
import sys, warnings; warnings.filterwarnings("ignore")
import numpy as np, lstar
v = np.asarray(lstar.read(sys.argv[1]).field("tns").values)
assert v.shape == (2, 3, 4), v.shape
# R col-major fill: array[i,j,k] = 1 + i + 2j + 6k (0-based) -> these exact cells
assert v[0,0,0]==1 and v[1,0,0]==2 and v[0,1,0]==3 and v[0,0,1]==7 and v[1,2,3]==24, v
print("  [py] read R-authored arity-3 tensor: shape + values exact (R-origin cross-read)")
PY
echo "arity-3 conformance PASSED."
