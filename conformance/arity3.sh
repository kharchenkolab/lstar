#!/usr/bin/env bash
# Arity-3 fields round-trip across languages. A 3-axis tensor (a cell-cell communication
# sender×receiver×lr_pair score) must survive Py -> store -> R -> store -> Py with its 3rd axis intact
# and values exact. The R dense reconstruction used to collapse arity-3 to 2-D (drop the 3rd axis); this
# guards the n-D fix in `.lstar_assemble` + the write payload.
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
echo "arity-3 conformance PASSED."
