#!/usr/bin/env bash
# Nullable / extension-dtype conformance: an explicit `uint8` validity mask (1 == missing) beside a
# field's values must round-trip across languages so nullable Int/boolean/string survive with their
# value-vs-missing distinction intact (distinct from categorical -1 and from float NaN). Python writes a
# masked integer + a masked string; the C++ core (via R) reconstructs NA-bearing vectors and re-writes;
# Python reads the mask + values back byte-identical.
set -e
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
RLIB="$ROOT/.Rlib"
PY=/tmp/nullable_py.lstar.zarr
R=/tmp/nullable_r.lstar.zarr

PYTHONPATH="$ROOT/python/src" python3 - "$PY" <<'PY'
import sys, warnings; warnings.filterwarnings("ignore")
import numpy as np, lstar
ds = lstar.Dataset(kind="sample")
ds.add_axis("cells", [f"c{i}" for i in range(5)])
ds.add_field("n_counts", np.array([10, 0, 7, 0, 3], dtype=np.int64), role="measure", span=["cells"],
             mask=np.array([0, 1, 0, 1, 0], dtype=np.uint8))               # positions 1,3 missing
ds.add_field("donor", np.array(["d1", "", "d2", "", "d1"], dtype=str), role="label", span=["cells"],
             mask=np.array([0, 1, 0, 1, 0], dtype=np.uint8))
assert not lstar.validate(ds)
lstar.write(ds, sys.argv[1])
print("  [py] wrote masked integer + masked string")
PY

Rscript -e '.libPaths(c("'"$RLIB"'", .libPaths())); suppressMessages(library(lstar))
ds <- lstar_read("'"$PY"'")                              # C++ core reconstructs NA-bearing vectors
n <- field_value(ds, "n_counts"); d <- field_value(ds, "donor")
stopifnot(is.na(n[2]), is.na(n[4]), !is.na(n[1]), n[1] == 10, n[3] == 7,
          is.na(d[2]), is.na(d[4]), identical(d[c(1,3,5)], c("d1","d2","d1")))
lstar_write(ds, "'"$R"'")                                # split NA back into values + mask on write
cat("  [R ] C++/R read mask as NA (integer + character); re-wrote with mask\n")' \
  2>&1 | grep -vE "^Warning|deprecat|masked|following object|Attaching|^$"

PYTHONPATH="$ROOT/python/src" python3 - "$R" <<'PY'
import sys, warnings; warnings.filterwarnings("ignore")
import numpy as np, lstar
ds = lstar.read(sys.argv[1])
n, d = ds.field("n_counts"), ds.field("donor")
assert n.mask is not None and list(n.mask) == [0, 1, 0, 1, 0]
assert list(np.asarray(n.values)[[0, 2, 4]]) == [10, 7, 3]                 # values (numerically) intact
assert d.mask is not None and list(d.mask) == [0, 1, 0, 1, 0]
assert list(np.asarray(d.values)[[0, 2, 4]]) == ["d1", "d2", "d1"]
assert not lstar.validate(ds)
print("  [py] read R-written masks back: integer + string validity masks byte-identical")
PY
