#!/usr/bin/env bash
# Zarr v3 format on the R surface (R + zarr-python; no cmake C++ core needed). lstar's R binding reads
# both v2 and v3 (read auto-probes) and writes v3 via lstar_write(format="v3"). On the maximal store:
#   1. R reads the v2 seed, and writes it back as BOTH v2 and v3 (gzip).
#   2. R reads its own v3 store back (self round-trip; field/axis counts preserved).
#   3. zarr-python 3 confirms R's v3 is genuine v3 (zarr_format 3, inline-consolidated) and that R's v2
#      and v3 outputs are value-identical -- so the on-disk format is faithful, independent of R's own
#      (format-agnostic) nullable/graph semantics.
# Self-contained; needs a zarr>=3.1 python and the installed R package.
set -e
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
RLIB="$ROOT/.Rlib"
export PYTHONPATH="$ROOT/python/src"

SEED=/tmp/v3r_seed.lstar.zarr
RV2=/tmp/v3r_r_v2.lstar.zarr
RV3=/tmp/v3r_r_v3.lstar.zarr

python3 "$ROOT/conformance/v3_gen.py" "$SEED" | sed 's/^/  [py] /'

Rscript -e '.libPaths(c("'"$RLIB"'", .libPaths()))
suppressMessages(library(lstar))
ds <- lstar_read("'"$SEED"'")                                  # 1. R reads the v2 seed
lstar_write(ds, "'"$RV2"'", compression = "gzip", format = "v2")
lstar_write(ds, "'"$RV3"'", compression = "gzip", format = "v3")
ds3 <- lstar_read("'"$RV3"'")                                  # 2. R reads its own v3 back
stopifnot(length(ds3$fields) == length(ds$fields), length(ds3$axes) == length(ds$axes))
cat(sprintf("  [R ] read v2 seed (%d fields); wrote + re-read R v3 (%d fields)\n",
            length(ds$fields), length(ds3$fields)))' \
  2>&1 | grep -vE "^Warning|deprecat|masked|following object|Attaching|^$"

python3 "$ROOT/conformance/v3_verify.py" compare "$RV2" "$RV3"   # 3. R v2 == R v3 (format faithful)
python3 - "$RV3" <<'PY'
import sys, json
rj = json.load(open(sys.argv[1] + "/zarr.json"))
assert rj["zarr_format"] == 3 and "consolidated_metadata" in rj, "R v3 not genuine v3"
print("  [py] R-written v3 is genuine v3 (zarr_format 3, inline-consolidated)")
PY
echo "  v3 format on the R surface: read both + write v3 conformant"
