#!/usr/bin/env bash
# Backwards-compatibility guard: the CURRENT lstar must still read the COMMITTED golden stores to their
# FROZEN values, on every surface. `golden_v2` holds the bytes a pre-flip lstar produced (the v2 write
# path is unchanged by the v2->v3 default flip); `golden_v3` freezes today's v3 layout. expected.json is
# ground truth from the fixture's input arrays. Regenerate with fixtures/backcompat_gen.py (then commit).
set -e
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
FIX="$ROOT/conformance/fixtures"; EXP="$FIX/expected.json"
V2="$FIX/golden_v2.lstar.zarr"; V3="$FIX/golden_v3.lstar.zarr"
PY="${LSTAR_PY:-python3}"; export PYTHONPATH="$ROOT/python/src${PYTHONPATH:+:$PYTHONPATH}"
TMP="$(mktemp -d)"; trap 'rm -rf "$TMP"' EXIT

# 1) Python reads both committed stores -> the frozen manifest (the anchor + Python's own reader).
echo "  -- Python: golden v2 + v3 -> frozen manifest --"
"$PY" "$FIX/backcompat_check.py" "$V2" "$EXP"
"$PY" "$FIX/backcompat_check.py" "$V3" "$EXP"

# 1b) No consolidated metadata -> per-node fallback still reads (v2: drop .zmetadata; v3: strip the inline
#     consolidated_metadata from the root zarr.json). Guards the "walk the tree" path on an old store.
echo "  -- no-consolidated fallback --"
cp -r "$V2" "$TMP/v2nc.lstar.zarr"; rm -f "$TMP/v2nc.lstar.zarr/.zmetadata"
"$PY" "$FIX/backcompat_check.py" "$TMP/v2nc.lstar.zarr" "$EXP"
cp -r "$V3" "$TMP/v3nc.lstar.zarr"
"$PY" - "$TMP/v3nc.lstar.zarr" <<'PY'
import sys, json
p = sys.argv[1] + "/zarr.json"; d = json.load(open(p)); d.pop("consolidated_metadata", None)
json.dump(d, open(p, "w"))
PY
"$PY" "$FIX/backcompat_check.py" "$TMP/v3nc.lstar.zarr" "$EXP"

# 2) C++ (libzarr core): reads both frozen stores to identical values; a v2->v3 read-write round-trip on
#    the frozen v2 store is value-identical.
BIN="$ROOT/core/build/test_v3"
if [ -x "$BIN" ]; then
  echo "  -- C++ (libstar core) --"
  "$BIN" compare "$V2" "$V3" && echo "    [c++] golden v2 == v3 (reads both)"
  "$BIN" write "$V2" "$TMP/cpp_v2_to_v3.lstar.zarr" && echo "    [c++] golden v2 -> v3 read-back identical"
else
  echo "  [skip] test_v3 not built — C++ leg skipped"
fi

# 3) JS/WASM (libstar reader): reads both frozen stores byte-for-byte == the zarr-python reference.
EMSDK="${EMSDK:-$HOME/emsdk}"; NODE="$(ls -d "$EMSDK"/node/*/bin/node 2>/dev/null | head -1)"; NODE="${NODE:-$(command -v node || true)}"
if [ -n "$NODE" ] && [ -f "$ROOT/js/dist/lstar_io.mjs" ]; then
  echo "  -- JS/WASM (libstar reader) --"
  for tag_store in "v2:$V2" "v3:$V3"; do
    tag="${tag_store%%:*}"; store="${tag_store#*:}"
    "$PY" "$ROOT/js/test/io_dump.py" "$store" "$TMP/$tag.json" >/dev/null
    "$NODE" "$ROOT/js/test/io_parity.mjs" "$store" "$TMP/$tag.json" && echo "    [js] golden $tag == zarr-python (arrays + manifest)"
  done
else
  echo "  [skip] node/WASM absent — JS leg skipped"
fi

# 4) R (cpp11 binding): reads both frozen stores -> key values match the manifest.
RLIB="${LSTAR_RLIB:-$ROOT/.Rlib}"
if command -v Rscript >/dev/null 2>&1 && Rscript -e '.libPaths(c("'"$RLIB"'",.libPaths())); quit(status=!requireNamespace("lstar",quietly=TRUE))' >/dev/null 2>&1; then
  echo "  -- R (cpp11 binding) --"
  if Rscript -e '.libPaths(c("'"$RLIB"'",.libPaths())); suppressMessages({library(lstar); library(Matrix)})
      exp <- jsonlite::fromJSON("'"$EXP"'")
      for (store in c("'"$V2"'","'"$V3"'")) {
        ds <- lstar_read(store)
        cnt <- as(field_value(ds,"counts"),"CsparseMatrix"); lei <- as.character(field_value(ds,"leiden"))
        stopifnot(length(cnt@x)==exp$counts$nnz, abs(sum(cnt@x)-exp$counts$sum)<1e-3, identical(lei, exp$leiden))
        cat(sprintf("    [R ] %s: counts nnz+sum + leiden match\n", basename(store)))
      }' >"$TMP/bc_r.log" 2>&1; then
    grep -vE "^Attaching|masked|following object" "$TMP/bc_r.log"
  else echo "  FAIL: R backcompat"; tail -8 "$TMP/bc_r.log"; exit 1; fi
else
  echo "  [skip] Rscript/lstar absent — R leg skipped"
fi
echo "backcompat conformance PASSED."
