#!/usr/bin/env bash
# Cross-surface parity of extend_for_viewer(primary=): the caller-declared "grouping the viewer opens on"
# must key the counts_cellmajor locality reorder IDENTICALLY on Python, R, and JS. Builds ONE multi-grouping
# base (default detected key = leiden; also louvain/annotation/phase), preps each surface with primary=phase
# (a NON-default grouping), and asserts (a) the preps agree on ALL viewer fields (equiv — incl. the reorder
# permutation) and (b) each surface keyed the reorder on `phase` (provgroup), while the DEFAULT keys
# elsewhere — so `primary` demonstrably took effect (not a vacuous pass). Skips a surface cleanly when its
# runtime is absent (R: no Rscript; JS: no built WASM dist / node / zarrita).
set -u
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
PY="${LSTAR_PY:-python3}"
export PYTHONPATH="$ROOT/python/src:$ROOT/python/tests${PYTHONPATH:+:$PYTHONPATH}"
export LSTAR_RLIB="${LSTAR_RLIB:-$ROOT/.Rlib}"
CHK="$ROOT/conformance/viewer_check.py"
TMP="$(mktemp -d)"; trap 'rm -rf "$TMP"' EXIT
fail=0
PRIMARY=phase

if ! "$PY" -c "import lstar" 2>/dev/null; then
  echo "  [skip] lstar not importable by '$PY' (set LSTAR_PY) — skipping primary= parity"; exit 0; fi

# ONE base store with competing groupings (default detected key = leiden; also louvain/annotation/phase)
"$PY" "$CHK" make-base "$TMP/base.lstar.zarr" csc multi >/dev/null || { echo "  FAIL: make-base"; exit 1; }

# Python: a DEFAULT prep (to prove primary changes the key) + the primary=phase prep (the reference)
"$PY" "$CHK" prep-lstar "$TMP/base.lstar.zarr" "$TMP/py_default.zarr" "" ""       >/dev/null || fail=1
"$PY" "$CHK" prep-lstar "$TMP/base.lstar.zarr" "$TMP/py.zarr"         "" "$PRIMARY" >/dev/null || fail=1
defkey="$("$PY" "$CHK" reorder-key "$TMP/py_default.zarr")"
echo "  default reorder key = '$defkey' ; primary = '$PRIMARY'"
[ "$defkey" = "$PRIMARY" ] && { echo "  FAIL: default key == primary — leg would be vacuous"; fail=1; }
"$PY" "$CHK" provgroup "$TMP/py.zarr" "$PRIMARY" || fail=1

# R: prep with primary=phase; must equal the Python primary prep on ALL viewer fields + key on phase
if command -v Rscript >/dev/null 2>&1; then
  if Rscript "$ROOT/conformance/viewer_extend_r.R" "$TMP/base.lstar.zarr" "$TMP/r.zarr" "" "$PRIMARY" >/tmp/lstar_vp_r.log 2>&1; then
    echo "  -- python == R (primary=$PRIMARY) --"
    "$PY" "$CHK" equiv "$TMP/py.zarr" "$TMP/r.zarr"     || fail=1
    "$PY" "$CHK" provgroup "$TMP/r.zarr" "$PRIMARY"     || fail=1
  else echo "  FAIL: R extend_for_viewer(primary) errored"; tail -12 /tmp/lstar_vp_r.log; fail=1; fi
else echo "  [skip] no Rscript — skipping R primary= parity"; fi

# JS: prep with --primary phase (in place on a copy of the base); must equal the Python primary prep + key on phase
NODE="$(ls -d "${EMSDK:-$HOME/emsdk}"/node/*/bin/node 2>/dev/null | head -1)"; NODE="${NODE:-node}"
DIST="$ROOT/js/dist/lstar_kernels.mjs"
# require a CURRENT dist (newer than its WASM + core sources) so a stale local build skips cleanly rather
# than erroring on a missing kernel; CI's js.sh builds a fresh one earlier in the js-wasm job.
if [ "$DIST" -nt "$ROOT/js/wasm/lstar_wasm.cpp" ] 2>/dev/null && [ "$DIST" -nt "$ROOT/core/include/lstar/lstar.hpp" ] 2>/dev/null \
   && [ -d "$ROOT/js/node_modules/zarrita" ] && command -v "$NODE" >/dev/null 2>&1; then
  cp -r "$TMP/base.lstar.zarr" "$TMP/js.zarr"
  if "$NODE" --experimental-strip-types "$ROOT/js/scripts/extend-viewer.ts" "$TMP/js.zarr" --primary "$PRIMARY" >/tmp/lstar_vp_js.log 2>&1; then
    echo "  -- python == JS (primary=$PRIMARY) --"
    "$PY" "$CHK" equiv "$TMP/py.zarr" "$TMP/js.zarr"    || fail=1
    "$PY" "$CHK" provgroup "$TMP/js.zarr" "$PRIMARY"    || fail=1
  else echo "  FAIL: JS extend-viewer.ts(--primary) errored"; tail -15 /tmp/lstar_vp_js.log; fail=1; fi
else echo "  [skip] no CURRENT WASM dist (build via js/build.sh) / node / zarrita — skipping JS primary= parity"; fi

[ $fail -eq 0 ] && echo "extend_for_viewer(primary=) cross-surface parity PASSED." || { echo "primary= parity FAILED."; exit 1; }
