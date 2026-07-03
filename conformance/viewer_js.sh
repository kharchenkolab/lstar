#!/usr/bin/env bash
# JS viewer@0.1 parity: lstar's OWN JS extend (js/scripts/extend-viewer.ts) must agree with the Python
# prep on every viewer field. This closes the gap where only the EXTERNAL pagoda3 prep.ts was ever
# cross-checked (viewer.sh leg c) -- lstar's own extend.ts had NO cross-surface coverage. Two inputs:
#   1. CSC counts  -> exercises the shared reorder (viewer_cell_order via WASM) + group ordering.
#   2. CSR counts  -> exercises A1: JS must NORMALIZE CSR (csrToCsc) like Python/R, not throw (Finding 1).
# Builds the WASM kernels from the current C++ core. Skips cleanly (exit 0) when emcc/node/zarrita are
# unavailable (e.g. locally without python>=3.10 for emcc); runs for real in CI (js-wasm job: emsdk 6 +
# node 22 + python 3.11). Set LSTAR_EMCC_PYTHON to a >=3.10 interpreter if the system python is older.
set -u
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
PY="${LSTAR_PY:-python3}"
export PYTHONPATH="$ROOT/python/src${PYTHONPATH:+:$PYTHONPATH}"
CHK="$ROOT/conformance/viewer_check.py"
EMSDK="${EMSDK:-$HOME/emsdk}"
TMP="$(mktemp -d)"; trap 'rm -rf "$TMP"' EXIT
fail=0

if ! "$PY" -c "import lstar" 2>/dev/null; then
  echo "  [skip] lstar not importable by '$PY' (set LSTAR_PY) — skipping JS viewer parity"; exit 0; fi
if [ ! -d "$ROOT/js/node_modules/zarrita" ]; then
  echo "  [skip] zarrita not installed (cd js && npm install) — skipping JS viewer parity"; exit 0; fi
if [ ! -f "$EMSDK/emsdk_env.sh" ] && ! command -v emcc >/dev/null 2>&1; then
  echo "  [skip] no emcc (emsdk at $EMSDK absent and emcc not on PATH) — skipping JS viewer parity"; exit 0; fi

# build the WASM kernels from the current core (picks up csrToCsc / viewerCellOrder)
if ! EMSDK="$EMSDK" LSTAR_EMCC_PYTHON="${LSTAR_EMCC_PYTHON:-}" bash "$ROOT/js/build.sh" >/tmp/lstar_wasm_vjs.log 2>&1; then
  # a too-old emcc python is an environment gap, not a parity failure -> skip cleanly (see header)
  if grep -q "python 3.10 or above" /tmp/lstar_wasm_vjs.log; then
    echo "  [skip] emcc needs python>=3.10 (set LSTAR_EMCC_PYTHON) — skipping JS viewer parity"; exit 0; fi
  echo "  FAIL: WASM build"; tail -15 /tmp/lstar_wasm_vjs.log; exit 1; fi
NODE="$(ls -d "$EMSDK"/node/*/bin/node 2>/dev/null | head -1)"; NODE="${NODE:-node}"
echo "  [wasm] built kernels; node=$("$NODE" --version 2>/dev/null)"

# reference: Python prep of the CSC base (all surfaces must match this)
"$PY" "$CHK" make-base "$TMP/base.lstar.zarr" csc || fail=1
"$PY" "$CHK" prep-lstar "$TMP/base.lstar.zarr" "$TMP/py.lstar.zarr" || fail=1

jsprep() {  # run lstar's own JS extend in place on a store copy
  "$NODE" --experimental-strip-types "$ROOT/js/scripts/extend-viewer.ts" "$1" >/tmp/lstar_jsprep.log 2>&1 \
    || { echo "  FAIL: extend-viewer.ts on $1"; tail -15 /tmp/lstar_jsprep.log; return 1; }
}

echo "  -- lstar JS extend-viewer.ts == python prep (CSC counts) --"
cp -r "$TMP/base.lstar.zarr" "$TMP/js.lstar.zarr"
if jsprep "$TMP/js.lstar.zarr"; then "$PY" "$CHK" equiv "$TMP/py.lstar.zarr" "$TMP/js.lstar.zarr" || fail=1; else fail=1; fi

echo "  -- lstar JS extend-viewer.ts on CSR counts == python prep (A1 encoding invariance) --"
"$PY" "$CHK" make-base "$TMP/base_csr.lstar.zarr" csr || fail=1
cp -r "$TMP/base_csr.lstar.zarr" "$TMP/js_csr.lstar.zarr"
if jsprep "$TMP/js_csr.lstar.zarr"; then "$PY" "$CHK" equiv "$TMP/py.lstar.zarr" "$TMP/js_csr.lstar.zarr" || fail=1; else fail=1; fi

echo "  -- lstar JS extend-viewer.ts on competing groupings == python prep (#4 detection + group order) --"
"$PY" "$CHK" make-base "$TMP/base_multi.lstar.zarr" csc multi || fail=1
"$PY" "$CHK" prep-lstar "$TMP/base_multi.lstar.zarr" "$TMP/py_multi.lstar.zarr" >/dev/null || fail=1
cp -r "$TMP/base_multi.lstar.zarr" "$TMP/js_multi.lstar.zarr"
if jsprep "$TMP/js_multi.lstar.zarr"; then "$PY" "$CHK" equiv "$TMP/py_multi.lstar.zarr" "$TMP/js_multi.lstar.zarr" || fail=1; else fail=1; fi

if [ "$fail" -eq 0 ]; then echo "  JS viewer parity OK"; else echo "  JS viewer parity FAIL"; fi
exit "$fail"
