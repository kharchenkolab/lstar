#!/usr/bin/env bash
# JS/WASM conformance: build the WASM kernels, then verify (in Node) the kernels, the zarrita-based
# L* reader, and the viewer query API against Python-written stores + references. Skips cleanly when
# emsdk is absent. Point LSTAR_EMCC_PYTHON at a >=3.10 interpreter if the system python is older.
set -e
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
EMSDK="${EMSDK:-$HOME/emsdk}"
if [ ! -f "$EMSDK/emsdk_env.sh" ]; then
  echo "  [skip] emsdk not found at $EMSDK — skipping JS/WASM conformance"
  exit 0
fi
if [ ! -d "$ROOT/js/node_modules/zarrita" ]; then
  echo "  [skip] zarrita not installed (run: cd js && npm install) — skipping JS conformance"
  exit 0
fi
NODE="$(ls -d "$EMSDK"/node/*/bin/node 2>/dev/null | head -1)"; NODE="${NODE:-node}"

# 1) build the WASM kernels
EMSDK="$EMSDK" bash "$ROOT/js/build.sh" >/tmp/lstar_wasm_build.log 2>&1 \
  || { echo "  FAIL: wasm build"; tail -15 /tmp/lstar_wasm_build.log; exit 1; }
echo "  [wasm] built kernels"

# 2) generate the test store (Python L* writer)
PYTHONPATH="$ROOT/python/src" python3 "$ROOT/js/test/make_store.py" >/tmp/lstar_js_store.log 2>&1 \
  || { echo "  FAIL: test store generation"; tail -15 /tmp/lstar_js_store.log; exit 1; }
echo "  [py ] generated test store + references"

# 3) run the Node tests (kernels: dense+Python conformance; reader: manifest+fields; view: queries)
echo "  -- kernels --"; "$NODE" "$ROOT/js/test/kernels.test.mjs"
echo "  -- reader  --"; "$NODE" --experimental-strip-types "$ROOT/js/test/reader.test.ts" 2>/dev/null
echo "  -- view    --"; "$NODE" --experimental-strip-types "$ROOT/js/test/view.test.ts" 2>/dev/null
echo "JS/WASM conformance PASSED."
