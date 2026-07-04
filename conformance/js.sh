#!/usr/bin/env bash
# JS/WASM conformance: build the WASM kernels, then verify (in Node) the kernels, the zarrita-based
# L* reader, and the viewer query API against Python-written stores + references. Skips cleanly when
# emsdk is absent. Point LSTAR_EMCC_PYTHON at a >=3.10 interpreter if the system python is older.
# Origin coverage: Py-authored ✓ (JS reads) | JS-authored ✓ (writer_make.ts emits every encoding, Python
# cross-reads via writer_crossread.py) — see conformance/README.md
set -e
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
EMSDK="${EMSDK:-$HOME/emsdk}"
if [ ! -f "$EMSDK/emsdk_env.sh" ] && ! command -v emcc >/dev/null 2>&1; then
  echo "  [skip] no emcc (emsdk at $EMSDK absent and emcc not on PATH) — skipping JS/WASM conformance"
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
echo "  -- extend  --"; "$NODE" --experimental-strip-types "$ROOT/js/test/extend_primary.test.ts" 2>/dev/null   # extend_for_viewer(primary=): hoist + compose + validation

# 4) the WRITE side: JS round-trip (writer.test.ts), then the cross-language gate -- JS writes a
# chunked + gzip-compressed store with every encoding (CSC/dense/categorical/mask/partial/aux, using the
# WASM zlib kernel) and the Python reader validates it clean + value-equal.
echo "  -- writer  --"; "$NODE" --experimental-strip-types "$ROOT/js/test/writer.test.ts" 2>/dev/null
echo "  -- writer x-read (JS-write chunked+gzip -> Python read+validate) --"
WX=/tmp/lstar_js_writer_cross.lstar.zarr; rm -rf "$WX"
"$NODE" --experimental-strip-types "$ROOT/js/test/writer_make.ts" "$WX" 2>/dev/null \
  || { echo "  FAIL: JS writer_make"; exit 1; }
PYTHONPATH="$ROOT/python/src" python3 "$ROOT/js/test/writer_crossread.py" "$WX" \
  || { echo "  FAIL: writer_crossread"; exit 1; }

# 5) reverse leg: JS addToStore appends a derived field to a PYTHON-written store -> Python re-reads it
echo "  -- writer extend (Python-write -> JS addToStore -> Python read) --"
WE=/tmp/lstar_js_writer_extend.lstar.zarr; rm -rf "$WE"; cp -r "$ROOT/js/test/data/sample.lstar.zarr" "$WE"
"$NODE" --experimental-strip-types "$ROOT/js/test/writer_extend.ts" "$WE" 2>/dev/null \
  || { echo "  FAIL: JS writer_extend"; exit 1; }
PYTHONPATH="$ROOT/python/src" python3 - "$WE" <<'PY' || { echo "  FAIL: extend re-read"; exit 1; }
import sys, lstar
ds = lstar.read(sys.argv[1])
assert "od_score" in ds.fields and "counts" in ds.fields, list(ds.fields)   # derived added, original kept
assert ds.axis("od_groups").origin == "derived" and "derived@0.1" in ds.profiles
assert not [e for e in lstar.validate(ds) if e.startswith("ERROR")]
print("  [py] JS addToStore on a Python store re-reads clean (derived field + axis + profile)")
PY
echo "JS/WASM conformance PASSED."
