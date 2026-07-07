#!/usr/bin/env bash
# Zstd read on the WASM/libzarr reader. zarr-python 3's DEFAULT v3 compressor is zstd, so a hosted v3 store
# may be zstd-encoded. Emscripten has no zstd port, so the reader is built with libzarr's
# LIBZARR_ZSTD_DECODE_ONLY guard + a vendored DECOMPRESS-ONLY zstd amalgamation
# (js/third_party/zstd/zstddeclib.c) — it decodes zstd without linking the compressor. Python writes a zstd
# v3 store; the WASM reader reads it == zarr-python. Needs a zarr>=3.1 python + emsdk/node + built js/dist.
set -e
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
export PYTHONPATH="$ROOT/python/src"
SEED=/tmp/zstdjs_seed.lstar.zarr
ZST=/tmp/zstdjs_zstd.lstar.zarr

python3 "$ROOT/conformance/v3_gen.py" "$SEED" | sed 's/^/  [py] /'
python3 - "$SEED" "$ZST" <<'PY'
import sys, json, numcodecs, lstar
lstar.write(lstar.read(sys.argv[1]), sys.argv[2], format="v3", compressor=numcodecs.Zstd(level=5))
codecs = [c.get("name") for c in json.load(open(sys.argv[2] + "/fields/counts/data/zarr.json"))["codecs"]]
assert "zstd" in codecs, f"expected zstd codec, got {codecs}"
print(f"  [py] wrote zstd v3 store (counts/data codecs: {codecs})")
PY

EMSDK="${EMSDK:-$HOME/emsdk}"
NODE="$(ls -d "$EMSDK"/node/*/bin/node 2>/dev/null | head -1)"; NODE="${NODE:-$(command -v node || true)}"
if [ -n "$NODE" ] && [ -f "$ROOT/js/dist/lstar_io.mjs" ]; then
  python3 "$ROOT/js/test/io_dump.py" "$ZST" /tmp/zstdjs_dump.json >/dev/null
  "$NODE" "$ROOT/js/test/io_parity.mjs" "$ZST" /tmp/zstdjs_dump.json
  echo "  zstd read on the WASM reader: decodes a zstd v3 store == zarr-python (decode-only build)"
else
  echo "  [skip] no node / WASM dist — WASM zstd-read leg skipped"
fi
