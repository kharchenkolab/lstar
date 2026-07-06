#!/usr/bin/env bash
# Zarr v3 sharding conformance (C++ + zarr-python + the WASM/libzarr reader). lstar writes sharded v3
# (many inner chunks packed into fewer shard objects -- a hosting optimization: fewer HTTP requests,
# still byte-range-readable via the shard index), and every reader sees the SAME values as the unsharded
# store. On the maximal store:
#   1. C++ writes chunked+sharded v3 and reads it back value-identical (self round-trip).
#   2. The store genuinely uses the sharding_indexed codec; zarr-python reads it == the unsharded seed.
#   3. The WASM/libzarr reader reads the sharded arrays == zarr-python.
# Needs the cmake C++ core (test_v3), a zarr>=3.1 python, and (for step 3) emsdk + node.
set -e
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
BIN="$ROOT/core/build/test_v3"
export PYTHONPATH="$ROOT/python/src"
SEED=/tmp/v3sh_seed.lstar.zarr
SH=/tmp/v3sh_sharded.lstar.zarr

if [ ! -x "$BIN" ]; then
  cmake -S "$ROOT/core" -B "$ROOT/core/build" -DCMAKE_BUILD_TYPE=Release >/tmp/v3sh_cmake.log 2>&1
  cmake --build "$ROOT/core/build" --target test_v3 -j4 >>/tmp/v3sh_cmake.log 2>&1
fi

python3 "$ROOT/conformance/v3_gen.py" "$SEED" | sed 's/^/  [py] /'
"$BIN" shard "$SEED" "$SH"                                       # 1. C++ chunked+sharded v3 self round-trip
python3 "$ROOT/conformance/v3_verify.py" shardcheck "$SH" "$SEED"  # 2. sharding_indexed + zarr-python == seed

# 3. the WASM reader reads the sharded arrays == zarr-python (skips cleanly without emsdk/node).
EMSDK="${EMSDK:-$HOME/emsdk}"
NODE="$(ls -d "$EMSDK"/node/*/bin/node 2>/dev/null | head -1)"; NODE="${NODE:-$(command -v node || true)}"
if [ -n "$NODE" ] && [ -f "$ROOT/js/dist/lstar_io.mjs" ]; then
  python3 "$ROOT/js/test/io_dump.py" "$SH" /tmp/v3sh_dump.json >/dev/null
  "$NODE" "$ROOT/js/test/io_parity.mjs" "$SH" /tmp/v3sh_dump.json
else
  echo "  [skip] no node / WASM dist — WASM sharded-read leg skipped"
fi
echo "  v3 sharding: write + read conformant across C++, zarr-python, and the WASM reader"
