#!/usr/bin/env bash
# Zarr v3 format conformance (C++ + zarr-python). lstar reads BOTH v2 and v3 (libzarr's open probes
# zarr.json before .zarray) and writes v3 on request (default stays v2). On a maximal store (every
# encoding + nullable/graph/partial/aux/viewer):
#   1. C++ re-emits v2 -> v3 and reads it back value-identical (self round-trip).
#   2. zarr-python 3 reads the C++-written v3: genuine v3, inline-consolidated, values == v2.
#   3. C++ reads a v3 store written by zarr-python (INDEPENDENT writer, zstd-by-default) == v2.
# Needs the cmake C++ core (test_v3), a zarr>=3.1 python, and libzstd (for step 3). The R surface is
# covered separately by v3_r.sh (R + zarr-python, no cmake). Self-contained.
set -e
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
BIN="$ROOT/core/build/test_v3"
export PYTHONPATH="$ROOT/python/src"

SEED=/tmp/v3_seed.lstar.zarr           # v2, maximal
V3=/tmp/v3_cpp.lstar.zarr              # C++-written v3
V3ZPY=/tmp/v3_zpy.lstar.zarr           # zarr-python-written v3 (zstd)

if [ ! -x "$BIN" ]; then               # standalone runs (run.sh builds the whole core first)
  cmake -S "$ROOT/core" -B "$ROOT/core/build" -DCMAKE_BUILD_TYPE=Release >/tmp/v3_cmake.log 2>&1
  cmake --build "$ROOT/core/build" --target test_v3 -j4 >>/tmp/v3_cmake.log 2>&1
fi

python3 "$ROOT/conformance/v3_gen.py" "$SEED" | sed 's/^/  [py] /'
"$BIN" write "$SEED" "$V3"                                     # 1. C++ v2 -> v3 -> read-back
python3 "$ROOT/conformance/v3_verify.py" check "$V3" "$SEED"   # 2. zarr-python reads C++ v3
python3 "$ROOT/conformance/v3_verify.py" reemit "$SEED" "$V3ZPY"  # 3a. zarr-python writes v3 (zstd)
"$BIN" compare "$SEED" "$V3ZPY"                                # 3b. C++ reads zarr-python v3 == v2
echo "  v3 format: read-both + write-v3 conformant across C++ and zarr-python (incl. zstd default)"
