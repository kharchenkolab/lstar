#!/usr/bin/env bash
# Build the lstar WASM kernels (Emscripten/embind) from the header-only C++ core.
#
# Needs `emcc` available — either an activated emsdk (EMSDK set, or ~/emsdk) or emcc on PATH
# (e.g. `brew install emscripten`) — and Python >= 3.10 for emcc. If the system python3 is older,
# point LSTAR_EMCC_PYTHON at a >=3.10 interpreter, e.g.:
#   LSTAR_EMCC_PYTHON=/path/to/python3.10 bash js/build.sh
set -e
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
JS="$ROOT/js"
: "${EMSDK:=$HOME/emsdk}"

# Capture our chosen interpreter BEFORE sourcing emsdk_env.sh (which overwrites EMSDK_PYTHON with
# whatever python emsdk found at install time). Source emsdk only if present; otherwise rely on a
# PATH emcc (Homebrew/system install).
PYBIN="${LSTAR_EMCC_PYTHON:-}"
# Building an emscripten PORT (e.g. -sUSE_ZLIB) spawns the `emcc` wrapper as a subprocess, whose
# `#!/usr/bin/env python3` shebang resolves python3 from PATH -- not from our chosen interpreter. So when
# LSTAR_EMCC_PYTHON is set (system python too old), put its dir first on PATH so the port build uses it too.
if [ -n "$PYBIN" ]; then export PATH="$(dirname "$PYBIN"):$PATH"; fi
if [ -f "$EMSDK/emsdk_env.sh" ]; then
  # shellcheck disable=SC1091
  source "$EMSDK/emsdk_env.sh" >/dev/null 2>&1 || true
fi

if [ -n "$PYBIN" ] && [ -f "$EMSDK/upstream/emscripten/emcc.py" ]; then
  EMCC=("$PYBIN" "$EMSDK/upstream/emscripten/emcc.py")
else
  EMCC=(emcc)
fi

mkdir -p "$JS/dist"
"${EMCC[@]}" "$JS/wasm/lstar_wasm.cpp" \
  -I"$ROOT/core/include" \
  -std=c++17 -O3 -lembind \
  -sUSE_ZLIB=1 -DLSTAR_HAVE_ZLIB \
  -sMODULARIZE=1 -sEXPORT_ES6=1 -sENVIRONMENT=node,web \
  -sALLOW_MEMORY_GROWTH=1 \
  -sEXPORT_NAME=createLstarKernels \
  -o "$JS/dist/lstar_kernels.mjs"
echo "built $JS/dist/lstar_kernels.mjs (+ lstar_kernels.wasm)"

# The I/O module: the libzarr-backed reader (retires the zarrita reimplementation). Needs libzarr's
# gzip codec (-DLIBZARR_HAS_ZLIB) atop the zlib port, AND its zstd codec (-DLIBZARR_HAS_ZSTD) — zarr-python
# 3's DEFAULT v3 compressor, so a hosted v3 store may be zstd-encoded. Emscripten has no zstd port, so we
# compile a vendored zstd DECODER: with -DLIBZARR_ZSTD_DECODE_ONLY libzarr omits the zstd compress side
# (the reader never encodes), so the small DECOMPRESS-ONLY amalgamation (js/third_party/zstd/zstddeclib.c)
# links — no ZSTD_compress* references. Both defines are required (decode-only is additive on HAS_ZSTD).
# A .c source can't share the .cpp's -std=c++17, so compile it to an object first, then link it in.
# Kept a separate module from the kernels for now (the kernels stay I/O-free); the viewer loads both.
"${EMCC[@]}" -c "$JS/third_party/zstd/zstddeclib.c" -O3 -o "$JS/dist/zstddeclib.o"
"${EMCC[@]}" "$JS/wasm/lstar_reader.cpp" "$JS/dist/zstddeclib.o" \
  -I"$ROOT/core/include" -I"$JS/third_party/zstd" \
  -std=c++17 -O3 -lembind \
  -sUSE_ZLIB=1 -DLSTAR_HAVE_ZLIB -DLIBZARR_HAS_ZLIB -DLIBZARR_HAS_ZSTD -DLIBZARR_ZSTD_DECODE_ONLY \
  -sMODULARIZE=1 -sEXPORT_ES6=1 -sENVIRONMENT=node,web \
  -sALLOW_MEMORY_GROWTH=1 \
  -sEXPORT_NAME=createLstarIO \
  -o "$JS/dist/lstar_io.mjs"
echo "built $JS/dist/lstar_io.mjs (+ lstar_io.wasm)"

# The WRITER module: the write-side pure functions from libzarr — chunk codec ENCODE (gzip/zstd) and
# shard-object assembly (shard::pack). Loaded only when WRITING (the pagoda3 prep / a store producer), so
# it's separate from the lean decode-only reader: it links the FULL zstd amalgamation (js/third_party/zstd/
# zstd.c, LIBZARR_HAS_ZSTD *without* DECODE_ONLY) so encode covers zstd. The drift-prone bytes (codec
# encode, shard index + crc32c) stay in libzarr; JS owns chunking + store writes + metadata.
"${EMCC[@]}" -c "$JS/third_party/zstd/zstd.c" -O3 -o "$JS/dist/zstd_full.o"
"${EMCC[@]}" "$JS/wasm/lstar_writer.cpp" "$JS/dist/zstd_full.o" \
  -I"$ROOT/core/include" -I"$JS/third_party/zstd" \
  -std=c++17 -O3 -lembind \
  -sUSE_ZLIB=1 -DLSTAR_HAVE_ZLIB -DLIBZARR_HAS_ZLIB -DLIBZARR_HAS_ZSTD \
  -sMODULARIZE=1 -sEXPORT_ES6=1 -sENVIRONMENT=node,web \
  -sALLOW_MEMORY_GROWTH=1 \
  -sEXPORT_NAME=createLstarWriter \
  -o "$JS/dist/lstar_writer.mjs"
echo "built $JS/dist/lstar_writer.mjs (+ lstar_writer.wasm)"

# --- resizable-heap UTF8 decode guard -----------------------------------------------------------------
# Emscripten's UTF8ArrayToString decodes the heap IN PLACE: `UTF8Decoder.decode(heapOrArray.subarray(idx,endPtr))`.
# With -sALLOW_MEMORY_GROWTH the heap is a *growable* buffer, and browsers that back it with a RESIZABLE
# ArrayBuffer throw `TextDecoder ... must not be resizable` on that decode — which fires for every embind
# std::string return (e.g. Reader.groupAttrs at store OPEN), crashing the viewer. (emsdk 5.x removed the
# `-sTEXTDECODER=0` manual-loop opt-out — "must be 1 or 2" — so the flag route is gone.) Copy off the
# growable/shared heap before decoding. Fail LOUD if the expected glue shape is absent, so a toolchain bump
# can't silently drop the guard and reintroduce the crash.
for m in "$JS/dist/lstar_kernels.mjs" "$JS/dist/lstar_io.mjs" "$JS/dist/lstar_writer.mjs"; do
  grep -qF 'UTF8Decoder.decode(heapOrArray.subarray(idx,endPtr))' "$m" \
    || { echo "ERROR: $m — Emscripten UTF8ArrayToString glue not found (toolchain changed?); resizable-heap guard NOT applied" >&2; exit 1; }
  perl -0pi -e 's/UTF8Decoder\.decode\(heapOrArray\.subarray\(idx,endPtr\)\)/UTF8Decoder.decode((heapOrArray.buffer.resizable||typeof SharedArrayBuffer!="undefined"&&heapOrArray.buffer instanceof SharedArrayBuffer)?heapOrArray.slice(idx,endPtr):heapOrArray.subarray(idx,endPtr))/g' "$m"
  echo "patched $m (resizable/shared-heap UTF8 decode guard)"
done
