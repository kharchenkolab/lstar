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
  -sMODULARIZE=1 -sEXPORT_ES6=1 -sENVIRONMENT=node,web \
  -sALLOW_MEMORY_GROWTH=1 \
  -sEXPORT_NAME=createLstarKernels \
  -o "$JS/dist/lstar_kernels.mjs"
echo "built $JS/dist/lstar_kernels.mjs (+ lstar_kernels.wasm)"
