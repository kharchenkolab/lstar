#!/usr/bin/env bash
# Real-data v3 conformance (LOCAL; needs the gitignored testdata/ corpus). For each cached real dataset:
# ingest it, write v2 + v3, assert the Python reference reads both to identical values (v3_corpus.py),
# and cross-read the two stores with the C++ libzarr core (test_v3 compare). Tolerant: missing datasets
# and unbuilt cores degrade to SKIP. This is the real-data analog of v3_format.sh (which uses a synthetic
# maximal store) -- it proves the v2->v3 migration on messy real AnnData/MuData, not just the fixture.
set -uo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"
export PYTHONPATH="$ROOT/python/src"
BIN="$ROOT/core/build/test_v3"
TD="$ROOT/testdata"

RLIB="$ROOT/.Rlib"
if [ ! -x "$BIN" ]; then
  cmake -S core -B core/build -DCMAKE_BUILD_TYPE=Release >/tmp/v3c_cmake.log 2>&1 || true
  cmake --build core/build --target test_v3 -j4 >>/tmp/v3c_cmake.log 2>&1 || true
fi

# optional JS/WASM reader (build once) + R package -- so the corpus check spans all four surfaces on real
# data when the toolchains are present; each degrades to a skip otherwise.
EMSDK="${EMSDK:-$HOME/emsdk}"
NODE="$(ls -d "$EMSDK"/node/*/bin/node 2>/dev/null | head -1)"; NODE="${NODE:-$(command -v node || true)}"
HAVE_JS=0
if [ -n "$NODE" ] && [ -d "$ROOT/js/node_modules/zarrita" ]; then
  if EMSDK="$EMSDK" bash "$ROOT/js/build.sh" >/tmp/v3c_wasm.log 2>&1; then HAVE_JS=1; else echo "  [skip] wasm build failed — JS surface skipped"; fi
fi
HAVE_R=0
if command -v Rscript >/dev/null 2>&1 && Rscript -e '.libPaths(c("'"$RLIB"'", .libPaths())); quit(status=!requireNamespace("lstar", quietly=TRUE))' >/dev/null 2>&1; then HAVE_R=1; fi

# cached real datasets to try (each optional). h5ad = AnnData, h5mu = MuData.
DATASETS=(
  "$TD/pbmc3k_processed.h5ad"
  "$TD/pancreas.h5ad"
  "$TD/pancreas_velocity_small.h5ad"
  "$TD/minipbcite.h5mu"
)

nok=0; nskip=0
for src in "${DATASETS[@]}"; do
  if [ ! -e "$src" ]; then echo "  SKIP $(basename "$src"): not cached"; nskip=$((nskip+1)); continue; fi
  v2="/tmp/v3corp_v2.lstar.zarr"; v3="/tmp/v3corp_v3.lstar.zarr"; rm -rf "$v2" "$v3"
  out="$(python3 "$ROOT/conformance/v3_corpus.py" "$src" "$v2" "$v3" 2>/tmp/v3corp_py.log)" || { echo "  FAIL $(basename "$src") (python)"; tail -8 /tmp/v3corp_py.log; exit 1; }
  echo "$out"
  case "$out" in *SKIP*) nskip=$((nskip+1)); continue;; esac
  if [ -x "$BIN" ]; then
    "$BIN" compare "$v2" "$v3" >/tmp/v3corp_cpp.log 2>&1 && echo "    [C++] libzarr reads both, v2==v3" \
      || { echo "  FAIL $(basename "$src") (C++ compare)"; tail -8 /tmp/v3corp_cpp.log; exit 1; }
  fi
  if [ "$HAVE_JS" = 1 ]; then
    # absolute: the WASM reader's arrays == zarr-python on real data, in BOTH formats; then format-invariance
    # across the full L* API (v2==v3). Zarr-python is the oracle here, not zarrita (which mis-reads bool).
    python3 "$ROOT/js/test/io_dump.py" "$v2" /tmp/v3corp_d2.json >/dev/null 2>&1
    python3 "$ROOT/js/test/io_dump.py" "$v3" /tmp/v3corp_d3.json >/dev/null 2>&1
    { "$NODE" "$ROOT/js/test/io_parity.mjs" "$v2" /tmp/v3corp_d2.json \
      && "$NODE" "$ROOT/js/test/io_parity.mjs" "$v3" /tmp/v3corp_d3.json \
      && "$NODE" --experimental-strip-types "$ROOT/js/test/wasm_corpus.mjs" "$v2" "$v3"; } >/tmp/v3corp_js.log 2>&1 \
      && echo "    [JS ] libzarr reads real v2+v3 == zarr-python; v2==v3 across the L* API" \
      || { echo "  FAIL $(basename "$src") (JS wasm reader)"; tail -10 /tmp/v3corp_js.log; exit 1; }
  fi
  if [ "$HAVE_R" = 1 ]; then
    Rscript -e '.libPaths(c("'"$RLIB"'", .libPaths())); suppressMessages(library(lstar))
      a <- lstar_read("'"$v2"'"); b <- lstar_read("'"$v3"'")
      stopifnot(length(a$fields) == length(b$fields), length(a$axes) == length(b$axes), length(b$fields) > 0)
      cat(sprintf("    [R  ] reads v2 + v3 (%d fields, %d axes)\n", length(b$fields), length(b$axes)))' \
      2>/tmp/v3corp_r.log || { echo "  FAIL $(basename "$src") (R read)"; tail -8 /tmp/v3corp_r.log; exit 1; }
  fi
  nok=$((nok+1))
done
echo "  v3 corpus: $nok datasets round-tripped v2<->v3 across C++/Python$([ "$HAVE_JS" = 1 ] && echo "/JS")$([ "$HAVE_R" = 1 ] && echo "/R"), $nskip skipped"
