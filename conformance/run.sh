#!/usr/bin/env bash
# Master conformance runner: build the C++ core, install the R package, run the Python
# tests, the Python<->C++ cross-impl test, and the cross-format (Seurat/SCE) chain.
set -e
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"
RLIB="$ROOT/.Rlib"
pass() { echo "  PASS  $1"; }

echo "== build C++ core (libstar) =="
cmake -S core -B core/build -DCMAKE_BUILD_TYPE=Release >/tmp/lstar_cmake.log 2>&1
cmake --build core/build -j4 >>/tmp/lstar_cmake.log 2>&1
pass "libstar build"

echo "== install R package =="
mkdir -p "$RLIB"
# --preclean: R's make has no header-dependency tracking, so a changed vendored header won't
# trigger recompilation of an existing .o; preclean forces a fresh build.
R CMD INSTALL --preclean --no-multiarch --library="$RLIB" R >/tmp/lstar_rinstall.log 2>&1 \
  && pass "R package install" || { echo "  FAIL R install"; tail -15 /tmp/lstar_rinstall.log; exit 1; }

echo "== Python tests =="
for t in test_roundtrip test_anndata_profile test_crossimpl test_validate test_versions test_lazy test_stream_write; do
  if PYTHONPATH=python/src python3 python/tests/$t.py >/tmp/lstar_$t.log 2>&1; then pass "$t"
  else echo "  FAIL  $t"; tail -15 /tmp/lstar_$t.log; exit 1; fi
done

echo "== cross-format conformance (R: Seurat + SCE) =="
bash conformance/cross_format.sh >/tmp/lstar_cf.log 2>&1 \
  && pass "AnnData<->Seurat<->SCE via L*" || { echo "  FAIL cross-format"; tail -15 /tmp/lstar_cf.log; exit 1; }

echo "== collection conformance (R collection -> L* -> Python) =="
bash conformance/collection.sh >/tmp/lstar_coll.log 2>&1 \
  && pass "collection of samples round-trips R->py" || { echo "  FAIL collection"; tail -15 /tmp/lstar_coll.log; exit 1; }

echo "== chunked+gzip conformance (Python -> C++ -> Python) =="
bash conformance/chunked.sh >/tmp/lstar_chunk.log 2>&1 \
  && pass "chunked+compressed cross-impl + transpose" || { echo "  FAIL chunked"; tail -15 /tmp/lstar_chunk.log; exit 1; }

echo "== blocked-reader conformance (R/C++ bounded col stats == full read) =="
bash conformance/stream_reduce.sh >/tmp/lstar_sr.log 2>&1 \
  && pass "blocked col-stats reducer matches full read" || { echo "  FAIL stream_reduce"; tail -15 /tmp/lstar_sr.log; exit 1; }

echo "== JS/WASM (Emscripten kernels + zarrita reader + viewer API; skips if emsdk absent) =="
LSTAR_EMCC_PYTHON="${LSTAR_EMCC_PYTHON:-}" bash conformance/js.sh 2>&1 | sed 's/^/  /'

echo "ALL CONFORMANCE TESTS PASSED"
