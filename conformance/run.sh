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
for t in test_roundtrip test_anndata_profile test_crossimpl test_validate test_versions test_lazy test_stream_write test_categorical test_induce test_nullable test_aux test_de; do
  if PYTHONPATH=python/src python3 python/tests/$t.py >/tmp/lstar_$t.log 2>&1; then pass "$t"
  else echo "  FAIL  $t"; tail -15 /tmp/lstar_$t.log; exit 1; fi
done

echo "== categorical-encoding conformance (codes/categories/ordered/-1 across Py/C++/R) =="
bash conformance/categorical.sh >/tmp/lstar_cat.log 2>&1 \
  && pass "categorical round-trips Py<->C++<->R" || { echo "  FAIL categorical"; tail -15 /tmp/lstar_cat.log; exit 1; }

echo "== induction conformance (factor axis + induced_by round-trip + checkable consistency across Py/C++/R) =="
bash conformance/induce.sh >/tmp/lstar_ind.log 2>&1 \
  && pass "induced factor axes + induced_by round-trip Py<->C++<->R" || { echo "  FAIL induce"; tail -15 /tmp/lstar_ind.log; exit 1; }

echo "== nullable conformance (validity mask: nullable Int/bool/string round-trip across Py/C++/R) =="
bash conformance/nullable.sh >/tmp/lstar_null.log 2>&1 \
  && pass "nullable validity masks round-trip Py<->C++<->R" || { echo "  FAIL nullable"; tail -15 /tmp/lstar_null.log; exit 1; }

echo "== aux passthrough conformance (uns/@misc subtree round-trips verbatim across Py/C++/R) =="
bash conformance/aux.sh >/tmp/lstar_aux.log 2>&1 \
  && pass "lossless passthrough subtree round-trips Py<->C++<->R" || { echo "  FAIL aux"; tail -15 /tmp/lstar_aux.log; exit 1; }

echo "== DE-bundle conformance (rank_genes_groups -> (factor,genes) bundle; lstar_markers in Py + R) =="
bash conformance/de.sh >/tmp/lstar_de.log 2>&1 \
  && pass "DE bundle round-trips + tidy markers (Py + R)" || { echo "  FAIL de"; tail -15 /tmp/lstar_de.log; exit 1; }

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

echo "== fused-view reducer conformance (depth+log1p mean/var + grouped sum == in-memory) =="
bash conformance/fused_view.sh >/tmp/lstar_fv.log 2>&1 \
  && pass "fused depth-view reducers match in-memory" || { echo "  FAIL fused_view"; tail -15 /tmp/lstar_fv.log; exit 1; }

echo "== block-reader conformance (lstar_read_block / read_genes == full read) =="
bash conformance/read_block.sh >/tmp/lstar_rb.log 2>&1 \
  && pass "block reader (contiguous + scattered) matches full read" || { echo "  FAIL read_block"; tail -15 /tmp/lstar_rb.log; exit 1; }

echo "== R writer chunking/compression (R-written chunked+gzip == default; cross-impl) =="
bash conformance/r_write_chunked.sh >/tmp/lstar_rwc.log 2>&1 \
  && pass "R writer emits chunked+gzip stores readable by Py/C++" || { echo "  FAIL r_write_chunked"; tail -15 /tmp/lstar_rwc.log; exit 1; }

echo "== disk-backed targets (L* -> h5ad -> backed AnnData / Seurat+BPCells / SCE+HDF5Array) =="
bash conformance/backed_targets.sh >/tmp/lstar_bt.log 2>&1 \
  && { pass "disk-backed conversion targets"; grep "SKIP" /tmp/lstar_bt.log | sed 's/^/      /'; } \
  || { echo "  FAIL backed_targets"; tail -15 /tmp/lstar_bt.log; exit 1; }

echo "== JS/WASM (Emscripten kernels + zarrita reader + viewer API; skips if emsdk absent) =="
LSTAR_EMCC_PYTHON="${LSTAR_EMCC_PYTHON:-}" bash conformance/js.sh 2>&1 | sed 's/^/  /'

echo "ALL CONFORMANCE TESTS PASSED"
