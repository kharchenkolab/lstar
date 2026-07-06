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
# The R package vendors the header-only core at R/inst/include (it must be self-contained for CRAN --
# core/ is a sibling dir, not shipped in the tarball). Sync it from the canonical core/ so a core change
# can't leave R compiling a stale copy. (Committed too, so CI/CRAN build the current header.)
cp core/include/lstar/lstar.hpp R/inst/include/lstar/lstar.hpp
# --preclean: R's make has no header-dependency tracking, so a changed vendored header won't
# trigger recompilation of an existing .o; preclean forces a fresh build.
R CMD INSTALL --preclean --no-multiarch --library="$RLIB" R >/tmp/lstar_rinstall.log 2>&1 \
  && pass "R package install" || { echo "  FAIL R install"; tail -15 /tmp/lstar_rinstall.log; exit 1; }

echo "== Python tests =="
for t in test_roundtrip test_anndata_profile test_crossimpl test_validate test_versions test_lazy test_stream_write test_categorical test_induce test_nullable test_aux test_de test_tier1_promote test_real_atlas test_mudata test_v3_format; do
  if PYTHONPATH=python/src python3 python/tests/$t.py >/tmp/lstar_$t.log 2>&1; then pass "$t"
  else echo "  FAIL  $t"; tail -15 /tmp/lstar_$t.log; exit 1; fi
done

# Two-tier contract gate (LOCAL only -- needs real data; not in CI): every synthetic CI fixture must
# structurally mirror the real dataset it stands in for. Catches synth drift when upstream libs change.
echo "== synthetic-faithfulness guard (synthetic structure ⊆ real corpus) =="
PYTHONPATH=python/src python3 python/tests/test_synth_faithful.py >/tmp/lstar_faithful.log 2>&1 \
  && pass "synthetic fixtures mirror the real corpus" || { echo "  FAIL faithfulness"; cat /tmp/lstar_faithful.log; exit 1; }

echo "== binding-parity tripwire (a shared kernel bound in one surface must be bound in all) =="
python3 conformance/binding_parity.py >/tmp/lstar_bindparity.log 2>&1 \
  && pass "compute kernels bound symmetrically across Python/R/WASM" \
  || { echo "  FAIL binding parity"; cat /tmp/lstar_bindparity.log; exit 1; }

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
bash conformance/passthrough.sh >/tmp/lstar_aux.log 2>&1 \
  && pass "lossless passthrough subtree round-trips Py<->C++<->R" || { echo "  FAIL aux"; tail -15 /tmp/lstar_aux.log; exit 1; }

echo "== R-path fidelity (float32 dtype + graph directed/weighted survive Py->R->Py) =="
bash conformance/r_fidelity.sh >/tmp/lstar_rfid.log 2>&1 \
  && pass "value dtype + directed/weighted preserved through the R path" || { echo "  FAIL r_fidelity"; tail -12 /tmp/lstar_rfid.log; exit 1; }

echo "== DE-bundle conformance (rank_genes_groups -> (factor,genes) bundle; lstar_markers in Py + R) =="
bash conformance/de.sh >/tmp/lstar_de.log 2>&1 \
  && pass "DE bundle round-trips + tidy markers (Py + R)" || { echo "  FAIL de"; tail -15 /tmp/lstar_de.log; exit 1; }

echo "== viewer@0.1 grouping-policy linter (Python/R/JS preferred-grouping lists single-sourced) =="
python3 conformance/policy_linter.py >/tmp/lstar_policy.log 2>&1 \
  && pass "grouping-detection policy single-sourced across Py/R/JS" \
  || { echo "  FAIL policy linter"; cat /tmp/lstar_policy.log; exit 1; }

echo "== viewer@0.1 profile conformance (python prep + R pagoda2 + native-R extend, vs the spec) =="
bash conformance/viewer.sh >/tmp/lstar_viewer.log 2>&1 \
  && pass "viewer@0.1 prep conformant + R==python (set PAGODA3 to also check pagoda3 prep.ts)" \
  || { echo "  FAIL viewer"; tail -20 /tmp/lstar_viewer.log; exit 1; }

echo "== viewer@0.1 corpus-driven parity (convert corpus -> prep -> Py==R[==JS] on all fields) =="
bash conformance/viewer_corpus.sh >/tmp/lstar_viewer_corpus.log 2>&1 \
  && pass "viewer@0.1 corpus-driven cross-surface parity" \
  || { echo "  FAIL viewer corpus"; tail -25 /tmp/lstar_viewer_corpus.log; exit 1; }

echo "== extend_for_viewer(primary=) cross-surface parity (Py==R[==JS] reorder keyed on a given primary) =="
bash conformance/viewer_primary.sh >/tmp/lstar_viewer_primary.log 2>&1 \
  && pass "extend_for_viewer(primary=) cross-surface parity" \
  || { echo "  FAIL viewer primary"; tail -20 /tmp/lstar_viewer_primary.log; exit 1; }

echo "== Seurat -> extend_for_viewer conformance (boolean QC excluded; active-idents not duplicated; real + synthetic) =="
bash conformance/viewer_seurat.sh >/tmp/lstar_viewer_seurat.log 2>&1 \
  && pass "Seurat -> viewer prep (real SeuratData + synthetic)" \
  || { echo "  FAIL viewer seurat"; tail -20 /tmp/lstar_viewer_seurat.log; exit 1; }

echo "== multimodal cross-format consistency (CITE-seq via Seurat == via MuData -> same feature axes) =="
bash conformance/multimodal_xformat.sh >/tmp/lstar_xfmt.log 2>&1 \
  && pass "Seurat/MuData multimodal land on the same canonical axes" \
  || { echo "  FAIL multimodal x-format"; tail -15 /tmp/lstar_xfmt.log; exit 1; }

echo "== cross-format conformance (R: Seurat + SCE) =="
bash conformance/cross_format.sh >/tmp/lstar_cf.log 2>&1 \
  && pass "AnnData<->Seurat<->SCE via L*" || { echo "  FAIL cross-format"; tail -15 /tmp/lstar_cf.log; exit 1; }

echo "== convert CLI (lstar convert: detect/route + fidelity report + native-acceptance of the target) =="
bash conformance/convert_cli.sh >/tmp/lstar_cli.log 2>&1 \
  && pass "lstar convert CLI (h5ad<->store<->Seurat, native-valid)" || { echo "  FAIL convert_cli"; tail -20 /tmp/lstar_cli.log; exit 1; }

echo "== Seurat Tier-1 extras (DimReduc stdev -> measure; active Idents captured + restored) =="
bash conformance/seurat_extras.sh >/tmp/lstar_se.log 2>&1 \
  && pass "Seurat stdev + active Idents round-trip" || { echo "  FAIL seurat_extras"; tail -15 /tmp/lstar_se.log; exit 1; }

echo "== Seurat version variety (v3 Assay / v5 Assay5 / v5 split / SCTAssay / multimodal round-trip) =="
bash conformance/seurat_versions.sh >/tmp/lstar_sv.log 2>&1 \
  && pass "Seurat version variants round-trip + validate" || { echo "  FAIL seurat_versions"; tail -20 /tmp/lstar_sv.log; exit 1; }

echo "== Seurat v2 (pre-Assay legacy 'seurat' class -> read + old-to-new conversion) =="
bash conformance/seurat_v2.sh >/tmp/lstar_sv2.log 2>&1 \
  && pass "Seurat v2 pre-Assay object read + validate" || { echo "  FAIL seurat_v2"; tail -20 /tmp/lstar_sv2.log; exit 1; }

echo "== SCE version variety (counts / +reducedDims / +altExps / +colData-rowData factors+metadata) =="
bash conformance/sce_versions.sh >/tmp/lstar_scev.log 2>&1 \
  && pass "SCE version variants round-trip + validate" || { echo "  FAIL sce_versions"; tail -20 /tmp/lstar_scev.log; exit 1; }

echo "== LOCAL real Seurat/SCE corpus (real published objects; skips if SeuratData/scRNAseq absent) =="
bash conformance/real_corpus_r.sh 2>&1 | sed "s/^/  /"

echo "== collection conformance (R collection -> L* -> Python) =="
bash conformance/collection.sh >/tmp/lstar_coll.log 2>&1 \
  && pass "collection of samples round-trips R->py" || { echo "  FAIL collection"; tail -15 /tmp/lstar_coll.log; exit 1; }

echo "== true-collection conformance (collection_from: divergent + disjoint/cross-species; Py<->R; pseudobulk) =="
bash conformance/collection_true.sh >/tmp/lstar_coltrue.log 2>&1 \
  && pass "heterogeneous collections (divergent/disjoint genes) round-trip Py<->R, never flattened" \
  || { echo "  FAIL collection_true"; tail -20 /tmp/lstar_coltrue.log; exit 1; }

echo "== conos conversions (Conos<->L*; collection -> Seurat v5 split + AnnData; graph-only, no corrected matrix) =="
bash conformance/conos.sh >/tmp/lstar_conos.log 2>&1 \
  && { pass "Conos collection round-trips + converts to Seurat v5 / AnnData"; grep "\[skip\]" /tmp/lstar_conos.log | sed 's/^/      /'; } \
  || { echo "  FAIL conos"; tail -25 /tmp/lstar_conos.log; exit 1; }

echo "== chunked+gzip conformance (Python -> C++ -> Python) =="
bash conformance/chunked.sh >/tmp/lstar_chunk.log 2>&1 \
  && pass "chunked+compressed cross-impl + transpose" || { echo "  FAIL chunked"; tail -15 /tmp/lstar_chunk.log; exit 1; }

echo "== Zarr v3 format (read both v2+v3; write v3; C++<->zarr-python<->R, incl. zstd-default) =="
bash conformance/v3_format.sh >/tmp/lstar_v3.log 2>&1 \
  && pass "v3 read+write conformant across C++/zarr-python/R" || { echo "  FAIL v3_format"; tail -20 /tmp/lstar_v3.log; exit 1; }

echo "== Zarr v3 sharding (lstar writes sharded v3; C++/zarr-python/WASM read == unsharded) =="
bash conformance/v3_shard.sh >/tmp/lstar_v3sh.log 2>&1 \
  && pass "v3 sharded write + read conformant across C++/zarr-python/WASM" || { echo "  FAIL v3_shard"; tail -20 /tmp/lstar_v3sh.log; exit 1; }

echo "== single-file .lstar.zarr.zip parity (Python + C++ + CLI + R + JS; STORED + ZIP64) =="
bash conformance/zip.sh >/tmp/lstar_zip.log 2>&1 \
  && pass "single-file .lstar.zarr.zip cross-surface parity" || { echo "  FAIL zip"; tail -25 /tmp/lstar_zip.log; exit 1; }

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

echo "== Zarr v3 on the R surface (R reads both v2+v3; writes v3; zarr-python confirms) =="
bash conformance/v3_r.sh >/tmp/lstar_v3r.log 2>&1 \
  && pass "v3 read+write on the R surface" || { echo "  FAIL v3_r"; tail -20 /tmp/lstar_v3r.log; exit 1; }

echo "== Zstd read parity on the R surface (R reads zarr-python's default zstd v3 store; needs libzstd) =="
bash conformance/zstd_r.sh >/tmp/lstar_zstdr.log 2>&1 \
  && pass "zstd read parity on the R surface" || { echo "  FAIL zstd_r"; tail -20 /tmp/lstar_zstdr.log; exit 1; }

echo "== disk-backed targets (L* -> h5ad -> backed AnnData / Seurat+BPCells / SCE+HDF5Array) =="
bash conformance/backed_targets.sh >/tmp/lstar_bt.log 2>&1 \
  && { pass "disk-backed conversion targets"; grep "SKIP" /tmp/lstar_bt.log | sed 's/^/      /'; } \
  || { echo "  FAIL backed_targets"; tail -15 /tmp/lstar_bt.log; exit 1; }

echo "== JS/WASM (Emscripten kernels + zarrita reader + viewer API; skips if emsdk absent) =="
LSTAR_EMCC_PYTHON="${LSTAR_EMCC_PYTHON:-}" bash conformance/js.sh 2>&1 | sed 's/^/  /'

echo "ALL CONFORMANCE TESTS PASSED"
