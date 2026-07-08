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
# The JS surface is the WASM/libzarr reader (zarrita retired) — needs only node + a successful build; the
# emcc build itself guards on emsdk/python. LSTAR_EMCC_PYTHON points emcc at a >=3.10 interpreter if needed.
if [ -n "$NODE" ]; then
  if EMSDK="$EMSDK" bash "$ROOT/js/build.sh" >/tmp/v3c_wasm.log 2>&1; then HAVE_JS=1; else echo "  [skip] wasm build failed — JS surface skipped"; fi
fi
HAVE_R=0
if command -v Rscript >/dev/null 2>&1 && Rscript -e '.libPaths(c("'"$RLIB"'", .libPaths())); quit(status=!requireNamespace("lstar", quietly=TRUE))' >/dev/null 2>&1; then HAVE_R=1; fi

# COMPREHENSIVE: every cached AnnData/MuData under testdata/ (the leads first, then the full spatial /
# perturbation / squidpy breadth). Auto-discovered so the sweep grows with the corpus; each is optional +
# tolerant (a dataset that won't load prints SKIP). LSTAR_V3_CORPUS_FAST=1 restricts to the 4 leads (CI/quick).
LEADS=("$TD/pbmc3k_processed.h5ad" "$TD/pancreas.h5ad" "$TD/pancreas_velocity_small.h5ad" "$TD/minipbcite.h5mu")
if [ "${LSTAR_V3_CORPUS_FAST:-0}" = 1 ]; then
  DATASETS=("${LEADS[@]}")
else
  mapfile -t MORE < <(find "$TD" -type f \( -name "*.h5ad" -o -name "*.h5mu" \) 2>/dev/null | sort)
  # Dedup by BASENAME: some datasets are cached under two paths (e.g. spatial/sq_*.h5ad and its
  # spatial/_sqcache/ copy) -- running the same file twice is pure waste, so first-seen basename wins.
  DATASETS=(); declare -A seen
  for f in "${LEADS[@]}" "${MORE[@]}"; do b="$(basename "$f")"; [ -n "${seen[$b]:-}" ] && continue; seen["$b"]=1; DATASETS+=("$f"); done
fi

nok=0; nskip=0
for src in "${DATASETS[@]}"; do
  if [ ! -e "$src" ]; then echo "  SKIP $(basename "$src"): not cached"; nskip=$((nskip+1)); continue; fi
  v2="/tmp/v3corp_v2.lstar.zarr"; v3="/tmp/v3corp_v3.lstar.zarr"; rm -rf "$v2" "$v3" "$v3".down_v2 "$v3".zstd_shard
  out="$(python3 "$ROOT/conformance/v3_corpus.py" "$src" "$v2" "$v3" 2>/tmp/v3corp_py.log)" || { echo "  FAIL $(basename "$src") (python)"; tail -8 /tmp/v3corp_py.log; exit 1; }
  decomp="$(printf '%s\n' "$out" | sed -n 's/^DECOMP_MB=//p')"      # decompressed field bytes (gates the JS leg)
  printf '%s\n' "$out" | grep -v '^DECOMP_MB='
  case "$out" in *SKIP*) nskip=$((nskip+1)); continue;; esac
  if [ -x "$BIN" ]; then
    "$BIN" compare "$v2" "$v3" >/tmp/v3corp_cpp.log 2>&1 && echo "    [C++] libzarr reads both, v2==v3" \
      || { echo "  FAIL $(basename "$src") (C++ compare)"; tail -8 /tmp/v3corp_cpp.log; exit 1; }
  fi
  # JS/WASM: the WASM reader's arrays == zarr-python on real data in BOTH formats, then L*-API v2==v3
  # (wasm_corpus). These paths materialize WHOLE arrays (io_parity dumps every array to one JSON;
  # wasm_corpus reads each into the WASM heap), which doesn't scale to very large real stores — the viewer
  # reads by BYTE RANGE, not whole-array, so this is the wrong tool above a size, not a real gap. Skip the
  # JS whole-array leg on big stores (C++/Py/R validate those); small/medium get the full byte-level check.
  if [ "$HAVE_JS" = 1 ]; then
    if [ "${decomp:-0}" -gt 300 ]; then           # DECOMPRESSED array bytes (io_dump base64s these into one JSON)
      echo "    [JS ] skip (whole-array read doesn't scale to ${decomp}MB decompressed; viewer byte-ranges — C++/Py/R cover this one)"
    else
      python3 "$ROOT/js/test/io_dump.py" "$v2" /tmp/v3corp_d2.json >/dev/null 2>&1
      python3 "$ROOT/js/test/io_dump.py" "$v3" /tmp/v3corp_d3.json >/dev/null 2>&1
      { "$NODE" "$ROOT/js/test/io_parity.mjs" "$v2" /tmp/v3corp_d2.json \
        && "$NODE" "$ROOT/js/test/io_parity.mjs" "$v3" /tmp/v3corp_d3.json \
        && "$NODE" --experimental-strip-types "$ROOT/js/test/wasm_corpus.mjs" "$v2" "$v3"; } >/tmp/v3corp_js.log 2>&1 \
        && echo "    [JS ] libzarr reads real v2+v3 == zarr-python; v2==v3 across the L* API" \
        || { echo "  FAIL $(basename "$src") (JS wasm reader)"; tail -10 /tmp/v3corp_js.log; exit 1; }
    fi
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
