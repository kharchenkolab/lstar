#!/usr/bin/env bash
# Corpus-driven viewer@0.1 parity: convert each curated corpus dataset to an L* store, prep it on every
# AVAILABLE surface (Python, R, and JS when a WASM dist is present), and assert cross-surface equivalence
# on ALL viewer fields. Exercises the prep on REALISTIC inputs -- real categoricals (non-alphabetical,
# several competing groupings), real embeddings, native CSR -- the structure the toy fixture can't reach.
# Real data locally; synthetic-but-faithful in CI (set LSTAR_SYNTHETIC_CORPUS=1; no real data committed).
# Skips a dataset cleanly when its loader is unavailable; skips a surface when its runtime is absent.
set -u
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
PY="${LSTAR_PY:-python3}"
export PYTHONPATH="$ROOT/python/src:$ROOT/python/tests${PYTHONPATH:+:$PYTHONPATH}"
export LSTAR_RLIB="${LSTAR_RLIB:-$ROOT/.Rlib}"
CHK="$ROOT/conformance/viewer_check.py"
COR="$ROOT/conformance/viewer_corpus.py"
TMP="$(mktemp -d)"; trap 'rm -rf "$TMP"' EXIT
fail=0; ran=0

if ! "$PY" -c "import lstar" 2>/dev/null; then echo "  [skip] lstar not importable — skipping corpus parity"; exit 0; fi
HAVE_R=0; command -v Rscript >/dev/null 2>&1 && HAVE_R=1
NODE="$(ls -d "${EMSDK:-$HOME/emsdk}"/node/*/bin/node 2>/dev/null | head -1)"; NODE="${NODE:-node}"
HAVE_JS=0; [ -f "$ROOT/js/dist/lstar_kernels.mjs" ] && [ -d "$ROOT/js/node_modules/zarrita" ] && command -v "$NODE" >/dev/null 2>&1 && HAVE_JS=1
echo "  surfaces: python$( [ "$HAVE_R" = 1 ] && echo ' + R' )$( [ "$HAVE_JS" = 1 ] && echo ' + JS' )"

for name in $("$PY" "$COR" datasets); do
  base="$TMP/$name.base.zarr"
  out="$("$PY" "$COR" base "$name" "$base" 2>&1)"; rc=$?
  echo "$out" | grep -v '^BASIS='
  [ $rc -eq 7 ] && continue                              # dataset unavailable -> clean skip
  [ $rc -ne 0 ] && { echo "  FAIL: corpus base $name"; fail=1; continue; }
  basis="$(printf '%s\n' "$out" | sed -n 's/^BASIS=//p')"
  ran=$((ran + 1))
  echo "  == $name (basis=${basis:-raw}) =="
  "$PY" "$CHK" prep-lstar "$base" "$TMP/$name.py.zarr" "$basis" >/dev/null || { fail=1; continue; }
  if [ "$ran" = 1 ]; then    # once: the per-field viewer compression default fires on this (real/faithful) data
    "$PY" - "$TMP/$name.py.zarr" <<'PY' || { echo "  FAIL: viewer compression layout on corpus data"; fail=1; }
import sys, os, json
z = json.load(open(os.path.join(sys.argv[1], "fields/counts_cellmajor/data/zarr.json")))["codecs"][0]
sharded = z["name"] == "sharding_indexed"
inner = z["configuration"]["codecs"] if sharded else json.load(open(os.path.join(sys.argv[1], "fields/counts_cellmajor/data/zarr.json")))["codecs"]
assert any(c["name"] == "zstd" for c in inner), ("counts_cellmajor not zstd", [c["name"] for c in inner])
print("  [py] corpus viewer store: counts_cellmajor is zstd%s (per-field compression default fires on real data)"
      % (" + sharded" if sharded else " (single-chunk at this scale)"))
PY
  fi

  if [ "$HAVE_R" = 1 ]; then
    if Rscript "$ROOT/conformance/viewer_extend_r.R" "$base" "$TMP/$name.r.zarr" "$basis" >/tmp/lstar_corpus_r.log 2>&1; then
      echo "  -- $name: python == R --"; "$PY" "$CHK" equiv "$TMP/$name.py.zarr" "$TMP/$name.r.zarr" || fail=1
    else echo "  FAIL: R prep $name"; tail -10 /tmp/lstar_corpus_r.log; fail=1; fi
  fi

  if [ "$HAVE_JS" = 1 ]; then
    cp -r "$base" "$TMP/$name.js.zarr"
    bflag=""; [ -n "$basis" ] && bflag="--basis $basis"
    if "$NODE" --experimental-strip-types "$ROOT/js/scripts/extend-viewer.ts" "$TMP/$name.js.zarr" $bflag >/tmp/lstar_corpus_js.log 2>&1; then
      echo "  -- $name: python == JS --"; "$PY" "$CHK" equiv "$TMP/$name.py.zarr" "$TMP/$name.js.zarr" || fail=1
    else echo "  FAIL: JS prep $name"; tail -12 /tmp/lstar_corpus_js.log; fail=1; fi
  fi
done

if [ "$ran" -eq 0 ]; then echo "  [skip] no corpus datasets available (real corpus is local-only; CI sets LSTAR_SYNTHETIC_CORPUS=1)"; exit 0; fi
if [ "$fail" -eq 0 ]; then echo "  corpus viewer parity OK ($ran dataset(s))"; else echo "  corpus viewer parity FAIL"; fi
exit "$fail"
