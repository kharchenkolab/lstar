#!/usr/bin/env bash
# viewer@0.1 cross-language conformance (docs/format.md "The viewer profile").
#   (a) lstar Python `extend_for_viewer` produces a spec-clean store.
#   (b) lstar R `profile_pagoda2` output validates against the spec.
#   (c) lstar-prep == pagoda3-prep on the viewer fields (added once both call the shared recipe).
# Env: LSTAR_PY (a python with lstar importable; default python3), LSTAR_RLIB (R lib with lstar),
#      PAGODA3 (path to the lstar-viewer checkout; enables leg (c)).
# Skips cleanly when a runtime is unavailable. Wired into run.sh only once all enabled legs are GREEN.
set -u
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
PY="${LSTAR_PY:-python3}"
export PYTHONPATH="$ROOT/python/src${PYTHONPATH:+:$PYTHONPATH}"
export LSTAR_RLIB="${LSTAR_RLIB:-$ROOT/.Rlib}"
CHK="$ROOT/conformance/viewer_check.py"
TMP="$(mktemp -d)"; trap 'rm -rf "$TMP"' EXIT
fail=0

if ! "$PY" -c "import lstar" 2>/dev/null; then
  echo "  [skip] lstar not importable by '$PY' (set LSTAR_PY) — skipping viewer conformance"; exit 0
fi

echo "  -- (a) lstar python extend_for_viewer is viewer@0.1-clean --"
"$PY" "$CHK" canonical || fail=1

if command -v Rscript >/dev/null 2>&1; then
  echo "  -- (b) lstar R profile_pagoda2 -> store -> python validate against the spec --"
  if Rscript "$ROOT/conformance/viewer_pagoda2.R" "$TMP/pagoda2.lstar.zarr" >/tmp/lstar_viewer_r.log 2>&1; then
    "$PY" "$CHK" validate "$TMP/pagoda2.lstar.zarr" || fail=1
  else
    echo "  FAIL: R profile_pagoda2 driver errored"; tail -8 /tmp/lstar_viewer_r.log; fail=1
  fi
else
  echo "  [skip] no Rscript — skipping leg (b)"
fi

# cross-prep equivalence: build ONE base store, lstar-prep it, and compare other producers against it
# on the semantic fields (stats/markers/od). (c) pagoda3 prep.ts (WASM); (d) native R extend_for_viewer.
HAVE_P3=0; [ -n "${PAGODA3:-}" ] && [ -f "${PAGODA3:-/nonexistent}/prep/prep.ts" ] && HAVE_P3=1
HAVE_R=0; command -v Rscript >/dev/null 2>&1 && HAVE_R=1
if [ "$HAVE_P3" = 1 ] || [ "$HAVE_R" = 1 ]; then
  "$PY" "$CHK" make-base "$TMP/base.lstar.zarr" || fail=1
  "$PY" "$CHK" prep-lstar "$TMP/base.lstar.zarr" "$TMP/lstar.lstar.zarr" || fail=1
fi

if [ "$HAVE_P3" = 1 ]; then
  NODE="$(ls -d "$HOME"/emsdk/node/*/bin/node 2>/dev/null | head -1)"; NODE="${NODE:-node}"
  if [ -f "$ROOT/js/dist/lstar_kernels.mjs" ] && command -v "$NODE" >/dev/null 2>&1; then
    echo "  -- (c) lstar-prep == pagoda3 prep.ts (equiv) --"
    cp -r "$TMP/base.lstar.zarr" "$TMP/p3.lstar.zarr"
    if "$NODE" --experimental-strip-types "$PAGODA3/prep/prep.ts" "$TMP/p3.lstar.zarr" leiden >/tmp/lstar_p3prep.log 2>&1; then
      "$PY" "$CHK" equiv "$TMP/lstar.lstar.zarr" "$TMP/p3.lstar.zarr" || fail=1
    else echo "  FAIL: pagoda3 prep.ts errored"; tail -10 /tmp/lstar_p3prep.log; fail=1; fi
  else echo "  [skip] no node / WASM dist — skipping leg (c)"; fi
else echo "  [skip] PAGODA3 unset or prep.ts absent — skipping leg (c)"; fi

if [ "$HAVE_R" = 1 ]; then
  echo "  -- (d) lstar-prep == native R extend_for_viewer (equiv) --"
  if Rscript "$ROOT/conformance/viewer_extend_r.R" "$TMP/base.lstar.zarr" "$TMP/r.lstar.zarr" >/tmp/lstar_rprep.log 2>&1; then
    "$PY" "$CHK" equiv "$TMP/lstar.lstar.zarr" "$TMP/r.lstar.zarr" || fail=1
  else echo "  FAIL: R extend_for_viewer errored"; tail -10 /tmp/lstar_rprep.log; fail=1; fi
else echo "  [skip] no Rscript — skipping leg (d)"; fi

if [ "$fail" -eq 0 ]; then echo "  viewer conformance OK"; else echo "  viewer conformance FAIL"; fi
exit "$fail"
