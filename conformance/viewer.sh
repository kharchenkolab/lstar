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

# (c) cross-prep equivalence: pagoda3 write_viewer vs lstar extend_for_viewer on one base store.
if [ -n "${PAGODA3:-}" ] && [ -d "${PAGODA3:-/nonexistent}" ]; then
  echo "  -- (c) lstar-prep == pagoda3-prep (equiv) --"
  "$PY" "$CHK" equiv-pagoda3 "$PAGODA3" "$TMP" || fail=1
else
  echo "  [skip] PAGODA3 unset — skipping leg (c) (cross-prep equivalence)"
fi

if [ "$fail" -eq 0 ]; then echo "  viewer conformance OK"; else echo "  viewer conformance FAIL"; fi
exit "$fail"
