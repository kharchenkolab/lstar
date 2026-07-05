#!/usr/bin/env bash
# Seurat -> extend_for_viewer conformance (the previously-untested seam): read a Seurat object into L*,
# run the viewer prep, and assert the store is sane -- a boolean QC flag is NOT a grouping, the viewer
# opens on a real clustering (not a QC/near-boolean field), the navigators are produced -- then have
# Python cross-validate the resulting store (viewer@0.1-clean + the same grouping checks). Synthetic
# Seurat in CI (no downloads); also a real SeuratData object locally. Skips cleanly without Rscript.
set -u
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
PY="${LSTAR_PY:-python3}"
export PYTHONPATH="$ROOT/python/src:$ROOT/python/tests${PYTHONPATH:+:$PYTHONPATH}"
export LSTAR_RLIB="${LSTAR_RLIB:-$ROOT/.Rlib}"
TMP="$(mktemp -d)"; trap 'rm -rf "$TMP"' EXIT
fail=0

if ! command -v Rscript >/dev/null 2>&1; then echo "  [skip] no Rscript — skipping Seurat->viewer parity"; exit 0; fi

OUT="$TMP/seurat_viewer.lstar.zarr"
Rscript "$ROOT/conformance/viewer_seurat.R" "$OUT" "pbmc3k" || fail=1   # real pbmc3k locally; skips in CI

if [ -d "$OUT" ]; then
  "$PY" "$ROOT/conformance/viewer_check.py" validate "$OUT" || fail=1   # viewer@0.1-clean (the canonical validator)
  "$PY" - "$OUT" <<'PYEOF' || fail=1
import sys
from lstar.zarr_io import read
from lstar.viewer import _detect_groupings
ds = read(sys.argv[1])
gp = _detect_groupings(ds)
assert "qc_kept" not in gp, f"Python detects qc_kept (a boolean QC flag) as a grouping on the Seurat store: {gp}"
assert "ident" not in gp, f"Python detects the active-idents mirror ('ident') as a grouping: {gp}"
prim = ds.field("counts_cellmajor_order").provenance.get("group")
assert prim and prim != "qc_kept", f"primary reorder key is {prim!r}"
print(f"  [py] cross-read the Seurat viewer store: primary={prim!r}; qc_kept excluded; groupings={gp}")
PYEOF
else
  echo "  [skip] no store written (Seurat/SeuratObject unavailable in R)"
fi

[ $fail -eq 0 ] && echo "Seurat -> viewer conformance PASSED." || { echo "Seurat -> viewer FAILED."; exit 1; }
