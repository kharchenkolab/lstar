#!/usr/bin/env bash
# R single-file `.lstar.zarr.zip` parity: R reads a Python-written STORED zip, writes its own zip (from
# both a zip and a directory source), and Python reads each back field-for-field, all entries STORED.
# R rides the C++ core's `.zip` dispatch (extract-to-read, pack-to-write) — no R-specific zip code.
# Needs the installed lstar R package (LSTAR_RLIB) + Python. Skips cleanly when R is unavailable.
set -u
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
PY="${LSTAR_PY:-python3}"
export PYTHONPATH="$ROOT/python/src:$ROOT/python/tests${PYTHONPATH:+:$PYTHONPATH}"
export LSTAR_RLIB="${LSTAR_RLIB:-$ROOT/.Rlib}"
TMP="$(mktemp -d)"; trap 'rm -rf "$TMP"' EXIT

if ! command -v Rscript >/dev/null 2>&1; then echo "  [skip] Rscript not found — skipping R zip parity"; exit 0; fi
if ! "$PY" -c "import lstar" 2>/dev/null; then echo "  [skip] lstar not importable — skipping R zip parity"; exit 0; fi

# 1) Python writes a dir + a STORED zip
"$PY" - "$TMP" <<'PY' || { echo "  FAIL: fixture"; exit 1; }
import sys, os, test_zip as T
from lstar.zarr_io import write
tmp = sys.argv[1]; ds = T._make_ds()
write(ds, os.path.join(tmp, "py.lstar.zarr")); write(ds, os.path.join(tmp, "py.lstar.zarr.zip"))
print("  [py ] wrote py.lstar.zarr(.zip)")
PY

# 2) R: read the Python zip, write R zips from both a zip source and a directory source
cat > "$TMP/zip_r.R" <<'RS'
a <- commandArgs(trailingOnly=TRUE); tmp <- a[1]
suppressMessages(library(lstar, lib.loc=Sys.getenv("LSTAR_RLIB")))
ds  <- lstar_read(file.path(tmp, "py.lstar.zarr.zip"))
cat("  [R ] read py.zip:", length(ds$axes), "axes,", length(ds$fields), "fields\n")
lstar_write(ds, file.path(tmp, "r_from_zip.lstar.zarr.zip"))
ds2 <- lstar_read(file.path(tmp, "py.lstar.zarr"))
lstar_write(ds2, file.path(tmp, "r_from_dir.lstar.zarr.zip"))
invisible(lstar_read(file.path(tmp, "r_from_zip.lstar.zarr.zip")))  # R reads its own zip
cat("  [R ] wrote r_from_zip.zip + r_from_dir.zip (and re-read its own)\n")
RS
Rscript "$TMP/zip_r.R" "$TMP" || { echo "  FAIL: R leg"; exit 1; }

# 3) Python cross-reads the R zips: all-STORED + field-identical to the original
"$PY" - "$TMP" <<'PY' || { echo "  FAIL: Python cross-read of R zip"; exit 1; }
import sys, os, test_zip as T
from lstar.zarr_io import read
tmp = sys.argv[1]; d = os.path.join(tmp, "py.lstar.zarr")
for name in ("r_from_zip", "r_from_dir"):
    z = os.path.join(tmp, f"{name}.lstar.zarr.zip")
    assert T._all_entries_stored(z), f"{name}: R zip not all-STORED"
    T._assert_ds_equal(read(d), read(z), ctx=f"py.dir vs {name}")
    print(f"  [py ] {name}: all-STORED + field-identical  OK")
PY

echo "R zip parity PASSED."
