#!/usr/bin/env bash
# Single-file `.lstar.zarr.zip` (STORED) cross-surface parity. A store round-trips through a zip
# identically to a directory, on every surface, and the artifact is always STORED (never DEFLATE):
#   Python  — write/read .zip (forces STORED; rejects DEFLATE; lazy+stream)      [test_zip.py]
#   C++/R   — extract-to-read + pack-to-write via the STORED-zip codec           [this script + zip_r.sh]
#   JS      — seek-into-zip ZipStore (HTTP Range) + zip writer                   [zip_js.sh]
#   CLI     — convert dir<->zip repackage + `--viewer` single-file viewer store
# This driver covers Python + C++ + the CLI + the guardrails (DEFLATE rejected, ZIP64 read), then runs
# the R and JS legs opportunistically when those runtimes are present (in CI each job runs the leg it
# has: python-cpp -> this; r-cross-format -> zip_r.sh; js-wasm -> zip_js.sh). Skips a surface cleanly.
set -u
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
PY="${LSTAR_PY:-python3}"
export PYTHONPATH="$ROOT/python/src:$ROOT/python/tests${PYTHONPATH:+:$PYTHONPATH}"
BIN="$ROOT/core/build/test_zip"
TMP="$(mktemp -d)"; trap 'rm -rf "$TMP"' EXIT
fail=0

if ! "$PY" -c "import lstar" 2>/dev/null; then
  echo "  [skip] lstar not importable by '$PY' (set LSTAR_PY) — skipping zip parity"; exit 0; fi

# ---- Python surface: the reference write/read + guardrails + lazy/stream ----
echo "  -- python (test_zip.py) --"
"$PY" "$ROOT/python/tests/test_zip.py" || { echo "  FAIL: python test_zip.py"; fail=1; }

# ---- C++ surface + CLI + guardrails (needs the built test_zip binary) ----
if [ -x "$BIN" ]; then
  echo "  -- C++ matrix + guardrails + CLI --"
  "$PY" - "$BIN" <<'PY' || fail=1
import sys, os, tempfile, subprocess, zipfile
import test_zip as T   # python/tests is on PYTHONPATH
from lstar.zarr_io import write, read
from lstar.cli import main
BIN = sys.argv[1]
tmp = tempfile.mkdtemp()
d   = os.path.join(tmp, "py.lstar.zarr");     write(T._make_ds(), d)
pyz = os.path.join(tmp, "py.lstar.zarr.zip"); write(T._make_ds(), pyz)
def cpp(i, o):
    r = subprocess.run([BIN, i, o], capture_output=True, text=True)
    if r.returncode != 0: raise SystemExit("C++ test_zip failed:\n" + r.stderr)
# (A) Py.zip -> C++ -> C++.zip -> Py ; (B) Py.dir -> C++ -> C++.zip -> Py ; (C) Py.zip -> C++ -> C++.dir -> Py
for tag, src in (("A", pyz), ("B", d)):
    o = os.path.join(tmp, f"cpp_{tag}.lstar.zarr.zip"); cpp(src, o)
    assert T._all_entries_stored(o), f"({tag}) C++ zip not all-STORED"
    T._assert_ds_equal(read(d), read(o), ctx=f"cpp-{tag}")
cdir = os.path.join(tmp, "cpp_C.lstar.zarr"); cpp(pyz, cdir)
T._assert_ds_equal(read(d), read(cdir), ctx="cpp-C")
print("  [c++] A/B/C matrix: Py.zip<->C++<->Py + all-STORED  OK")
# (D) DEFLATE rejected by C++ ; (E) ZIP64 read by C++
def pack(zp, comp, limit=None):
    old = zipfile.ZIP64_LIMIT
    if limit is not None: zipfile.ZIP64_LIMIT = limit
    try:
        with zipfile.ZipFile(zp, "w", comp, allowZip64=True) as zf:
            for r,_,fs in os.walk(d):
                for fn in fs: fp=os.path.join(r,fn); zf.write(fp, arcname=os.path.relpath(fp,d))
    finally: zipfile.ZIP64_LIMIT = old
bad = os.path.join(tmp, "bad.lstar.zarr.zip"); pack(bad, zipfile.ZIP_DEFLATED)
r = subprocess.run([BIN, bad, os.path.join(tmp,"x.lstar.zarr")], capture_output=True, text=True)
assert r.returncode != 0 and ("stored" in (r.stderr+r.stdout).lower() or "deflate" in (r.stderr+r.stdout).lower()), \
    "C++ must reject DEFLATE with an actionable message"
z64 = os.path.join(tmp, "z64.lstar.zarr.zip"); pack(z64, zipfile.ZIP_STORED, limit=4)
o64 = os.path.join(tmp, "z64_out.lstar.zarr.zip"); cpp(z64, o64)
T._assert_ds_equal(read(d), read(o64), ctx="cpp-zip64")
print("  [c++] guardrails: DEFLATE rejected + ZIP64 read  OK")
# CLI: repackage dir<->zip + --viewer single-file
z = os.path.join(tmp, "cli.lstar.zarr.zip"); assert main(["convert", d, z]) == 0
assert T._all_entries_stored(z); T._assert_ds_equal(read(d), read(z), ctx="cli-repack")
d2 = os.path.join(tmp, "cli_out.lstar.zarr"); assert main(["convert", z, d2]) == 0
T._assert_ds_equal(read(d), read(d2), ctx="cli-unpack")
vz = os.path.join(tmp, "cli_viewer.lstar.zarr.zip"); assert main(["convert", d, vz, "--viewer"]) == 0
assert T._all_entries_stored(vz)
assert any(n.startswith("stats_") or "counts_cellmajor" in n for n in read(vz).fields), "viewer fields missing"
print("  [cli] repackage dir<->zip + --viewer single-file store  OK")
PY
else
  echo "  [skip] C++ test_zip not built ($BIN) — build with: cmake --build core/build --target test_zip"
fi

# ---- opportunistic R + JS legs (each CI job otherwise runs its own) ----
if command -v Rscript >/dev/null 2>&1 && [ -d "${LSTAR_RLIB:-$ROOT/.Rlib}" ]; then
  echo "  -- R (zip_r.sh) --"; LSTAR_RLIB="${LSTAR_RLIB:-$ROOT/.Rlib}" bash "$ROOT/conformance/zip_r.sh" || fail=1
fi
NODE="$(ls -d "${EMSDK:-$HOME/emsdk}"/node/*/bin/node 2>/dev/null | head -1)"; NODE="${NODE:-node}"
if [ -d "$ROOT/js/node_modules/zarrita" ] && command -v "$NODE" >/dev/null 2>&1; then
  echo "  -- JS (zip_js.sh) --"; bash "$ROOT/conformance/zip_js.sh" || fail=1
fi

[ $fail -eq 0 ] && echo "zip cross-surface parity PASSED." || { echo "zip parity FAILED."; exit 1; }
