#!/usr/bin/env bash
# JS single-file `.lstar.zarr.zip` parity (no WASM needed — pure TS reader/writer over zarrita + node fs):
#   1. ZipStore reads a Python-written STORED zip byte-for-byte like the directory store (seek-into-zip),
#   2. JS writes a STORED zip that Python reads back field-for-field, all entries STORED,
#   3. guardrails: a DEFLATE-packed store is rejected; a ZIP64 STORED archive reads.
# Skips cleanly (exit 0) when node or zarrita is unavailable. Node comes from emsdk if present, else PATH.
set -u
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
PY="${LSTAR_PY:-python3}"
export PYTHONPATH="$ROOT/python/src:$ROOT/python/tests${PYTHONPATH:+:$PYTHONPATH}"
TMP="$(mktemp -d)"; trap 'rm -rf "$TMP"' EXIT

if ! "$PY" -c "import lstar" 2>/dev/null; then
  echo "  [skip] lstar not importable by '$PY' (set LSTAR_PY) — skipping JS zip parity"; exit 0; fi
if [ ! -d "$ROOT/js/node_modules/zarrita" ]; then
  echo "  [skip] zarrita not installed (cd js && npm install) — skipping JS zip parity"; exit 0; fi
NODE="$(ls -d "${EMSDK:-$HOME/emsdk}"/node/*/bin/node 2>/dev/null | head -1)"; NODE="${NODE:-node}"
if ! command -v "$NODE" >/dev/null 2>&1; then
  echo "  [skip] node not found — skipping JS zip parity"; exit 0; fi

DIR="$TMP/py.lstar.zarr"; PYZIP="$TMP/py.lstar.zarr.zip"
DEFZIP="$TMP/deflate.lstar.zarr.zip"; Z64ZIP="$TMP/z64.lstar.zarr.zip"; OUTZIP="$TMP/js.lstar.zarr.zip"

# 1) Python builds the fixtures: a dir store, a STORED zip, a DEFLATE zip, and a (forced) ZIP64 STORED zip
"$PY" - "$DIR" "$PYZIP" "$DEFZIP" "$Z64ZIP" <<'PY' || { echo "  FAIL: fixture build"; exit 1; }
import sys, os, zipfile, test_zip as T
from lstar.zarr_io import write
d, pyzip, defzip, z64zip = sys.argv[1:5]
ds = T._make_ds()
write(ds, d); write(ds, pyzip)
def pack(zp, comp, limit=None):
    old = zipfile.ZIP64_LIMIT
    if limit is not None: zipfile.ZIP64_LIMIT = limit
    try:
        with zipfile.ZipFile(zp, "w", comp, allowZip64=True) as zf:
            for root,_,files in os.walk(d):
                for fn in files:
                    fp=os.path.join(root,fn); zf.write(fp, arcname=os.path.relpath(fp,d))
    finally:
        zipfile.ZIP64_LIMIT = old
pack(defzip, zipfile.ZIP_DEFLATED)
pack(z64zip, zipfile.ZIP_STORED, limit=4)   # force ZIP64 structures on small entries
print("  [py ] fixtures: dir + STORED zip + DEFLATE zip + ZIP64 zip")
PY

# 2) Node: read parity + write outzip + guardrails
"$NODE" --experimental-strip-types "$ROOT/js/test/zip.test.ts" "$DIR" "$PYZIP" "$OUTZIP" "$DEFZIP" "$Z64ZIP" \
  || { echo "  FAIL: JS zip test"; exit 1; }

# 2b) store-backend VALUE parity: read real field values through openLstar across FS-dir / FS-zip /
# HTTP-dir / HTTP-zip and assert they are EQUAL (not just that each opens). This is the consumer layer —
# it catches a backend that returns empty chunks (opens fine, data silently collapses to zeros).
"$NODE" --experimental-strip-types "$ROOT/js/test/store_backends.test.ts" "$DIR" "$PYZIP" \
  || { echo "  FAIL: store-backend value parity"; exit 1; }

# 2c) throughput parity: the zip must read a first-screen of fields with the same concurrency and ~the
# same round-trips as the directory store — no per-read local-header hop (which would ~double requests
# and serialize the browser's cold open). Self-contained (builds its own multi-field fixture).
"$NODE" --experimental-strip-types "$ROOT/js/test/zip_concurrency.test.ts" \
  || { echo "  FAIL: zip throughput/concurrency parity"; exit 1; }

# 2d) httpZipSource must opt out of the browser same-URL cache lock (cache:"no-store"), or a hosted zip's
# concurrent reads serialize behind one cache entry. Browser-only effect; guard the intent via mock fetch.
"$NODE" --experimental-strip-types "$ROOT/js/test/httpzip_nostore.test.ts" \
  || { echo "  FAIL: httpZipSource no-store cache"; exit 1; }

# 3) Python cross-reads the JS-written zip: all-STORED + expected fields/values
"$PY" - "$OUTZIP" <<'PY' || { echo "  FAIL: Python cross-read of JS zip"; exit 1; }
import sys, numpy as np, test_zip as T
from lstar.zarr_io import read
z = sys.argv[1]
assert T._all_entries_stored(z), "JS-written zip is not all-STORED"
ds = read(z)
assert set(ds.fields) == {"counts","umap","leiden"}, f"fields={set(ds.fields)}"
assert list(map(str, ds.axis("cells").labels)) == ["c0","c1","c2","c3"]
assert list(map(str, ds.field("leiden").values)) == ["a","b","a","b"]
umap = np.asarray(ds.field("umap").values)
assert umap.shape == (4,3) and np.allclose(umap.ravel(), np.arange(12)), "umap mismatch"
c = ds.field("counts").values.tocsc()
assert c.shape == (4,3) and c.nnz == 4, f"counts {c.shape} nnz={c.nnz}"
print("  [py ] cross-read JS zip: all-STORED + fields/values OK")
PY

echo "JS zip parity PASSED."
