#!/usr/bin/env bash
# Viewer-prep per-field compression layout (pagoda3-measured, single-sourced in viewer_policy.json).
# A viewer store is compressed-by-default, per field: the gene-major raw-counts basis stays RAW single-chunk
# (hero gene-color = exact-byte reads) unless --compress-primary; counts_cellmajor is zstd chunked+sharded
# (chunk-granular subset reads via per-chunk decompress); every other array is zstd single-chunk (read
# whole). Asserts the on-disk codecs/chunking on Python (write viewer=True, 3 modes) and JS (extendForViewer
# with the injected WASM writer codec). The JS leg needs node + a WASM build; it skips otherwise.
set -e
ROOT="$(cd "$(dirname "$0")/.." && pwd)"

echo "  -- Python: write(viewer=True) per-field layout (latency / size / none) --"
PYTHONPATH="$ROOT/python/src" python3 - <<'PY'
import warnings; warnings.filterwarnings("ignore")
import numpy as np, scipy.sparse as sp, lstar, json, os, shutil

CHUNK = 16384
def mk():
    # counts nnz > 16384 so counts_cellmajor's data is multi-chunk (-> sharded), exercising the full layout
    ds = lstar.Dataset(kind="sample")
    n, g = 4000, 100
    ds.add_axis("cells", [f"c{i}" for i in range(n)], role="observation")
    ds.add_axis("genes", [f"g{j}" for j in range(g)], role="feature")
    X = sp.random(n, g, density=0.06, format="csc", random_state=1)
    X.data = np.ceil(X.data * 9).astype("f4")
    ds.add_field("counts", X, role="measure", span=["cells", "genes"], state="raw")
    ds.add_field("leiden", [f"cl{i%5}" for i in range(n)], role="label", span=["cells"])
    assert X.nnz > CHUNK, X.nnz
    return ds

def lay(out, field, arr="data"):
    z = json.load(open(f"{out}/fields/{field}/{arr}/zarr.json"))
    c0 = z["codecs"][0]
    sharded = c0["name"] == "sharding_indexed"
    inner = [c["name"] for c in (c0["configuration"]["codecs"] if sharded else z["codecs"])]
    chunk = (c0["configuration"]["chunk_shape"] if sharded else z["chunk_grid"]["configuration"]["chunk_shape"])
    return {"zstd": "zstd" in inner, "sharded": sharded, "chunk0": chunk[0], "shape0": z["shape"][0]}

out = "/tmp/vcomp_py.lstar.zarr"

# latency (default): gene-major raw single-chunk; cell-major zstd+chunked+sharded; dense zstd single-chunk
shutil.rmtree(out, ignore_errors=True); lstar.write(mk(), out, viewer=True)
gm, cm, od = lay(out, "counts"), lay(out, "counts_cellmajor"), lay(out, "od_score", "values")
assert not gm["zstd"] and gm["chunk0"] == gm["shape0"], ("latency gene-major must be raw single-chunk", gm)
assert cm["zstd"] and cm["sharded"] and cm["chunk0"] == 16384, ("cell-major must be zstd+sharded@16384", cm)
assert od["zstd"] and not od["sharded"] and od["chunk0"] == od["shape0"], ("od_score must be zstd single-chunk", od)
assert not [e for e in lstar.validate(lstar.read(out)) if e.startswith("ERROR")]
print("    [py] latency: gene-major raw single-chunk | counts_cellmajor zstd+sharded@16384 | od_score zstd single-chunk; validate clean")

# size (compress_primary): gene-major ALSO zstd (chunked), still range-readable via per-chunk decompress
shutil.rmtree(out, ignore_errors=True); lstar.write(mk(), out, viewer=True, compress_primary=True)
gm = lay(out, "counts")
assert gm["zstd"] and gm["chunk0"] == 16384, ("size gene-major must be zstd@16384", gm)
print("    [py] size (compress_primary): gene-major zstd@16384 (chunked, range-readable)")

# none (compress=False): everything raw single-chunk
shutil.rmtree(out, ignore_errors=True)
ds = mk(); lstar.extend_for_viewer(ds, compress=False); lstar.write(ds, out)
gm, cm = lay(out, "counts"), lay(out, "counts_cellmajor")
assert not gm["zstd"] and not cm["zstd"], ("none must be all-raw", gm, cm)
print("    [py] none (compress=False): all fields raw")
PY

# JS leg: extendForViewer with the injected WASM writer codec produces the same per-field layout.
EMSDK="${EMSDK:-$HOME/emsdk}"
NODE="$(ls -d "$EMSDK"/node/*/bin/node 2>/dev/null | head -1)"; NODE="${NODE:-$(command -v node || true)}"
if [ -n "$NODE" ] && [ -f "$ROOT/js/dist/lstar_writer.mjs" ]; then
  echo "  -- JS: extendForViewer(codec) per-field layout --"
  "$NODE" --experimental-strip-types "$ROOT/js/test/viewer_compression.test.ts" 2>/dev/null \
    || { echo "  FAIL: JS viewer_compression"; exit 1; }
else
  echo "  [skip] node/WASM absent — JS viewer-compression leg skipped"
fi
echo "viewer-compression conformance PASSED."
