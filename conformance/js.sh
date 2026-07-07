#!/usr/bin/env bash
# JS/WASM conformance: build the WASM modules (compute kernels + the libzarr I/O reader), then verify (in
# Node) the kernels, the libzarr-backed L* reader, and the viewer query API against Python-written stores +
# references. The reader has NO JS Zarr dependency — it reads v2 and v3 through the same libzarr core as
# R/Python. Skips cleanly when emsdk is absent. Point LSTAR_EMCC_PYTHON at a >=3.10 interpreter if the
# system python is older.
# Origin coverage: Py-authored ✓ (JS reads) | JS-authored ✓ (writer_make.ts emits every encoding, Python
# cross-reads via writer_crossread.py) — see conformance/README.md
set -e
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
EMSDK="${EMSDK:-$HOME/emsdk}"
if [ ! -f "$EMSDK/emsdk_env.sh" ] && ! command -v emcc >/dev/null 2>&1; then
  echo "  [skip] no emcc (emsdk at $EMSDK absent and emcc not on PATH) — skipping JS/WASM conformance"
  exit 0
fi
NODE="$(ls -d "$EMSDK"/node/*/bin/node 2>/dev/null | head -1)"; NODE="${NODE:-node}"

# 1) build the WASM kernels
EMSDK="$EMSDK" bash "$ROOT/js/build.sh" >/tmp/lstar_wasm_build.log 2>&1 \
  || { echo "  FAIL: wasm build"; tail -15 /tmp/lstar_wasm_build.log; exit 1; }
echo "  [wasm] built kernels"

# 2) generate the test store (Python L* writer)
PYTHONPATH="$ROOT/python/src" python3 "$ROOT/js/test/make_store.py" >/tmp/lstar_js_store.log 2>&1 \
  || { echo "  FAIL: test store generation"; tail -15 /tmp/lstar_js_store.log; exit 1; }
echo "  [py ] generated test store + references"

# 3) run the Node tests (kernels: dense+Python conformance; reader: manifest+fields; view: queries)
echo "  -- kernels --"; "$NODE" "$ROOT/js/test/kernels.test.mjs"
echo "  -- reader  --"; "$NODE" --experimental-strip-types "$ROOT/js/test/reader.test.ts" 2>/dev/null

# libzarr reader (the sole reader): on a maximal store, whole-array reads == zarr-python for BOTH v2 and
# v3 (gzip), and the full L* reader API (openLstar) reads v2 and v3 to identical values across every
# encoding. The v3 store is written by zarr-python (an independent writer).
echo "  -- libzarr io (v2+v3) --"
JV2=/tmp/lstar_js_v2.lstar.zarr; JV3=/tmp/lstar_js_v3.lstar.zarr
PYTHONPATH="$ROOT/python/src" python3 "$ROOT/conformance/v3_gen.py" "$JV2" >/dev/null
PYTHONPATH="$ROOT/python/src" python3 "$ROOT/conformance/v3_verify.py" reemit_gzip "$JV2" "$JV3" >/dev/null
PYTHONPATH="$ROOT/python/src" python3 "$ROOT/js/test/io_dump.py" "$JV2" /tmp/lstar_js_v2.json >/dev/null
PYTHONPATH="$ROOT/python/src" python3 "$ROOT/js/test/io_dump.py" "$JV3" /tmp/lstar_js_v3.json >/dev/null
"$NODE" "$ROOT/js/test/io_parity.mjs" "$JV2" /tmp/lstar_js_v2.json               # arrays == zarr-python (v2)
"$NODE" "$ROOT/js/test/io_parity.mjs" "$JV3" /tmp/lstar_js_v3.json               # arrays == zarr-python (v3)
"$NODE" --experimental-strip-types "$ROOT/js/test/wasm_corpus.mjs" "$JV2" "$JV3" 2>/dev/null   # L* API: v2 == v3

# v3-aware byte-range fast path (the viewer's hot path): on UNCOMPRESSED v2 + v3 stores, cscColumn/csrRow
# issue exact byte-range reads on BOTH formats (v3 chunk keys c/0) and return identical values.
BV2=/tmp/lstar_js_brv2.lstar.zarr; BV3=/tmp/lstar_js_brv3.lstar.zarr
PYTHONPATH="$ROOT/python/src" python3 - "$BV2" "$BV3" <<'PY'
import sys, warnings; warnings.filterwarnings("ignore")
import lstar, numpy as np, scipy.sparse as sp
n, g = 120, 50
ds = lstar.Dataset(kind="sample")
ds.add_axis("cells", [f"c{i}" for i in range(n)]); ds.add_axis("genes", [f"g{j}" for j in range(g)])
C = sp.random(n, g, density=0.2, format="csc", random_state=1); C.data = np.ceil(C.data * 9).astype("float32")
ds.add_field("counts", C, role="measure", span=["cells", "genes"], state="raw")
R = sp.random(n, g, density=0.2, format="csr", random_state=2).astype("float32")
ds.add_field("counts_cellmajor", R, role="measure", span=["cells", "genes"], state="raw")
lstar.write(ds, sys.argv[1], format="v2")            # uncompressed (byte-range fast path needs raw chunks)
lstar.write(ds, sys.argv[2], format="v3")
PY
"$NODE" --experimental-strip-types "$ROOT/js/test/v3_range.mjs" "$BV2" "$BV3"

echo "  -- parallel--"; "$NODE" --experimental-strip-types "$ROOT/js/test/reader_parallel.test.ts" 2>/dev/null   # independent component arrays read concurrently
echo "  -- view    --"; "$NODE" --experimental-strip-types "$ROOT/js/test/view.test.ts" 2>/dev/null
echo "  -- enc-inv --"; "$NODE" --experimental-strip-types "$ROOT/js/test/encoding_invariance.test.ts" 2>/dev/null   # live-viewer compute is encoding-invariant (dense/csc/csr) — guards the sparse-hardcoded measure-read regression
echo "  -- extend  --"; "$NODE" --experimental-strip-types "$ROOT/js/test/extend_primary.test.ts" 2>/dev/null   # extend_for_viewer(primary=): hoist + compose + validation

# 4) the WRITE side: JS round-trip (writer.test.ts), then the cross-language gate -- JS writes a
# chunked + gzip-compressed store with every encoding (CSC/dense/categorical/mask/partial/aux, using the
# WASM zlib kernel) and the Python reader validates it clean + value-equal.
echo "  -- writer  --"; "$NODE" --experimental-strip-types "$ROOT/js/test/writer.test.ts" 2>/dev/null
echo "  -- writer x-read (JS-write chunked+gzip -> Python read+validate) --"
WX=/tmp/lstar_js_writer_cross.lstar.zarr; rm -rf "$WX"
"$NODE" --experimental-strip-types "$ROOT/js/test/writer_make.ts" "$WX" 2>/dev/null \
  || { echo "  FAIL: JS writer_make"; exit 1; }
PYTHONPATH="$ROOT/python/src" python3 "$ROOT/js/test/writer_crossread.py" "$WX" \
  || { echo "  FAIL: writer_crossread"; exit 1; }
# same, but the JS writer emits Zarr v3 (per-node zarr.json + inline consolidated + c/ chunk keys) --
# proves the JS v3 emitter is interoperable with zarr-python (the flip prerequisite).
echo "  -- writer x-read v3 (JS-write v3 -> Python read+validate) --"
WX3=/tmp/lstar_js_writer_cross_v3.lstar.zarr; rm -rf "$WX3"
"$NODE" --experimental-strip-types "$ROOT/js/test/writer_make.ts" "$WX3" v3 2>/dev/null \
  || { echo "  FAIL: JS writer_make v3"; exit 1; }
PYTHONPATH="$ROOT/python/src" python3 "$ROOT/js/test/writer_crossread.py" "$WX3" \
  || { echo "  FAIL: writer_crossread v3"; exit 1; }
# full JS writer parity: SHARDED v3 with a ZSTD inner codec (both via the libzarr WASM writer -- shard::pack
# + CodecPipeline encode). Python reads every encoding back value-equal, and the store is genuinely
# sharding_indexed with zstd inside the shard.
echo "  -- writer x-read v3 SHARDED+zstd (JS-write -> Python read+validate) --"
WXS=/tmp/lstar_js_writer_cross_v3shard.lstar.zarr; rm -rf "$WXS"
"$NODE" --experimental-strip-types "$ROOT/js/test/writer_make.ts" "$WXS" v3 zstd 16 2>/dev/null \
  || { echo "  FAIL: JS writer_make v3 sharded+zstd"; exit 1; }
PYTHONPATH="$ROOT/python/src" python3 "$ROOT/js/test/writer_crossread.py" "$WXS" \
  || { echo "  FAIL: writer_crossread v3 sharded"; exit 1; }
python3 - "$WXS" <<'PY' || { echo "  FAIL: JS sharded store not sharding_indexed+zstd"; exit 1; }
import sys, json
c = json.load(open(sys.argv[1] + "/fields/counts/data/zarr.json"))["codecs"][0]
assert c["name"] == "sharding_indexed", c
inner = [x["name"] for x in c["configuration"]["codecs"]]
assert "zstd" in inner, inner
print("  [py] JS-written v3 store: sharding_indexed with a zstd inner codec, reads == unsharded")
PY

# 5) reverse leg: JS addToStore appends a derived field to a PYTHON-written store -> Python re-reads it
echo "  -- writer extend (Python-write -> JS addToStore -> Python read) --"
WE=/tmp/lstar_js_writer_extend.lstar.zarr; rm -rf "$WE"; cp -r "$ROOT/js/test/data/sample.lstar.zarr" "$WE"
"$NODE" --experimental-strip-types "$ROOT/js/test/writer_extend.ts" "$WE" 2>/dev/null \
  || { echo "  FAIL: JS writer_extend"; exit 1; }
PYTHONPATH="$ROOT/python/src" python3 - "$WE" <<'PY' || { echo "  FAIL: extend re-read"; exit 1; }
import sys, lstar
ds = lstar.read(sys.argv[1])
assert "od_score" in ds.fields and "counts" in ds.fields, list(ds.fields)   # derived added, original kept
assert ds.axis("od_groups").origin == "derived" and "derived@0.1" in ds.profiles
assert not [e for e in lstar.validate(ds) if e.startswith("ERROR")]
print("  [py] JS addToStore on a Python store re-reads clean (derived field + axis + profile)")
PY
# same, but the base store is Zarr v3 -- addToStore must AUTO-DETECT v3 and append v3 nodes (never mix
# formats). This is the exact pagoda3-prep scenario once the default flips to v3.
echo "  -- writer extend v3 (Python-write v3 -> JS addToStore auto-detect -> Python read) --"
WE3=/tmp/lstar_js_writer_extend_v3.lstar.zarr; rm -rf "$WE3"
PYTHONPATH="$ROOT/python/src" python3 -c "import lstar; lstar.write(lstar.read('$ROOT/js/test/data/sample.lstar.zarr'), '$WE3', format='v3')" \
  || { echo "  FAIL: v3 base write"; exit 1; }
"$NODE" --experimental-strip-types "$ROOT/js/test/writer_extend.ts" "$WE3" 2>/dev/null \
  || { echo "  FAIL: JS writer_extend v3"; exit 1; }
[ -f "$WE3/.zmetadata" ] && { echo "  FAIL: addToStore mixed formats (v3 base gained .zmetadata)"; exit 1; }
PYTHONPATH="$ROOT/python/src" python3 - "$WE3" <<'PY' || { echo "  FAIL: extend v3 re-read"; exit 1; }
import sys, lstar
ds = lstar.read(sys.argv[1])
assert "od_score" in ds.fields and "counts" in ds.fields, list(ds.fields)
assert ds.axis("od_groups").origin == "derived" and "derived@0.1" in ds.profiles
assert not [e for e in lstar.validate(ds) if e.startswith("ERROR")]
print("  [py] JS addToStore on a v3 Python store re-reads clean (v3 preserved, no format mixing)")
PY
echo "JS/WASM conformance PASSED."
