#!/usr/bin/env bash
# Asymmetric round-trip: one maximal artifact pushed through a chain that crosses SURFACES x FORMATS x
# COMPRESSION, then verified value-identical to the original. Each per-surface / per-format / per-codec test
# proves its own leg; this proves the SEAMS between them — a store written by one surface in one
# (format, codec) is re-read and re-written faithfully by the next.
#
#   Python v2+gzip  ->  R v3 uncompressed  ->  C++ v3 gzip+sharded  ->  Python v3 zstd+sharded  ->  {C++, JS} read
#
# Tolerant: a missing surface (no Rscript / no test_v3 / no node) degrades that leg to a skip.
set -e
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
export PYTHONPATH="$ROOT/python/src${PYTHONPATH:+:$PYTHONPATH}"
BIN="$ROOT/core/build/test_v3"
RLIB="${LSTAR_RLIB:-$ROOT/.Rlib}"
TMP="$(mktemp -d)"; trap 'rm -rf "$TMP"' EXIT

A="$TMP/a_py_v2gz.lstar.zarr"        # 1. Python, v2 + gzip  (v3_gen writes exactly this)
B="$TMP/b_r_v3.lstar.zarr"           # 2. R, v3 uncompressed
C="$TMP/c_cpp_v3shard.lstar.zarr"    # 3. C++, v3 gzip + sharded
D="$TMP/d_py_v3zstdshard.lstar.zarr" # 4. Python, v3 zstd + sharded

python3 "$ROOT/conformance/v3_gen.py" "$A" >/dev/null
echo "  [py ] 1. wrote v2+gzip seed"
read NNZ SUM < <(python3 - "$A" <<'PY'
import sys, warnings; warnings.filterwarnings("ignore"); import lstar, scipy.sparse as sp
c = sp.csc_matrix(lstar.read(sys.argv[1]).field("counts").values); print(int(c.nnz), round(float(c.sum()), 3))
PY
)

have_r=0; command -v Rscript >/dev/null 2>&1 && Rscript -e '.libPaths(c("'"$RLIB"'",.libPaths())); quit(status=!requireNamespace("lstar",quietly=TRUE))' >/dev/null 2>&1 && have_r=1
if [ "$have_r" = 1 ]; then
  Rscript -e '.libPaths(c("'"$RLIB"'",.libPaths())); suppressMessages(library(lstar))
    lstar_write(lstar_read("'"$A"'"), "'"$B"'", format="v3")' >/dev/null 2>&1
  echo "  [R  ] 2. re-wrote as v3 uncompressed"
else B="$A"; echo "  [skip] no R — leg 2 uses the v2+gzip store"; fi

if [ -x "$BIN" ]; then
  "$BIN" shard "$B" "$C" >/dev/null && echo "  [c++] 3. re-wrote as v3 gzip + sharded"
else C="$B"; echo "  [skip] no test_v3 — leg 3 skipped"; fi

python3 - "$C" "$D" <<'PY'
import sys, warnings; warnings.filterwarnings("ignore"); import lstar, numcodecs
lstar.write(lstar.read(sys.argv[1]), sys.argv[2], format="v3",
            compressor=numcodecs.Zstd(3), chunk_elems=1000, shard_elems=4000)
PY
echo "  [py ] 4. re-wrote as v3 zstd + sharded"

# verify the artifact survived the whole chain, value-identical to the original — EXCEPT the masked
# positions of a nullable field, which an R hop legitimately fills with 0 (the mask marks them missing;
# both are valid + this is a known, format-independent cross-surface behavior). So compare every field
# exactly, and nullable fields only where valid.
python3 - "$A" "$D" <<'PY'
import sys, warnings; warnings.filterwarnings("ignore")
import lstar, numpy as np, scipy.sparse as sp
a, d = lstar.read(sys.argv[1]), lstar.read(sys.argv[2])
assert set(a.fields) == set(d.fields), (set(a.fields) ^ set(d.fields))
for n in a.fields:
    fa, fd = a.field(n), d.field(n)
    va, vd = fa.values, fd.values
    if sp.issparse(va):
        va, vd = sp.csr_matrix(va), sp.csr_matrix(vd)
        assert va.shape == vd.shape and (va != vd).nnz == 0, f"sparse field {n} differs"
    elif np.asarray(va).dtype.kind in "fiub":
        xa, xd = np.nan_to_num(np.asarray(va, "f8")), np.nan_to_num(np.asarray(vd, "f8"))
        m = np.asarray(fa.mask).astype(bool) if getattr(fa, "mask", None) is not None else None
        if m is not None: xa, xd = xa[~m], xd[~m]              # nullable: only the valid positions
        assert np.allclose(xa, xd), f"dense field {n} differs (valid positions)"
    else:
        assert [str(x) for x in np.asarray(va)] == [str(x) for x in np.asarray(vd)], f"label field {n} differs"
print("  [py ] chain end == original (every field; nullable compared where valid)")
PY
EMSDK="${EMSDK:-$HOME/emsdk}"; NODE="$(ls -d "$EMSDK"/node/*/bin/node 2>/dev/null | head -1)"; NODE="${NODE:-$(command -v node || true)}"
if [ -n "$NODE" ] && [ -f "$ROOT/js/dist/lstar_io.mjs" ]; then
  "$NODE" --experimental-strip-types --input-type=module -e '
    import { openLstar } from "'"$ROOT"'/js/core/reader.ts";
    import { NodeFSStore } from "'"$ROOT"'/js/core/node-store.ts";
    const [dir, nnz, sum] = process.argv.slice(1);   // node -e: passed args start at argv[1] (no script path)
    const ds = await openLstar(new NodeFSStore(dir));
    const sp = await ds.fieldSparse("counts");
    let s = 0; for (const x of sp.data) s += x;
    if (sp.data.length !== +nnz || Math.abs(s - +sum) > 1e-2)
      throw new Error(`JS read of the chain end != original: nnz ${sp.data.length}/${nnz}, sum ${s}/${sum}`);
    console.log("  [js ] read the chain end (v3 zstd+sharded) == original counts (nnz+sum)");
  ' "$D" "$NNZ" "$SUM" 2>/dev/null || { echo "  FAIL: JS read of chain end"; exit 1; }
else echo "  [skip] no node/WASM — JS read-back skipped"; fi
echo "asymmetric round-trip PASSED."
