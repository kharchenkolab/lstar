#!/usr/bin/env bash
# Format round-trips BOTH directions preserve values: v2->v3->v2 and v3->v2->v3 on the maximal synthetic
# store (every encoding). Complements v3_format.sh, which only tests the v2->v3 (upgrade) direction — this
# adds the v3->v2 DOWNGRADE. Covers the two distinct v2-write paths: the C++ core (also R's, via cpp11) and
# Python (zarr-python). The JS writer's v2+v3 output is covered by js.sh (writer.test.ts writes both).
set -e
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
export PYTHONPATH="$ROOT/python/src${PYTHONPATH:+:$PYTHONPATH}"
BIN="$ROOT/core/build/test_v3"
TMP="$(mktemp -d)"; trap 'rm -rf "$TMP"' EXIT

SEED="$TMP/seed_v2.lstar.zarr"                                   # v2 maximal seed (every encoding)
python3 "$ROOT/conformance/v3_gen.py" "$SEED" >/dev/null
SEED3="$TMP/seed_v3.lstar.zarr"                                  # same, as v3
python3 - "$SEED" "$SEED3" <<'PY'
import sys, lstar; lstar.write(lstar.read(sys.argv[1]), sys.argv[2], format="v3")
PY

if [ -x "$BIN" ]; then
  echo "  -- C++ (libstar core; also R's write path) --"
  "$BIN" write   "$SEED"  "$TMP/c_v3"  >/dev/null                # v2 -> v3
  "$BIN" writev2 "$TMP/c_v3" "$TMP/c_v2b" >/dev/null             # v3 -> v2
  "$BIN" compare "$SEED" "$TMP/c_v2b"  && echo "    [c++] v2->v3->v2 == original"
  "$BIN" writev2 "$SEED3" "$TMP/c_v2"  >/dev/null                # v3 -> v2
  "$BIN" write   "$TMP/c_v2" "$TMP/c_v3b" >/dev/null             # v2 -> v3
  "$BIN" compare "$SEED3" "$TMP/c_v3b" && echo "    [c++] v3->v2->v3 == original"
else
  echo "  [skip] test_v3 not built — C++ leg skipped"
fi

echo "  -- Python (zarr-python writer) --"
python3 - "$SEED" "$SEED3" "$TMP" <<'PY'
import sys, warnings; warnings.filterwarnings("ignore")
import lstar, numpy as np, scipy.sparse as sp
seed, seed3, tmp = sys.argv[1:4]
def vals(ds):                                                   # value signature over every field
    out = {}
    for n, f in ds.fields.items():
        v = f.values
        if sp.issparse(v):
            out[n] = (int(v.nnz), round(float(v.sum()), 3))
        else:
            a = np.asarray(v)
            out[n] = round(float(np.nan_to_num(a.astype("f8")).sum()), 3) if a.dtype.kind in "fiub" \
                     else [str(x) for x in a.tolist()]            # labels/categorical: compare the strings
    return out
ref = vals(lstar.read(seed))
# v2 -> v3 -> v2
lstar.write(lstar.read(seed), tmp + "/p_v3", format="v3")
lstar.write(lstar.read(tmp + "/p_v3"), tmp + "/p_v2b", format="v2")
got = vals(lstar.read(tmp + "/p_v2b"))
assert got == ref, ("v2->v3->v2 mismatch", {k: (got[k], ref[k]) for k in ref if got.get(k) != ref[k]})
# v3 -> v2 -> v3
lstar.write(lstar.read(seed3), tmp + "/p_v2", format="v2")
lstar.write(lstar.read(tmp + "/p_v2"), tmp + "/p_v3b", format="v3")
got = vals(lstar.read(tmp + "/p_v3b"))
assert got == ref, ("v3->v2->v3 mismatch", {k: (got[k], ref[k]) for k in ref if got.get(k) != ref[k]})
print("    [py] v2->v3->v2 and v3->v2->v3 == original")
PY
echo "format round-trip (both directions) PASSED."
