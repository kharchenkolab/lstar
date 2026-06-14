#!/usr/bin/env bash
# CLI conversion conformance. `lstar convert` detects formats from paths, routes the data through the L*
# store, and the produced object round-trips. L0: in-process Python h5ad <-> store <-> h5ad (detection +
# routing). Self-contained — a tiny synthetic AnnData through the real anndata pipeline; SKIPs cleanly if
# anndata is absent.
set -eo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
export PYTHONPATH="$ROOT/python/src"
python3 -c 'import anndata' 2>/dev/null || { echo "  [skip] anndata not installed — skipping convert-cli conformance"; exit 0; }

TMP="$(mktemp -d)"; trap 'rm -rf "$TMP"' EXIT
H5AD="$TMP/in.h5ad"; STORE="$TMP/mid.lstar.zarr"; OUT="$TMP/out.h5ad"

python3 - "$H5AD" <<'PY'
import sys, warnings; warnings.filterwarnings("ignore")
import numpy as np, scipy.sparse as sp, pandas as pd, anndata as ad
np.random.seed(0)
X = sp.random(20, 8, density=0.4, format="csr").astype("float32")
obs = pd.DataFrame({"leiden": pd.Categorical(np.random.choice(list("AB"), 20))},
                   index=[f"c{i}" for i in range(20)])
var = pd.DataFrame(index=[f"g{i}" for i in range(8)])
a = ad.AnnData(X=X, obs=obs, var=var)
a.obsm["X_pca"] = np.random.randn(20, 3).astype("float32")
a.write_h5ad(sys.argv[1])
print(f"  [py] wrote synthetic AnnData {a.shape} -> in.h5ad")
PY

# detection by extension: .h5ad -> anndata, .lstar.zarr -> store
python3 -m lstar convert "$H5AD" "$STORE" | sed 's/^/  /'
python3 -m lstar convert "$STORE" "$OUT"   | sed 's/^/  /'

python3 - "$H5AD" "$OUT" <<'PY'
import sys, warnings; warnings.filterwarnings("ignore")
import numpy as np, scipy.sparse as sp, anndata as ad
a = ad.read_h5ad(sys.argv[1]); b = ad.read_h5ad(sys.argv[2])
assert a.shape == b.shape, (a.shape, b.shape)
assert list(a.obs_names) == list(b.obs_names) and list(a.var_names) == list(b.var_names)
den = lambda M: M.toarray() if sp.issparse(M) else np.asarray(M)
assert np.allclose(den(a.X), den(b.X)), "X values changed across h5ad -> store -> h5ad"
assert "X_pca" in b.obsm and b.obsm["X_pca"].shape == (20, 3), "pca embedding lost"
print(f"  [py] h5ad -> store -> h5ad: shape {b.shape}, X intact, pca preserved")
PY

# L1: the fidelity report (full text + machine-readable JSON)
RJSON="$TMP/report.json"
python3 -m lstar convert "$H5AD" "$STORE" --report --report-json "$RJSON" >/dev/null
python3 -m lstar inspect "$STORE" >/dev/null                      # inspect renders the report, no write
python3 - "$RJSON" <<'PY'
import sys, json
r = json.load(open(sys.argv[1]))
assert r["source"]["format"] == "anndata" and r["target"]["format"] == "store", r["source"]
assert isinstance(r["dropped"], list)
names = {f["name"] for f in r["fields"]}
assert {"X", "leiden", "pca"} <= names, names
assert any(f["name"] == "X" and f["role"] == "measure" for f in r["fields"])
assert any(a["name"] == "leiden" and a["role"] == "factor" for a in r["axes"])      # induced factor axis
assert all({"role", "span", "encoding", "coverage", "nullable"} <= set(f) for f in r["fields"])
print(f"  [py] report JSON: {len(r['axes'])} axes, {len(r['fields'])} fields, "
      f"dropped={len(r['dropped'])}; roles + structure correct")
PY

# L3: native-acceptance — open the produced object in its native library + a canonical-ops smoke. For an
# anndata target this runs scanpy (normalize/log1p/pca/rank_genes_groups) when scanpy is present, else
# falls back to open + structural invariants. A valid conversion must never FAIL the check (--strict).
RJC="$TMP/report_check.json"
python3 -m lstar convert "$H5AD" "$OUT" --strict --report-json "$RJC" >/dev/null
python3 - "$RJC" <<'PY'
import sys, json
nc = json.load(open(sys.argv[1]))["native_check"]
assert nc["format"] == "anndata" and nc["status"] in ("pass", "skip"), nc
print(f"  [py] native-acceptance (anndata): {nc['status']} — {nc['detail'][:64]}")
PY

# Tier-A / P1: the package-free h5py reader (--backend direct) must produce a store value-equal to the
# native (anndata) read on the core. Forcing --backend direct exercises the fallback even though anndata
# IS present, so CI covers the path that would otherwise run only on package-absent machines.
if python3 -c 'import h5py' 2>/dev/null; then
  SN="$TMP/ta_native.lstar.zarr"; SD="$TMP/ta_direct.lstar.zarr"
  python3 -m lstar convert "$H5AD" "$SN" --backend native --no-check -q
  python3 -m lstar convert "$H5AD" "$SD" --backend direct --no-check -q
  python3 - "$SN" "$SD" <<'PY'
import sys, warnings; warnings.filterwarnings("ignore")
import numpy as np, scipy.sparse as sp, lstar
from lstar import Categorical
n = lstar.read(sys.argv[1]); d = lstar.read(sys.argv[2])
dense = lambda f: (f.values.toarray() if sp.issparse(f.values) else np.asarray(f.values))
assert list(np.asarray(n.axis("cells").labels)) == list(np.asarray(d.axis("cells").labels))
assert list(np.asarray(n.axis("genes").labels)) == list(np.asarray(d.axis("genes").labels))
assert np.allclose(dense(n.field("X")), dense(d.field("X"))), "X differs (direct vs native)"
for nm in d.fields:                                  # every core field present in both must match
    if nm not in n.fields:
        continue
    fn, fd = n.field(nm), d.field(nm)
    if isinstance(fd.values, Categorical):
        assert list(fd.values.codes) == list(fn.values.codes), nm
    else:
        assert np.allclose(dense(fd), dense(fn), equal_nan=True), nm
assert d.axis("leiden").role == "factor", "induction missing on the direct path"
assert not [e for e in lstar.validate(d) if e.startswith("ERROR")]
print("  [py] Tier-A direct READ: package-free h5py read == native read on the core (X, obs/var, obsm, induction)")
PY
  # P2: the package-free h5py WRITER (--backend direct) must emit an h5ad that native anndata reads and
  # ACCEPTS (--strict native-acceptance), value-equal to the native writer on the core.
  HWD="$TMP/ta_w_direct.h5ad"; HWN="$TMP/ta_w_native.h5ad"
  python3 -m lstar convert "$STORE" "$HWD" --backend direct --strict >/dev/null   # write + native-acceptance
  python3 -m lstar convert "$STORE" "$HWN" --backend native --no-check -q
  python3 - "$HWD" "$HWN" <<'PY'
import sys, warnings; warnings.filterwarnings("ignore")
import numpy as np, scipy.sparse as sp, anndata as ad
d = ad.read_h5ad(sys.argv[1]); n = ad.read_h5ad(sys.argv[2])
den = lambda M: (M.toarray() if sp.issparse(M) else np.asarray(M))
assert d.shape == n.shape, (d.shape, n.shape)
assert list(d.obs_names) == list(n.obs_names) and list(d.var_names) == list(n.var_names)
Xof = lambda a: (den(a.X) if a.X is not None else den(next(iter(a.layers.values()))))
assert np.allclose(Xof(d), Xof(n)), "X differs (direct vs native write)"
assert np.allclose(np.asarray(d.obsm["X_pca"]), np.asarray(n.obsm["X_pca"])), "pca differs"
assert str(d.obs["leiden"].dtype) == "category" and list(d.obs["leiden"]) == list(n.obs["leiden"])
print("  [py] Tier-A direct WRITE: package-free h5py write read + accepted by native anndata, == native write")
PY
else
  echo "  [skip] h5py absent — skipping the package-free (direct) AnnData read/write checks"
fi

# L2: cross-language routing (h5ad <-> Seurat .rds, bridged by the store + an Rscript driver). Needs R
# with the lstar package (.Rlib) + SeuratObject; SKIPs cleanly otherwise.
if command -v Rscript >/dev/null 2>&1 && [ -d "$ROOT/.Rlib/lstar" ] && \
   Rscript -e '.libPaths(c("'"$ROOT"'/.Rlib", .libPaths())); quit(status=!requireNamespace("SeuratObject", quietly=TRUE))' </dev/null >/dev/null 2>&1; then
  export LSTAR_RLIB="$ROOT/.Rlib"
  RDS="$TMP/mid.rds"; BACK="$TMP/back.h5ad"; RJR="$TMP/rds_check.json"
  python3 -m lstar convert "$H5AD" "$RDS" --strict --report-json "$RJR" | sed 's/^/  /'   # h5ad -> Seurat .rds
  python3 - "$RJR" <<'PY'
import sys, json
nc = json.load(open(sys.argv[1]))["native_check"]
assert nc["format"] == "seurat" and nc["status"] == "pass", nc      # must OPEN + pass canonical Seurat ops
print(f"  [py] native-acceptance (seurat): pass — {nc['detail'][:64]}")
PY
  python3 -m lstar convert "$RDS" "$BACK"  | sed 's/^/  /'       # Seurat .rds -> h5ad (via R)
  python3 - "$H5AD" "$BACK" <<'PY'
import sys, warnings; warnings.filterwarnings("ignore")
import numpy as np, scipy.sparse as sp, anndata as ad
a = ad.read_h5ad(sys.argv[1]); b = ad.read_h5ad(sys.argv[2])
assert a.shape == b.shape, (a.shape, b.shape)
assert set(a.obs_names) == set(b.obs_names) and set(a.var_names) == set(b.var_names)
den = lambda M: M.toarray() if sp.issparse(M) else np.asarray(M)
bi = [list(b.obs_names).index(n) for n in a.obs_names]          # align order (Seurat may reorder)
bj = [list(b.var_names).index(n) for n in a.var_names]
ax = den(a.X)
# the Seurat path reshapes the raw measure into an Assay5 `counts` layer, so the matrix comes back in
# X *or* layers["counts"] -- accept either; we're checking the values survived, not their slot.
cands = ([den(b.X)] if b.X is not None else []) + [den(b.layers[L]) for L in b.layers]
ok = any(c.ndim == 2 and c.shape == ax.shape and np.allclose(ax, c[np.ix_(bi, bj)]) for c in cands)
assert ok, f"expression values not recovered (checked X + layers {list(b.layers)})"
print(f"  [py] h5ad -> Seurat(.rds) -> h5ad: shape {b.shape}, expression intact (cross-language round-trip)")
PY

  # Tier-A / P3: the package-free Seurat reader (--backend direct, base R + Matrix, NO SeuratObject
  # accessors — it walks S4 slots via attr()) must produce a store value-equal to the native read of the
  # SAME .rds. Forcing --backend direct exercises the slot-walk even with SeuratObject present.
  SRD="$TMP/ta_seu_direct.lstar.zarr"; SRN="$TMP/ta_seu_native.lstar.zarr"
  python3 -m lstar convert "$RDS" "$SRD" --backend direct --no-check -q
  python3 -m lstar convert "$RDS" "$SRN" --backend native --no-check -q
  python3 - "$SRD" "$SRN" <<'PY'
import sys, warnings; warnings.filterwarnings("ignore")
import numpy as np, scipy.sparse as sp, lstar
from lstar import Categorical
d = lstar.read(sys.argv[1]); n = lstar.read(sys.argv[2])
dense = lambda f: (f.values.toarray() if sp.issparse(f.values) else np.asarray(f.values))
assert list(np.asarray(d.axis("cells").labels)) == list(np.asarray(n.axis("cells").labels))
assert list(np.asarray(d.axis("genes").labels)) == list(np.asarray(n.axis("genes").labels))
assert set(d.fields) == set(n.fields), (sorted(d.fields), sorted(n.fields))
for nm in d.fields:
    fd, fn = d.field(nm), n.field(nm)
    if isinstance(fd.values, Categorical):
        assert list(fd.values.codes) == list(fn.values.codes), nm
    else:
        assert np.allclose(dense(fd), dense(fn), equal_nan=True), nm
assert not [e for e in lstar.validate(d) if e.startswith("ERROR")]
print(f"  [py] Tier-A direct READ (Seurat): base-R slot-walk == native read ({len(d.fields)} fields + axes)")
PY

  # Tier-A / P4: the package-free Seurat WRITER (--backend direct, base R + a PINNED Assay5 setClass schema
  # with the S4 class identity forged to SeuratObject, NO SeuratObject installed-or-used) must emit a .rds
  # that native SeuratObject reads + ACCEPTS (--strict: real NormalizeData/ScaleData/RunPCA), re-reading
  # value-equal to the native writer's output.
  RWD="$TMP/ta_w_direct.rds"; RWN="$TMP/ta_w_native.rds"
  SWD="$TMP/ta_w_d.lstar.zarr"; SWN="$TMP/ta_w_n.lstar.zarr"
  python3 -m lstar convert "$STORE" "$RWD" --backend direct --strict >/dev/null    # write + native-acceptance
  python3 -m lstar convert "$STORE" "$RWN" --backend native --no-check -q
  python3 -m lstar convert "$RWD" "$SWD" --backend native --no-check -q            # native re-reads both
  python3 -m lstar convert "$RWN" "$SWN" --backend native --no-check -q
  python3 - "$SWD" "$SWN" <<'PY'
import sys, warnings; warnings.filterwarnings("ignore")
import numpy as np, scipy.sparse as sp, lstar
from lstar import Categorical
d = lstar.read(sys.argv[1]); n = lstar.read(sys.argv[2])
dense = lambda f: (f.values.toarray() if sp.issparse(f.values) else np.asarray(f.values))
assert set(d.fields) == set(n.fields), (sorted(d.fields), sorted(n.fields))
for nm in d.fields:
    fd, fn = d.field(nm), n.field(nm)
    if isinstance(fd.values, Categorical):
        assert list(fd.values.codes) == list(fn.values.codes), nm
    else:
        assert np.allclose(dense(fd), dense(fn), equal_nan=True), nm
print(f"  [py] Tier-A direct WRITE (Seurat): pinned-schema build accepted by native, == native write ({len(d.fields)} fields)")
PY
else
  echo "  [skip] R / .Rlib / SeuratObject unavailable — skipping the cross-language (Seurat) leg"
fi
echo "convert-cli (L0-L3) PASSED."
