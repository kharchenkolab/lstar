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
else
  echo "  [skip] R / .Rlib / SeuratObject unavailable — skipping the cross-language (Seurat) leg"
fi
echo "convert-cli (L0-L3) PASSED."
