#!/usr/bin/env bash
# Cross-format multimodal consistency: a CITE-seq dataset (RNA + protein) converted from a Seurat
# object and from a MuData object must land on the SAME canonical feature axes (genes, proteins) -- the
# shared modality vocabulary (read_seurat .modality_axis == read_mudata _MODALITY_AXIS). So a modality
# is the same lstar feature space regardless of source format. Skips cleanly when a toolchain is absent.
set -u
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
PY="${LSTAR_PY:-python3}"
export PYTHONPATH="$ROOT/python/src${PYTHONPATH:+:$PYTHONPATH}"
export LSTAR_RLIB="${LSTAR_RLIB:-$ROOT/.Rlib}"
TMP="$(mktemp -d)"; trap 'rm -rf "$TMP"' EXIT

"$PY" -c "import lstar, mudata, anndata" 2>/dev/null || { echo "  [skip] python lstar/mudata/anndata missing"; exit 0; }
command -v Rscript >/dev/null 2>&1 || { echo "  [skip] no Rscript"; exit 0; }
Rscript -e '.libPaths(c(Sys.getenv("LSTAR_RLIB"),.libPaths())); if(!requireNamespace("SeuratObject",quietly=TRUE)) quit(status=1)' 2>/dev/null \
  || { echo "  [skip] R SeuratObject missing"; exit 0; }

# (1) Seurat CITE-seq (RNA + ADT) -> lstar store
Rscript - "$TMP/se.lstar.zarr" >/tmp/lstar_xfmt_r.log 2>&1 <<'EOF'
.libPaths(c(Sys.getenv("LSTAR_RLIB"), .libPaths())); suppressMessages({library(lstar); library(SeuratObject); library(Matrix)})
out <- commandArgs(trailingOnly = TRUE)[1]
rna <- matrix(rpois(20*30, 2), 20, 30, dimnames = list(paste0("g", 1:20), paste0("c", 1:30)))
adt <- matrix(rpois(6*30, 5), 6, 30, dimnames = list(paste0("p", 1:6), paste0("c", 1:30)))
so <- CreateSeuratObject(as(rna, "CsparseMatrix"), assay = "RNA"); so[["ADT"]] <- CreateAssayObject(counts = as(adt, "CsparseMatrix"))
if (dir.exists(out)) unlink(out, recursive = TRUE); lstar_write(read_seurat(so), out)
EOF
[ -d "$TMP/se.lstar.zarr" ] || { echo "  FAIL: Seurat->store"; tail -8 /tmp/lstar_xfmt_r.log; exit 1; }

# (2) MuData CITE-seq (rna + adt) -> lstar store, then compare feature axes
"$PY" - "$TMP/se.lstar.zarr" "$TMP/mu.lstar.zarr" <<'EOF'
import sys, warnings; warnings.filterwarnings("ignore")
import numpy as np, anndata as ad, mudata as mu, scipy.sparse as sp, lstar
se, muo = sys.argv[1], sys.argv[2]
rna = ad.AnnData(sp.csr_matrix(np.random.default_rng(2).poisson(2, (30, 20)).astype("f4")))
rna.obs_names = [f"c{i+1}" for i in range(30)]; rna.var_names = [f"g{j+1}" for j in range(20)]
adt = ad.AnnData(sp.csr_matrix(np.random.default_rng(3).poisson(5, (30, 6)).astype("f4")))
adt.obs_names = rna.obs_names; adt.var_names = [f"p{j+1}" for j in range(6)]
lstar.write(lstar.read_mudata(mu.MuData({"rna": rna, "adt": adt})), muo)
feat = lambda p: sorted(n for n, a in lstar.read(p).axes.items() if a.role == "feature")
a, b = feat(se), feat(muo)
print(f"  Seurat feature axes: {a}   MuData feature axes: {b}")
if a == b:
    print("  OK: CITE-seq lands on the same feature axes regardless of source")
else:
    print("  FAIL: cross-format feature-axis mismatch"); sys.exit(1)
EOF
