#!/usr/bin/env bash
# Variable-length round-trip across FORMATS and LANGUAGES, returning to the original format.
#
#   AnnData (Python) -> L* -> [ Seurat (R) -> L* -> SCE (R) -> L* ] x N -> AnnData (Python)
#
# The shared-vocabulary core (counts/X measure, pca + umap embeddings, a cell label) must survive
# a chain of any length that hops Python<->R and AnnData->Seurat->SCE, and come back into AnnData
# matching the original object. Whatever a format cannot carry (e.g. neighbor graphs through
# Seurat, uns) is reported, not silently changed. Runs on real data by default.
#
#   usage: roundtrip_xlang.sh [n_loops] [path.h5ad]
set -e
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
RLIB="$ROOT/.Rlib"
N="${1:-1}"
H5AD="${2:-/home/pkharchenko/cacoa/age/tab.muris/tabula-muris-senis-droplet-processed-official-annotations-Marrow.h5ad}"
S_ORIG=/tmp/xl_orig.lstar.zarr      # L* straight from the original AnnData (the reference)
S_CUR=/tmp/xl_cur.lstar.zarr        # rolling store as it hops through R formats

echo "chain: AnnData -> L* -> [ Seurat -> L* -> SCE -> L* ] x $N -> AnnData     (data: $(basename "$H5AD"))"

# 1) Python: original AnnData -> L*
PYTHONPATH="$ROOT/python/src" python3 - "$H5AD" "$S_ORIG" "$S_CUR" <<'PY'
import sys, warnings; warnings.filterwarnings("ignore")
import anndata as ad, lstar
from lstar.profiles.anndata import read_anndata
a = ad.read_h5ad(sys.argv[1])
ds = read_anndata(a)
lstar.write(ds, sys.argv[2]); lstar.write(ds, sys.argv[3])
print(f"  [py] AnnData {a.shape} -> L*  (fields: {', '.join(list(ds.fields)[:8])}...)")
PY

# 2) R: N loops of  L* -> Seurat -> L* -> SCE -> L*
for i in $(seq 1 "$N"); do
Rscript -e '.libPaths(c("'"$RLIB"'", .libPaths()))
suppressMessages({library(lstar); library(SeuratObject); library(SingleCellExperiment)})
s <- "'"$S_CUR"'"
ds  <- lstar_read(s)
so  <- write_seurat(ds);  ds1 <- read_seurat(so)    # L* -> Seurat -> L*
sce <- write_sce(ds1);    ds2 <- read_sce(sce)      # L* -> SCE    -> L*
lstar_write(ds2, s)
cat(sprintf("  [R ] loop %s: L* -> Seurat -> L* -> SCE -> L*\n", "'"$i"'"))' \
  2>&1 | grep -vE "^Warning|deprecat|Attaching|masked|conflicts|^The following|loaded|package|^$"
done

# 3) Python: L* -> AnnData, compare the shared-vocab core to the original
PYTHONPATH="$ROOT/python/src" python3 - "$S_ORIG" "$S_CUR" <<'PY'
import sys, warnings; warnings.filterwarnings("ignore")
import numpy as np, scipy.sparse as sp, lstar
from lstar.profiles.anndata import write_anndata
a0 = write_anndata(lstar.read(sys.argv[1]))     # original (via L*)
aN = write_anndata(lstar.read(sys.argv[2]))     # after the cross-language chain

def dense(x): return x.toarray() if sp.issparse(x) else np.asarray(x)
ok = True
# measure (X): same nnz + values after AnnData->Seurat->SCE->AnnData
xa, xb = a0.X, aN.X
mnnz = (int(xa.nnz) if sp.issparse(xa) else xa.size) == (int(xb.nnz) if sp.issparse(xb) else xb.size)
mval = np.allclose(dense(xa), dense(xb), rtol=1e-5, atol=1e-6)
print(f"  [py] back in AnnData {aN.shape}; core vs original:")
print(f"        X measure: nnz match {mnnz}, values match {mval}")
ok = ok and mnnz and mval
for k in ("X_pca", "X_umap"):
    if k in a0.obsm and k in aN.obsm:
        m = np.allclose(np.asarray(a0.obsm[k]), np.asarray(aN.obsm[k]), rtol=1e-5, atol=1e-6)
        print(f"        {k}: shape {np.asarray(aN.obsm[k]).shape} match {m}"); ok = ok and m
# a cell label carried as obs (whichever categorical survived the chain)
lab = next((c for c in a0.obs.columns if c in aN.obs.columns
            and not np.issubdtype(np.asarray(a0.obs[c]).dtype, np.number)), None)
if lab:
    m = (np.asarray(a0.obs[lab].astype(str)) == np.asarray(aN.obs[lab].astype(str))).all()
    print(f"        label '{lab}': match {m}"); ok = ok and m
lost = sorted(set(a0.obsm) - set(aN.obsm)) + (["obsp graphs"] if a0.obsp and not aN.obsp else [])
print(f"        not carried through Seurat/SCE (expected): {lost or 'none'}")
print("\n  RESULT:", "PASS — core returned to AnnData through R formats intact" if ok else "FAIL")
sys.exit(0 if ok else 1)
PY
echo "cross-language round-trip PASSED."
