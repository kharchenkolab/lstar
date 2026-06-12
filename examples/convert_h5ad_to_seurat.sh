#!/usr/bin/env bash
#
# Convert an AnnData (.h5ad) into a Seurat object (.rds), using lstar as the glue.
# ---------------------------------------------------------------------------------------------
# AnnData lives in Python and Seurat in R, so the conversion crosses a language boundary. lstar
# bridges them with an on-disk L* store: Python writes it, R reads it. There is no shared memory
# and no re-implementation of either format in the other language -- each side just uses its own
# profile (read_anndata / write_seurat) against the common L* model.
#
#   AnnData (Python)  --read_anndata-->  L* store on disk  --write_seurat-->  Seurat (R, .rds)
#
# Usage:  bash convert_h5ad_to_seurat.sh [input.h5ad] [output.rds]
#   With no input, a small synthetic AnnData is generated so the script runs anywhere.
set -e
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
RLIB="$ROOT/.Rlib"                       # the lstar R package, installed locally (see conformance/run.sh)
H5AD="${1:-}"                            # optional input .h5ad
OUT_RDS="${2:-/tmp/converted_seurat.rds}"
STORE=/tmp/convert_h5ad.lstar.zarr       # the L* store that bridges Python and R

# ---- Step 1 (Python): AnnData -> L* -------------------------------------------------------------
# read_anndata applies the `anndata` profile: it walks the AnnData's slots and turns each into an
# L* field, naming them in the shared vocabulary (X -> a measure; obsm['X_pca'] -> a `pca`
# embedding; varm['PCs'] -> `pca_loadings`; obsp -> relations; obs/var columns -> per-cell/gene
# fields). Anything with no L* representation (e.g. miscellaneous `uns` entries) is recorded in
# `ds.dropped`, so the loss is visible rather than silent. We then write the store to disk.
PYTHONPATH="$ROOT/python/src:$ROOT/python/tests" python3 - "$H5AD" "$STORE" <<'PY'
import sys, warnings; warnings.filterwarnings("ignore")
import lstar
from lstar.profiles.anndata import read_anndata

h5ad, store = sys.argv[1], sys.argv[2]
if h5ad:
    import anndata as ad
    adata = ad.read_h5ad(h5ad)
else:
    # No input given: build a small AnnData with the test helper, so the demo is self-contained.
    from test_anndata_profile import make_adata
    adata = make_adata()

ds = read_anndata(adata)                                  # AnnData object -> L* dataset (axes + fields)
lstar.write(ds, store)                                    # L* dataset -> on-disk store (the bridge)
print(f"  [py] AnnData {adata.shape} -> L*: {len(ds.fields)} fields "
      f"({', '.join(list(ds.fields)[:6])}{', ...' if len(ds.fields) > 6 else ''})")
print(f"  [py] recorded as not representable in L* (dropped): {ds.dropped or 'nothing'}")
print(f"  [py] wrote bridge store -> {store}")
PY

# ---- Step 2 (R): L* -> Seurat ------------------------------------------------------------------
# write_seurat applies the `seurat` profile in reverse: L* measures over (cells, genes) become
# Seurat assay layers (transposed to Seurat's genes x cells orientation); `pca`/`umap` embeddings
# become DimReduc objects, INCLUDING their gene loadings (which a direct h5ad->Seurat conversion
# usually drops -- here they survive because L* keeps scores and loadings on one shared axis);
# per-cell fields become meta.data columns. We save the result as a standard .rds.
Rscript -e '
.libPaths(c("'"$RLIB"'", .libPaths()))
suppressMessages({ library(lstar); library(SeuratObject) })
so <- write_seurat(lstar_read("'"$STORE"'"))              # L* dataset -> Seurat object
saveRDS(so, "'"$OUT_RDS"'")
# Report what came across, so you can see the conversion actually preserved the analysis:
red <- Reductions(so)
load_dim <- if ("pca" %in% red) paste(dim(Loadings(so[["pca"]])), collapse = "x") else "n/a"
cat(sprintf("  [R ] L* -> Seurat: %s (genes x cells); layers={%s}; reductions={%s}; pca loadings=%s\n",
            paste(dim(so[["RNA"]]), collapse = " x "),
            paste(Layers(so, assay = "RNA"), collapse = ","),
            paste(red, collapse = ","), load_dim))
cat(sprintf("  [R ] saved Seurat object -> %s\n", "'"$OUT_RDS"'"))
' 2>&1 | grep -vE "^Warning|Attaching|masked|^The following|loaded|package|^$"

echo "done: $( [ -n "$H5AD" ] && echo "$H5AD" || echo "(synthetic AnnData)" ) -> $OUT_RDS"
