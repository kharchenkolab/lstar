#!/usr/bin/env bash
# Cross-format conformance, grounded in REAL data (pbmc68k_reduced): chain a real dataset through every
# profile and confirm the shared-vocabulary core survives end to end:
#
#   AnnData (py) -> L* -> Seurat (R) -> L* -> SCE (R) -> L* -> Python (verify)
set -e
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
RLIB="$ROOT/.Rlib"
S0=/tmp/cf0.lstar.zarr
S1=/tmp/cf1.lstar.zarr

if ! PYTHONPATH="$ROOT/python/src:$ROOT/python/tests" python3 - "$S0" <<'PY'
import sys, warnings; warnings.filterwarnings("ignore")
import corpus
from lstar import read_anndata, write
a = corpus.pbmc68k_reduced()
if a is None:
    print("  SKIP cross_format (corpus unavailable)"); sys.exit(7)
write(read_anndata(a), sys.argv[1]); print("  [py] real AnnData (pbmc68k) -> L*")
PY
then [ $? -eq 7 ] && exit 0 || exit 1; fi

Rscript -e '.libPaths(c("'"$RLIB"'", .libPaths()))
suppressMessages({library(lstar); library(SeuratObject); library(SingleCellExperiment)})
ds  <- lstar_read("'"$S0"'")
so  <- write_seurat(ds)      # L* -> Seurat
ds2 <- read_seurat(so)       # Seurat -> L*
sce <- write_sce(ds2)        # L* -> SCE
ds3 <- read_sce(sce)         # SCE -> L*
lstar_write(ds3, "'"$S1"'")
cat("  [R ] L* -> Seurat -> L* -> SCE -> L*\n")' 2>&1 | grep -vE "^Warning|deprecat|Attaching|masked|conflicts|^The following|loaded|package|^$"

PYTHONPATH="$ROOT/python/src" python3 - "$S0" "$S1" <<'PY'
import sys, warnings; warnings.filterwarnings("ignore")
import numpy as np, scipy.sparse as sp
from lstar.zarr_io import read
a, b = read(sys.argv[1]), read(sys.argv[2])
def dense(ds, k):
    v = ds.field(k).values
    return v.toarray() if sp.issparse(v) else np.asarray(v)
# the cross-format-robust shared vocabulary survives end to end. The count measure is named by each
# format's convention -- anndata `.raw` -> L* 'raw'; Seurat counts -> L* 'counts' -- so compare across.
assert np.allclose(dense(a, "raw"), dense(b, "counts"), rtol=1e-4), "counts mismatch through the chain"
for k in ("pca", "umap"):                                          # embeddings preserved
    assert np.allclose(dense(a, k), dense(b, k), rtol=1e-4, equal_nan=True), "mismatch in " + k
assert (np.asarray(a.field("louvain").values) == np.asarray(b.field("louvain").values)).all()
print("  [py] shared-vocab core preserved through Seurat<->SCE: counts, pca, umap, louvain  OK")
PY
echo "cross-format conformance PASSED."
