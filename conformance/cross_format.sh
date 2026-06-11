#!/usr/bin/env bash
# Cross-format conformance: chain a dataset through every profile and confirm the
# shared-vocabulary core (counts, X, pca, umap, leiden) is preserved end to end:
#
#   AnnData (py) -> L* -> Seurat (R) -> L* -> SCE (R) -> L* -> Python (verify)
set -e
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
RLIB="$ROOT/.Rlib"
S0=/tmp/cf0.lstar.zarr
S1=/tmp/cf1.lstar.zarr

PYTHONPATH="$ROOT/python/src:$ROOT/python/tests" python3 - "$S0" <<'PY'
import sys, warnings; warnings.filterwarnings("ignore")
from test_anndata_profile import make_adata
from lstar import read_anndata, write
write(read_anndata(make_adata()), sys.argv[1]); print("  [py] AnnData -> L*")
PY

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
for k in ("counts", "X", "pca", "umap"):
    assert np.allclose(dense(a, k), dense(b, k)), "mismatch in " + k
assert (np.asarray(a.field("leiden").values) == np.asarray(b.field("leiden").values)).all()
print("  [py] shared-vocab core preserved: counts, X, pca, umap, leiden  OK")
PY
echo "cross-format conformance PASSED."
