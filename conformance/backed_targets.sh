#!/usr/bin/env bash
# Disk-backed conversion targets: an L* store streamed to .h5ad must open as a *disk-backed* native
# object -- a Seurat v5 assay backed by BPCells, a SingleCellExperiment assay backed by HDF5Array,
# and a backed AnnData -- with the expression matrix staying on disk (genes x cells, values intact).
# This closes the bounded-memory loop end to end. The R targets need BPCells / HDF5Array; each is
# skipped (not failed) when its package is absent. Self-contained (synthetic data).
set -e
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
RLIB="$ROOT/.Rlib"
STORE=/tmp/backed_conf.lstar.zarr
H5=/tmp/backed_conf.h5ad

# (1) Python: build a small store, stream it out to .h5ad, and confirm the backed AnnData target.
PYTHONPATH="$ROOT/python/src" python3 - "$STORE" "$H5" <<'PY'
import sys, warnings; warnings.filterwarnings("ignore")
import numpy as np, scipy.sparse as sp, anndata as ad
import lstar
store, h5 = sys.argv[1], sys.argv[2]
ncell, ngene = 300, 120
X = sp.random(ncell, ngene, density=0.1, format="csr", random_state=0)
X.data = np.round(X.data * 9 + 1).astype("f4")
ds = lstar.Dataset(kind="sample")
ds.add_axis("cells", [f"cell{i}" for i in range(ncell)])
ds.add_axis("genes", [f"gene{i}" for i in range(ngene)])
ds.add_field("X", X, role="measure", span=["cells", "genes"], state="lognorm")
lstar.write(ds, store, chunk_elems=5000)
lstar.convert_to_h5ad(store, h5, chunk_elems=5000)          # streamed L* -> h5ad

a = ad.read_h5ad(h5, backed="r")                            # AnnData disk-backed target
# X must be a disk-backed proxy, not an in-memory matrix. (Don't assert the class *name*: anndata
# renamed it SparseDataset -> CSRDataset/CSCDataset at 0.10, which is what broke this across versions.)
assert a.isbacked and not isinstance(a.X, (np.ndarray, sp.spmatrix))
assert (sp.csr_matrix(a.X[10:20]) != X[10:20]).nnz == 0     # a slice reads from disk, matches
a.file.close()
print("  [py] streamed L* -> h5ad; AnnData backed='r' target: X on disk (backed proxy), slice OK")
PY

# (2) R: open the same .h5ad as disk-backed Seurat v5 (BPCells) and SCE (HDF5Array); skip if absent.
# Each target's matrix must be on disk, genes x cells, with correct dimnames; and the two targets
# must read the identical block off disk (cross-checking orientation and values without a separate
# ground truth).
Rscript -e '.libPaths(c("'"$RLIB"'", .libPaths()))
suppressMessages(library(lstar))
h5 <- "'"$H5"'"

ok_block <- function(M) {
  d <- dim(M)
  stopifnot(d[1] == 120, d[2] == 300)                       # genes x cells
  stopifnot(identical(rownames(M)[1:2], c("gene0","gene1")))
  stopifnot(identical(colnames(M)[1:2], c("cell0","cell1")))
  as.matrix(M[1:8, 1:6])
}

if (requireNamespace("BPCells", quietly = TRUE) && requireNamespace("SeuratObject", quietly = TRUE)) {
  so  <- read_seurat_backed(h5)
  cnt <- SeuratObject::LayerData(so, "counts")
  stopifnot(inherits(cnt, "IterableMatrix"))                # on disk
  b1 <- ok_block(cnt)
  cat("  [R ] Seurat v5 / BPCells: disk-backed, genes x cells = 120 x 300, dimnames + block OK\n")
} else {
  cat("  [R ] Seurat v5 / BPCells: SKIP (BPCells not installed)\n"); b1 <- NULL
}

if (requireNamespace("HDF5Array", quietly = TRUE) && requireNamespace("SingleCellExperiment", quietly = TRUE)) {
  sce <- read_sce_backed(h5)
  X   <- SummarizedExperiment::assay(sce, "counts")
  stopifnot(is(X, "DelayedMatrix"))                         # on disk
  b2 <- ok_block(X)
  cat("  [R ] SCE / HDF5Array: disk-backed, genes x cells = 120 x 300, dimnames + block OK\n")
} else {
  cat("  [R ] SCE / HDF5Array: SKIP (HDF5Array not installed)\n"); b2 <- NULL
}

if (!is.null(b1) && !is.null(b2))
  stopifnot(max(abs(b1 - b2)) < 1e-5)                       # the two disk-backed targets agree exactly
' 2>&1 | grep -vE "deprecat|^Warning|masked|following object|Attaching|once every|^This message"
