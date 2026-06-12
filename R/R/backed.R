# Disk-backed conversion targets.
#
# The streamed L* -> .h5ad writer (Python `lstar.convert_to_h5ad`) produces an .h5ad whose big
# measures live on disk. These R readers open that .h5ad as a *disk-backed* object -- a Seurat v5
# assay backed by BPCells, or a SingleCellExperiment assay backed by HDF5Array -- so the expression
# matrix is never loaded into memory. Together with the Python `read_h5ad(path, backed="r")` path,
# this closes the bounded-memory conversion loop end to end: a multi-gigabyte atlas converts from an
# L* store into a usable native object on a laptop.
#
# BPCells and HDF5Array are optional (Suggests); each reader errors with a clear message if its
# package is absent, and is skipped by tests when unavailable.

.need_pkg <- function(pkg, fn) {
  if (!requireNamespace(pkg, quietly = TRUE))
    stop(sprintf("%s() requires the '%s' package (install it to use this disk-backed target).",
                 fn, pkg), call. = FALSE)
}

#' Open an .h5ad's expression matrix as a disk-backed Seurat v5 assay (via BPCells).
#'
#' Reads the matrix from `h5ad` with BPCells (`open_matrix_anndata_hdf5`) so it stays on disk as a
#' streaming `IterableMatrix` -- already oriented genes x cells (rownames genes, colnames cells), the
#' Seurat convention -- and wraps it in a Seurat v5 `Assay5`. The matrix is never materialized: peak
#' memory is a few megabytes regardless of atlas size. Typical use is the bounded end of an L*
#' conversion, after `lstar.convert_to_h5ad(store, h5ad)` in Python:
#'
#' \preformatted{
#'   so <- read_seurat_backed("atlas.h5ad")          # counts live on disk (BPCells)
#'   so <- Seurat::NormalizeData(so)                  # Seurat v5 ops stream off disk
#' }
#'
#' @param h5ad path to an .h5ad file (its `X`/layer stays on disk).
#' @param group h5ad group to read as the matrix (default `"X"`; e.g. `"layers/counts"`).
#' @param assay name for the Seurat assay (default `"RNA"`).
#' @param project Seurat project name.
#' @return a Seurat object whose assay counts are a disk-backed BPCells matrix.
#' @export
read_seurat_backed <- function(h5ad, group = "X", assay = "RNA", project = "lstar") {
  .need_pkg("BPCells", "read_seurat_backed")
  .need_pkg("SeuratObject", "read_seurat_backed")
  m <- BPCells::open_matrix_anndata_hdf5(h5ad, group = group)   # genes x cells, on disk (Seurat orientation)
  assay5 <- SeuratObject::CreateAssay5Object(counts = m)
  SeuratObject::CreateSeuratObject(assay5, assay = assay, project = project)
}

#' Open an .h5ad's expression matrix as a disk-backed SingleCellExperiment assay (via HDF5Array).
#'
#' Reads the matrix from `h5ad` as an `HDF5Array::H5ADMatrix` -- a `DelayedMatrix` that stays on disk
#' (genes x cells, the Bioconductor convention) -- and wraps it in a `SingleCellExperiment`. The
#' matrix is never materialized. Pairs with Python `lstar.convert_to_h5ad(store, h5ad)`.
#'
#' @param h5ad path to an .h5ad file.
#' @param layer h5ad layer to read; `NULL` (default) reads `X`.
#' @param assay_name name for the SCE assay (default inferred: `"counts"`).
#' @return a SingleCellExperiment whose assay is a disk-backed DelayedMatrix.
#' @export
read_sce_backed <- function(h5ad, layer = NULL, assay_name = "counts") {
  .need_pkg("HDF5Array", "read_sce_backed")
  .need_pkg("SingleCellExperiment", "read_sce_backed")
  X <- HDF5Array::H5ADMatrix(h5ad, layer = layer)               # genes x cells DelayedMatrix, on disk
  assays <- stats::setNames(list(X), assay_name)
  SingleCellExperiment::SingleCellExperiment(assays = assays)
}
