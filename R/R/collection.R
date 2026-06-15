# Assemble a collection of heterogeneous samples from per-sample objects (the non-Conos path). Mirrors
# the canonical structure the Conos / Seurat-v5-split profiles produce, so every collection has one shape.

#' Assemble a collection of heterogeneous samples from per-sample objects.
#'
#' A *collection* keeps each sample's own data --- per-sample `cells.<s>`/`genes.<s>` axes and
#' `<field>.<s>` fields, over gene sets that may overlap, differ, or be entirely **disjoint** across
#' samples --- alongside a `samples` axis and a *union* `cells` axis carrying the joint analysis (a shared
#' embedding, clusters, a graph). This is lstar's "a collection is not one aligned tensor" model, built
#' from any list of separately-processed samples rather than hand-assembled.
#'
#' @param samples a named list of per-sample `lstar_dataset`s (or `Seurat` / `SingleCellExperiment`
#'   objects, read via the profiles). Each must have `cells` and `genes` axes.
#' @param joint optional named list of fields over the **union** cells (the integration outputs): a matrix
#'   becomes a joint embedding; a factor/character vector a clustering (inducing a factor axis); a
#'   `(cells x cells)` sparse matrix a graph relation.
#' @param sample_field name of the design label recording each union cell's sample (default `"sample"`).
#' @param prefix_cells prefix each cell label with its sample name so union labels are unique (default `TRUE`).
#' @return an `lstar_dataset` of kind `"collection"`.
#' @examples
#' mk <- function(s, nc, genes) {
#'   ds <- list(kind = "sample", axes = list(), fields = list())
#'   ds$axes$cells <- list(labels = paste0("c", seq_len(nc)), origin = "observed", role = "observation")
#'   ds$axes$genes <- list(labels = genes, origin = "observed", role = "feature")
#'   m <- as(Matrix::Matrix(matrix(rpois(nc * length(genes), 2), nc), sparse = TRUE), "CsparseMatrix")
#'   ds$fields$counts <- list(role = "measure", span = c("cells", "genes"), state = "raw", values = m)
#'   class(ds) <- "lstar_dataset"; ds
#' }
#' col <- collection_from(list(A = mk("A", 5, paste0("g", 1:8)),
#'                             B = mk("B", 7, paste0("g", 3:12))))   # divergent gene sets
#' col$kind
#' @seealso [write_conos()], [read_seurat()]
#' @export
collection_from <- function(samples, joint = NULL, sample_field = "sample", prefix_cells = TRUE) {
  if (!length(samples)) stop("collection_from: no samples given")
  if (is.null(names(samples))) names(samples) <- paste0("s", seq_along(samples))
  as_ds <- function(x) {
    if (inherits(x, "lstar_dataset")) return(x)
    if (methods::is(x, "Seurat")) return(read_seurat(x))
    if (methods::is(x, "SingleCellExperiment")) return(read_sce(x))
    stop("collection_from: each sample must be an lstar_dataset or a Seurat / SingleCellExperiment object")
  }
  ds <- list(kind = "collection", spec_version = "0.1", profiles = "collection@0.1",
             dropped = character(0), axes = list(), fields = list())
  union_cells <- character(0); sample_of <- character(0); n_per <- integer(0)
  for (s in names(samples)) {
    dss <- as_ds(samples[[s]])
    if (is.null(dss$axes$cells) || is.null(dss$axes$genes))
      stop(sprintf("collection_from: sample '%s' has no cells/genes axes", s))
    clabs <- as.character(dss$axes$cells$labels)
    if (prefix_cells) clabs <- paste0(s, "_", clabs)
    rn <- list()
    for (ax in names(dss$axes)) {
      a <- dss$axes[[ax]]; rn[[ax]] <- paste0(ax, ".", s)
      if (identical(a$role, "factor") && !is.null(a$induced_by)) next   # re-induced by its (renamed) field
      labs <- if (identical(ax, "cells")) clabs else as.character(a$labels)
      ds$axes[[rn[[ax]]]] <- list(labels = labs, origin = a$origin %||% "observed", role = a$role)
    }
    for (fn in names(dss$fields)) {
      f <- dss$fields[[fn]]
      f$span <- vapply(as.character(f$span), function(a) rn[[a]], character(1), USE.NAMES = FALSE)
      if (!is.null(f$index_axis)) f$index_axis <- rn[[f$index_axis]]
      ds$fields[[paste0(fn, ".", s)]] <- f
    }
    union_cells <- c(union_cells, clabs); sample_of <- c(sample_of, rep(s, length(clabs)))
    n_per <- c(n_per, length(clabs))
  }
  ds$axes$samples <- list(labels = names(samples), origin = "observed", role = "sample")
  ds$fields$n_cells <- list(role = "measure", span = "samples", state = "", values = as.numeric(n_per))
  if (anyDuplicated(union_cells)) stop("collection_from: union cell labels are not unique (use prefix_cells = TRUE)")
  ds$axes$cells <- list(labels = union_cells, origin = "derived", role = "observation")
  ds$fields[[sample_field]] <- list(role = "label", span = "cells", subtype = "design", values = sample_of)

  for (nm in names(joint)) {
    v <- joint[[nm]]
    if (methods::is(v, "sparseMatrix")) {
      ds$fields[[nm]] <- list(role = "relation", span = c("cells", "cells"), subtype = "knn",
                              values = as(v, "CsparseMatrix"))
    } else if (is.matrix(v) && !is.factor(v) && ncol(v) > 1) {
      ds$axes[[nm]] <- list(labels = paste0(nm, seq_len(ncol(v))), origin = "derived", role = "coordinate")
      ds$fields[[nm]] <- list(role = "embedding", span = c("cells", nm), values = unname(v))
    } else {                                              # clustering -> categorical label + factor axis
      vf <- as.factor(v)
      ds$fields[[nm]] <- list(role = "label", span = "cells", encoding = "categorical", values = vf)
      ds$axes[[nm]] <- list(labels = levels(vf), origin = "derived", role = "factor", induced_by = nm)
    }
  }
  class(ds) <- "lstar_dataset"
  ds
}
