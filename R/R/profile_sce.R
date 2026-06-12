# SingleCellExperiment / SummarizedExperiment profile.
# L* measures over (cells, genes) <-> SCE assays (genes x cells); embeddings <-> reducedDims;
# arity-1 fields over cells/genes <-> colData/rowData; loadings <-> the reducedDim "rotation" attr.

.sce_assay_name <- function(state, nm) {
  if (identical(state, "raw") || nm == "counts") "counts"
  else if (identical(state, "lognorm") || nm %in% c("X", "data", "logcounts")) "logcounts"
  else nm
}

#' Build a SingleCellExperiment from an L* dataset.
#'
#' Measures become assays (transposed to genes x cells), embeddings become `reducedDims`, and
#' arity-1 cell fields become `colData`.
#'
#' @param ds an `lstar_dataset`
#' @return a `SingleCellExperiment` object.
#' @seealso [read_sce()]
#' @export
write_sce <- function(ds) {
  if (!requireNamespace("SingleCellExperiment", quietly = TRUE) ||
      !requireNamespace("S4Vectors", quietly = TRUE))
    stop("SingleCellExperiment and S4Vectors are required")
  cells <- as.character(ds$axes$cells$labels)
  genes <- as.character(ds$axes$genes$labels)
  gxc <- function(nm) { m <- Matrix::t(ds$fields[[nm]]$values); dimnames(m) <- list(genes, cells); m }

  assays <- list()
  for (nm in .fields_over(ds, c("cells", "genes"), role = "measure")) {
    assays[[.sce_assay_name(ds$fields[[nm]]$state, nm)]] <- gxc(nm)
  }

  cd <- list(); rd <- list()
  for (nm in names(ds$fields)) {
    f <- ds$fields[[nm]]
    if (identical(as.character(f$span), "cells") && length(f$values) == length(cells)) cd[[nm]] <- f$values
    if (identical(as.character(f$span), "genes") && length(f$values) == length(genes)) rd[[nm]] <- f$values
  }

  reduced <- list()
  for (nm in names(ds$fields)) {
    f <- ds$fields[[nm]]
    if (identical(f$role, "embedding") && length(f$span) == 2 && f$span[[1]] == "cells") {
      emb <- as.matrix(f$values); rownames(emb) <- cells
      load_nm <- paste0(f$span[[2]], "_loadings")
      if (!is.null(ds$fields[[load_nm]])) {
        L <- as.matrix(ds$fields[[load_nm]]$values); rownames(L) <- genes
        attr(emb, "rotation") <- L
      }
      reduced[[toupper(nm)]] <- emb
    }
  }

  SingleCellExperiment::SingleCellExperiment(
    assays = assays,
    colData = if (length(cd)) S4Vectors::DataFrame(cd, row.names = cells) else S4Vectors::DataFrame(row.names = cells),
    rowData = if (length(rd)) S4Vectors::DataFrame(rd, row.names = genes) else S4Vectors::DataFrame(row.names = genes),
    reducedDims = reduced)
}

#' Read a SingleCellExperiment into an L* dataset.
#'
#' @param sce a `SingleCellExperiment`
#' @return an `lstar_dataset` of kind `"sample"`.
#' @seealso [write_sce()]
#' @export
read_sce <- function(sce) {
  cells <- colnames(sce); genes <- rownames(sce)
  scev <- tryCatch(as.character(utils::packageVersion("SingleCellExperiment")),
                   error = function(e) "?")
  ds <- list(kind = "sample", spec_version = "0.1",
             profiles = c("singlecellexperiment@0.1", paste0("SingleCellExperiment@", scev)),
             dropped = character(0), axes = list(), fields = list())
  ds$axes$cells <- list(labels = cells, origin = "observed", role = "observation")
  ds$axes$genes <- list(labels = genes, origin = "observed", role = "feature")
  add <- function(nm, v, role, span, state = "")
    ds$fields[[nm]] <<- list(role = role, span = span, state = state, subtype = "", values = v)
  add_factor <- function(nm, v, span) {                   # factor col/rowData -> categorical + factor axis
    v <- as.factor(v)
    ds$fields[[nm]] <<- list(role = "label", span = span, state = "", subtype = "",
                             encoding = "categorical", values = v)
    if (is.null(ds$axes[[nm]]))
      ds$axes[[nm]] <<- list(labels = levels(v), origin = "derived", role = "factor", induced_by = nm)
  }

  state_of <- function(a) switch(a, counts = "raw", logcounts = "lognorm", "")
  name_of <- function(a) switch(a, counts = "counts", logcounts = "X", a)
  for (a in SummarizedExperiment::assayNames(sce)) {
    add(name_of(a), Matrix::t(SummarizedExperiment::assay(sce, a)), "measure", c("cells", "genes"), state = state_of(a))
  }
  cdt <- SummarizedExperiment::colData(sce)
  for (col in colnames(cdt)) {
    v <- cdt[[col]]
    if (is.numeric(v)) add(col, as.numeric(v), "measure", "cells")
    else if (is.factor(v)) add_factor(col, v, "cells")
    else add(col, as.character(v), "label", "cells")
  }
  rdt <- SummarizedExperiment::rowData(sce)
  for (col in colnames(rdt)) {
    v <- rdt[[col]]
    if (is.numeric(v)) add(col, as.numeric(v), "measure", "genes")
    else if (is.factor(v)) add_factor(col, v, "genes")
    else add(col, as.character(v), "label", "genes")
  }
  for (rn in SingleCellExperiment::reducedDimNames(sce)) {
    emb <- SingleCellExperiment::reducedDim(sce, rn)
    coord <- tolower(rn)
    labs <- if (!is.null(colnames(emb))) colnames(emb) else paste0(coord, seq_len(ncol(emb)))
    ds$axes[[coord]] <- list(labels = labs, origin = "derived", role = "coordinate")
    add(coord, unname(as.matrix(emb)), "embedding", c("cells", coord))
    rot <- attr(emb, "rotation")
    if (!is.null(rot)) add(paste0(coord, "_loadings"), unname(as.matrix(rot)), "loading", c("genes", coord))
  }
  class(ds) <- "lstar_dataset"
  ds
}
