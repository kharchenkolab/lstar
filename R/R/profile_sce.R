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

  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = assays,
    colData = if (length(cd)) S4Vectors::DataFrame(cd, row.names = cells) else S4Vectors::DataFrame(row.names = cells),
    rowData = if (length(rd)) S4Vectors::DataFrame(rd, row.names = genes) else S4Vectors::DataFrame(row.names = genes),
    reducedDims = reduced)

  # multimodal: rebuild each non-`genes` feature axis as an altExp (the SCE analogue of a Seurat assay)
  for (fax in names(ds$axes)) {
    ax <- ds$axes[[fax]]
    if (!identical(ax$role, "feature") || fax == "genes") next
    feats <- as.character(ax$labels); aas <- list()
    for (nm in names(ds$fields)) {
      f <- ds$fields[[nm]]; sp <- as.character(f$span)
      if (identical(f$role, "measure") && length(sp) == 2 && sp[1] == "cells" && sp[2] == fax) {
        mat <- Matrix::t(as(f$values, "CsparseMatrix")); dimnames(mat) <- list(feats, cells)
        aas[[sub(paste0("^", fax, "\\."), "", nm)]] <- mat
      }
    }
    if (length(aas)) SingleCellExperiment::altExp(sce, fax) <-
      SummarizedExperiment::SummarizedExperiment(assays = aas)
  }

  # colPairs/rowPairs (cell-cell / gene-gene graphs) reconstructed from relation fields
  for (nm in names(ds$fields)) {
    f <- ds$fields[[nm]]; sp <- as.character(f$span)
    if (!identical(f$role, "relation") || length(sp) != 2) next
    if (sp[1] == "cells" && sp[2] == "cells" && startsWith(nm, "colpair_"))
      SingleCellExperiment::colPair(sce, sub("^colpair_", "", nm)) <- as(f$values, "CsparseMatrix")
    else if (sp[1] == "genes" && sp[2] == "genes" && startsWith(nm, "rowpair_"))
      SingleCellExperiment::rowPair(sce, sub("^rowpair_", "", nm)) <- as(f$values, "CsparseMatrix")
  }
  sce
}

#' Read a SingleCellExperiment into an L* dataset.
#'
#' @param sce a `SingleCellExperiment`
#' @return an `lstar_dataset` of kind `"sample"`.
#' @seealso [write_sce()]
#' @export
read_sce <- function(sce) {
  # Real SCEs often have NULL dimnames (cells keyed by a `Barcode` colData column, not colnames) --
  # synthesize stable labels so axes aren't empty (which crashes assay reconstruction on write-back).
  cells <- colnames(sce); genes <- rownames(sce)
  if (is.null(cells)) cells <- if ("Barcode" %in% colnames(SummarizedExperiment::colData(sce)))
    as.character(SummarizedExperiment::colData(sce)$Barcode) else paste0("cell", seq_len(ncol(sce)))
  if (is.null(genes)) genes <- paste0("gene", seq_len(nrow(sce)))
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
  # Real col/rowData columns aren't always plain vectors: S4Vectors `Rle` (run-length) needs unpacking,
  # and a nested `DataFrame`/`GRanges`/list column can't be coerced at all -> record it, don't crash.
  sce_col <- function(col, v, span, n) {
    if (methods::is(v, "Rle")) v <- tryCatch(as.vector(v), error = function(e) v)
    if (is.numeric(v) && length(v) == n) add(col, as.numeric(v), "measure", span)
    else if (is.factor(v) && length(v) == n) add_factor(col, v, span)
    else {
      vv <- tryCatch(if (is.atomic(v)) as.character(v) else NULL, error = function(e) NULL)
      if (!is.null(vv) && length(vv) == n) add(col, vv, "label", span)
      else ds$dropped <<- c(ds$dropped, sprintf("%sData/%s (%s)",
                                                if (span == "cells") "col" else "row", col, class(v)[1]))
    }
  }
  cdt <- SummarizedExperiment::colData(sce)
  for (col in colnames(cdt)) sce_col(col, cdt[[col]], "cells", length(cells))
  rdt <- SummarizedExperiment::rowData(sce)
  for (col in colnames(rdt)) sce_col(col, rdt[[col]], "genes", length(genes))
  # reducedDims/altExps are SingleCellExperiment-only; a plain SummarizedExperiment (or an old/odd SCE
  # subclass, e.g. ReprocessedFluidigmData) lacks the method -> guard, don't crash.
  rdn <- if (methods::is(sce, "SingleCellExperiment"))
    tryCatch(SingleCellExperiment::reducedDimNames(sce), error = function(e) character(0)) else character(0)
  for (rn in rdn) {
    emb <- SingleCellExperiment::reducedDim(sce, rn)
    coord <- tolower(rn)
    labs <- if (!is.null(colnames(emb))) colnames(emb) else paste0(coord, seq_len(ncol(emb)))
    ds$axes[[coord]] <- list(labels = labs, origin = "derived", role = "coordinate")
    add(coord, unname(as.matrix(emb)), "embedding", c("cells", coord))
    rot <- attr(emb, "rotation")
    if (!is.null(rot)) add(paste0(coord, "_loadings"), unname(as.matrix(rot)), "loading", c("genes", coord))
  }
  # altExps are a *second feature space* (ADT / spike-ins) over the same cells -- the SCE analogue of a
  # Seurat assay. Capture each as a feature axis named after the altExp + measures `<altExp>.<assay>`.
  for (ae in tryCatch(SingleCellExperiment::altExpNames(sce), error = function(e) character(0))) {
    asce <- SingleCellExperiment::altExp(sce, ae)
    feats <- rownames(asce); if (is.null(feats)) feats <- paste0(ae, seq_len(nrow(asce)))
    if (ae %in% names(ds$axes)) { ds$dropped <- c(ds$dropped, paste0("altExp/", ae, " (axis-name clash)")); next }
    ds$axes[[ae]] <- list(labels = as.character(feats), origin = "observed", role = "feature")
    for (an in SummarizedExperiment::assayNames(asce))
      add(paste0(ae, ".", an), Matrix::t(as(SummarizedExperiment::assay(asce, an), "CsparseMatrix")),
          "measure", c("cells", ae), state = state_of(an))
  }
  # colPairs/rowPairs (cell-cell SNN/kNN graph, gene-gene graph -- e.g. scran::buildSNNGraph) are the SCE
  # analogue of AnnData obsp/varp -> type as relations over (cells,cells)/(genes,genes), not dropped.
  for (cp in tryCatch(SingleCellExperiment::colPairNames(sce), error = function(e) character(0))) {
    m <- tryCatch(SingleCellExperiment::colPair(sce, cp, asSparse = TRUE), error = function(e) NULL)
    if (!is.null(m) && all(dim(m) == length(cells)))
      add(paste0("colpair_", cp), as(m, "CsparseMatrix"), "relation", c("cells", "cells"))
    else ds$dropped <- c(ds$dropped, paste0("colPair/", cp))
  }
  for (rp in tryCatch(SingleCellExperiment::rowPairNames(sce), error = function(e) character(0))) {
    m <- tryCatch(SingleCellExperiment::rowPair(sce, rp, asSparse = TRUE), error = function(e) NULL)
    if (!is.null(m) && all(dim(m) == length(genes)))
      add(paste0("rowpair_", rp), as(m, "CsparseMatrix"), "relation", c("genes", "genes"))
    else ds$dropped <- c(ds$dropped, paste0("rowPair/", rp))
  }
  md <- names(S4Vectors::metadata(sce))                    # free-form study-level list: not typed -> record
  if (length(md)) ds$dropped <- c(ds$dropped, paste0("metadata/", md))
  class(ds) <- "lstar_dataset"
  ds
}

# ---- SCE v? PACKAGE-FREE read (base R, no SingleCellExperiment) ---------------------------------------
# The `--backend direct` fallback: read an SCE .rds by walking its S4 slots via attr() -- the SE/SCE class
# hierarchy stores everything in plain slots (assays = Assays -> SimpleList -> listData; colData /
# elementMetadata = DataFrames whose columns are in `listData` and whose names are in `rownames`;
# reducedDims live in `int_colData$reducedDims`). Only base R is needed (sparse matrices are dgCMatrix
# from Matrix). Produces the same core L* dataset `read_sce` builds; verified value-equal by convert_cli.sh.
.df_cols <- function(df) { ld <- tryCatch(attr(df, "listData"), error = function(e) NULL); if (is.null(ld)) list() else ld }

.read_sce_direct <- function(so) {
  S <- function(obj, nm) tryCatch(attr(obj, nm), error = function(e) NULL)
  amats <- S(S(S(so, "assays"), "data"), "listData")       # Assays -> SimpleList -> named matrices (genes x cells)
  if (is.null(amats) || !length(amats)) stop("read_sce (direct): object has no assays -- not a SingleCellExperiment?")
  m1 <- amats[[1]]; cd <- S(so, "colData")           # names: the assay matrix's own dimnames are most reliable;
  cells <- colnames(m1); if (is.null(cells)) cells <- S(cd, "rownames")          # else colData rownames
  genes <- rownames(m1); if (is.null(genes)) genes <- S(so, "NAMES")             # else NAMES / rowData rownames
  if (is.null(genes)) genes <- attr(S(so, "elementMetadata"), "rownames")
  if (is.null(cells)) cells <- paste0("cell", seq_len(ncol(m1)))
  if (is.null(genes)) genes <- paste0("g", seq_len(nrow(m1)))
  cells <- as.character(cells); genes <- as.character(genes)

  ds <- list(kind = "sample", spec_version = "0.1", profiles = c("sce@0.1", "object@sce"),
             dropped = character(0), axes = list(), fields = list())
  ds$axes$cells <- list(labels = cells, origin = "observed", role = "observation")
  ds$axes$genes <- list(labels = genes, origin = "observed", role = "feature")
  add <- function(nm, v, role, span, state = "") {
    ds$fields[[nm]] <<- list(role = role, span = span, state = state, subtype = "", values = v)
  }
  add_factor <- function(nm, v, span) {
    v <- as.factor(v)
    ds$fields[[nm]] <<- list(role = "label", span = span, state = "", subtype = "", encoding = "categorical", values = v)
    if (is.null(ds$axes[[nm]])) ds$axes[[nm]] <<- list(labels = levels(v), origin = "derived", role = "factor", induced_by = nm)
  }
  name_of <- function(a) if (identical(a, "counts")) "counts" else if (identical(a, "logcounts")) "X" else a
  state_of <- function(a) if (identical(a, "counts")) "raw" else if (identical(a, "logcounts")) "lognorm" else ""
  for (a in names(amats))                            # genes x cells -> (cells, genes); axes carry the names
    add(name_of(a), Matrix::t(as(amats[[a]], "CsparseMatrix")), "measure", c("cells", "genes"), state = state_of(a))
  sce_col <- function(nm, v, span, n) {                    # decode S4 Rle/factor/vector columns
    v <- tryCatch(if (methods::is(v, "Rle")) as.vector(v) else v, error = function(e) v)
    if (is.numeric(v) && length(v) == n) add(nm, as.numeric(v), "measure", span)
    else if (is.factor(v) && length(v) == n) add_factor(nm, v, span)
    else { vv <- tryCatch(as.character(v), error = function(e) NULL); if (!is.null(vv) && length(vv) == n) add(nm, vv, "label", span) }
  }
  for (col in names(.df_cols(cd))) sce_col(col, .df_cols(cd)[[col]], "cells", length(cells))
  for (col in names(.df_cols(S(so, "elementMetadata")))) sce_col(col, .df_cols(S(so, "elementMetadata"))[[col]], "genes", length(genes))

  rdims <- .df_cols(S(so, "int_colData"))$reducedDims       # int_colData -> reducedDims DFrame -> named embeddings
  if (!is.null(rdims)) for (rn in names(.df_cols(rdims))) {
    emb <- as.matrix(.df_cols(rdims)[[rn]]); if (!nrow(emb)) next
    coord <- tolower(rn)
    labs <- if (!is.null(colnames(emb))) colnames(emb) else paste0(coord, seq_len(ncol(emb)))
    ds$axes[[coord]] <- list(labels = labs, origin = "derived", role = "coordinate")
    add(coord, unname(emb), "embedding", c("cells", coord))
  }
  class(ds) <- "lstar_dataset"
  ds
}
