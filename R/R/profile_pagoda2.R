# pagoda2 profile: a Pagoda2 (pagoda2.1 R6) object -> an L* store with the canonical viewer
# schema. Mirrors profile_conos / profile_sce; reads through the documented public accessors
# (getRawCounts, embeddings, cellMeta) so it survives slot churn, routing what it can't find to
# `dropped`. The cluster sufficient stats are computed via the bound libstar kernel
# (lstar_cpp_col_sum_by_group) — the same C++ core compiled to WASM and bound to Python.

.p2_raw_counts <- function(p2) {
  # pagoda2.1: getRawCounts() returns a cell x gene dgCMatrix; $counts is a hard error.
  if (is.function(p2$getRawCounts)) return(p2$getRawCounts())
  stop("pagoda2 object exposes no getRawCounts(); is this pagoda2.1?")
}

.first_embedding <- function(emb) {
  if (is.matrix(emb) || methods::is(emb, "Matrix")) {
    m <- as.matrix(emb); if (ncol(m) >= 2) return(m[, 1:2, drop = FALSE])
  }
  if (is.list(emb)) for (e in emb) { r <- .first_embedding(e); if (!is.null(r)) return(r) }
  NULL
}

.p2_feature_axis <- function(ft) switch(if (is.null(ft)) "gene" else ft,
  gene = "genes", protein = "proteins", peak = "peaks", paste0(ft, "_features"))

#' Read a Pagoda2 object into an L* dataset.
#'
#' The pagoda2 counterpart of [read_seurat()]: extracts **every facet** (RNA/ADT/ATAC… raw counts over
#' their feature axes), **all embeddings** (`embeddings[[reduction]][[name]]`), and **all `cellMeta`
#' columns** (categorical → a `label` field inducing a factor axis; numeric → a `measure`), each carrying
#' the provenance (`pagoda2::Pagoda2$fromLstar()` uses it to reconstruct the object).
#' Duck-typed: a facet-aware pagoda2.1 object uses `listFacets()`/`getFacet()`; a simpler object falls
#' back to `getRawCounts()` (single RNA facet). Viewer navigators are NOT added here — call
#' [extend_for_viewer()] for that.
#'
#' @param p2 a Pagoda2 (pagoda2.1) object (or a `getRawCounts()`/`embeddings`/`cellMeta`-shaped object).
#' @return an `lstar_dataset` (kind `"sample"`).
#' @seealso [extend_for_viewer()], [read_seurat()]
#' @export
read_pagoda2 <- function(p2) {
  facet_names <- if (is.function(p2$listFacets)) p2$listFacets() else "RNA"
  dfac <- if (!is.null(p2$defaultFacet)) p2$defaultFacet else "RNA"
  axes <- list(); fields <- list(); dropped <- character(); cells <- NULL

  for (fn in facet_names) {
    if (is.function(p2$getFacet)) {
      f <- p2$getFacet(fn); raw <- as(f$rawCounts, "CsparseMatrix")
      ft <- f$featureType; model <- f$modelType; dred <- f$defaultReduction
    } else {
      raw <- as(.p2_raw_counts(p2), "CsparseMatrix"); ft <- "gene"; model <- "plain"; dred <- "PCA"
    }
    if (is.null(cells)) cells <- rownames(raw) %||% paste0("cell", seq_len(nrow(raw)))
    fax <- .p2_feature_axis(ft)
    axes[[fax]] <- list(labels = as.character(colnames(raw) %||% paste0(fax, seq_len(ncol(raw)))),
                        origin = "observed", role = "feature")
    fname <- if (identical(fn, dfac)) "counts" else paste0(fn, ".counts")
    fields[[fname]] <- list(values = raw, role = "measure", span = c("cells", fax), state = "raw",
                            encoding = "csc", provenance = list(facet = fn, feature_axis = fax,
                            featureType = if (is.null(ft)) "gene" else ft,
                            model = if (is.null(model)) "plain" else model,
                            defaultReduction = if (is.null(dred)) "PCA" else dred))
  }
  axes <- c(list(cells = list(labels = as.character(cells), origin = "observed", role = "observation")), axes)

  # embeddings: one field per embeddings[[reduction]][[name]] (role=embedding), provenance keeps both.
  embs <- p2$embeddings
  if (is.list(embs)) for (red in names(embs)) for (nm in names(embs[[red]])) {
    em <- as.matrix(embs[[red]][[nm]]); if (length(dim(em)) != 2 || ncol(em) < 1) next
    field <- make.names(tolower(if (length(embs) > 1) paste0(red, "_", nm) else nm))
    eax <- paste0(field, "_dim")
    axes[[eax]] <- list(labels = as.character(colnames(em) %||% paste0(tolower(nm), seq_len(ncol(em)))),
                        origin = "derived", role = "coordinate")
    fields[[field]] <- list(values = em, role = "embedding", span = c("cells", eax),
                            provenance = list(reduction = red, embedding = nm))
  }

  # cellMeta: categorical -> label (induces a factor axis); numeric -> measure. provenance routes it back.
  meta <- p2$cellMeta
  if (is.data.frame(meta)) for (col in colnames(meta)) {
    v <- meta[[col]]
    if (is.factor(v) || is.character(v)) {
      cv <- as.character(v); lv <- sort(unique(cv[!is.na(cv)])); if (!length(lv)) next
      fields[[col]] <- list(values = factor(cv, levels = lv), role = "label", span = "cells",
                            encoding = "categorical", provenance = list(cellmeta = col))
      axes[[col]] <- list(labels = lv, origin = "derived", role = "factor", induced_by = col)
    } else if (is.numeric(v)) {
      fields[[col]] <- list(values = as.numeric(v), role = "measure", span = "cells",
                            provenance = list(cellmeta = col))
    } else dropped <- c(dropped, paste0("cellMeta/", col))
  }

  ds <- list(kind = "sample", spec_version = "0.1", profiles = "pagoda2@0.1", dropped = dropped,
             axes = axes, fields = fields)
  class(ds) <- "lstar_dataset"
  ds
}
