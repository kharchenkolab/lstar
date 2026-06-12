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

#' Export a Pagoda2 object to an L* (`*.lstar.zarr`) store.
#'
#' Writes counts (raw, cell x gene), the embedding, cluster/cell-type/QC labels, and the
#' viewer profile's cluster sufficient stats + marker tables (`viewer@0.1`). Computes the
#' cluster stats with the shared libstar kernel.
#'
#' @param p2 a Pagoda2 (pagoda2.1) object.
#' @param path output store path (`*.lstar.zarr`); if `NULL`, only the `lstar_dataset` is returned.
#' @param grouping a `cellMeta` column to use as the primary clustering (default `"leiden"`).
#' @return an `lstar_dataset` (invisibly if written).
#' @export
write_pagoda2 <- function(p2, path = NULL, grouping = "leiden") {
  cnt <- as(.p2_raw_counts(p2), "CsparseMatrix")          # cells x genes (CSC gene-major)
  nc <- nrow(cnt); ng <- ncol(cnt)
  cells <- rownames(cnt) %||% paste0("cell", seq_len(nc))
  genes <- colnames(cnt) %||% paste0("g", seq_len(ng))
  dropped <- character()

  axes <- list(
    cells = list(labels = as.character(cells), origin = "observed", role = "observation"),
    genes = list(labels = as.character(genes), origin = "observed", role = "feature"))
  fields <- list(
    counts = list(values = cnt, role = "measure", span = c("cells", "genes"), state = "raw", encoding = "csc"))

  # embedding
  emb <- .first_embedding(p2$embeddings)
  if (!is.null(emb)) {
    axes$umap <- list(labels = c("umap0", "umap1"), origin = "derived", role = "coordinate")
    fields$umap <- list(values = emb, role = "embedding", span = c("cells", "umap"))
  } else dropped <- c(dropped, "embedding")

  # cell metadata: grouping, cell_type, qc
  meta <- p2$cellMeta
  lei <- NULL
  if (!is.null(meta) && grouping %in% names(meta)) {
    lei <- as.character(meta[[grouping]])
    fields[[grouping]] <- list(values = lei, role = "label", span = "cells", encoding = "utf8")
  }
  for (col in c("cell_type", "sample", "condition")) if (!is.null(meta) && col %in% names(meta))
    fields[[col]] <- list(values = as.character(meta[[col]]), role = "label", span = "cells", encoding = "utf8")
  for (col in c("mito", "percent_mito", "n_molecules", "n_genes"))
    if (!is.null(meta) && col %in% names(meta))
      fields[[col]] <- list(values = as.numeric(meta[[col]]), role = "measure", span = "cells")

  # viewer profile: cluster sufficient stats (libstar kernel) + marker tables
  profiles <- "pagoda2@0.1"
  if (!is.null(lei)) {
    groups <- sort(unique(lei))
    code <- as.integer(factor(lei, levels = groups)) - 1L
    K <- length(groups)
    gs <- lstar_cpp_col_sum_by_group(as.double(cnt@x), cnt@p, cnt@i, nc, ng, code, K, TRUE)
    S <- matrix(gs$sum, nrow = K, byrow = TRUE)
    SS <- matrix(gs$sumsq, nrow = K, byrow = TRUE)
    NE <- matrix(gs$n_expr, nrow = K, byrow = TRUE)
    Xl <- cnt; Xl@x <- log1p(Xl@x); grand <- Matrix::colSums(Xl)
    nper <- as.integer(table(factor(lei, levels = groups)))
    lfc <- vapply(seq_len(K), function(g) S[g, ] / max(nper[g], 1) - (grand - S[g, ]) / max(nc - nper[g], 1), numeric(ng))
    padj <- pmin(pmax(exp(-abs(lfc * sqrt(t(NE) + 1))), 1e-12), 1)
    axes[[paste0("groups_", grouping)]] <- list(labels = groups, origin = "derived", role = "feature")
    sg <- c(paste0("groups_", grouping), "genes")
    fields[[paste0("stats_", grouping, "_sum")]]   <- list(values = S,  role = "measure", span = sg)
    fields[[paste0("stats_", grouping, "_sumsq")]] <- list(values = SS, role = "measure", span = sg)
    fields[[paste0("stats_", grouping, "_nexpr")]] <- list(values = NE, role = "measure", span = sg)
    fields[[paste0("markers_", grouping, "_lfc")]]  <- list(values = lfc,  role = "measure", span = c("genes", paste0("groups_", grouping)))
    fields[[paste0("markers_", grouping, "_padj")]] <- list(values = padj, role = "measure", span = c("genes", paste0("groups_", grouping)))
    profiles <- c(profiles, "viewer@0.1")
  } else dropped <- c(dropped, "clustering")

  ds <- list(kind = "sample", spec_version = "0.1", profiles = profiles, dropped = dropped,
             axes = axes, fields = fields)
  class(ds) <- "lstar_dataset"
  if (!is.null(path)) { lstar_write(ds, path); return(invisible(ds)) }
  ds
}
