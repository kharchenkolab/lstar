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

  # cell metadata: grouping, cell_type, qc. A categorical cell label is stored as a `categorical`
  # field (factor: codes + ordered category set) and *induces* a bare-named `factor` axis (role=factor,
  # induced_by) whose labels ARE its categories -- so per-group results (below) are fields over it.
  meta <- p2$cellMeta
  lei <- NULL; groups <- NULL; K <- 0L
  if (!is.null(meta) && grouping %in% names(meta)) {
    lei <- as.character(meta[[grouping]])
    groups <- sort(unique(lei[!is.na(lei)])); K <- length(groups)
    fields[[grouping]] <- list(values = factor(lei, levels = groups), role = "label",
                               span = "cells", encoding = "categorical")
    axes[[grouping]] <- list(labels = groups, origin = "derived", role = "factor", induced_by = grouping)
  }
  for (col in c("cell_type", "sample", "condition")) if (!is.null(meta) && col %in% names(meta)) {
    cv <- as.character(meta[[col]]); lv <- sort(unique(cv[!is.na(cv)]))
    fields[[col]] <- list(values = factor(cv, levels = lv), role = "label", span = "cells", encoding = "categorical")
    axes[[col]] <- list(labels = lv, origin = "derived", role = "factor", induced_by = col)
  }
  for (col in c("mito", "percent_mito", "n_molecules", "n_genes"))
    if (!is.null(meta) && col %in% names(meta))
      fields[[col]] <- list(values = as.numeric(meta[[col]]), role = "measure", span = "cells")

  # viewer@0.1 profile (docs/format.md): cluster sufficient stats + marker tables + a whole-dataset
  # overdispersion score + a cell-major counts copy, all via the shared libstar kernels so the store
  # matches what Python/JS produce. stats are group-major (K x genes); markers are gene-major
  # (genes x K); od_score is per-gene.
  profiles <- "pagoda2@0.1"
  if (!is.null(lei)) {
    code <- as.integer(factor(lei, levels = groups)) - 1L
    gs <- lstar_cpp_col_sum_by_group(as.double(cnt@x), cnt@p, cnt@i, nc, ng, code, K, TRUE)
    S <- matrix(gs$sum, nrow = K, byrow = TRUE)
    SS <- matrix(gs$sumsq, nrow = K, byrow = TRUE)
    NE <- matrix(gs$n_expr, nrow = K, byrow = TRUE)
    nper <- as.integer(table(factor(lei, levels = groups)))
    mk <- lstar_cpp_markers_one_vs_rest(gs$sum, gs$n_expr, nper, K, ng, as.double(nc))   # shared kernel
    lfc  <- matrix(mk$lfc,  nrow = ng, ncol = K, byrow = TRUE)        # genes x K (gene-major)
    padj <- matrix(mk$padj, nrow = ng, ncol = K, byrow = TRUE)
    sg <- c(grouping, "genes")                                       # stats: group-major (factor, genes)
    mg <- c("genes", grouping)                                       # markers: gene-major (genes, factor)
    fields[[paste0("stats_", grouping, "_sum")]]   <- list(values = S,  role = "measure", span = sg)
    fields[[paste0("stats_", grouping, "_sumsq")]] <- list(values = SS, role = "measure", span = sg)
    fields[[paste0("stats_", grouping, "_nexpr")]] <- list(values = NE, role = "measure", span = sg)
    fields[[paste0("markers_", grouping, "_lfc")]]  <- list(values = lfc,  role = "measure", span = mg)
    fields[[paste0("markers_", grouping, "_padj")]] <- list(values = padj, role = "measure", span = mg)

    # whole-dataset overdispersion (pagoda2 lowess + F-test, shared kernel): mean/var/nobs over log1p.
    g0 <- lstar_cpp_col_sum_by_group(as.double(cnt@x), cnt@p, cnt@i, nc, ng, integer(nc), 1L, TRUE)
    om <- g0$sum / nc; ov <- pmax(g0$sumsq / nc - om^2, 0)
    od <- lstar_cpp_overdispersion(om, ov, as.integer(g0$n_expr))
    fields[["od_score"]] <- list(values = od, role = "measure", span = "genes")

    # cell-major (CSR) copy of counts -- the per-cell-read substrate.
    fields[["counts_cellmajor"]] <- list(values = methods::as(cnt, "RsparseMatrix"),
                                          role = "measure", span = c("cells", "genes"),
                                          state = "raw", encoding = "csr")
    profiles <- c(profiles, "viewer@0.1")
  } else dropped <- c(dropped, "clustering")

  ds <- list(kind = "sample", spec_version = "0.1", profiles = profiles, dropped = dropped,
             axes = axes, fields = fields)
  class(ds) <- "lstar_dataset"
  if (!is.null(path)) { lstar_write(ds, path); return(invisible(ds)) }
  ds
}
