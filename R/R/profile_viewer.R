# viewer@0.1 profile (R) — precompute the reductions an interactive viewer needs so common
# queries touch no matrix. Mirrors python/src/lstar/profiles/viewer.py field-for-field (same
# numbers; cross-checked in tests/testthat/test-viewer.R). The heavy per-(group, gene) reduction
# runs on the shared libstar kernel (lstar_cpp_col_sum_by_group, the same C++ bound to WASM and
# Python); the light glue (markers, od-gene selection, the cell-major panel) uses Matrix. This
# keeps R self-sufficient: read_*() -> write_viewer() -> lstar_write(), no other runtime in the loop.

#' Add the `viewer@0.1` profile to an L* dataset.
#'
#' For a single-cell dataset with a `counts` measure over `(cells, genes)` and a categorical
#' grouping label over `(cells)`, precomputes the grouping's cluster sufficient stats and ranked
#' marker tables (a *global* question over a fixed clustering), a cluster-coherent `cell_order`, and
#' `counts_cellmajor` — the counts in **cell-major (CSR)** orientation over ALL genes. Scope-dependent work
#' (selection DE, overdispersed-gene/HVG selection) is intentionally NOT precomputed: it is computed
#' on the fly by subsampling the cells in the current scope and reducing over all genes, because a
#' globally-chosen gene subset is wrong for any local question. The profile is **recomputed from the
#' current `(counts, grouping)` and overwrites** same-named fields, so it can never go stale.
#'
#' @param ds an `lstar_dataset` (with a counts measure and a grouping label over cells).
#' @param grouping the cell label to summarize by (default `"leiden"`).
#' @param counts name of the raw counts measure (default `"counts"`).
#' @param n_od retained for signature compatibility; unused (the gene scope is on-the-fly, not baked).
#' @return `ds` with the viewer profile added (`viewer@0.1` in `ds$profiles`).
#' @export
write_viewer <- function(ds, grouping = "leiden", counts = "counts", n_od = 300L) {
  if (is.null(ds$fields[[counts]])) stop("write_viewer: no counts measure '", counts, "'")
  if (is.null(ds$fields[[grouping]])) stop("write_viewer: no grouping field '", grouping, "'")
  cnt <- methods::as(ds$fields[[counts]]$values, "CsparseMatrix")     # cells x genes (CSC)
  nc <- nrow(cnt); ng <- ncol(cnt)
  genes <- as.character(ds$axes$genes$labels)

  lab <- as.character(ds$fields[[grouping]]$values)
  groups <- sort(unique(lab))
  K <- length(groups)
  code <- as.integer(factor(lab, levels = groups)) - 1L              # 0-based group code per cell

  # cluster sufficient stats over log1p — the shared libstar kernel (returns flat row-major K x ng)
  gs <- lstar_cpp_col_sum_by_group(as.double(cnt@x), cnt@p, cnt@i, nc, ng, code, K, TRUE)
  S  <- matrix(gs$sum,    nrow = K, byrow = TRUE)
  SS <- matrix(gs$sumsq,  nrow = K, byrow = TRUE)
  NE <- matrix(gs$n_expr, nrow = K, byrow = TRUE)

  # log1p counts for the grand totals the marker tables need
  Xl <- cnt; Xl@x <- log1p(Xl@x)
  grand <- Matrix::colSums(Xl)
  nper <- as.integer(table(factor(lab, levels = groups)))

  # marker tables: lfc = group mean(log1p) - rest mean(log1p); padj a monotone proxy  (genes x groups)
  lfc <- vapply(seq_len(K), function(g) S[g, ] / max(nper[g], 1) -
                  (grand - S[g, ]) / max(nc - nper[g], 1), numeric(ng))
  padj <- pmin(pmax(exp(-abs(lfc * sqrt(t(NE) + 1))), 1e-12), 1)

  order_cells <- order(code, method = "radix") - 1L                       # 0-based stable cell permutation

  # counts in cell-major (CSR) orientation, all genes — the substrate for on-the-fly, scope-correct
  # selection DE / overdispersion (subsample cells, reduce over ALL genes). Never a baked gene subset.
  panel <- methods::as(cnt, "RsparseMatrix")

  gax <- paste0("groups_", grouping)
  ds$axes[[gax]] <- list(labels = groups, origin = "derived", role = "feature")
  sg <- c(gax, "genes")
  ds$fields[[paste0("stats_", grouping, "_sum")]]   <- list(role = "measure", span = sg, values = S)
  ds$fields[[paste0("stats_", grouping, "_sumsq")]] <- list(role = "measure", span = sg, values = SS)
  ds$fields[[paste0("stats_", grouping, "_nexpr")]] <- list(role = "measure", span = sg, values = NE)
  ds$fields[[paste0("markers_", grouping, "_lfc")]]  <- list(role = "measure", span = c("genes", gax), values = lfc)
  ds$fields[[paste0("markers_", grouping, "_padj")]] <- list(role = "measure", span = c("genes", gax), values = padj)
  ds$fields[["cell_order"]] <- list(role = "measure", span = "cells", state = "permutation", values = order_cells)
  ds$fields[["counts_cellmajor"]] <- list(role = "measure", span = c("cells", "genes"),
                                           state = "raw", encoding = "csr", values = panel)
  if (!("viewer@0.1" %in% ds$profiles)) ds$profiles <- c(ds$profiles, "viewer@0.1")
  if (!methods::is(ds, "lstar_dataset")) class(ds) <- "lstar_dataset"
  ds
}
