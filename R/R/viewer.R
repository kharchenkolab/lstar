# Native viewer@0.1 prep (R) ------------------------------------------------------------------------
#
# Adds the viewer@0.1 navigator fields (docs/format.md) to an L* dataset using the SHARED libstar
# kernels (the same C++ core bound to Python/WASM), so an R-prepped store is byte-equivalent to a
# Python/JS-prepped one and to what the viewer computes live. No Python or shell-out.

# Pick a categorical cell label to summarize by: a `label`-role field over the cell axis with 2..60
# distinct values, preferring clustering / cell-type names.
.detect_grouping <- function(ds, cell_axis) {
  pref <- c("leiden", "cluster", "clusters", "cell_type", "celltype", "cell_types",
            "louvain", "seurat_clusters", "annotation")
  cand <- character(0)
  for (nm in names(ds$fields)) {
    f <- ds$fields[[nm]]
    if (is.null(f$role) || f$role != "label") next
    if (is.null(f$span) || length(f$span) != 1L || f$span[1] != cell_axis) next
    v <- f$values
    lv <- if (is.factor(v)) levels(v) else unique(as.character(v[!is.na(v)]))
    if (length(lv) >= 2L && length(lv) <= 60L) cand <- c(cand, nm)
  }
  if (!length(cand))
    stop("extend_for_viewer: no categorical grouping (2..60 levels) over '", cell_axis,
         "'; pass grouping=", call. = FALSE)
  isp <- vapply(cand, function(nm) any(vapply(pref, function(p) grepl(p, tolower(nm), fixed = TRUE),
                                              logical(1))), logical(1))
  cand[order(!isp, cand)][1]
}

#' Extend an L* dataset with the viewer@0.1 navigator fields (native R).
#'
#' Precomputes, via the shared libstar kernels: per-grouping cluster sufficient stats
#' (`stats_<g>_{sum,sumsq,nexpr}`, group-major), 1-vs-rest marker tables (`markers_<g>_{lfc,padj}`,
#' gene-major), a whole-dataset `od_score` (pagoda2 lowess + F-test), a cell-major `counts_cellmajor`
#' physically reordered cluster-contiguous, and its `counts_cellmajor_order` permutation. Stamps the
#' `viewer@0.1` profile. The store this writes is interchangeable with a Python/JS-prepped one.
#'
#' @param ds an `lstar_dataset` (a `counts` measure + at least one categorical cell label).
#' @param grouping primary grouping label (default: auto-detect; clustering/cell-type names preferred).
#' @param also additional grouping labels to also summarize (e.g. `"cell_type"`).
#' @param counts name of the raw counts measure (default `"counts"`).
#' @return `ds` with the navigator fields added and `viewer@0.1` in `ds$profiles`.
#' @seealso [viewer_extend()], [lstar_write_viewer()]
#' @export
extend_for_viewer <- function(ds, grouping = NULL, also = character(0), counts = "counts") {
  cf <- ds$fields[[counts]]
  if (is.null(cf)) stop("extend_for_viewer: no counts measure '", counts, "'", call. = FALSE)
  cnt <- methods::as(cf$values, "CsparseMatrix")                 # cells x genes (CSC)
  span <- if (!is.null(cf$span)) cf$span else c("cells", "genes")
  cell_axis <- span[1]; gene_axis <- span[2]
  nc <- nrow(cnt); ng <- ncol(cnt)

  if (is.null(grouping)) grouping <- .detect_grouping(ds, cell_axis)
  groupings <- unique(c(grouping, also))
  groupings <- groupings[vapply(groupings, function(g) !is.null(ds$fields[[g]]), logical(1))]
  if (!length(groupings))
    stop("extend_for_viewer: no categorical grouping found (pass grouping=)", call. = FALSE)

  # whole-dataset overdispersion (pagoda2 lowess + F-test): mean/var/nobs over log1p, shared kernel.
  g0 <- lstar_cpp_col_sum_by_group(as.double(cnt@x), cnt@p, cnt@i, nc, ng, integer(nc), 1L, TRUE)
  om <- g0$sum / nc; ov <- pmax(g0$sumsq / nc - om^2, 0)
  ds$fields[["od_score"]] <- list(values = lstar_cpp_overdispersion(om, ov, as.integer(g0$n_expr)),
                                  role = "measure", span = gene_axis)

  primary_code <- NULL
  for (gp in groupings) {
    lab <- as.character(ds$fields[[gp]]$values)
    groups <- sort(unique(lab[!is.na(lab)])); K <- length(groups)
    code <- as.integer(factor(lab, levels = groups)) - 1L
    if (gp == grouping) primary_code <- code
    gs <- lstar_cpp_col_sum_by_group(as.double(cnt@x), cnt@p, cnt@i, nc, ng, code, K, TRUE)
    S <- matrix(gs$sum, nrow = K, byrow = TRUE)
    SS <- matrix(gs$sumsq, nrow = K, byrow = TRUE)
    NE <- matrix(gs$n_expr, nrow = K, byrow = TRUE)
    nper <- as.integer(table(factor(lab, levels = groups)))
    mk <- lstar_cpp_markers_one_vs_rest(gs$sum, gs$n_expr, nper, K, ng, as.double(nc))
    lfc  <- matrix(mk$lfc,  nrow = ng, ncol = K, byrow = TRUE)      # genes x K (gene-major)
    padj <- matrix(mk$padj, nrow = ng, ncol = K, byrow = TRUE)
    gax <- paste0("groups_", gp)
    ds$axes[[gax]] <- list(labels = groups, origin = "derived", role = "feature")
    sg <- c(gax, gene_axis); mg <- c(gene_axis, gax)
    ds$fields[[paste0("stats_", gp, "_sum")]]   <- list(values = S,  role = "measure", span = sg)
    ds$fields[[paste0("stats_", gp, "_sumsq")]] <- list(values = SS, role = "measure", span = sg)
    ds$fields[[paste0("stats_", gp, "_nexpr")]] <- list(values = NE, role = "measure", span = sg)
    ds$fields[[paste0("markers_", gp, "_lfc")]]  <- list(values = lfc,  role = "measure", span = mg)
    ds$fields[[paste0("markers_", gp, "_padj")]] <- list(values = padj, role = "measure", span = mg)
  }
  if (is.null(primary_code)) primary_code <- integer(nc)

  # cell-major counts, physically reordered cluster-contiguous; counts_cellmajor_order = pos_of (cell
  # -> physical row), so the reader's `<field>_order` sibling coalesces a cluster/lasso into ~1 read.
  perm <- order(primary_code, method = "radix")                  # 1-based: perm[p] = cell at row p
  pos_of <- integer(nc); pos_of[perm] <- seq_len(nc) - 1L
  ds$fields[["counts_cellmajor"]] <- list(values = methods::as(cnt[perm, , drop = FALSE], "RsparseMatrix"),
                                          role = "measure", span = c(cell_axis, gene_axis),
                                          state = "raw", encoding = "csr")
  ds$fields[["counts_cellmajor_order"]] <- list(values = as.double(pos_of), role = "measure",
                                                span = cell_axis, state = "permutation")

  if (!("viewer@0.1" %in% ds$profiles)) ds$profiles <- c(ds$profiles, "viewer@0.1")
  if (!methods::is(ds, "lstar_dataset")) class(ds) <- "lstar_dataset"
  ds
}

#' Extend an existing L* store in place with the viewer@0.1 navigator fields.
#'
#' Reads the store, runs [extend_for_viewer()] (native R, shared libstar kernels), and writes it back.
#' Mirrors the `lstar viewer <store>` CLI verb.
#'
#' @param path path to an existing `*.lstar.zarr` store (must contain `counts` + a categorical label).
#' @param grouping,also passed to [extend_for_viewer()].
#' @return the store `path`, invisibly.
#' @seealso [extend_for_viewer()], [lstar_write_viewer()]
#' @export
viewer_extend <- function(path, grouping = NULL, also = character(0)) {
  path <- path.expand(path)
  if (!dir.exists(path)) stop(sprintf("viewer_extend: store not found: %s", path), call. = FALSE)
  ds <- lstar_read(path)
  ds <- extend_for_viewer(ds, grouping = grouping, also = also)
  lstar_write(ds, path)
  invisible(path)
}

#' Write an R dataset to an L* store, optionally extending it for the viewer first.
#'
#' Convenience wrapper: when `viewer = TRUE` (default) runs [extend_for_viewer()] before writing.
#'
#' @param ds an `lstar_dataset`.
#' @param path output store path (a `*.lstar.zarr` directory).
#' @param viewer if `TRUE` (default), extend via [extend_for_viewer()] before writing.
#' @param ... further arguments forwarded to [lstar_write()] (e.g. `chunk_elems`, `compression`).
#' @return the output `path`, invisibly.
#' @seealso [extend_for_viewer()], [lstar_write()]
#' @export
lstar_write_viewer <- function(ds, path, viewer = TRUE, ...) {
  if (isTRUE(viewer)) ds <- extend_for_viewer(ds)
  lstar_write(ds, path, ...)
  invisible(path)
}
