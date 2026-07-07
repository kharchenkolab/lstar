# Native viewer@0.1 prep (R) ------------------------------------------------------------------------
#
# Adds the viewer@0.1 navigator fields (docs/format.md) to an L* dataset using the SHARED libstar
# kernels (the same C++ core bound to Python/WASM), so an R-prepped store is byte-equivalent to a
# Python/JS-prepped one and to what the viewer computes live. No Python or shell-out.

# Canonical viewer@0.1 grouping-detection policy (single source: mirrors Python `_PREFERRED_GROUPINGS`
# and js/core/policy.ts; enforced against conformance/viewer_policy.json by conformance/policy_linter.py).
.VIEWER_PREFERRED_GROUPINGS <- c("leiden", "cluster", "clusters", "cell_type", "celltype", "cell_types",
                                 "louvain", "seurat_clusters", "annotation", "cluster_label")
.VIEWER_MIN_GROUPS <- 2L; .VIEWER_MAX_GROUPS <- 60L    # single-sourced with viewer_policy.json (policy_linter)
.VIEWER_LOGNORM_NAMES <- c("X", "data", "logcounts")   # lognorm measure-name fallback
.VIEWER_PREFERRED_EMBEDDINGS <- c("umap")
.VIEWER_HILBERT_GRID <- 1024L

# Preference rank of an embedding name (mirrors .viewer_grouping_rank): first preferred term that is a
# substring of the lowercased name; non-matches rank last.
.viewer_embedding_rank <- function(nm) {
  low <- tolower(nm)
  for (i in seq_along(.VIEWER_PREFERRED_EMBEDDINGS))
    if (grepl(.VIEWER_PREFERRED_EMBEDDINGS[i], low, fixed = TRUE)) return(i)
  length(.VIEWER_PREFERRED_EMBEDDINGS) + 1L
}

# Preference rank of a label name: index of the first preferred term that is a substring of the lowercased
# name (matches sort first, by list position); non-matches rank last. Mirrors Python `_rank`.
.viewer_grouping_rank <- function(nm) {
  low <- tolower(nm)
  for (i in seq_along(.VIEWER_PREFERRED_GROUPINGS))
    if (grepl(.VIEWER_PREFERRED_GROUPINGS[i], low, fixed = TRUE)) return(i)
  length(.VIEWER_PREFERRED_GROUPINGS) + 1L
}

# ALL usable cell groupings (not just one): `label`-role fields over the cell axis with 2..60 distinct
# values, preferred clustering/cell-type names first (by list position), then alphabetical -- identical to
# Python `_detect_groupings`. (Was `.detect_grouping`, which returned a single primary and diverged.)
.detect_groupings <- function(ds, cell_axis) {
  cand <- character(0)
  for (nm in names(ds$fields)) {
    f <- ds$fields[[nm]]
    if (is.null(f$role) || f$role != "label") next
    if (is.null(f$span) || length(f$span) != 1L || f$span[1] != cell_axis) next
    if (identical(f$subtype, "active_ident")) next         # Seurat active-idents mirror (== a clustering) -> not a grouping
    v <- f$values
    if (!is.factor(v) && !is.character(v)) next           # string-like labels only (match Python: skip numeric/logical)
    lv <- if (is.factor(v)) levels(v) else unique(v[!is.na(v)])
    if (length(lv) >= .VIEWER_MIN_GROUPS && length(lv) <= .VIEWER_MAX_GROUPS) cand <- c(cand, nm)
  }
  if (!length(cand))
    stop("extend_for_viewer: no categorical grouping (2..60 levels) over '", cell_axis,
         "'; pass grouping=", call. = FALSE)
  ranks <- vapply(cand, .viewer_grouping_rank, integer(1))
  cand[order(ranks, cand)]
}

# Pick the primary embedding that keys the within-cluster (Hilbert) locality order: an `embedding`-role
# field over the cell axis with >=2 dims, preferring umap. Mirrors Python's `_detect_embedding` so the
# shared core reorder gets the same secondary key on every surface. NULL when no embedding is present.
.detect_embedding <- function(ds, cell_axis) {
  cand <- character(0)
  for (nm in names(ds$fields)) {
    f <- ds$fields[[nm]]
    if (is.null(f$role) || f$role != "embedding") next
    if (is.null(f$span) || length(f$span) < 1L || f$span[[1]] != cell_axis) next
    v <- f$values
    if (is.null(dim(v)) || ncol(v) < 2L) next
    cand <- c(cand, nm)
  }
  if (!length(cand)) return(NULL)
  cand[order(vapply(cand, .viewer_embedding_rank, integer(1)), cand)][1]
}

#' Extend an L* dataset with the viewer@0.1 navigator fields (native R).
#'
#' Precomputes, via the shared libstar kernels: per-grouping cluster sufficient stats
#' (`stats_<g>_{sum,sumsq,nexpr}`, group-major), 1-vs-rest marker tables (`markers_<g>_{lfc,padj}`,
#' gene-major), a whole-dataset `od_score` (pagoda2 lowess + F-test), a cell-major `counts_cellmajor`
#' physically reordered cluster-contiguous, and its `counts_cellmajor_order` permutation. Stamps the
#' `viewer@0.1` profile. The store this writes is interchangeable with a Python/JS-prepped one.
#'
#' @param ds an `lstar_dataset` (a raw-counts measure + at least one categorical cell label).
#' @param grouping primary grouping label (default: auto-detect; clustering/cell-type names preferred).
#' @param also additional grouping labels to also summarize (e.g. `"cell_type"`).
#' @param primary the grouping the *viewer opens on*. Hoisted to the front of the prepared groupings, so it
#'   keys the `counts_cellmajor` locality reorder AND is summarized first -- the eager-prepare a fast launch
#'   waits on. Unlike ordering the groupings by hand, `primary` composes with auto-detect: `primary="cell_type"`
#'   with `grouping=NULL` preps *every* detected grouping but keys the reorder on `cell_type` (the auto-detect
#'   policy prefers clusterings, but the viewer may open on a cell-type annotation). `NULL` = current behavior.
#' @param counts force the count measure to build from (`log1p` unless it is already log-normalized).
#'   `NULL` (default) lets `basis` choose.
#' @param basis how to choose the measure by state (`NULL` == `"auto"`): `"auto"` prefers a raw measure
#'   (`log1p`-transformed) and, when none is present, falls back to a log-normalized one (used as-is;
#'   stats are var-of-lognorm, so HVG/markers are approximate -- a warning is emitted). `"raw"`/
#'   `"lognorm"` force that state. A scaled/z-scored measure is never used.
#' @param order `"hybrid"` (default) physically reorders `counts_cellmajor` (cluster-contiguous, then a
#'   Hilbert curve over the embedding) and adds `counts_cellmajor_order`; `"none"` keeps rows in cell
#'   order and omits the permutation.
#' @param markers if `TRUE` (default), also compute the 1-vs-rest `markers_<g>_{lfc,padj}` tables.
#' @return `ds` with the navigator fields added and `viewer@0.1` in `ds$profiles`.
#' @seealso [viewer_extend()], [lstar_write_viewer()]
#' @export
extend_for_viewer <- function(ds, grouping = NULL, also = character(0), counts = NULL, basis = NULL,
                              order = "hybrid", markers = TRUE, primary = NULL) {
  sel <- .viewer_counts_basis(ds, counts, basis)   # pick basis by state, not the literal name "counts"
  cf <- ds$fields[[sel$name]]; use_lognorm <- sel$log1p
  if (!sel$log1p && is.null(counts) && (is.null(basis) || identical(basis, "auto")))
    warning(sprintf(paste0("extend_for_viewer: no raw counts found; prepped from the log-normalized measure '%s'. ",
                           "od_score (HVG) and markers are approximate (var-of-lognorm, not var-of-log1p(counts)). ",
                           "Pass counts=<field> or basis='raw' if a raw measure is available."), sel$name), call. = FALSE)
  cnt <- methods::as(cf$values, "CsparseMatrix")                 # cells x genes (CSC)
  span <- if (!is.null(cf$span)) cf$span else c("cells", "genes")
  cell_axis <- span[1]; gene_axis <- span[2]
  nc <- nrow(cnt); ng <- ncol(cnt)

  if (is.null(grouping)) groupings <- .detect_groupings(ds, cell_axis)   # ALL detected (parity with Python)
  else groupings <- unique(c(grouping, also))
  groupings <- groupings[vapply(groupings, function(g) !is.null(ds$fields[[g]]), logical(1))]
  # Hoist the viewer's primary grouping to the front (guaranteed present): it keys the reorder + is summarized
  # first. Composes with auto-detect above -- the rest of the groupings are still prepped. (Parity: Python/JS.)
  if (!is.null(primary)) {
    if (is.null(ds$fields[[primary]]))
      stop(sprintf("extend_for_viewer: primary='%s' is not a field in the dataset", primary), call. = FALSE)
    # must be a 1-D grouping over the CELL axis (else a cryptic reorder crash); span==cell_axis is the check
    # identical across Py/R/JS (their detection predicates differ, but this structural one does not).
    if (!identical(as.character(ds$fields[[primary]]$span), as.character(cell_axis)))
      stop(sprintf("extend_for_viewer: primary='%s' must be a grouping over the cell axis '%s' (a 1-D label)",
                   primary, cell_axis), call. = FALSE)
    groupings <- unique(c(primary, groupings[groupings != primary]))
  }
  if (!length(groupings))
    stop("extend_for_viewer: no categorical grouping found (pass grouping=)", call. = FALSE)
  primary_grouping <- groupings[1]                                       # reorder keys on the first grouping

  # All fields below are viewer@0.1 caches (regenerable from counts; non-viewer converters drop+record).
  cache <- list(cache = "viewer@0.1")

  # whole-dataset overdispersion (pagoda2 lowess + F-test): mean/var/nobs over log1p, shared kernel.
  g0 <- lstar_cpp_col_sum_by_group(as.double(cnt@x), cnt@p, cnt@i, nc, ng, integer(nc), 1L, use_lognorm)
  om <- g0$sum / nc; ov <- pmax(g0$sumsq / nc - om^2, 0)
  ds$fields[["od_score"]] <- list(values = lstar_cpp_overdispersion(om, ov, as.integer(g0$n_expr)),
                                  role = "measure", span = gene_axis, provenance = cache)

  primary_code <- NULL
  for (gp in groupings) {
    lab <- as.character(ds$fields[[gp]]$values)
    groups <- sort(unique(lab[!is.na(lab)])); K <- length(groups)
    code <- as.integer(factor(lab, levels = groups)) - 1L
    if (identical(gp, primary_grouping)) primary_code <- code
    gs <- lstar_cpp_col_sum_by_group(as.double(cnt@x), cnt@p, cnt@i, nc, ng, code, K, use_lognorm)
    S <- matrix(gs$sum, nrow = K, byrow = TRUE)
    SS <- matrix(gs$sumsq, nrow = K, byrow = TRUE)
    NE <- matrix(gs$n_expr, nrow = K, byrow = TRUE)
    gax <- paste0("groups_", gp)
    ds$axes[[gax]] <- list(labels = groups, origin = "derived", role = "feature")
    sg <- c(gax, gene_axis); mg <- c(gene_axis, gax)
    ds$fields[[paste0("stats_", gp, "_sum")]]   <- list(values = S,  role = "measure", span = sg, provenance = cache)
    ds$fields[[paste0("stats_", gp, "_sumsq")]] <- list(values = SS, role = "measure", span = sg, provenance = cache)
    ds$fields[[paste0("stats_", gp, "_nexpr")]] <- list(values = NE, role = "measure", span = sg, provenance = cache)
    if (markers) {                                                 # optional 1-vs-rest marker tables (parity with Python markers=)
      nper <- as.integer(table(factor(lab, levels = groups)))
      mk <- lstar_cpp_markers_one_vs_rest(gs$sum, gs$n_expr, nper, K, ng, as.double(nc))
      ds$fields[[paste0("markers_", gp, "_lfc")]]  <- list(values = matrix(mk$lfc,  nrow = ng, ncol = K, byrow = TRUE), role = "measure", span = mg, provenance = cache)
      ds$fields[[paste0("markers_", gp, "_padj")]] <- list(values = matrix(mk$padj, nrow = ng, ncol = K, byrow = TRUE), role = "measure", span = mg, provenance = cache)
    }
  }
  if (is.null(primary_code)) primary_code <- integer(nc)

  # cell-major counts, physically reordered via the SHARED core reorder (cluster code, then Hilbert of
  # the embedding when present) -- byte-identical to the Python/JS surfaces. counts_cellmajor_order =
  # pos_of (cell -> physical row), so the reader's `<field>_order` sibling coalesces a cluster/lasso read.
  if (identical(order, "hybrid")) {                              # locality reorder + the _order permutation
    emb_nm <- .detect_embedding(ds, cell_axis)
    emb <- if (!is.null(emb_nm)) as.matrix(ds$fields[[emb_nm]]$values)[, 1:2, drop = FALSE] else numeric(0)
    pos_of <- lstar_cpp_viewer_cell_order(as.integer(primary_code), as.double(emb), nc, .VIEWER_HILBERT_GRID)  # 0-based
    perm <- integer(nc); perm[pos_of + 1L] <- seq_len(nc)        # inverse: perm[p] = cell at 1-based row p
    ds$fields[["counts_cellmajor_order"]] <- list(values = as.double(pos_of), role = "measure",
                                                  span = cell_axis, state = "permutation",
                                                  provenance = c(cache, list(group = primary_grouping)))
  } else {
    perm <- seq_len(nc)                                          # order="none": rows stay in cell order, no _order
  }
  ds$fields[["counts_cellmajor"]] <- list(values = methods::as(cnt[perm, , drop = FALSE], "RsparseMatrix"),
                                          role = "measure", span = c(cell_axis, gene_axis),
                                          state = if (!is.null(cf$state)) cf$state else "raw",
                                          encoding = "csr", provenance = cache)

  if (!("viewer@0.1" %in% ds$profiles)) ds$profiles <- c(ds$profiles, "viewer@0.1")
  if (!methods::is(ds, "lstar_dataset")) class(ds) <- "lstar_dataset"
  ds
}

# Select the viewer count basis by content/state (not the literal name "counts"), so a converter that
# named its raw matrix `X` or a modality is still preppable. A raw measure is preferred and log1p'd;
# `counts=` forces a field; `basis="lognorm"` preps (approximately) from a log-normalized measure.
.viewer_counts_basis <- function(ds, counts = NULL, basis = NULL) {
  # Contract identical to Python's _select_counts_basis and JS's selectCountsBasis: basis 'auto'
  # (default) prefers a raw measure (log1p) and falls back to a log-normalized one (as-is); 'raw'/
  # 'lognorm' force one; a scaled/z-scored measure is never used. counts=<field> forces a measure.
  is_cxg <- function(f) identical(f$role, "measure") && !is.null(f$span) && length(f$span) == 2 &&
    grepl("^cells", as.character(f$span[[1]]))
  twod <- Filter(function(nm) is_cxg(ds$fields[[nm]]), names(ds$fields))
  st <- function(nm) { s <- ds$fields[[nm]]$state; if (is.null(s)) NA_character_ else s }
  present <- if (length(twod)) paste(vapply(twod, function(nm) sprintf("%s[%s]", nm, st(nm)), ""), collapse = ", ") else "(none)"
  raw_pick <- function() {
    # the literal "counts" shortcut EXCLUDES a scaled/z-scored measure (symmetric with lognorm_pick): a
    # field named "counts" that is actually scaled must NOT be picked as raw + log1p'd -- fall through.
    if ("counts" %in% twod && !identical(st("counts"), "scaled")) return("counts")
    p <- twod[vapply(twod, function(nm) identical(st(nm), "raw"), logical(1))]
    if (length(p)) p[1] else NULL
  }
  lognorm_pick <- function() {
    p <- twod[vapply(twod, function(nm) identical(st(nm), "lognorm"), logical(1))]
    if (!length(p))    # name fallback (FIELD order, match Py/JS), EXCLUDING a scaled measure -- a scaled X is not lognorm
      p <- twod[vapply(twod, function(nm) (nm %in% .VIEWER_LOGNORM_NAMES) && !identical(st(nm), "scaled"), logical(1))]
    if (length(p)) p[1] else NULL
  }
  if (!is.null(counts)) {
    if (is.null(ds$fields[[counts]]))
      stop(sprintf("extend_for_viewer: counts='%s' is not a measure (present: %s)", counts, present), call. = FALSE)
    return(list(name = counts, log1p = !identical(st(counts), "lognorm")))
  }
  b <- if (is.null(basis)) "auto" else basis
  if (identical(b, "raw")) {
    p <- raw_pick()
    if (is.null(p)) stop(sprintf("extend_for_viewer: basis='raw' but no raw counts measure found (present cells x genes measures: %s).", present), call. = FALSE)
    return(list(name = p, log1p = TRUE))
  }
  if (identical(b, "lognorm")) {
    p <- lognorm_pick()
    if (is.null(p)) stop(sprintf("extend_for_viewer: basis='lognorm' but no log-normalized measure found (present cells x genes measures: %s).", present), call. = FALSE)
    return(list(name = p, log1p = FALSE))
  }
  if (identical(b, "auto")) {
    p <- raw_pick();     if (!is.null(p)) return(list(name = p, log1p = TRUE))
    p <- lognorm_pick(); if (!is.null(p)) return(list(name = p, log1p = FALSE))
    stop(sprintf(paste0("extend_for_viewer: no raw or log-normalized measure found (present cells x genes measures: %s). ",
                        "Viewer prep needs raw counts or log-normalized values; pass counts=<field> to force a measure. ",
                        "(A scaled/z-scored measure cannot be used as a basis.)"), present), call. = FALSE)
  }
  stop(sprintf("extend_for_viewer: basis must be 'auto', 'raw', or 'lognorm' (or pass counts=<field>); got '%s'", b), call. = FALSE)
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
