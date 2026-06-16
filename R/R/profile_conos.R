# Conos profile: map a Conos R6 collection to an L* dataset.
#
# A Conos object is a *collection of samples*, not one aligned matrix, and L* represents it
# as such:
#   - a `samples` axis (the collection itself)
#   - per-sample `cells.<s>` and `genes.<s>` axes -- different sizes AND different gene sets
#     per sample; this heterogeneity is exactly what an aligned cells x genes tensor erases
#   - per-sample `counts.<s>` measures and `pca.<s>` embeddings, each over its own sample's axes
#   - a union `cells` axis for the joint analysis layer, with a `sample` label (getDatasetPerCell)
#   - the joint conos `embedding`, the joint clustering label(s), and the joint `graph` as a
#     `relation` over (cells x cells)
#
# write_conos() is the import direction (Conos -> L*); read_conos() is the inverse (L* -> Conos),
# reconstituting the per-sample Pagoda2 objects (raw counts + PCA reduction) and restoring the joint
# graph / embedding / clustering layer so a stored collection can be re-opened for plotting, markers and
# label transfer. (Re-running runGraph() recomputes the per-sample variance model, which is not stored.)

.conos_sample_counts <- function(p) {
  # Recognize the pagoda2 version gracefully: pagoda2.1 (devel) serves raw counts through the
  # getRawCounts() accessor (the `$counts` slot was removed); older pagoda2 stored `$counts`
  # directly. Try the accessor first, fall back to the slot. Result is a cells x genes
  # dgCMatrix (already L* orientation).
  rc <- tryCatch(p$getRawCounts(), error = function(e) NULL)
  if (is.null(rc)) rc <- tryCatch(p$counts, error = function(e) NULL)
  rc
}

.conos_versions <- function(co) {
  cv <- tryCatch(as.character(utils::packageVersion("conos")), error = function(e) "?")
  pv <- tryCatch(as.character(utils::packageVersion("pagoda2")), error = function(e) "?")
  # also note which per-sample counts API is live (accessor vs legacy slot)
  api <- if (length(co$samples) &&
             !is.null(tryCatch(co$samples[[1]]$getRawCounts, error = function(e) NULL)))
           "pagoda2-accessor" else "pagoda2-slot"
  c("conos@0.1", paste0("conos@", cv), paste0("pagoda2@", pv), api)
}

#' Build an L* dataset from a Conos object (a collection of samples).
#'
#' @param co a `Conos` R6 object
#' @param clustering optional name of a joint clustering in `co$clusters` (default: all)
#' @return an `lstar_dataset` of kind `collection`
#' @export
write_conos <- function(co, clustering = NULL) {
  if (!inherits(co, "Conos")) stop("write_conos expects a Conos object")
  `%||%` <- function(a, b) if (is.null(a)) b else a
  ds <- list(kind = "collection", spec_version = "0.1", profiles = .conos_versions(co),
             dropped = character(0), axes = list(), fields = list())
  add <- function(nm, values, role, span, state = "", subtype = "") {
    ds$fields[[nm]] <<- list(role = role, span = span, state = state, subtype = subtype,
                             values = values)
  }
  add_factor <- function(nm, v, span) {                   # a clustering factor -> categorical field +
    v <- droplevels(as.factor(v))                         # induced bare-named `factor` axis (its levels)
    ds$fields[[nm]] <<- list(role = "label", span = span, state = "", subtype = "",
                             encoding = "categorical", values = v)
    if (is.null(ds$axes[[nm]]))
      ds$axes[[nm]] <<- list(labels = levels(v), origin = "derived", role = "factor", induced_by = nm)
  }

  sample_names <- names(co$samples)
  n_per <- integer(length(sample_names)); names(n_per) <- sample_names

  # ---- per-sample axes + fields (the collection members) -------------------------------
  for (sn in sample_names) {
    p <- co$samples[[sn]]
    rc <- .conos_sample_counts(p)
    if (is.null(rc)) { ds$dropped <- c(ds$dropped, sprintf("samples/%s/counts", sn)); next }
    scells <- rownames(rc) %||% paste0(sn, "_", seq_len(nrow(rc)))
    sgenes <- colnames(rc) %||% paste0("g", seq_len(ncol(rc)))
    n_per[sn] <- nrow(rc)
    ca <- paste0("cells.", sn); ga <- paste0("genes.", sn)
    ds$axes[[ca]] <- list(labels = scells, origin = "observed", role = "observation")
    ds$axes[[ga]] <- list(labels = sgenes, origin = "observed", role = "feature")
    add(paste0("counts.", sn), as(rc, "CsparseMatrix"), "measure", c(ca, ga), state = "raw")

    pca <- tryCatch(p$reductions$PCA, error = function(e) NULL)
    if (!is.null(pca) && nrow(pca) == nrow(rc)) {
      pa <- paste0("pca.", sn)
      ds$axes[[pa]] <- list(labels = colnames(pca) %||% paste0("PC", seq_len(ncol(pca))),
                            origin = "derived", role = "coordinate")
      add(paste0("pca.", sn), unname(as.matrix(pca)), "embedding", c(ca, pa))
    }
  }

  # ---- the samples axis itself + a per-sample size measure -----------------------------
  ds$axes$samples <- list(labels = sample_names, origin = "observed", role = "sample")
  add("n_cells", as.numeric(n_per), "measure", "samples")

  # ---- the joint analysis layer over a union `cells` axis ------------------------------
  dpc <- tryCatch(co$getDatasetPerCell(), error = function(e) NULL)
  emb <- tryCatch(co$embedding, error = function(e) NULL)
  ucells <- if (!is.null(dpc)) names(dpc) else if (!is.null(emb)) rownames(emb) else NULL
  if (!is.null(emb)) ucells <- intersect(ucells, rownames(emb))
  g <- tryCatch(co$graph, error = function(e) NULL)
  if (!is.null(g) && requireNamespace("igraph", quietly = TRUE)) {
    vn <- igraph::V(g)$name
    if (!is.null(vn)) ucells <- intersect(ucells, vn)
  }

  if (!is.null(ucells) && length(ucells)) {
    ds$axes$cells <- list(labels = ucells, origin = "derived", role = "observation")
    if (!is.null(dpc)) add("sample", as.character(dpc[ucells]), "label", "cells", subtype = "design")

    if (!is.null(emb)) {
      e <- emb[ucells, , drop = FALSE]
      ds$axes$embedding <- list(labels = colnames(e) %||% paste0("E", seq_len(ncol(e))),
                                origin = "derived", role = "coordinate")
      add("embedding", unname(as.matrix(e)), "embedding", c("cells", "embedding"))
    }

    cl_names <- if (!is.null(clustering)) clustering else names(co$clusters)
    for (cn in cl_names) {
      grp <- tryCatch(co$clusters[[cn]]$groups, error = function(e) NULL)
      if (!is.null(grp)) add_factor(cn, grp[ucells], "cells")     # conos clusterings are factors -> factor axis
    }

    if (!is.null(g) && requireNamespace("igraph", quietly = TRUE)) {
      A <- tryCatch(igraph::as_adjacency_matrix(g, attr = "weight", sparse = TRUE),
                    error = function(e) igraph::as_adjacency_matrix(g, sparse = TRUE))
      A <- A[ucells, ucells]
      add("graph", as(A, "CsparseMatrix"), "relation", c("cells", "cells"), subtype = "knn")
    }
  }

  class(ds) <- "lstar_dataset"
  ds
}

#' Reconstruct a Conos object from an L* collection
#'
#' The inverse of \code{write_conos()}: rebuilds the per-sample \code{Pagoda2} objects (raw counts plus the
#' stored PCA reduction) and restores the joint graph, embedding and clustering(s), returning a live
#' \code{conos::Conos} object ready for plotting, marker detection and label transfer. Re-running
#' \code{runGraph()} will recompute the per-sample variance model (not stored).
#'
#' @param ds an \code{lstar_dataset} of kind \code{"collection"} (e.g. from \code{lstar_read()}), or a path
#'   to an lstar store holding one.
#' @return a \code{conos::Conos} object.
#' @export
read_conos <- function(ds) {
  if (is.character(ds) && length(ds) == 1) ds <- lstar_read(ds)
  if (!inherits(ds, "lstar_dataset") || !identical(ds$kind, "collection"))
    stop("read_conos expects an L* 'collection' dataset (or a path to one)")
  if (!requireNamespace("conos", quietly = TRUE) || !requireNamespace("pagoda2", quietly = TRUE))
    stop("read_conos requires the conos and pagoda2 packages")
  `%||%` <- function(a, b) if (is.null(a)) b else a
  axl <- function(nm) { a <- ds$axes[[nm]]; if (is.null(a)) NULL else as.character(a$labels) }
  val <- function(nm) ds$fields[[nm]]$values

  sample_names <- axl("samples")
  if (is.null(sample_names)) stop("collection has no `samples` axis")

  samples <- list()
  for (sn in sample_names) {
    cnm <- paste0("counts.", sn)
    if (is.null(ds$fields[[cnm]])) next
    counts <- methods::as(val(cnm), "CsparseMatrix")               # cells x genes
    rownames(counts) <- axl(paste0("cells.", sn))
    colnames(counts) <- axl(paste0("genes.", sn))
    p <- pagoda2::Pagoda2$new(Matrix::t(counts), min.transcripts.per.cell = 0, min.cells.per.gene = 0,
                              log.scale = TRUE, n.cores = 1, verbose = FALSE)   # genes x cells
    pnm <- paste0("pca.", sn)
    if (!is.null(ds$fields[[pnm]])) {
      pca <- as.matrix(val(pnm))
      rownames(pca) <- axl(paste0("cells.", sn))
      colnames(pca) <- axl(paste0("pca.", sn)) %||% paste0("PC", seq_len(ncol(pca)))
      p$reductions$PCA <- pca
    }
    samples[[sn]] <- p
  }
  if (length(samples) == 0) stop("collection has no per-sample counts to reconstruct")

  con <- conos::Conos$new(samples, n.cores = 1)

  cells <- axl("cells")
  if (!is.null(cells) && length(cells)) {
    if (!is.null(ds$fields$graph) && requireNamespace("igraph", quietly = TRUE)) {
      A <- methods::as(val("graph"), "CsparseMatrix"); dimnames(A) <- list(cells, cells)
      con$graph <- igraph::graph_from_adjacency_matrix(A, mode = "undirected", weighted = TRUE)
    }
    if (!is.null(ds$fields$embedding)) {
      emb <- as.matrix(val("embedding")); rownames(emb) <- cells
      colnames(emb) <- axl("embedding") %||% paste0("E", seq_len(ncol(emb)))
      con$embedding <- emb
    }
    ## clustering layers: categorical labels over the joint `cells` axis (the `sample` design label is not
    ## categorical, so it is skipped here)
    for (nm in names(ds$fields)) {
      f <- ds$fields[[nm]]
      if (identical(f$role, "label") && identical(f$encoding %||% "", "categorical") &&
          identical(as.character(f$span), "cells")) {
        grp <- f$values
        lv <- axl(nm)
        grp <- if (!is.null(lv)) factor(as.character(grp), levels = lv) else as.factor(grp)
        names(grp) <- cells
        con$clusters[[nm]] <- list(groups = grp)
      }
    }
  }
  con
}
