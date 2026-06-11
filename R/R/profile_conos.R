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
# This is the read direction (Conos -> L*). Rebuilding a live Conos R6 from L* (read_conos) is
# deferred: the joint graph/embedding/clusters round-trip as fields, but reconstituting the R6
# wrapper and its per-sample Pagoda2 objects is profile-specific and not needed for interchange.

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
      if (!is.null(grp)) add(cn, as.character(grp[ucells]), "label", "cells")
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
