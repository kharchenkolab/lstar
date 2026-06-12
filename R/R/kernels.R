# Public compute kernels — the shared libstar primitives, exposed so downstream tools (e.g. the
# pagoda3 viewer's store prep) build on lstar's fast path rather than reimplementing reductions.

#' Per-(group, gene) sufficient stats over a CSC measure (the shared libstar kernel).
#'
#' For a `cells x genes` sparse matrix and a per-cell group assignment, returns each group's
#' sum, sum-of-squares, and number of expressing cells per gene, computed over `log1p(m)` (or raw).
#' This is the reduction cluster stats and marker tables are built from; the same C++ core backs
#' the WASM and Python bindings.
#'
#' @param m a `cells x genes` sparse matrix (coerced to CSC).
#' @param code length-`nrow(m)` integer group assignment in `[0, ngroups)` (0-based).
#' @param ngroups number of groups.
#' @param lognorm compute over `log1p(m)` (default `TRUE`).
#' @return a list with `sum`, `sumsq`, `n_expr`, each a flat row-major `(group, gene)` numeric vector
#'   of length `ngroups * ncol(m)`.
#' @export
col_sum_by_group <- function(m, code, ngroups, lognorm = TRUE) {
  m <- methods::as(m, "CsparseMatrix")
  lstar_cpp_col_sum_by_group(as.double(m@x), m@p, m@i, nrow(m), ncol(m),
                             as.integer(code), as.integer(ngroups), isTRUE(lognorm))
}

#' Per-gene mean/variance of a CSC measure in a store, read with bounded memory.
#'
#' Computes the zero-aware per-gene mean, variance, and number of expressing cells of a `cells x
#' genes` measure directly from an L\* store, reading it **block-by-block** so the whole matrix never
#' lands in memory -- the C++/R counterpart of Python's `stream_col_stats`. Bounded memory requires a
#' chunked store (one written with `chunk_elems` set, e.g. by a streamed conversion); on an unchunked
#' store it still works but reads the `data` array whole. The field must be CSC (gene-major); use it
#' for HVG selection / variance modeling over an atlas too large to load.
#'
#' @param path path to an L\* store (`.lstar.zarr`).
#' @param field measure field name (e.g. `"counts"`), stored CSC.
#' @param block number of gene columns per streamed block (default 4096).
#' @param n_threads threading policy for each block's reduction: 1 = serial (default), N = N threads,
#'   `<=0` = all cores. Results are thread-count invariant.
#' @param lognorm reduce over `log1p(x)` per nonzero, on the fly, without building the normalized
#'   matrix (default `FALSE`).
#' @return a list with `mean`, `var` (length `ngenes` numeric) and `nnz` (length `ngenes` integer).
#' @export
stream_col_stats <- function(path, field, block = 4096L, n_threads = 1L, lognorm = FALSE) {
  lstar_cpp_stream_col_stats(path, field, as.integer(block), as.integer(n_threads), isTRUE(lognorm))
}

#' Read a contiguous gene (column) range of a CSC measure from an L* store, bounded-memory.
#'
#' Reads genes `[g_lo, g_hi)` (0-based, half-open) of a CSC measure as a `dgCMatrix` (cells x genes),
#' decoding only the store chunks that overlap the range. The general block-read primitive a consumer
#' drives to build out-of-core reductions over an L* store without implementing them in lstar.
#' @export
lstar_read_block <- function(path, field, g_lo, g_hi, cell_names = NULL, gene_names = NULL) {
  b <- lstar_cpp_read_csc_block(path, field, as.integer(g_lo), as.integer(g_hi))
  m <- new("dgCMatrix", i = as.integer(b$indices), p = as.integer(b$indptr),
           x = as.numeric(b$data), Dim = as.integer(c(b$nrows, b$ncols)))
  if (!is.null(cell_names)) rownames(m) <- cell_names
  if (!is.null(gene_names)) colnames(m) <- gene_names[(g_lo + 1L):g_hi]
  m
}

#' Read an arbitrary set of gene columns of a CSC measure, returning cells x genes.
#'
#' Gathers the requested gene columns from a chunked CSC store, decoding each touched chunk **at most
#' once** (an ascending sweep over sorted-unique columns), then restores the caller's order. Efficient
#' for scattered subsets (e.g. overdispersed genes for PCA) -- unlike a per-column read it does not
#' re-decode a chunk once per gene it contains.
#' @export
lstar_read_genes <- function(path, field, genes, all_genes, cell_names = NULL) {
  idx <- if (is.character(genes)) match(genes, all_genes) else as.integer(genes)
  if (anyNA(idx)) stop("some requested genes are not in the store")
  u <- sort(unique(idx))                                    # 1-based, sorted, unique
  b <- lstar_cpp_read_csc_cols(path, field, as.integer(u - 1L))      # decode each chunk once
  m <- new("dgCMatrix", i = as.integer(b$indices), p = as.integer(b$indptr),
           x = as.numeric(b$data), Dim = as.integer(c(b$nrows, length(u))))
  m <- m[, match(idx, u), drop = FALSE]                     # restore caller order (and any duplicates)
  if (!is.null(cell_names)) rownames(m) <- cell_names
  colnames(m) <- if (is.character(genes)) genes else all_genes[genes]
  m
}
