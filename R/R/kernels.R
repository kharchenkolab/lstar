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
