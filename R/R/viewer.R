# Viewer extension hook (R) -------------------------------------------------------------------------
#
# Adds the lstar-viewer precomputed fields (counts_cellmajor, per-group stats + marker tables,
# od_score, and a hybrid cell order) to a `.lstar.zarr` store so the lstar-viewer web app can browse
# it with fast differential-expression / variable-gene / dotplot views and latency-cheap reads.
#
# This is a THIN WRAPPER over the Python implementation (`lstar.extend_for_viewer`, exposed on the CLI
# as `lstar convert <store> <store> --viewer`): native-R parity is a deliberate follow-up. The R side
# writes the store (or accepts an existing one) and shells to the Python CLI to extend it in place.

#' Extend an L* store with the lstar-viewer precomputed fields.
#'
#' Shells out to the lstar Python CLI (`lstar convert <store> <store> --viewer`) to add the viewer's
#' precomputed fields to an existing `.lstar.zarr` store: a cell-major (CSR) `counts_cellmajor`, per
#' categorical grouping its `stats_<g>_{sum,sumsq,nexpr}` sufficient statistics and `markers_<g>_{lfc,padj}`
#' 1-vs-rest marker tables, a global per-gene `od_score`, and (for the default hybrid order) a
#' `counts_cellmajor_order` permutation that makes cluster/lasso selections coalesce into a few
#' byte-range reads. The store is extended IN PLACE.
#'
#' Native-R parity is a follow-up; this wrapper requires a working `lstar` Python install on the PATH
#' (or via the `LSTAR_PYTHON` environment variable / the `python` argument).
#'
#' @param path path to an existing `*.lstar.zarr` store (must already contain `counts`, an embedding,
#'   and at least one categorical cell label).
#' @param python the Python interpreter to invoke (default: `$LSTAR_PYTHON` or `"python3"`).
#' @param check if `TRUE`, run the CLI's native-acceptance check on the result (default `FALSE`).
#' @return the store `path`, invisibly.
#' @examples
#' \dontrun{
#'   p <- "pbmc.lstar.zarr"
#'   lstar_write(ds, p)
#'   viewer_extend(p)            # adds counts_cellmajor, stats_*, markers_*, od_score, order
#' }
#' @seealso [lstar_write()]
#' @export
viewer_extend <- function(path, python = Sys.getenv("LSTAR_PYTHON", "python3"), check = FALSE) {
  path <- path.expand(path)
  if (!dir.exists(path))
    stop(sprintf("viewer_extend: store not found: %s", path), call. = FALSE)
  args <- c("-m", "lstar.cli", "convert", path, path, "--viewer")
  if (!isTRUE(check)) args <- c(args, "--no-check")
  status <- suppressWarnings(system2(python, args, stdout = TRUE, stderr = TRUE))
  rc <- attr(status, "status")
  if (!is.null(rc) && rc != 0L)
    stop(sprintf("viewer_extend: `%s %s` failed (exit %s):\n%s",
                 python, paste(args, collapse = " "), rc, paste(status, collapse = "\n")),
         call. = FALSE)
  invisible(path)
}

#' Write an R dataset to an L* store, optionally extending it for the viewer.
#'
#' A convenience wrapper that calls [lstar_write()] and then, when `viewer = TRUE`, [viewer_extend()]
#' to add the lstar-viewer precomputed fields. Equivalent to `lstar_write(ds, path); viewer_extend(path)`.
#'
#' @param ds an `lstar_dataset` (as returned by [lstar_read()] or a profile reader).
#' @param path output store path (a `*.lstar.zarr` directory).
#' @param viewer if `TRUE` (default), extend the written store via [viewer_extend()].
#' @param python the Python interpreter for the viewer extension (passed to [viewer_extend()]).
#' @param ... further arguments forwarded to [lstar_write()] (e.g. `chunk_elems`, `compression`).
#' @return the output `path`, invisibly.
#' @seealso [lstar_write()], [viewer_extend()]
#' @export
lstar_write_viewer <- function(ds, path, viewer = TRUE, python = Sys.getenv("LSTAR_PYTHON", "python3"), ...) {
  lstar_write(ds, path, ...)
  if (isTRUE(viewer)) viewer_extend(path, python = python)
  invisible(path)
}
