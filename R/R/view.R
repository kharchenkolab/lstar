# Soft delegate to the pagoda3 viewer ---------------------------------------------------------------
#
# lstar is the substrate; the interactive viewer is the *separate* `pagoda3` package (which Imports
# lstar). We keep the dependency one-way: lstar only SUGGESTS pagoda3 and forwards to it when it is
# installed, so lstar installs and works with no viewer present. `pagoda3::view()` owns the real
# logic (coerce -> prep -> serve -> open); this is just the convenience entry from an lstar session.
# The helper is a named internal so tests can mock it (local_mocked_bindings) without touching disk.

.pagoda3_installed <- function() requireNamespace("pagoda3", quietly = TRUE)

#' Open an object or L* store in the pagoda3 viewer (delegates to the 'pagoda3' package)
#'
#' A convenience shim: interactive viewing lives in the separate \pkg{pagoda3} package. When it is
#' installed this forwards to \code{pagoda3::view()}; otherwise it errors with an install hint.
#' lstar only \emph{Suggests} pagoda3, so the dependency stays one-way (pagoda3 imports lstar, not
#' the reverse) and a plain lstar install carries no viewer weight.
#'
#' @param obj a \code{*.lstar.zarr} store path, an \code{lstar_dataset}, or a Seurat/SCE object.
#' @param ... passed through to \code{pagoda3::view()} (e.g. \code{prepare}, \code{local},
#'   \code{host}, \code{port}, \code{open}).
#' @return the viewer URL (invisibly), from \code{pagoda3::view()}.
#' @seealso \code{\link{extend_for_viewer}} (the prep this shares), the \pkg{pagoda3} package.
#' @export
view <- function(obj, ...) {
  if (!.pagoda3_installed())
    stop("Interactive viewing needs the 'pagoda3' package:\n  install.packages(\"pagoda3\")",
         call. = FALSE)
  pagoda3::view(obj, ...)
}
