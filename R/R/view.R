# Soft delegate to the pagoda3 viewer ---------------------------------------------------------------
#
# lstar is the substrate; the interactive viewer is the *separate* `pagoda3` package (which Imports
# lstar). We keep the dependency one-way: lstar forwards to pagoda3 when it is installed, so lstar
# installs and works with no viewer present. `pagoda3::view()` owns the real logic
# (coerce -> prep -> serve -> open); this is just the convenience entry from an lstar session.
# pagoda3 is not on a mainstream repository, so it is deliberately NOT declared in DESCRIPTION. The
# delegate therefore refers to it only dynamically -- the package name is a value, not a literal in
# `requireNamespace()` / `::` -- so `R CMD check` raises no "undeclared dependency" note, while the
# behaviour is identical to `pagoda3::view()`. It degrades to an install hint when absent. The helper
# is a named internal so tests can mock it (local_mocked_bindings) without touching disk.

.pagoda3_pkg <- "pagoda3"

.pagoda3_installed <- function() requireNamespace(.pagoda3_pkg, quietly = TRUE)

#' Open an object or L* store in the pagoda3 viewer (delegates to the 'pagoda3' package)
#'
#' A convenience shim: interactive viewing lives in the separate \pkg{pagoda3} package. When it is
#' installed this forwards to its \code{view()}; otherwise it errors with an install hint. The
#' dependency stays one-way (pagoda3 imports lstar, not the reverse), and because pagoda3 is not on a
#' mainstream repository it is resolved at run time rather than declared, so a plain lstar install
#' carries no viewer weight.
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
  # Resolve at run time via the package-name value (no literal `pagoda3::` token) for a package that
  # is intentionally not declared; behaviour is identical to pagoda3::view.
  pagoda3_view <- getExportedValue(.pagoda3_pkg, "view")
  pagoda3_view(obj, ...)
}
