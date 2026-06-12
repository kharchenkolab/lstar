#' @useDynLib lstar, .registration = TRUE
#' @importFrom Matrix sparseMatrix t
NULL

`%||%` <- function(a, b) if (is.null(a)) b else a

#' Read an L* Zarr store into an R dataset.
#'
#' @param path path to a `*.lstar.zarr` store
#' @return an `lstar_dataset`: a list with `axes` and `fields`, each field's `values`
#'   assembled as a base vector, matrix, or `Matrix` sparse matrix.
#' @examples
#' p <- tempfile(fileext = ".lstar.zarr")
#' ds <- list(kind = "sample", axes = list(), fields = list())
#' ds$axes$cells <- list(labels = paste0("c", 1:3), origin = "observed", role = "observation")
#' ds$fields$depth <- list(role = "measure", span = "cells", values = c(1, 2, 3))
#' class(ds) <- "lstar_dataset"
#' lstar_write(ds, p)
#' ds2 <- lstar_read(p)
#' field_value(ds2, "depth")
#' @export
lstar_read <- function(path) {
  ds <- lstar_cpp_read(path.expand(path))
  for (nm in names(ds$fields)) {
    ds$fields[[nm]]$values <- .lstar_assemble(ds$fields[[nm]])
  }
  class(ds) <- "lstar_dataset"
  ds
}

.lstar_assemble <- function(f) {
  enc <- f$encoding
  v <- if (enc %in% c("csc", "csr")) {
    dims <- as.integer(f$shape)
    if (enc == "csc") {
      Matrix::sparseMatrix(i = f$indices, p = f$indptr, x = f$data,
                           dims = dims, index1 = FALSE, repr = "C")
    } else {
      Matrix::sparseMatrix(j = f$indices, p = f$indptr, x = f$data,
                           dims = dims, index1 = FALSE, repr = "R")
    }
  } else if (enc == "utf8") {
    f$strings
  } else if (enc == "categorical") {
    codes <- as.integer(f$codes)                 # 0-based; -1 = missing
    iv <- codes + 1L
    iv[codes < 0] <- NA
    structure(iv, levels = as.character(f$categories),
              class = if (isTRUE(f$ordered)) c("ordered", "factor") else "factor")
  } else {
    shp <- as.integer(f$shape)
    if (length(shp) <= 1) f$dense else t(matrix(f$dense, nrow = shp[2], ncol = shp[1]))
  }
  if (!is.null(f$mask) && is.null(dim(v))) {       # nullable: 1 == missing -> NA in the R vector
    miss <- as.integer(f$mask) == 1L
    if (any(miss)) v[miss] <- NA
  }
  v
}

.infer_encoding <- function(v) {
  if (is.factor(v)) "categorical"
  else if (is.character(v)) "utf8"
  else if (methods::is(v, "RsparseMatrix")) "csr"
  else if (methods::is(v, "sparseMatrix")) "csc"
  else "dense"
}

#' Write an R dataset to an L* Zarr store.
#'
#' @param ds an `lstar_dataset` (as returned by [lstar_read()] or a profile reader)
#' @param path output store path (a `*.lstar.zarr` directory)
#' @param chunk_elems if non-NULL, chunk each array along its first axis so each chunk holds about
#'   this many elements (e.g. `1e6`). This is what lets a reader stream/block-read only the touched
#'   chunks (e.g. [lstar_read_block()], `stream_col_stats()`); the default (NULL) writes each array as
#'   a single chunk -- the portable, byte-identical-to-before default.
#' @param compression chunk codec: `"none"` (default), `"gzip"`, or `"zlib"` (numcodecs-compatible;
#'   readable by the C++ core and zarr-python).
#' @param level compression level 1-9 (default 5), used when `compression` is `"gzip"`/`"zlib"`.
#' @return the output `path`, invisibly.
#' @seealso [lstar_read()], [lstar_read_block()]
#' @export
lstar_write <- function(ds, path, chunk_elems = NULL, compression = c("none", "gzip", "zlib"),
                        level = 5L) {
  compression <- match.arg(compression)
  axes <- lapply(names(ds$axes), function(nm) {
    a <- ds$axes[[nm]]
    list(labels = as.character(a$labels), origin = a$origin %||% "observed",
         role = a$role %||% "", induced_by = a$induced_by %||% "")
  })
  names(axes) <- names(ds$axes)

  fields <- lapply(names(ds$fields), function(nm) {
    f <- ds$fields[[nm]]
    v <- f$values
    enc <- f$encoding %||% .infer_encoding(v)
    out <- list(role = f$role %||% "", span = as.character(f$span), encoding = enc,
                state = f$state %||% "", subtype = f$subtype %||% "")
    if (enc %in% c("csc", "csr")) {
      m <- if (enc == "csc") as(v, "CsparseMatrix") else as(v, "RsparseMatrix")
      out$data <- as.numeric(m@x)
      out$indices <- if (enc == "csc") as.numeric(m@i) else as.numeric(m@j)
      out$indptr <- as.numeric(m@p)
      out$shape <- as.integer(dim(m))
    } else if (enc == "utf8") {
      sv <- as.character(v)
      if (any(is.na(sv))) { out$mask <- as.integer(is.na(sv)); sv[is.na(sv)] <- "" }  # nullable string
      out$strings <- sv
    } else if (enc == "categorical") {
      v <- if (is.factor(v)) v else factor(v)
      code <- as.integer(v) - 1L                   # R factor: 1-based, NA -> -1, 0-based
      code[is.na(code)] <- -1L
      out$codes <- as.integer(code)
      out$categories <- as.character(levels(v))
      out$ordered <- is.ordered(v)
    } else {
      if (is.null(dim(v))) {                       # a plain vector -> arity-1 field
        if (any(is.na(v))) { out$mask <- as.integer(is.na(v)); v[is.na(v)] <- 0 }  # nullable int/num
        out$dense <- as.numeric(v); out$shape <- as.integer(length(v))
      } else {                                     # a matrix -> C-order + 2-D shape
        m <- as.matrix(v)
        out$dense <- as.numeric(t(m)); out$shape <- as.integer(dim(m))
      }
    }
    out
  })
  names(fields) <- names(ds$fields)

  payload <- list(kind = ds$kind %||% "sample", spec_version = ds$spec_version %||% "0.1",
                  profiles = as.character(ds$profiles %||% character(0)),
                  dropped = as.character(ds$dropped %||% character(0)),
                  axes = axes, fields = fields, aux = ds$aux %||% list())  # passthrough, round-tripped verbatim
  lstar_cpp_write(payload, path.expand(path),
                  as.integer(if (is.null(chunk_elems)) 0L else chunk_elems),
                  if (compression == "none") "" else compression, as.integer(level))
  invisible(path)
}

#' Print an L* dataset
#'
#' @param x an `lstar_dataset`
#' @param ... ignored, for S3 compatibility
#' @return `x`, invisibly (called for the side effect of printing a summary of axes and fields).
#' @export
print.lstar_dataset <- function(x, ...) {
  cat(sprintf("lstar_dataset (%s): %d axes, %d fields\n",
              x$kind %||% "?", length(x$axes), length(x$fields)))
  for (nm in names(x$axes)) cat(sprintf("  axis  %-10s %d\n", nm, length(x$axes[[nm]]$labels)))
  for (nm in names(x$fields)) {
    f <- x$fields[[nm]]
    cat(sprintf("  field %-14s %-10s [%s]\n", nm, f$role %||% "", paste(f$span, collapse = " x ")))
  }
  invisible(x)
}

#' Accessor: a field's value by name.
#' @param ds an `lstar_dataset`
#' @param name field name
#' @return the field's `values` (a vector, matrix, or sparse `Matrix`), or `NULL` if absent.
#' @export
field_value <- function(ds, name) ds$fields[[name]]$values
