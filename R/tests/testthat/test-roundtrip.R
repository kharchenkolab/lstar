test_that("a dataset round-trips through an L* store", {
  ds <- list(kind = "sample", spec_version = "0.1", axes = list(), fields = list())
  ds$axes$cells <- list(labels = paste0("c", 1:4), origin = "observed", role = "observation")
  ds$axes$genes <- list(labels = paste0("g", 1:3), origin = "observed", role = "feature")
  m <- as(matrix(c(0, 1, 0, 2, 0, 0, 3, 0, 0, 0, 4, 0), 4, 3), "CsparseMatrix")
  ds$fields$counts <- list(role = "measure", span = c("cells", "genes"), state = "raw", values = m)
  ds$fields$cluster <- list(role = "label", span = "cells", values = c("a", "a", "b", "b"))
  class(ds) <- "lstar_dataset"

  p <- tempfile(fileext = ".lstar.zarr")
  lstar_write(ds, p)
  ds2 <- lstar_read(p)

  expect_s3_class(ds2, "lstar_dataset")
  expect_setequal(names(ds2$fields), c("counts", "cluster"))
  expect_equal(as.numeric(sum(field_value(ds2, "counts"))), sum(m))
  expect_equal(dim(field_value(ds2, "counts")), dim(m))
  expect_equal(field_value(ds2, "cluster"), c("a", "a", "b", "b"))
})

test_that("only state=raw integer data narrows to i4; non-raw float layers stay f8", {
  skip_if_not_installed("Matrix")
  m <- as(Matrix::Matrix(c(0, 3, 0, 5, 2, 0, 0, 4, 1, 0, 0, 6), 4, 3, sparse = TRUE), "CsparseMatrix")
  ds <- structure(list(kind = "sample", spec_version = "0.1", profiles = character(0), dropped = character(0),
    axes = list(cells = list(labels = paste0("c", 1:4), origin = "observed", role = "observation"),
                genes = list(labels = paste0("g", 1:3), origin = "observed", role = "feature")),
    fields = list(
      counts  = list(role = "measure", span = c("cells", "genes"), state = "raw",     encoding = "csc", values = m),
      lognorm = list(role = "measure", span = c("cells", "genes"), state = "lognorm", encoding = "csc", values = m))),
    class = "lstar_dataset")
  p <- tempfile(fileext = ".lstar.zarr"); lstar_write(ds, p)
  zarray <- function(field) paste(readLines(file.path(p, "fields", field, "data", ".zarray"), warn = FALSE), collapse = "")
  expect_match(zarray("counts"), "<i4", fixed = TRUE)    # raw integer counts -> i4 (de-widened)
  expect_match(zarray("lognorm"), "<f8", fixed = TRUE)   # SAME integer values but state!=raw -> stays f8
  ds2 <- lstar_read(p)                                   # values round-trip regardless of dtype
  expect_equal(as.numeric(sum(field_value(ds2, "counts"))), sum(m))
  expect_equal(as.numeric(sum(field_value(ds2, "lognorm"))), sum(m))
})

test_that("print returns the dataset invisibly", {
  ds <- structure(list(kind = "sample", axes = list(), fields = list()),
                  class = "lstar_dataset")
  expect_output(print(ds), "lstar_dataset")
  expect_invisible(print(ds))
})
