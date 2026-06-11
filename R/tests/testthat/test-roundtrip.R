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

test_that("print returns the dataset invisibly", {
  ds <- structure(list(kind = "sample", axes = list(), fields = list()),
                  class = "lstar_dataset")
  expect_output(print(ds), "lstar_dataset")
  expect_invisible(print(ds))
})
