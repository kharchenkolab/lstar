# Profiles (Seurat / SCE) + the package-free direct backends, exercised by the package's own suite (the
# heavy cross-language coverage lives in conformance/, which R CMD check does not run). Each test skips
# cleanly when its native package is absent.

.mini_ds <- function() {
  cells <- paste0("c", 1:6); genes <- paste0("g", 1:4)
  m <- as(matrix(as.numeric(1:24), 6, 4, dimnames = list(cells, genes)), "CsparseMatrix")  # cells x genes
  emb <- matrix(as.numeric(1:12), 6, 2, dimnames = list(cells, c("PC_1", "PC_2")))
  ds <- list(
    kind = "sample", spec_version = "0.1", profiles = character(0), dropped = character(0),
    axes = list(
      cells = list(labels = cells, origin = "observed", role = "observation"),
      genes = list(labels = genes, origin = "observed", role = "feature"),
      pca   = list(labels = c("PC_1", "PC_2"), origin = "derived", role = "coordinate")),
    fields = list(
      counts = list(role = "measure", span = c("cells", "genes"), state = "raw", values = m),
      leiden = list(role = "label", span = "cells", encoding = "categorical",
                    values = factor(c("a", "a", "b", "b", "a", "b"))),
      pca    = list(role = "embedding", span = c("cells", "pca"), values = emb)))
  class(ds) <- "lstar_dataset"
  ds
}

test_that("Seurat profile round-trips the core", {
  skip_if_not_installed("SeuratObject")
  ds <- .mini_ds()
  so <- write_seurat(ds)
  expect_true(methods::is(so, "Seurat"))
  ds2 <- read_seurat(so)
  expect_true(all(c("counts", "leiden", "pca") %in% names(ds2$fields)))
  expect_equal(as.numeric(sum(field_value(ds2, "counts"))), sum(ds$fields$counts$values))
  expect_equal(dim(field_value(ds2, "pca")), c(6L, 2L))
})

test_that("SingleCellExperiment profile round-trips the core", {
  skip_if_not_installed("SingleCellExperiment")
  ds <- .mini_ds()
  sce <- write_sce(ds)
  expect_true(methods::is(sce, "SingleCellExperiment"))
  ds2 <- read_sce(sce)
  expect_true("counts" %in% names(ds2$fields))
  expect_equal(as.numeric(sum(field_value(ds2, "counts"))), sum(ds$fields$counts$values))
})

test_that("package-free Seurat read (S4 slot-walk) matches the native read", {
  skip_if_not_installed("SeuratObject")
  so <- write_seurat(.mini_ds())
  direct <- lstar:::.read_seurat_direct(so)   # attr()-based, no SeuratObject accessors
  native <- read_seurat(so)                   # SeuratObject accessors
  expect_setequal(names(direct$fields), names(native$fields))
  expect_equal(as.numeric(sum(field_value(direct, "counts"))),
               as.numeric(sum(field_value(native, "counts"))))
})

test_that("package-free Seurat write builds a native-valid object", {
  skip_if_not_installed("SeuratObject")
  obj <- lstar:::.build_seurat_direct(.mini_ds())   # pinned-schema build (uses the real class when present)
  expect_true(methods::is(obj, "Seurat"))
  expect_silent(methods::validObject(obj))
  expect_true("counts" %in% names(read_seurat(obj)$fields))   # native reads the forged object
})

test_that("package-free SCE read (S4 slot-walk) matches the native read", {
  skip_if_not_installed("SingleCellExperiment")
  sce <- write_sce(.mini_ds())
  direct <- lstar:::.read_sce_direct(sce)
  native <- read_sce(sce)
  expect_setequal(names(direct$fields), names(native$fields))
  expect_equal(as.numeric(sum(field_value(direct, "counts"))),
               as.numeric(sum(field_value(native, "counts"))))
})
