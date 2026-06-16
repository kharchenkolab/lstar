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


.mini_collection_ds <- function() {
  ax <- list(); fl <- list(); joint_cells <- character(0); joint_sample <- character(0)
  for (s in c("A", "B")) {
    set.seed(match(s, c("A", "B")))
    cells <- paste0(s, "_c", 1:8); genes <- paste0("g", 1:20)
    m <- as(Matrix::Matrix(matrix(rpois(8 * 20, 2) + 1L, 8, 20, dimnames = list(cells, genes)), sparse = TRUE), "CsparseMatrix")
    pca <- matrix(rnorm(8 * 5), 8, 5, dimnames = list(cells, paste0("PC", 1:5)))
    ax[[paste0("cells.", s)]] <- list(labels = cells, origin = "observed", role = "observation")
    ax[[paste0("genes.", s)]] <- list(labels = genes, origin = "observed", role = "feature")
    ax[[paste0("pca.", s)]]   <- list(labels = colnames(pca), origin = "derived", role = "coordinate")
    fl[[paste0("counts.", s)]] <- list(role = "measure", span = c(paste0("cells.", s), paste0("genes.", s)), state = "raw", values = m)
    fl[[paste0("pca.", s)]]    <- list(role = "embedding", span = c(paste0("cells.", s), paste0("pca.", s)), values = unname(pca))
    joint_cells <- c(joint_cells, cells); joint_sample <- c(joint_sample, rep(s, length(cells)))
  }
  n <- length(joint_cells)
  A <- as(Matrix::Matrix(0, n, n, dimnames = list(joint_cells, joint_cells)), "CsparseMatrix"); A[1, 2] <- A[2, 1] <- 1
  ax$samples   <- list(labels = c("A", "B"), origin = "observed", role = "sample")
  ax$cells     <- list(labels = joint_cells, origin = "derived", role = "observation")
  ax$embedding <- list(labels = c("E1", "E2"), origin = "derived", role = "coordinate")
  ax$leiden    <- list(labels = c("1", "2"), origin = "derived", role = "factor", induced_by = "leiden")
  fl$sample    <- list(role = "label", span = "cells", subtype = "design", values = joint_sample)
  fl$embedding <- list(role = "embedding", span = c("cells", "embedding"), values = matrix(rnorm(n * 2), n, 2))
  fl$leiden    <- list(role = "label", span = "cells", encoding = "categorical", values = factor(rep(c("1", "2"), length.out = n)))
  fl$graph     <- list(role = "relation", span = c("cells", "cells"), subtype = "knn", values = A)
  structure(list(kind = "collection", spec_version = "0.1", profiles = character(0), dropped = character(0),
                 axes = ax, fields = fl), class = "lstar_dataset")
}

test_that("Conos profile reconstructs a collection (read_conos)", {
  skip_if_not_installed("conos")
  skip_if_not_installed("pagoda2")
  skip_if_not_installed("igraph")
  con <- read_conos(.mini_collection_ds())
  expect_s3_class(con, "Conos")
  expect_equal(length(con$samples), 2L)
  expect_true(all(vapply(con$samples, function(s) inherits(s, "Pagoda2"), logical(1))))
  expect_equal(nrow(getPca(con$samples[["A"]])), 8L)        # per-sample PCA restored
  expect_false(is.null(con$graph))                          # joint graph restored
  expect_false(is.null(con$embedding))                      # joint embedding restored
  expect_equal(nlevels(con$clusters$leiden$groups), 2L)     # joint clustering restored
})
