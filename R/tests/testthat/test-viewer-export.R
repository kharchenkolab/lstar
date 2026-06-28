# A non-viewer converter (write_seurat/write_sce) drops the viewer@0.1 `cache` navigators (regenerable;
# records them in `dropped`) rather than carrying a redundant/mis-aligned copy. Primaries are kept.
test_that("write_seurat / .lstar_drop_cache drop viewer@0.1 cache navigators, keep primaries", {
  skip_if_not_installed("Matrix")
  set.seed(1); nc <- 40L; ng <- 15L
  cnt <- as(Matrix::Matrix(rpois(nc * ng, 0.8), nc, ng, sparse = TRUE), "CsparseMatrix")
  rownames(cnt) <- paste0("c", 1:nc); colnames(cnt) <- paste0("g", 1:ng)
  ds <- structure(list(kind = "sample", spec_version = "0.1", profiles = character(0), dropped = character(0),
    axes = list(cells = list(labels = rownames(cnt), origin = "observed", role = "observation"),
                genes = list(labels = colnames(cnt), origin = "observed", role = "feature"),
                umap = list(labels = c("u1", "u2"), origin = "derived", role = "coordinate")),
    fields = list(counts = list(values = cnt, role = "measure", span = c("cells", "genes"), state = "raw", encoding = "csc"),
                  umap = list(values = matrix(rnorm(nc * 2), nc, 2), role = "embedding", span = c("cells", "umap")),
                  leiden = list(values = factor(paste0("k", (0:(nc - 1)) %% 3)), role = "label", span = "cells", encoding = "categorical"))),
    class = "lstar_dataset")
  ds <- extend_for_viewer(ds, grouping = "leiden")

  # the shared helper: caches recorded in `dropped` + removed; primaries (counts/leiden/umap) kept
  d2 <- lstar:::.lstar_drop_cache(ds)
  expect_true(all(c("counts_cellmajor", "counts_cellmajor_order", "od_score",
                    "stats_leiden_sum", "markers_leiden_lfc") %in% d2$dropped))
  expect_false(any(grepl("cellmajor|^stats_|^markers_|^od_score$", names(d2$fields))))
  expect_true(all(c("counts", "leiden", "umap") %in% names(d2$fields)))

  # the full converter doesn't carry navigators into the assay; counts/clusters/embedding survive
  skip_if_not_installed("SeuratObject")
  so <- write_seurat(ds)
  expect_false(any(grepl("cellmajor|stats_|markers_|od_score", SeuratObject::Layers(so))))
  expect_true("RNA" %in% SeuratObject::Assays(so))
  expect_true("umap" %in% SeuratObject::Reductions(so))
  expect_true("leiden" %in% colnames(so[[]]))
})
