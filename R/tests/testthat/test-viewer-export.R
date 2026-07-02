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

# extend_for_viewer selects the count basis by STATE, not the literal name "counts": a raw measure
# named otherwise is auto-picked; no raw basis gives a clear error; basis="lognorm" preps from a log
# measure. (Regression for the anndata-class bug where converters named the raw matrix `X`/modality.)
test_that("extend_for_viewer selects basis by state (not the name 'counts')", {
  skip_if_not_installed("Matrix")
  set.seed(2); nc <- 60L; ng <- 20L
  mk <- function(state) {
    cnt <- as(Matrix::Matrix(rpois(nc * ng, 1.2) + 1L, nc, ng, sparse = TRUE), "CsparseMatrix")
    rownames(cnt) <- paste0("c", 1:nc); colnames(cnt) <- paste0("g", 1:ng)
    structure(list(kind = "sample", spec_version = "0.1", profiles = character(0), dropped = character(0),
      axes = list(cells = list(labels = rownames(cnt), origin = "observed", role = "observation"),
                  genes = list(labels = colnames(cnt), origin = "observed", role = "feature")),
      fields = list(X = list(values = cnt, role = "measure", span = c("cells", "genes"), state = state),
                    leiden = list(values = factor(paste0("k", (0:(nc - 1)) %% 3)), role = "label", span = "cells", encoding = "categorical"))),
      class = "lstar_dataset")
  }
  # raw counts under the name "X" -> auto-picked by state
  d <- extend_for_viewer(mk("raw"), grouping = "leiden")
  expect_false(is.null(d$fields[["od_score"]]))
  expect_false(is.null(d$fields[["counts_cellmajor"]]))
  # no raw basis -> clear, actionable error (not a bare "no counts measure")
  expect_error(extend_for_viewer(mk("scaled"), grouping = "leiden"), "no raw counts measure")
  # opt-in lognorm basis
  d3 <- extend_for_viewer(mk("lognorm"), grouping = "leiden", basis = "lognorm")
  expect_equal(d3$fields[["counts_cellmajor"]]$state, "lognorm")
})
