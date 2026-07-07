# `primary=` (the grouping the viewer opens on) hoists that grouping to the front: it keys the
# counts_cellmajor locality reorder + is summarized first, and COMPOSES with auto-detect (the other
# groupings are still prepped). Mirrors python/tests/test_viewer.py::test_primary_* and js/test/extend_primary.
test_that("extend_for_viewer primary= hoists the reorder-key grouping and composes with auto-detect", {
  skip_if_not_installed("Matrix")
  set.seed(1); nc <- 30L; ng <- 8L
  cnt <- as(Matrix::rsparsematrix(nc, ng, density = 0.4, rand.x = function(n) rpois(n, 3) + 1), "CsparseMatrix")
  rownames(cnt) <- paste0("c", 0:(nc - 1)); colnames(cnt) <- paste0("g", 0:(ng - 1))
  mkds <- function() structure(list(kind = "sample", spec_version = "0.1", profiles = character(0), dropped = character(0),
    axes = list(cells = list(labels = rownames(cnt), origin = "observed", role = "observation"),
                genes = list(labels = colnames(cnt), origin = "observed", role = "feature"),
                umap = list(labels = c("u1", "u2"), origin = "derived", role = "coordinate")),
    fields = list(counts = list(values = cnt, role = "measure", span = c("cells", "genes"), state = "raw", encoding = "csc"),
                  umap = list(values = matrix(rnorm(nc * 2), nc, 2), role = "embedding", span = c("cells", "umap")),
                  leiden = list(values = factor(paste0("k", (0:(nc - 1)) %% 3)), role = "label", span = "cells", encoding = "categorical"),
                  cell_type = list(values = factor(c("T", "B", "NK")[(0:(nc - 1)) %% 3 + 1]), role = "label", span = "cells", encoding = "categorical"))),
    class = "lstar_dataset")
  # primary hoists cell_type but still preps the auto-detected leiden; the reorder is keyed on cell_type.
  d1 <- extend_for_viewer(mkds(), primary = "cell_type")
  expect_true(all(c("stats_cell_type_sum", "stats_leiden_sum", "markers_cell_type_lfc") %in% names(d1$fields)))
  expect_identical(d1$fields[["counts_cellmajor_order"]]$provenance$group, "cell_type")
  # default: detection prefers leiden over cell_type, so the default reorder key is leiden.
  d2 <- extend_for_viewer(mkds())
  expect_identical(d2$fields[["counts_cellmajor_order"]]$provenance$group, "leiden")
  expect_error(extend_for_viewer(mkds(), primary = "not_a_field"), "primary")
  expect_error(extend_for_viewer(mkds(), primary = "umap"), "cell axis")   # a non-grouping field (2-D embedding)
})

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
# named otherwise is auto-picked; with no raw, basis="auto" FALLS BACK to a log-normalized measure (with
# a warning); a scaled-only store gives a clear error (scaled is never a basis). Identical contract to
# Python/JS. (Regression for the anndata-class bug where converters named the raw matrix `X`/modality.)
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
  # a scaled-only measure cannot be a viewer basis -> clear error; auto NEVER falls into it
  expect_error(extend_for_viewer(mk("scaled"), grouping = "leiden"), "no raw or log-normalized measure")
  # no raw, but a lognorm measure present -> auto FALLS BACK to lognorm (with a warning), does not error
  expect_warning(d2 <- extend_for_viewer(mk("lognorm"), grouping = "leiden"), "log-normalized")
  expect_equal(d2$fields[["counts_cellmajor"]]$state, "lognorm")
  # opt-in explicit lognorm basis (no warning)
  d3 <- extend_for_viewer(mk("lognorm"), grouping = "leiden", basis = "lognorm")
  expect_equal(d3$fields[["counts_cellmajor"]]$state, "lognorm")
})

test_that(".viewer_counts_basis: a measure NAMED 'counts' that is scaled is not picked as raw", {
  # regression (the name-shortcut hole): the literal-"counts" fast path must exclude a scaled measure,
  # symmetric with the lognorm name fallback -- else a scaled "counts" would be picked as raw and log1p'd.
  mkds <- function(fl) structure(list(fields = fl), class = "lstar_dataset")
  mkf <- function(state) list(role = "measure", span = c("cells", "genes"), state = state)
  # scaled "counts" + a real raw measure -> skip the scaled "counts", pick the raw one (log1p)
  b <- lstar:::.viewer_counts_basis(mkds(list(counts = mkf("scaled"), X = mkf("raw"))))
  expect_equal(b$name, "X"); expect_true(b$log1p)
  # only a scaled "counts" -> clear error (a scaled measure is never a basis), not a silent log1p
  expect_error(lstar:::.viewer_counts_basis(mkds(list(counts = mkf("scaled")))), "no raw or log-normalized measure")
})

# A1 contract: extend_for_viewer output is identical whether counts arrive CSC (CsparseMatrix) or CSR
# (RsparseMatrix) -- R normalizes on input like Python, so encoding must not change any navigator field.
test_that("extend_for_viewer is invariant to the counts encoding (CSC vs CSR)", {
  skip_if_not_installed("Matrix")
  mk <- function(fmt) {
    set.seed(7); nc <- 80L; ng <- 25L                       # same seed -> identical counts values, umap, labels
    m <- Matrix::Matrix(rpois(nc * ng, 1.0) + 1L, nc, ng, sparse = TRUE)
    cnt <- as(m, if (fmt == "csr") "RsparseMatrix" else "CsparseMatrix")
    rownames(cnt) <- paste0("c", 1:nc); colnames(cnt) <- paste0("g", 1:ng)
    structure(list(kind = "sample", spec_version = "0.1", profiles = character(0), dropped = character(0),
      axes = list(cells = list(labels = rownames(cnt), origin = "observed", role = "observation"),
                  genes = list(labels = colnames(cnt), origin = "observed", role = "feature"),
                  umap = list(labels = c("u1", "u2"), origin = "derived", role = "coordinate")),
      fields = list(counts = list(values = cnt, role = "measure", span = c("cells", "genes"), state = "raw",
                                  encoding = if (fmt == "csr") "csr" else "csc"),
                    umap = list(values = matrix(rnorm(nc * 2), nc, 2), role = "embedding", span = c("cells", "umap")),
                    leiden = list(values = factor(paste0("k", (0:(nc - 1)) %% 4)), role = "label", span = "cells", encoding = "categorical"))),
      class = "lstar_dataset")
  }
  a <- extend_for_viewer(mk("csc"), grouping = "leiden")
  b <- extend_for_viewer(mk("csr"), grouping = "leiden")
  expect_equal(a$fields[["counts_cellmajor_order"]]$values, b$fields[["counts_cellmajor_order"]]$values)
  expect_equal(as.matrix(a$fields[["counts_cellmajor"]]$values), as.matrix(b$fields[["counts_cellmajor"]]$values))
  expect_equal(a$fields[["od_score"]]$values, b$fields[["od_score"]]$values)
  expect_equal(a$fields[["stats_leiden_sum"]]$values, b$fields[["stats_leiden_sum"]]$values)
  expect_equal(a$fields[["markers_leiden_lfc"]]$values, b$fields[["markers_leiden_lfc"]]$values)
})
