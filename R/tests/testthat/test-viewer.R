# write_viewer (R) builds the viewer profile via the bound kernel; mirrors the Python test so the
# two stay in lockstep (the bridge: R and Python produce the same fields with the same numbers).
test_that("write_viewer adds the viewer profile; stats kernel-exact, panel = counts re-oriented", {
  skip_if_not_installed("Matrix")
  set.seed(7); nc <- 140L; ng <- 22L
  cnt <- as(Matrix::Matrix(rpois(nc * ng, 0.7), nc, ng, sparse = TRUE), "CsparseMatrix")
  rownames(cnt) <- paste0("cell", 1:nc); colnames(cnt) <- paste0("g", 1:ng)
  lab <- paste0("k", (0:(nc - 1)) %% 5)
  ds <- list(kind = "sample", spec_version = "0.1", profiles = character(0), dropped = character(0),
             axes = list(cells = list(labels = rownames(cnt), origin = "observed", role = "observation"),
                         genes = list(labels = colnames(cnt), origin = "observed", role = "feature")),
             fields = list(counts = list(role = "measure", span = c("cells", "genes"), state = "raw", values = cnt),
                           leiden = list(role = "label", span = "cells", values = lab)))
  class(ds) <- "lstar_dataset"

  ds <- write_viewer(ds, grouping = "leiden")
  expect_true("viewer@0.1" %in% ds$profiles)
  expect_true(all(c("stats_leiden_sum", "markers_leiden_lfc", "cell_order", "counts_cellmajor")
                  %in% names(ds$fields)))
  expect_false("od_genes" %in% names(ds$axes))          # gene scope is on-the-fly, never precomputed

  # cluster sufficient stats == Matrix per-group colSums(log1p)
  Xl <- cnt; Xl@x <- log1p(Xl@x); groups <- sort(unique(lab))
  ref <- t(sapply(groups, function(g) Matrix::colSums(Xl[lab == g, , drop = FALSE])))
  expect_lt(max(abs(ds$fields$stats_leiden_sum$values - ref)), 1e-6)

  # counts_cellmajor is the same counts, cell-major (CSR), raw
  expect_equal(ds$fields$counts_cellmajor$encoding, "csr")
  expect_equal(max(abs(as.matrix(ds$fields$counts_cellmajor$values) - as.matrix(cnt))), 0)

  # round-trips through the store
  p <- file.path(tempdir(), "viewer.lstar.zarr"); if (dir.exists(p)) unlink(p, recursive = TRUE)
  lstar_write(ds, p); ds2 <- lstar_read(p)
  expect_true(all(c("counts_cellmajor", "stats_leiden_sum") %in% names(ds2$fields)))
})
