# write_pagoda2: a mock Pagoda2-shaped object -> the canonical viewer store. (A real Pagoda2
# object needs the pagoda2 package; this mocks the public accessors the exporter reads.)
test_that("write_pagoda2 builds the canonical store + viewer profile", {
  skip_if_not_installed("Matrix")
  set.seed(3); nc <- 90L; ng <- 18L
  cnt <- as(Matrix::Matrix(rpois(nc * ng, 0.6), nc, ng, sparse = TRUE), "CsparseMatrix")
  rownames(cnt) <- paste0("cell", 1:nc); colnames(cnt) <- paste0("g", 1:ng)
  emb <- matrix(rnorm(nc * 2), nc, 2); rownames(emb) <- rownames(cnt)
  meta <- data.frame(leiden = paste0("c", (0:(nc - 1)) %% 4), cell_type = paste0("t", (0:(nc - 1)) %% 4),
                     mito = runif(nc, 1, 10), row.names = rownames(cnt), stringsAsFactors = FALSE)
  p2 <- list(getRawCounts = function(...) cnt, embeddings = list(PCA = list(UMAP = emb)), cellMeta = meta, misc = list())

  ds <- write_pagoda2(p2, grouping = "leiden")
  expect_true("viewer@0.1" %in% ds$profiles)
  expect_true(all(c("counts", "umap", "leiden", "cell_type", "mito",
                    "stats_leiden_sum", "markers_leiden_lfc") %in% names(ds$fields)))

  # cluster stats (via the bound libstar kernel) == Matrix reference
  Xl <- cnt; Xl@x <- log1p(Xl@x); groups <- sort(unique(meta$leiden))
  ref <- t(sapply(groups, function(g) Matrix::colSums(Xl[meta$leiden == g, , drop = FALSE])))
  expect_lt(max(abs(ds$fields$stats_leiden_sum$values - ref)), 1e-5)

  # round-trips through the L* store
  p <- file.path(tempdir(), "p2.lstar.zarr"); if (dir.exists(p)) unlink(p, recursive = TRUE)
  write_pagoda2(p2, p)
  ds2 <- lstar_read(p)
  expect_true("viewer@0.1" %in% ds2$profiles)
  expect_true("stats_leiden_sum" %in% names(ds2$fields))
})
