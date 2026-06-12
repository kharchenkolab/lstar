# Conformance: the cpp11-bound libstar kernels (the same C++ core compiled to WASM and
# bound to Python) must match Matrix-based references. "One kernel, every runtime."
test_that("col_sum_by_group matches per-group colSums(log1p)", {
  skip_if_not_installed("Matrix")
  set.seed(1); ncells <- 60L; ngenes <- 10L
  X <- as(Matrix::Matrix(rpois(ncells * ngenes, 0.6), ncells, ngenes, sparse = TRUE), "CsparseMatrix")
  Xl <- X; Xl@x <- log1p(Xl@x)
  group <- as.integer((seq_len(ncells) - 1L) %% 3L)
  r <- lstar:::lstar_cpp_col_sum_by_group(as.double(X@x), X@p, X@i, nrow(X), ncol(X), group, 3L, TRUE)
  sumM <- matrix(r$sum, nrow = 3, byrow = TRUE)
  ref <- t(sapply(0:2, function(g) Matrix::colSums(Xl[group == g, , drop = FALSE])))
  expect_equal(r$ngenes, ngenes)
  expect_lt(max(abs(sumM - ref)), 1e-9)
})

test_that("subsample_de_rank matches group colMeans(log1p) lfc", {
  skip_if_not_installed("Matrix")
  set.seed(2); ncells <- 60L; ngenes <- 10L
  X <- as(Matrix::Matrix(rpois(ncells * ngenes, 0.6), ncells, ngenes, sparse = TRUE), "CsparseMatrix")
  Xl <- X; Xl@x <- log1p(Xl@x)
  Xr <- as(X, "RsparseMatrix")
  mem <- as.integer(ifelse(seq_len(ncells) <= ncells / 2, 0L, 1L))
  d <- lstar:::lstar_cpp_subsample_de_rank(as.double(Xr@x), Xr@p, Xr@j, nrow(Xr), ncol(Xr), mem, TRUE)
  mA <- Matrix::colMeans(Xl[mem == 0, , drop = FALSE]); mB <- Matrix::colMeans(Xl[mem == 1, , drop = FALSE])
  expect_equal(d$nA, sum(mem == 0)); expect_equal(d$nB, sum(mem == 1))
  expect_lt(max(abs(d$meanA - mA)), 1e-9)
  expect_lt(max(abs(d$lfc - (mA - mB))), 1e-9)
})
