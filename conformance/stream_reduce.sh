#!/usr/bin/env bash
# Blocked-reader conformance: a per-gene reduction over a measure must be computable in *bounded
# memory* -- reading the store block-by-block -- and agree exactly with the whole-matrix result.
# Python writes a chunked CSC store; the R/C++ blocked reducer (lstar::stream_col_stats) reads only
# each column block's chunks and must match a full-read reduction, for both raw and log1p, and be
# invariant to the thread count. This exercises read_array_range + stream_csc_col_mean_var in the
# C++ core through the R binding. Self-contained (synthetic data).
set -e
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
RLIB="$ROOT/.Rlib"
STORE=/tmp/stream_reduce_conf.lstar.zarr

PYTHONPATH="$ROOT/python/src" python3 - "$STORE" <<'PY'
import sys, warnings; warnings.filterwarnings("ignore")
import numpy as np, scipy.sparse as sp
import lstar
np.random.seed(0)
X = sp.random(600, 400, density=0.07, format="csc", random_state=0)
X.data = np.round(X.data * 9 + 1).astype("f4")
ds = lstar.Dataset(kind="sample")
ds.add_axis("cells", [f"c{i}" for i in range(600)])
ds.add_axis("genes", [f"g{i}" for i in range(400)])
ds.add_field("counts", X, role="measure", span=["cells", "genes"], state="raw")
lstar.write(ds, sys.argv[1], chunk_elems=1500)        # chunked so the blocked reader truly streams
# Python's own streamed reduction must equal an eager one (the reference side).
lazy = lstar.read(sys.argv[1], lazy=True).field("counts").values
m, v, n = lstar.stream_col_stats(lazy, lognorm=True)
full = sp.csc_matrix(lstar.read(sys.argv[1]).field("counts").values)
fm, fv, fn = lstar.stream_col_stats(full, lognorm=True)
assert np.allclose(m, fm) and np.allclose(v, fv) and (n == fn).all()
print("  [py] wrote chunked CSC store; python streamed reduction == eager")
PY

Rscript -e '.libPaths(c("'"$RLIB"'", .libPaths()))
suppressMessages({library(lstar); library(Matrix)})
store <- "'"$STORE"'"

# blocked C++ reduction (small block forces multi-chunk streaming), raw + log1p
s_raw <- stream_col_stats(store, "counts", block = 29L, n_threads = 1L, lognorm = FALSE)
s_ln  <- stream_col_stats(store, "counts", block = 29L, n_threads = 1L, lognorm = TRUE)

# ground truth: read the whole measure and reduce it with Matrix (zero-aware)
m  <- as(field_value(lstar_read(store), "counts"), "CsparseMatrix")
nr <- nrow(m)
cmv <- function(M, lognorm) {
  X <- M; if (lognorm) X@x <- log1p(X@x)
  cs <- Matrix::colSums(X); mu <- cs / nr
  ss <- Matrix::colSums(X * X) - cs * cs / nr
  list(mean = as.numeric(mu), var = as.numeric(ss / (nr - 1)), nnz = diff(M@p))
}
g_raw <- cmv(m, FALSE); g_ln <- cmv(m, TRUE)

# DETERMINISM CONTRACT: the multithreaded C++ kernel is BIT-identical across thread counts (float64
# accumulation, column-parallel, no cross-thread reduction). mean/var/nnz at 2/4/8 threads must equal
# the 1-thread result exactly (==0), not just within tolerance.
for (nt in c(2L, 4L, 8L)) {
  sN <- stream_col_stats(store, "counts", block = 29L, n_threads = nt, lognorm = TRUE)
  stopifnot(max(abs(s_ln$mean - sN$mean)) == 0, max(abs(s_ln$var - sN$var)) == 0, all(s_ln$nnz == sN$nnz))
}

tol <- 1e-9
stopifnot(max(abs(s_raw$mean - g_raw$mean)) < tol,
          max(abs(s_raw$var  - g_raw$var )) < tol,
          all(s_raw$nnz == g_raw$nnz),
          max(abs(s_ln$mean - g_ln$mean)) < tol,
          max(abs(s_ln$var  - g_ln$var )) < tol)
cat(sprintf("  [R ] blocked reduction == full read (%d genes; raw + log1p; bit-identical 1==2==4==8 threads)\n",
            length(s_raw$mean)))' 2>&1 | grep -vE "^Warning|deprecat|Attaching|masked|following object|^$"
