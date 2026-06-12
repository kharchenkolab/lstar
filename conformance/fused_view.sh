#!/usr/bin/env bash
# Fused-view reducer conformance: the depth-normalized log1p variants of stream_col_stats and
# lstar_stream_col_sum_by_group (computed in one streamed pass off a chunked CSC store) must match an
# independent in-memory reference -- per-gene mean/var (population) and per-(group,gene) sums of the
# "plain" view log1p(x * depthScale / depth[cell]). Self-contained.
set -e
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
RLIB="$ROOT/.Rlib"
STORE=/tmp/fused_view_conf.lstar.zarr

PYTHONPATH="$ROOT/python/src" python3 - "$STORE" <<'PY'
import sys, warnings; warnings.filterwarnings("ignore")
import numpy as np, scipy.sparse as sp, lstar
np.random.seed(0)
X = sp.csc_matrix(sp.random(500, 300, density=0.1, format="csc", random_state=0))
X.data = np.round(X.data * 9 + 1).astype("f4")
ds = lstar.Dataset(kind="sample")
ds.add_axis("cells", [f"c{i}" for i in range(500)]); ds.add_axis("genes", [f"g{i}" for i in range(300)])
ds.add_field("counts", X, role="measure", span=["cells", "genes"], state="raw")
lstar.write(ds, sys.argv[1], chunk_elems=1500)
print("  [py] wrote chunked CSC store (500x300)")
PY

Rscript -e '.libPaths(c("'"$RLIB"'", .libPaths()))
suppressMessages({library(lstar); library(Matrix)})
store <- "'"$STORE"'"
m  <- as(field_value(lstar_read(store), "counts"), "CsparseMatrix")   # 500 cells x 300 genes
nr <- nrow(m); depth <- Matrix::rowSums(m); depthScale <- 1000

## reference "plain" view: log1p(x * depthScale / depth[cell]); per-gene population mean/var
V <- m; V@x <- log1p(V@x * depthScale / depth[V@i + 1L])
cs <- Matrix::colSums(V); ref_m <- as.numeric(cs / nr)
ref_v <- as.numeric(Matrix::colSums(V * V) / nr - ref_m^2)            # population (/n)
s <- stream_col_stats(store, "counts", block = 64L, n_threads = 4L, lognorm = TRUE,
                      depth = as.numeric(depth), depthScale = depthScale, population = TRUE)
stopifnot(max(abs(s$mean - ref_m)) < 1e-9, max(abs(s$var - ref_v)) < 1e-9)

## reference grouped sums of the view; groups 1..G with some NA -> bucket 0
g <- sample(c(NA, 1:4), nr, replace = TRUE); codes <- ifelse(is.na(g), 0L, g)
ref_sum <- matrix(0, 5, ncol(m))                                     # rows: 0(NA),1..4
for (cell in seq_len(nr)) { b <- codes[cell] + 1L
  nz <- (m@p[1]+1):0; }                                              # (unused; vectorized below)
Vt <- as(V, "TsparseMatrix")
ref_sum <- as.matrix(Matrix::sparseMatrix(i = codes[Vt@i + 1L] + 1L, j = Vt@j + 1L, x = Vt@x,
                                          dims = c(5, ncol(m))))
M <- lstar_stream_col_sum_by_group(store, "counts", codes, ngroups = 5L, lognorm = TRUE,
                                   depth = as.numeric(depth), depthScale = depthScale, block = 64L, n_threads = 4L)
stopifnot(max(abs(M - ref_sum)) < 1e-9)
cat("  [R ] fused stream_col_stats (depth+population) & grouped-sum == in-memory reference (exact)\n")' \
  2>&1 | grep -vE "^Warning|deprecat|masked|following object|Attaching|^$"
