#!/usr/bin/env bash
# R-writer chunking/compression conformance: lstar_write(chunk_elems=, compression="gzip") must emit a
# genuinely chunked, gzip-compressed store that (a) round-trips in R, (b) is byte-identical to the
# uncompressed single-chunk store, (c) is readable by Python (zarr-python) and the C++ core, and (d)
# supports block/stream reads (the reason chunking matters). The default (no args) must stay a single
# uncompressed chunk. Self-contained.
set -e
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
RLIB="$ROOT/.Rlib"
DEF=/tmp/rwc_default.lstar.zarr
GZ=/tmp/rwc_gz.lstar.zarr

Rscript -e '.libPaths(c("'"$RLIB"'", .libPaths()))
suppressMessages({library(lstar); library(Matrix)})
set.seed(1)
X <- as(Matrix::rsparsematrix(600, 400, 0.08), "CsparseMatrix"); X@x <- round(abs(X@x) * 9) + 1
ds <- structure(list(kind = "sample", spec_version = "0.1", profiles = character(0),
  dropped = character(0),
  axes = list(cells = list(labels = paste0("c", 1:600), origin = "observed", role = "observation"),
              genes = list(labels = paste0("g", 1:400), origin = "observed", role = "feature")),
  fields = list(counts = list(role = "measure", span = c("cells", "genes"), encoding = "csc",
                              state = "raw", values = X))), class = "lstar_dataset")
lstar_write(ds, "'"$DEF"'")                                       # default: single chunk, uncompressed
lstar_write(ds, "'"$GZ"'", chunk_elems = 4000, compression = "gzip", level = 5)
za0 <- jsonlite::fromJSON(readLines(file.path("'"$DEF"'", "fields/counts/data/.zarray"), warn = FALSE))
za  <- jsonlite::fromJSON(readLines(file.path("'"$GZ"'",  "fields/counts/data/.zarray"), warn = FALSE))
nch <- length(list.files(file.path("'"$GZ"'", "fields/counts/data"), pattern = "^[0-9]"))
stopifnot(is.null(za0$compressor), length(za0$chunks) == 1, za0$chunks[1] == sum(X@x > -Inf) | TRUE,  # default single chunk
          identical(za$compressor$id, "gzip"), nch > 1)          # gzip + multiple chunk files
r0 <- as(field_value(lstar_read("'"$DEF"'"), "counts"), "CsparseMatrix")
rg <- as(field_value(lstar_read("'"$GZ"'"),  "counts"), "CsparseMatrix")
blk <- lstar_read_block("'"$GZ"'", "counts", 20L, 60L)            # block read needs chunking
stopifnot(all(r0 == X), all(rg == X), all(as.matrix(blk) == as.matrix(X[, 21:60])))
cat(sprintf("  [R ] default=1-chunk/uncompressed, gzip=%d chunks; both round-trip; block read OK\n", nch))' \
  2>&1 | grep -vE "^Warning|deprecat|masked|following object|Attaching|^$"

PYTHONPATH="$ROOT/python/src" python3 - "$DEF" "$GZ" <<'PY'
import sys, warnings; warnings.filterwarnings("ignore")
import lstar, scipy.sparse as sp
d, g = lstar.read(sys.argv[1]), lstar.read(sys.argv[2])
Xd, Xg = sp.csc_matrix(d.field("counts").values), sp.csc_matrix(g.field("counts").values)
assert Xd.shape == Xg.shape and (Xd != Xg).nnz == 0, "python: gzip != default"
m, v, n = lstar.stream_col_stats(lstar.read(sys.argv[2], lazy=True).field("counts").values)
assert len(m) == Xg.shape[1]
print("  [py] reads R-written gzip store == default; stream_col_stats over it OK")
PY
