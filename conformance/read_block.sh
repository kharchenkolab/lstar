#!/usr/bin/env bash
# Block-reader conformance: lstar_read_block / lstar_read_genes must return exactly the same values
# as a full read of the measure -- for a contiguous gene range AND a scattered gene subset -- while
# touching only the overlapping chunks of a chunked CSC store. This is the general bounded block-read
# primitive consumers (e.g. pagoda2's disk backend) drive to build out-of-core ops. Self-contained.
set -e
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
RLIB="$ROOT/.Rlib"
STORE=/tmp/read_block_conf.lstar.zarr

PYTHONPATH="$ROOT/python/src" python3 - "$STORE" <<'PY'
import sys, warnings; warnings.filterwarnings("ignore")
import numpy as np, scipy.sparse as sp, lstar
np.random.seed(0)
X = sp.csc_matrix(sp.random(400, 250, density=0.08, format="csc", random_state=0))
X.data = np.round(X.data * 9 + 1).astype("f4")
ds = lstar.Dataset(kind="sample")
ds.add_axis("cells", [f"c{i}" for i in range(400)])
ds.add_axis("genes", [f"g{i}" for i in range(250)])
ds.add_field("counts", X, role="measure", span=["cells", "genes"], state="raw")
lstar.write(ds, sys.argv[1], chunk_elems=1200)        # chunked so block reads touch few chunks
print("  [py] wrote chunked CSC store (400x250)")
PY

Rscript -e '.libPaths(c("'"$RLIB"'", .libPaths()))
suppressMessages({library(lstar); library(Matrix)})
store <- "'"$STORE"'"
full <- as(field_value(lstar_read(store), "counts"), "CsparseMatrix")   # 400 x 250 (cells x genes)
gn <- paste0("g", 0:249); cn <- paste0("c", 0:399)

## contiguous gene range [40, 90)
b1 <- lstar_read_block(store, "counts", 40L, 90L, cell_names = cn, gene_names = gn)
stopifnot(identical(dim(b1), c(400L, 50L)),
          max(abs(as.matrix(b1) - as.matrix(full[, 41:90]))) == 0)

## scattered gene subset (run-coalesced), arbitrary order
sel <- c("g100","g3","g3","g248","g7","g199","g8","g9")       # incl. a duplicate + reverse-ish order
sel <- unique(sel)
b2 <- lstar_read_genes(store, "counts", sel, gn, cell_names = cn)
stopifnot(identical(colnames(b2), sel),
          max(abs(as.matrix(b2) - as.matrix(full[, match(sel, gn)]))) == 0)

cat("  [R ] read_block (contiguous) + read_genes (scattered) == full read, exact\n")' 2>&1 \
  | grep -vE "^Warning|deprecat|masked|following object|Attaching|^$"
