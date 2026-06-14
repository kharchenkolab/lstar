#!/usr/bin/env bash
# Partial coverage round-trips through the R path, including when `index_axis` is a **derived union**
# axis of a collection (the pagoda2/Conos case: a facet measured on a subset of the joint union cells).
# Py -> store -> R -> store -> Py preserves coverage/index/index_axis and the partial values.
# Origin coverage: Py-authored ✓ | R-authored ✓ (each cross-read by the other language) — see conformance/README.md
set -eo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
RLIB="$ROOT/.Rlib"
S0=/tmp/partc0.lstar.zarr
S1=/tmp/partc1.lstar.zarr

PYTHONPATH="$ROOT/python/src" python3 - "$S0" <<'PY'
import sys, warnings; warnings.filterwarnings("ignore")
import numpy as np, scipy.sparse as sp, lstar
ds = lstar.Dataset(kind="collection")
ds.add_axis("cells", [f"u{i}" for i in range(10)], origin="derived", role="observation")  # DERIVED union
ds.add_axis("genes", [f"g{i}" for i in range(4)])
idx = np.array([0, 2, 4, 6, 8])                              # a facet on 5 of the 10 joint cells
ds.add_field("adt", sp.random(5, 4, density=0.5, format="csc").astype("float32"),
             role="measure", span=["cells", "genes"], state="raw", index=idx, index_axis="cells")
assert not lstar.validate(ds)
lstar.write(ds, sys.argv[1])
print("  [py] wrote a partial measure over a DERIVED union cells axis (5 of 10)")
PY

Rscript -e '.libPaths(c("'"$RLIB"'", .libPaths())); suppressMessages(library(lstar))
ds <- lstar_read("'"$S0"'"); f <- ds$fields[["adt"]]
stopifnot(identical(f$coverage, "partial"), identical(f$index_axis, "cells"),
          length(f$index) == 5, identical(ds$axes$cells$origin, "derived"), nrow(f$values) == 5)
lstar_write(ds, "'"$S1"'")
cat("  [R ] read partial coverage over the derived axis + rewrote\n")' 2>&1 | grep -E "^  \[R"

PYTHONPATH="$ROOT/python/src" python3 - "$S1" <<'PY'
import sys, warnings; warnings.filterwarnings("ignore")
import numpy as np, lstar
ds = lstar.read(sys.argv[1]); f = ds.field("adt")
assert f.coverage == "partial" and f.index_axis == "cells"
assert list(np.asarray(f.index)) == [0, 2, 4, 6, 8] and f.values.shape == (5, 4)
assert ds.axes["cells"].origin == "derived" and not lstar.validate(ds)
print("  [py] partial coverage over a derived union axis survived Py -> R -> Py")
PY
# R-AUTHORED partial coverage (origin coverage): R builds a partial CSC facet over a DERIVED union axis
# from scratch -> Python reads coverage/index/values. The Py-origin flow exercises R only as a rewriter.
S2=/tmp/partc2.lstar.zarr
Rscript -e '.libPaths(c("'"$RLIB"'", .libPaths())); suppressMessages({library(lstar); library(Matrix)})
M <- matrix(0, 5, 4); M[1,1]<-10; M[2,2]<-20; M[3,3]<-30; M[5,4]<-40      # deterministic, value-checkable
vals <- as(Matrix::Matrix(M, sparse=TRUE), "dgCMatrix")
ds <- structure(list(kind="collection", spec_version="0.1", profiles=character(0), dropped=character(0),
  axes=list(cells=list(labels=paste0("u",0:9), origin="derived", role="observation"),
            genes=list(labels=paste0("g",0:3), origin="observed", role="feature")),
  fields=list(adt=list(values=vals, role="measure", span=c("cells","genes"), state="raw",
                       encoding="csc", index=c(0L,2L,4L,6L,8L), index_axis="cells"))),
  class="lstar_dataset")
lstar_write(ds, "'"$S2"'")
f <- lstar_read("'"$S2"'")$fields[["adt"]]
stopifnot(identical(f$coverage,"partial"), identical(f$index_axis,"cells"), length(f$index)==5, nrow(f$values)==5)
cat("  [R ] authored a partial facet (5 of 10 derived cells) from scratch, rewrote\n")' 2>&1 | grep -E "^  \[R"

PYTHONPATH="$ROOT/python/src" python3 - "$S2" <<'PY'
import sys, warnings; warnings.filterwarnings("ignore")
import numpy as np, scipy.sparse as sp, lstar
ds = lstar.read(sys.argv[1]); f = ds.field("adt")
assert f.coverage == "partial" and f.index_axis == "cells"
assert list(np.asarray(f.index)) == [0, 2, 4, 6, 8] and f.values.shape == (5, 4)
m = sp.csc_matrix((f.values.data, f.values.indices, f.values.indptr), shape=(5, 4)).toarray()
assert m[0,0]==10 and m[1,1]==20 and m[2,2]==30 and m[4,3]==40, m
assert ds.axes["cells"].origin == "derived" and not lstar.validate(ds)
print("  [py] read R-authored partial facet: coverage/index/values exact (R-origin cross-read)")
PY
echo "partial-coverage conformance PASSED."
