#!/usr/bin/env bash
# Categorical-encoding conformance: a label stored as integer codes + an ordered category set + a `-1`
# missing sentinel must round-trip across languages with codes, category order, and missingness intact.
# Python writes -> Python/C++(via R) read it back as the same factor; R writes a factor -> Python reads
# the same codes/categories/ordered. This is the Tier-1 dtype-fidelity gate and the factor-axis substrate.
# Origin coverage: Py-authored ✓ | R-authored ✓ (each cross-read by the other language) — see conformance/README.md
set -e
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
RLIB="$ROOT/.Rlib"
PY=/tmp/cat_conf_py.lstar.zarr
R=/tmp/cat_conf_r.lstar.zarr

PYTHONPATH="$ROOT/python/src" python3 - "$PY" <<'PY'
import sys, warnings; warnings.filterwarnings("ignore")
import numpy as np, lstar
from lstar import Categorical
ds = lstar.Dataset(kind="sample")
ds.add_axis("cells", [f"c{i}" for i in range(7)])
ds.add_field("leiden", Categorical(np.array([0, 2, -1, 1, 1, -1, 0]), np.array(["A", "B", "C"]),
                                   ordered=True), span=["cells"])
assert not lstar.validate(ds)
lstar.write(ds, sys.argv[1])
v = lstar.read(sys.argv[1]).field("leiden").values            # Python read-back
assert (v.codes == [0, 2, -1, 1, 1, -1, 0]).all() and list(v.categories) == ["A", "B", "C"] and v.ordered
print("  [py] wrote + read categorical (codes/categories/ordered/-1 exact)")
PY

Rscript -e '.libPaths(c("'"$RLIB"'", .libPaths())); suppressMessages(library(lstar))
# C++ core (via the R reader) reads the Python store as an ordered factor with NA for -1
v <- field_value(lstar_read("'"$PY"'"), "leiden")
stopifnot(is.factor(v), is.ordered(v), identical(levels(v), c("A","B","C")),
          identical(as.character(v), c("A","C",NA,"B","B",NA,"A")))
# R writes a factor -> read back identical
ds <- structure(list(kind="sample", spec_version="0.1", profiles=character(0), dropped=character(0),
  axes=list(cells=list(labels=paste0("c",1:5), origin="observed", role="observation")),
  fields=list(ct=list(values=factor(c("x","y","x",NA,"z"), levels=c("x","y","z"), ordered=TRUE),
                      role="label", span="cells"))), class="lstar_dataset")
lstar_write(ds, "'"$R"'")
rb <- field_value(lstar_read("'"$R"'"), "ct")
stopifnot(is.ordered(rb), identical(levels(rb), c("x","y","z")),
          identical(as.character(rb), c("x","y","x",NA,"z")))
cat("  [R ] C++/R read Python categorical as ordered factor; R-written factor round-trips\n")' \
  2>&1 | grep -vE "^Warning|deprecat|masked|following object|Attaching|^$"

PYTHONPATH="$ROOT/python/src" python3 - "$R" <<'PY'
import sys, warnings; warnings.filterwarnings("ignore")
import numpy as np, lstar
v = lstar.read(sys.argv[1]).field("ct").values                # Python reads the R-written categorical
assert list(v.categories) == ["x", "y", "z"] and v.ordered and list(v.codes) == [0, 1, 0, -1, 2]
print("  [py] read R-written categorical (codes/categories/ordered match)")
PY
