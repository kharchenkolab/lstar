#!/usr/bin/env bash
# Induction conformance: a categorical `label` induces a bare-named `factor` axis whose labels ARE its
# categories, with an `induced_by` back link. That axis + link must round-trip across languages, and a
# reader must be able to *check* the link (validate: axis labels == inducing field's categories). This
# is the Tier-2 factor-axis gate (per induction_design.md §4): independent per-group results land on one
# axis and align, and drift between an induced axis and its field is caught, never silent.
# Origin coverage: Py-authored ✓ | R-authored ✓ (each cross-read by the other language) — see conformance/README.md
set -e
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
RLIB="$ROOT/.Rlib"
PY=/tmp/induce_py.lstar.zarr
R=/tmp/induce_r.lstar.zarr

# (1) Python: a categorical add auto-induces its factor axis; induced_by survives a Python round-trip.
PYTHONPATH="$ROOT/python/src" python3 - "$PY" <<'PY'
import sys, warnings; warnings.filterwarnings("ignore")
import numpy as np, lstar
from lstar import Categorical
ds = lstar.Dataset(kind="sample")
ds.add_axis("cells", [f"c{i}" for i in range(6)])
ds.add_field("leiden", Categorical(np.array([0, 1, 2, -1, 1, 0]), np.array(["A", "B", "C"]),
                                   ordered=True), span=["cells"])
assert not [i for i in lstar.validate(ds) if i.startswith("ERROR")]
ax = ds.axis("leiden")                                  # auto-induced
assert ax.role == "factor" and ax.origin == "derived" and ax.induced_by == "leiden"
lstar.write(ds, sys.argv[1])
ax = lstar.read(sys.argv[1]).axis("leiden")             # survives the store
assert ax.role == "factor" and ax.induced_by == "leiden" and list(ax.labels) == ["A", "B", "C"]
print("  [py] categorical label auto-induced a factor axis; induced_by round-trips in Python")
PY

# (2) C++ core (via the R reader) reads the Python store's factor axis with role + induced_by intact.
Rscript -e '.libPaths(c("'"$RLIB"'", .libPaths())); suppressMessages(library(lstar))
ax <- lstar_read("'"$PY"'")$axes$leiden
stopifnot(identical(ax$role, "factor"), identical(ax$induced_by, "leiden"),
          identical(ax$labels, c("A","B","C")))
# (3) R writes a factor axis induced by an R factor field; read it back unchanged.
ds <- structure(list(kind="sample", spec_version="0.1", profiles=character(0), dropped=character(0),
  axes=list(cells=list(labels=paste0("c",1:4), origin="observed", role="observation"),
            ct=list(labels=c("x","y"), origin="derived", role="factor", induced_by="ct")),
  fields=list(ct=list(values=factor(c("x","y","x","y"), levels=c("x","y")),
                      role="label", span="cells", encoding="categorical"))), class="lstar_dataset")
lstar_write(ds, "'"$R"'")
ax2 <- lstar_read("'"$R"'")$axes$ct
stopifnot(identical(ax2$role, "factor"), identical(ax2$induced_by, "ct"))
cat("  [R ] C++ core read Python factor axis (role+induced_by); R-written factor axis round-trips\n")' \
  2>&1 | grep -vE "^Warning|deprecat|masked|following object|Attaching|^$"

# (4) Python reads the R-written factor axis: induced_by present AND validate confirms consistency
#     (the axis labels equal the inducing field's categories -- induction is checkable, not conventional).
PYTHONPATH="$ROOT/python/src" python3 - "$R" <<'PY'
import sys, warnings; warnings.filterwarnings("ignore")
import lstar
ds = lstar.read(sys.argv[1])
ax = ds.axis("ct")
assert ax.role == "factor" and ax.induced_by == "ct", (ax.role, ax.induced_by)
assert not [i for i in lstar.validate(ds) if i.startswith("ERROR")]
print("  [py] read R-written factor axis; induced_by + axis<->field consistency confirmed")
PY
