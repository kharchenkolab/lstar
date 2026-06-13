#!/usr/bin/env bash
# Provenance round-trips through the R path. A field's `provenance` (method params, native location, and
# -- for pagoda2 -- a normalization *recipe*: model + params) must survive Py -> store -> R -> store ->
# Py. The R cpp11 binding used to drop it; this guards the fix (provenance is carried as an opaque JSON
# string at the R boundary, preserving arbitrary nesting).
set -eo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
RLIB="$ROOT/.Rlib"
S0=/tmp/prov0.lstar.zarr
S1=/tmp/prov1.lstar.zarr

PYTHONPATH="$ROOT/python/src" python3 - "$S0" <<'PY'
import sys, warnings; warnings.filterwarnings("ignore")
import lstar, scipy.sparse as sp
ds = lstar.Dataset(kind="sample")
ds.add_axis("cells", [f"c{i}" for i in range(6)]); ds.add_axis("genes", [f"g{i}" for i in range(3)])
prov = {"pagoda2": "rawCounts",
        "recipe": {"model": "clr", "depthScale": 1e4, "log_base": None, "winsor_caps": [0.0, 0.99]}}
ds.add_field("counts", sp.random(6, 3, density=0.5, format="csc").astype("float32"),
             role="measure", span=["cells", "genes"], state="raw", provenance=prov)
lstar.write(ds, sys.argv[1])
print("  [py] wrote store with a recipe-bearing provenance")
PY

Rscript -e '.libPaths(c("'"$RLIB"'", .libPaths())); suppressMessages(library(lstar))
ds <- lstar_read("'"$S0"'")
p <- ds$fields[["counts"]]$provenance
stopifnot(is.character(p), grepl("\"model\":\"clr\"", gsub("\\s","",p)))   # R sees the recipe (JSON string)
lstar_write(ds, "'"$S1"'")
cat("  [R ] read provenance, rewrote store (recipe seen + preserved)\n")' 2>&1 | grep -E "^  \[R"

PYTHONPATH="$ROOT/python/src" python3 - "$S1" <<'PY'
import sys, warnings; warnings.filterwarnings("ignore")
import lstar
f = lstar.read(sys.argv[1]).field("counts")
r = (f.provenance or {}).get("recipe", {})
assert r.get("model") == "clr" and r.get("depthScale") == 1e4 and r.get("winsor_caps") == [0.0, 0.99], f.provenance
assert (f.provenance or {}).get("pagoda2") == "rawCounts"
print("  [py] provenance (incl. the normalization recipe) survived Py -> R -> Py intact")
PY
echo "provenance conformance PASSED."
