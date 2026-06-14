#!/usr/bin/env bash
# Provenance round-trips through the R path. A field's `provenance` (method params, native location, and
# -- for pagoda2 -- a normalization *recipe*: model + params) must survive Py -> store -> R -> store ->
# Py. The R cpp11 binding used to drop it; this guards the fix. Provenance round-trips at the R boundary
# as a native **named list** <-> JSON object (matching Python's dict), incl. provenance ORIGINATING as an
# R list (pagoda2's facet/recipe + a joint product's input_axes -- the case the axis->facet fallback can't
# recover). Also guards that a malformed dataset raises a clean R error instead of a core dump.
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
stopifnot(is.list(p), identical(p$pagoda2, "rawCounts"),                   # R sees the recipe as a NAMED LIST
          identical(p$recipe$model, "clr"), isTRUE(all.equal(p$recipe$depthScale, 1e4)),
          identical(as.numeric(p$recipe$winsor_caps), c(0, 0.99)))
lstar_write(ds, "'"$S1"'")
cat("  [R ] read provenance as a named list, rewrote store (recipe preserved)\n")' 2>&1 | grep -E "^  \[R"

PYTHONPATH="$ROOT/python/src" python3 - "$S1" <<'PY'
import sys, warnings; warnings.filterwarnings("ignore")
import lstar
f = lstar.read(sys.argv[1]).field("counts")
r = (f.provenance or {}).get("recipe", {})
assert r.get("model") == "clr" and r.get("depthScale") == 1e4 and r.get("winsor_caps") == [0.0, 0.99], f.provenance
assert (f.provenance or {}).get("pagoda2") == "rawCounts"
print("  [py] provenance (incl. the normalization recipe) survived Py -> R -> Py intact")
PY

# Case 6 (pagoda2 regression): provenance ORIGINATING as an R named list -- a NON-STANDARD facet (not
# inferable from an axis name) + a joint product's `input_axes` (S5). This is exactly the case the
# axis->facet fallback can't recover, so it fails loudly if the R writer drops the payload.
S2=/tmp/prov2.lstar.zarr
Rscript -e '.libPaths(c("'"$RLIB"'", .libPaths())); suppressMessages({library(lstar); library(Matrix)})
nc <- 8; ng <- 12; nk <- 3
cx <- as(Matrix::Matrix(matrix(rpois(nc*ng,2),nc,ng),sparse=TRUE),"dgCMatrix")
dimnames(cx) <- list(paste0("c",1:nc), paste0("g",1:ng))
wnn <- matrix(rnorm(nc*nk), nc, nk)
ds <- list(kind="sample", spec_version="0.1", profiles=character(0), dropped=character(0),
  axes=list(cells=list(labels=rownames(cx),origin="observed",role="observation"),
            genes=list(labels=colnames(cx),origin="observed",role="feature"),
            wnn=list(labels=paste0("WNN",1:nk),origin="derived",role="coordinate")),
  fields=list(
    counts=list(role="measure",span=c("cells","genes"),state="raw",encoding="csc",values=cx,
                provenance=list(facet="custom2", model="foo", defaultReduction="CCA")),
    WNN=list(role="embedding",span=c("cells","wnn"),state="",encoding="dense",values=wnn,
             provenance=list(method="WNN", input_axes=c("genes","proteins")))))
class(ds) <- "lstar_dataset"
lstar_write(ds, "'"$S2"'")
rd <- lstar_read("'"$S2"'")
pc <- rd$fields$counts$provenance; pw <- rd$fields$WNN$provenance
stopifnot(is.list(pc), identical(pc$facet,"custom2"), identical(pc$model,"foo"), identical(pc$defaultReduction,"CCA"))
stopifnot(is.list(pw), identical(pw$method,"WNN"), identical(as.character(pw$input_axes), c("genes","proteins")))
cat("  [R ] non-standard facet + joint input_axes provenance round-trips (R origin -> store -> R)\n")' 2>&1 | grep -E "^  \[R"

PYTHONPATH="$ROOT/python/src" python3 - "$S2" <<'PY'
import sys, warnings; warnings.filterwarnings("ignore")
import lstar
ds = lstar.read(sys.argv[1])
pc = ds.field("counts").provenance or {}; pw = ds.field("WNN").provenance or {}
assert pc.get("facet") == "custom2" and pc.get("defaultReduction") == "CCA", pc
assert pw.get("method") == "WNN" and list(pw.get("input_axes", [])) == ["genes", "proteins"], pw
print("  [py] R-origin non-standard facet + input_axes provenance readable in Python (the pagoda2 case)")
PY

# malformed dataset -> a clean R error, not a core dump (pagoda2 robustness ask)
Rscript -e '.libPaths(c("'"$RLIB"'", .libPaths())); suppressMessages({library(lstar); library(Matrix)})
cx <- as(Matrix::Matrix(matrix(rpois(200,2),10,20),sparse=TRUE),"dgCMatrix")
dimnames(cx) <- list(paste0("c",1:10), paste0("g",1:20))
ds <- list(axes=list(cells=list(values=rownames(cx)), genes=list(values=colnames(cx))),   # values= not labels=
           fields=list(counts=list(values=Matrix::t(cx), role="measure", span=c("genes","cells"),
                       state="raw", encoding="csc")))
class(ds) <- "lstar_dataset"
r <- tryCatch({ lstar_write(ds, tempfile(fileext=".lstar.zarr")); "NO-ERROR" },
             error=function(e) if (grepl("no .labels.", conditionMessage(e))) "GUARDED" else "WRONG-ERROR")
stopifnot(identical(r, "GUARDED"))
cat("  [R ] malformed dataset raises a clean error (no segfault)\n")' </dev/null 2>&1 | grep -E "^  \[R"
echo "provenance conformance PASSED."
