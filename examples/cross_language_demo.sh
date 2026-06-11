#!/usr/bin/env bash
# Cross-language demo: AnnData (Python) -> L* Zarr -> Seurat (R) -> L* -> Python.
# Requires: the Python package (python/src), the R package installed in .Rlib,
# and the C++ core built (for the Python<->C++ leg of test_crossimpl, run separately).
set -e
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
STORE=/tmp/lstar_demo.lstar.zarr
RT=/tmp/lstar_demo_seurat_rt.lstar.zarr
RLIB="$ROOT/.Rlib"

echo "== 1. AnnData (Python) -> L* Zarr =="
PYTHONPATH="$ROOT/python/src:$ROOT/python/tests" python3 - "$STORE" <<'PY'
import sys, warnings; warnings.filterwarnings("ignore")
from test_anndata_profile import make_adata
from lstar import read_anndata, write
ds = read_anndata(make_adata()); write(ds, sys.argv[1])
print("   [py] wrote %s : %d fields, %d axes" % (sys.argv[1], len(ds.fields), len(ds.axes)))
PY

echo "== 2. L* Zarr -> Seurat (R), then Seurat -> L* =="
Rscript -e '.libPaths(c("'"$RLIB"'", .libPaths())); suppressMessages(library(lstar)); suppressMessages(library(SeuratObject))
so <- write_seurat(lstar_read("'"$STORE"'"))
cat(sprintf("   [R ] Seurat: %s genes x cells; layers={%s}; reductions={%s}; loadings %s\n",
    paste(dim(so[["RNA"]]),collapse=" x "),
    paste(Layers(so,assay="RNA"),collapse=","), paste(Reductions(so),collapse=","),
    paste(dim(Loadings(so[["pca"]])),collapse="x")))
lstar_write(read_seurat(so), "'"$RT"'")'

echo "== 3. Seurat-derived L* -> Python (validate) =="
PYTHONPATH="$ROOT/python/src" python3 - "$RT" <<'PY'
import sys, warnings; warnings.filterwarnings("ignore")
import lstar; from lstar.zarr_io import read
ds = read(sys.argv[1]); errs=[i for i in lstar.validate(ds) if i.startswith("ERROR")]
print("   [py] read back: %d fields; validate: %s" % (len(ds.fields), "OK" if not errs else errs))
PY
echo "cross-language demo complete (AnnData <-> L* <-> Seurat)."
