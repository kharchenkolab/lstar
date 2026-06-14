#!/usr/bin/env bash
# Collection-model conformance: a *collection of samples* (not an aligned matrix) must
# round-trip through L* and read back identically in another language. This builds a small
# synthetic two-sample collection in R -- per-sample cells/genes axes + counts, a samples axis,
# a union cells axis with a sample label, and a joint graph as a (cells x cells) relation --
# writes it, and verifies Python reads & validates it with the structure and heavy fields intact.
# Self-contained: no external datasets (the real-data version is examples/conos_collection_demo.R).
# Origin coverage: R-authored ✓ (builds the collection) | Python cross-reads — see conformance/README.md
set -e
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
RLIB="$ROOT/.Rlib"
STORE=/tmp/collection_conf.lstar.zarr

Rscript -e '.libPaths(c("'"$RLIB"'", .libPaths()))
suppressMessages({library(lstar); library(Matrix)})
set.seed(1)
mk <- function(prefix, nc, ng) {
  m <- as(matrix(rpois(nc*ng, 0.3), nc, ng), "CsparseMatrix")
  rownames(m) <- paste0(prefix, "_c", seq_len(nc)); colnames(m) <- paste0("g", seq_len(ng)); m
}
A <- mk("S1", 40, 25); B <- mk("S2", 60, 25)          # two samples, different cell counts
ucells <- c(rownames(A), rownames(B))                  # union cells
ds <- list(kind="collection", spec_version="0.1", profiles="synthetic@0.1",
           dropped=character(0), axes=list(), fields=list())
ds$axes[["cells.S1"]] <- list(labels=rownames(A), origin="observed", role="observation")
ds$axes[["genes.S1"]] <- list(labels=colnames(A), origin="observed", role="feature")
ds$axes[["cells.S2"]] <- list(labels=rownames(B), origin="observed", role="observation")
ds$axes[["genes.S2"]] <- list(labels=colnames(B), origin="observed", role="feature")
ds$axes[["samples"]]  <- list(labels=c("S1","S2"), origin="observed", role="sample")
ds$axes[["cells"]]    <- list(labels=ucells, origin="derived", role="observation")
ds$fields[["counts.S1"]] <- list(role="measure", span=c("cells.S1","genes.S1"), state="raw", values=A)
ds$fields[["counts.S2"]] <- list(role="measure", span=c("cells.S2","genes.S2"), state="raw", values=B)
ds$fields[["sample"]]    <- list(role="label", span="cells", values=rep(c("S1","S2"), c(40,60)))
# a joint knn-ish graph over the union cells, as a relation
G <- as(Matrix::rsparsematrix(100, 100, 0.05, symmetric=TRUE), "CsparseMatrix")
dimnames(G) <- list(ucells, ucells)
ds$fields[["graph"]] <- list(role="relation", span=c("cells","cells"), values=G)
class(ds) <- "lstar_dataset"
lstar_write(ds, "'"$STORE"'")
cat("  [R ] synthetic collection (2 samples) -> L*\n")' 2>&1 | grep -vE "^Warning|deprecat|Attaching|masked|^$"

PYTHONPATH="$ROOT/python/src" python3 - "$STORE" <<'PY'
import sys, warnings; warnings.filterwarnings("ignore")
import lstar
ds = lstar.read(sys.argv[1])
assert ds.kind == "collection", ds.kind
errs = lstar.validate(ds)
assert not errs, errs
# structure: a samples axis, per-sample measures over their own axes, a relation over cells x cells
assert "samples" in ds.axes and len(ds.axes["samples"]) == 2
a, b = ds.fields["counts.S1"].values, ds.fields["counts.S2"].values
assert a.shape == (40, 25) and b.shape == (60, 25), (a.shape, b.shape)
g = ds.fields["graph"]
assert g.role == "relation" and list(g.span) == ["cells", "cells"], (g.role, g.span)
assert g.values.shape == (100, 100)
assert ds.fields["sample"].role == "label"
print("  [py] read collection: kind=collection, samples axis, per-sample measures, "
      "graph relation -- validate clean")
PY

# Seurat v5 split-assay collection (skips gracefully if SeuratObject is unavailable).
SSTORE=/tmp/collection_seurat_conf.lstar.zarr
Rscript -e '.libPaths(c("'"$RLIB"'", .libPaths()))
if (!requireNamespace("SeuratObject", quietly=TRUE)) { cat("  [R ] SeuratObject absent - skip\n"); quit(status=0) }
suppressMessages({library(lstar); library(SeuratObject); library(Matrix)})
set.seed(2)
m <- matrix(rpois(80*150,0.5),80,150,dimnames=list(paste0("g",1:80),paste0("c",1:150)))
obj <- CreateSeuratObject(counts=as(m,"dgCMatrix")); obj$sample <- rep(c("A","B"),length.out=150)
obj[["RNA"]] <- split(obj[["RNA"]], f=obj$sample)
ds <- read_seurat(obj)
stopifnot(ds$kind=="collection", "samples" %in% names(ds$axes), "counts.A" %in% names(ds$fields))
lstar_write(ds, "'"$SSTORE"'"); cat("  [R ] Seurat v5 split assay -> L* collection\n")' \
  2>&1 | grep -vE "^Warning|deprecat|Attaching|masked|^$"
if [ -d "$SSTORE" ]; then
PYTHONPATH="$ROOT/python/src" python3 - "$SSTORE" <<'PY'
import sys, warnings; warnings.filterwarnings("ignore")
import lstar
ds = lstar.read(sys.argv[1])
assert ds.kind == "collection" and not lstar.validate(ds)
assert ds.fields["counts.A"].values.shape[1] == 80  # shared genes; per-sample cells
print("  [py] Seurat collection read & validated clean")
PY
fi
echo "collection conformance PASSED."
