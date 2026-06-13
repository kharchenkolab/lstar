#!/usr/bin/env bash
# SingleCellExperiment variety: counts-only; +logcounts + reducedDims (PCA/UMAP); +altExps (ADT, a
# second feature space -> Tier-3, must be *recorded* in dropped); +colData/rowData factors + free-form
# metadata (recorded). Constructed via SCE's own constructors (real class, deterministic), grounded in
# an analysis of real SCEs; each round-trips through the profile + validates.
set -eo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
RLIB="$ROOT/.Rlib"

# Synthetic CITE-seq via the shared Python generator (same RNA+ADT as the Seurat/MuData tests); the real
# breadth is local-only (sweep / real_corpus_r.sh). Nothing committed.
CDIR="$ROOT/testdata/citeseq"
python3 "$ROOT/python/tests/synth.py" citeseq "$CDIR" >/dev/null

Rscript -e '.libPaths(c("'"$RLIB"'", .libPaths())); suppressMessages({library(SingleCellExperiment); library(SummarizedExperiment); library(S4Vectors); library(lstar)})
set.seed(1)
# Synthetic CITE-seq counts: synthetic RNA as the base assay, synthetic ADT as the altExp. reducedDims /
# factors / metadata are constructed structures on top.
cdir <- "'"$CDIR"'"
ccells <- readLines(file.path(cdir, "cells.txt"))
counts <- as.matrix(Matrix::t(Matrix::readMM(file.path(cdir, "rna.mtx"))))   # genes x cells
dimnames(counts) <- list(readLines(file.path(cdir, "genes.txt")), ccells)
adtm <- as.matrix(Matrix::t(Matrix::readMM(file.path(cdir, "adt.mtx"))))
dimnames(adtm) <- list(readLines(file.path(cdir, "proteins.txt")), ccells)
ng <- nrow(counts); nc <- ncol(counts)                                       # 27 genes x 80 cells
rt <- function(tag, sce, checks) {
  ds <- read_sce(sce); sce2 <- write_sce(ds); stopifnot(checks(ds, sce2))
  cat(sprintf("  [R] %-30s OK\n", tag))
  p <- file.path("/tmp", paste0("sce_", gsub("[^a-z0-9]","",tolower(tag)), ".lstar.zarr"))
  if (dir.exists(p)) unlink(p, recursive = TRUE); lstar_write(ds, p); p
}
stores <- c()
stores <- c(stores, rt("counts only", SingleCellExperiment(assays = list(counts = counts)),
  function(ds, s2) identical(assayNames(s2), "counts")))

sce_b <- SingleCellExperiment(assays = list(counts = counts, logcounts = log1p(counts)))
reducedDims(sce_b) <- list(PCA = matrix(rnorm(nc*5), nc, 5, dimnames = list(colnames(counts), paste0("PC",1:5))),
                           UMAP = matrix(rnorm(nc*2), nc, 2, dimnames = list(colnames(counts), c("UMAP1","UMAP2"))))
stores <- c(stores, rt("+logcounts +reducedDims", sce_b, function(ds, s2)
  setequal(assayNames(s2), c("counts","logcounts")) && setequal(reducedDimNames(s2), c("PCA","UMAP"))))

sce_c <- sce_b
altExp(sce_c, "ADT") <- SummarizedExperiment(assays = list(counts = adtm))   # synthetic ADT (29 proteins)
stores <- c(stores, rt("+altExps(ADT synth)", sce_c, function(ds, s2)   # altExp captured as 2nd feature space
  "ADT" %in% names(ds$axes) && identical(ds$axes$ADT$role, "feature") && "ADT" %in% altExpNames(s2) &&
  length(ds$axes$ADT$labels) == 29))

sce_d <- sce_b
colData(sce_d)$cluster <- factor(rep(c("c1","c2"), length.out = nc))
rowData(sce_d)$type <- factor(rep(c("g1","g2"), length.out = ng))
metadata(sce_d)$study <- "demo"
stores <- c(stores, rt("+colData/rowData factors +metadata", sce_d, function(ds, s2)
  identical(ds$axes$cluster$role, "factor") && identical(ds$axes$type$role, "factor") &&
  "metadata/study" %in% ds$dropped))

# colPairs/rowPairs (cell-cell SNN / gene-gene graph) -> relations over (cells,cells)/(genes,genes)
sce_e <- sce_b
SingleCellExperiment::colPair(sce_e, "knn") <- Matrix::sparseMatrix(
  i = rep(1:nc, each = 2), j = ((seq_len(2*nc)) %% nc) + 1, x = 1, dims = c(nc, nc))
SingleCellExperiment::rowPair(sce_e, "corr") <- Matrix::sparseMatrix(
  i = rep(1:ng, each = 2), j = ((seq_len(2*ng)) %% ng) + 1, x = 1, dims = c(ng, ng))
stores <- c(stores, rt("+colPairs/rowPairs (graphs)", sce_e, function(ds, s2)
  identical(as.character(ds$fields[["colpair_knn"]]$span), c("cells","cells")) &&
  identical(ds$fields[["colpair_knn"]]$role, "relation") &&
  setequal(SingleCellExperiment::colPairNames(s2), "knn") &&
  setequal(SingleCellExperiment::rowPairNames(s2), "corr")))

# --- sweep-caught real structures, now guarded synthetically (CI parity with the local scRNAseq sweep) ---
# NULL dimnames: cells keyed by a `Barcode` colData column, no colnames (BachMammary/Ernst) -> labels synthesized
sce_nd <- SingleCellExperiment(assays = list(counts = unname(counts)))
colData(sce_nd)$Barcode <- paste0("BC", seq_len(nc))
stores <- c(stores, rt("NULL dimnames (Barcode colData)", sce_nd, function(ds, s2)
  length(ds$axes$cells$labels) == nc && all(grepl("^BC", ds$axes$cells$labels))))

# S4Vectors `Rle` run-length colData -> unpacked; a nested DataFrame column -> recorded, not crashed
# (Buettner/Bunis/Darmanis). Blind as.character() on these was the sweep-caught bug.
sce_s4 <- sce_b
colData(sce_s4)$rle_label <- S4Vectors::Rle(rep(c("a", "b"), length.out = nc))
colData(sce_s4)$nested <- S4Vectors::DataFrame(x = seq_len(nc), y = seq_len(nc))
stores <- c(stores, rt("S4 Rle + nested colData", sce_s4, function(ds, s2)
  "rle_label" %in% names(ds$fields) && any(grepl("^colData/nested", ds$dropped))))

# a plain `SummarizedExperiment` (NOT an SCE -- ReprocessedFluidigm): SCE-only accessors are guarded,
# degrades to assays + colData/rowData only.
se <- SummarizedExperiment(assays = list(counts = counts))
colData(se)$grp <- factor(rep(c("x", "y"), length.out = nc))
stores <- c(stores, rt("plain SummarizedExperiment", se, function(ds, s2)
  identical(assayNames(s2), "counts") && identical(ds$axes$grp$role, "factor")))
cat(paste(stores, collapse="\n"), "\n", file = "/tmp/scev_stores.txt")
' 2>&1 | grep -E "^  \[R\]|Error|Execution halted|cannot|unable|no method"

PYTHONPATH="$ROOT/python/src" python3 - <<'PY'
import warnings; warnings.filterwarnings("ignore")
import lstar
stores = [s for s in open("/tmp/scev_stores.txt").read().split() if s]
for p in stores:
    errs = [i for i in lstar.validate(lstar.read(p)) if i.startswith("ERROR")]
    assert not errs, (p, errs)
print("  [py] all %d SCE version-variant stores validate clean" % len(stores))
PY
