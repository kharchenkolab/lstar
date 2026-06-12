#!/usr/bin/env bash
# SingleCellExperiment variety: counts-only; +logcounts + reducedDims (PCA/UMAP); +altExps (ADT, a
# second feature space -> Tier-3, must be *recorded* in dropped); +colData/rowData factors + free-form
# metadata (recorded). Constructed via SCE's own constructors (real class, deterministic), grounded in
# an analysis of real SCEs; each round-trips through the profile + validates.
set -e
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
RLIB="$ROOT/.Rlib"

Rscript -e '.libPaths(c("'"$RLIB"'", .libPaths())); suppressMessages({library(SingleCellExperiment); library(SummarizedExperiment); library(S4Vectors); library(lstar)})
set.seed(1); ng <- 20; nc <- 12
counts <- matrix(rpois(ng*nc, 2), ng, nc, dimnames = list(paste0("Gene",1:ng), paste0("Cell",1:nc)))
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
altExp(sce_c, "ADT") <- SummarizedExperiment(assays = list(counts =
  matrix(rpois(5*nc,10), 5, nc, dimnames = list(paste0("ADT",1:5), colnames(counts)))))
stores <- c(stores, rt("+altExps(ADT)", sce_c, function(ds, s2)   # altExp captured as a 2nd feature space
  "ADT" %in% names(ds$axes) && identical(ds$axes$ADT$role, "feature") && "ADT" %in% altExpNames(s2)))

sce_d <- sce_b
colData(sce_d)$cluster <- factor(rep(c("c1","c2"), 6)); rowData(sce_d)$type <- factor(rep(c("g1","g2"), 10))
metadata(sce_d)$study <- "demo"
stores <- c(stores, rt("+colData/rowData factors +metadata", sce_d, function(ds, s2)
  identical(ds$axes$cluster$role, "factor") && identical(ds$axes$type$role, "factor") &&
  "metadata/study" %in% ds$dropped))
cat(paste(stores, collapse="\n"), "\n", file = "/tmp/scev_stores.txt")
' 2>&1 | grep -E "^  \[R\]"

PYTHONPATH="$ROOT/python/src" python3 - <<'PY'
import warnings; warnings.filterwarnings("ignore")
import lstar
stores = [s for s in open("/tmp/scev_stores.txt").read().split() if s]
for p in stores:
    errs = [i for i in lstar.validate(lstar.read(p)) if i.startswith("ERROR")]
    assert not errs, (p, errs)
print("  [py] all %d SCE version-variant stores validate clean" % len(stores))
PY
