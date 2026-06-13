#!/usr/bin/env bash
# LOCAL real Seurat/SCE corpus (not CI — these objects are 100s of MB). Runs the profiles against real
# *published* objects of different versions, so the synthetic version-variety fixtures (seurat_versions.sh
# / sce_versions.sh) stay honest. Skips cleanly where SeuratData / scRNAseq (+ the datasets) aren't
# installed. These objects are exactly what surfaced the real bugs (v5-split write-back; scale.data /
# PCA loadings over a variable-feature subset; ADT / altExp / metadata silent loss).
set -e
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
RLIB="$ROOT/.Rlib"

Rscript -e '.libPaths(c("'"$RLIB"'", .libPaths()))
have <- function(p) requireNamespace(p, quietly = TRUE)
if (!have("SeuratData") && !have("scRNAseq")) { cat("  [skip] real corpus: SeuratData + scRNAseq absent (local-only)\n"); quit(status = 0) }
suppressMessages({library(Seurat); library(SeuratObject); library(lstar)})
n <- 0
chk <- function(tag, ds, so2_or_sce2) { stopifnot(length(ds$fields) > 0); n <<- n + 1
  vers <- grep("^(object@|assay@|SingleCellExperiment@)", ds$profiles, value = TRUE)  # track versions
  cat(sprintf("  [R] %-20s {%s} %d fld/%d ax (%s) | dropped: %s\n",
      tag, paste(vers, collapse=" "), length(ds$fields), length(ds$axes), ds$kind, paste(ds$dropped, collapse="; "))) }

if (have("SeuratData")) {
  suppressMessages(library(SeuratData)); inst <- tryCatch(InstalledData()$Dataset, error = function(e) character(0))
  if ("pbmc3k" %in% inst) {                      # real v3/v4 Assay object (triggers UpdateSeuratObject)
    suppressWarnings(suppressMessages(data("pbmc3k.final"))); so <- UpdateSeuratObject(pbmc3k.final)
    ds <- read_seurat(so); so2 <- write_seurat(ds)
    stopifnot("pca_features" %in% names(ds$axes),          # HVG-subset PCA loadings now TYPED over a
              identical(as.character(ds$fields[["pca_loadings"]]$span), c("pca_features", "pca")))  # subset axis (was dropped)
    chk("pbmc3k.final (v4)", ds, so2)
  }
  if ("cbmc" %in% inst) {                        # real CITE-seq: RNA + ADT (multimodal)
    suppressWarnings(suppressMessages(data("cbmc"))); so <- UpdateSeuratObject(cbmc)
    ds <- read_seurat(so); so2 <- write_seurat(ds)
    stopifnot("ADT" %in% names(ds$axes), identical(ds$axes$ADT$role, "feature"),  # ADT captured as a
              setequal(Assays(so2), c("RNA","ADT")))                              # 2nd feature space
    chk("cbmc (CITE-seq)", ds, so2)
  }
}
if (have("scRNAseq")) {
  suppressMessages({library(scRNAseq); library(SingleCellExperiment)})
  sce <- tryCatch(suppressMessages(ZeiselBrainData()), error = function(e) NULL)
  if (!is.null(sce)) {                           # real SCE with ERCC/repeat spike-in altExps (multimodal)
    ds <- read_sce(sce); sce2 <- write_sce(ds)
    stopifnot("ERCC" %in% names(ds$axes), identical(ds$axes$ERCC$role, "feature"),  # altExps captured as
              "ERCC" %in% SingleCellExperiment::altExpNames(sce2))                   # feature spaces
    chk("ZeiselBrain (SCE)", ds, sce2)
  }
}
cat(sprintf("  [R] real corpus: %d real published object(s) round-tripped + losses recorded\n", n))
' 2>&1 | grep -E "^  \[(R|skip)\]|Error|stop" | grep -vE "deprecat"
