#!/usr/bin/env bash
# Seurat Tier-1 extras: a DimReduc's per-dim standard deviation (`Stdev`) is captured as a measure over
# the reduction's coordinate axis (the Seurat analogue of anndata's PCA variance), and the active
# identity (`Idents`) is captured as a flagged `ident` factor field -- both restored on write-back.
set -e
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
RLIB="$ROOT/.Rlib"

Rscript -e '.libPaths(c("'"$RLIB"'", .libPaths())); suppressMessages({library(lstar); library(SeuratObject)})
set.seed(1); n <- 40; g <- 20
m <- matrix(rpois(n * g, 2), g, n); rownames(m) <- paste0("g", 1:g); colnames(m) <- paste0("c", 1:n)
so <- CreateSeuratObject(counts = m)
emb <- matrix(rnorm(n * 5), n, 5); rownames(emb) <- colnames(m); colnames(emb) <- paste0("PC_", 1:5)
so[["pca"]] <- CreateDimReducObject(embeddings = emb, key = "PC_", assay = "RNA", stdev = sqrt(5:1))
Idents(so) <- factor(rep(c("A", "B", "C", "D"), length.out = n), levels = c("A", "B", "C", "D"))

ds <- read_seurat(so)
stopifnot("pca_stdev" %in% names(ds$fields), identical(as.character(ds$fields$pca_stdev$span), "pca"),
          isTRUE(all.equal(as.numeric(ds$fields$pca_stdev$values), sqrt(5:1))),
          "ident" %in% names(ds$fields), ds$fields$ident$subtype == "active_ident",
          identical(ds$axes$ident$role, "factor"))
so2 <- write_seurat(ds)
stopifnot(isTRUE(all.equal(as.numeric(Stdev(so2[["pca"]])), sqrt(5:1))),
          identical(as.character(Idents(so2)), as.character(Idents(so))))
cat("  [R] DimReduc stdev (measure over pca axis) + active Idents captured and restored\n")' \
  2>&1 | grep -vE "^Warning|deprecat|masked|following object|Attaching|^$|^package|loaded|conflicts|^The following"
