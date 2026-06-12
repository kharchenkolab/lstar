#!/usr/bin/env bash
# Seurat version/class variety: "a Seurat object" is several structurally-different classes. These
# small fixtures are constructed by Seurat's own constructors (real classes, deterministic) -- grounded
# in an analysis of real objects -- and each must round-trip through the profile + produce a store that
# validates. Covers: v3/v4 `Assay`, v5 `Assay5`, a v5 *split* (integration) object (reads as a
# collection, re-splits on write-back), an `SCTAssay`, and a multimodal RNA+ADT object (the second
# assay is a Tier-3 gap -> must be *recorded* as dropped, never silently lost).
set -e
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
RLIB="$ROOT/.Rlib"

Rscript -e '.libPaths(c("'"$RLIB"'", .libPaths())); suppressMessages({library(Seurat); library(SeuratObject); library(lstar)})
set.seed(1); ng <- 120; nc <- 60
m <- matrix(rpois(ng*nc, 3), ng, nc, dimnames = list(paste0("Gene",1:ng), paste0("Cell",1:nc)))
emb <- matrix(rnorm(nc*5), nc, 5, dimnames = list(colnames(m), paste0("PC_",1:5)))
withpca <- function(so) { so[["pca"]] <- CreateDimReducObject(emb, key="PC_", assay=DefaultAssay(so), stdev=sqrt(5:1)); Idents(so) <- factor(rep(c("A","B"), length.out=nc)); so }

rt <- function(tag, so, checks) {
  ds <- read_seurat(so); so2 <- write_seurat(ds); stopifnot(checks(ds, so2)); cat(sprintf("  [R] %-22s OK\n", tag))
  p <- file.path("/tmp", paste0("sv_", gsub("[^a-z0-9]","", tolower(tag)), ".lstar.zarr"))
  if (dir.exists(p)) unlink(p, recursive=TRUE); lstar_write(ds, p); p
}
stores <- c()

options(Seurat.object.assay.version = "v3")
stores <- c(stores, rt("v3 Assay", withpca(CreateSeuratObject(m)), function(ds, so2)
  inherits(so2[["RNA"]], "Assay") && identical(ds$kind, "sample") && "pca" %in% Reductions(so2) &&
  any(grepl("assay@RNA:Assay$", ds$profiles))))         # version tracked: v3 Assay class recorded
options(Seurat.object.assay.version = "v5")
stores <- c(stores, rt("v5 Assay5", withpca(CreateSeuratObject(m)), function(ds, so2)
  inherits(so2[["RNA"]], "Assay5") && "pca" %in% Reductions(so2) &&
  any(grepl("assay@RNA:Assay5", ds$profiles))))         # version tracked: v5 Assay5 class recorded

so5 <- CreateSeuratObject(m); so5$samp <- rep(c("s1","s2"), each = nc/2)
so5[["RNA"]] <- split(so5[["RNA"]], f = so5$samp)
stores <- c(stores, rt("v5 split (collection)", so5, function(ds, so2)
  identical(ds$kind, "collection") && setequal(Layers(so2[["RNA"]]), c("counts.s1","counts.s2"))))

sct <- tryCatch(suppressWarnings(suppressMessages(SCTransform(CreateSeuratObject(m), verbose = FALSE))),
                error = function(e) NULL)
if (!is.null(sct)) stores <- c(stores, rt("SCTAssay", sct, function(ds, so2) length(ds$fields) >= 3)) else
  cat("  [R] SCTAssay               SKIP (SCTransform unavailable)\n")

adt <- matrix(rpois(8*nc, 10), 8, nc, dimnames = list(paste0("ADT",1:8), colnames(m)))
mm <- CreateSeuratObject(m); mm[["ADT"]] <- CreateAssay5Object(counts = adt)
stores <- c(stores, rt("multimodal RNA+ADT", mm, function(ds, so2)   # ADT captured as a 2nd feature space
  "ADT" %in% names(ds$axes) && identical(ds$axes$ADT$role, "feature") &&
  identical(as.character(ds$fields[["ADT.counts"]]$span), c("cells","ADT")) &&
  setequal(Assays(so2), c("RNA","ADT")) &&
  isTRUE(all.equal(as.matrix(LayerData(so2[["ADT"]],"counts")), adt, check.attributes = FALSE))))
cat(paste(stores, collapse="\n"), "\n", file = "/tmp/sv_stores.txt")
' 2>&1 | grep -E "^  \[R\]"

# every version fixture produced a store that validates in Python (cross-language)
PYTHONPATH="$ROOT/python/src" python3 - <<'PY'
import warnings; warnings.filterwarnings("ignore")
import os, tempfile, lstar
stores = [s for s in open("/tmp/sv_stores.txt").read().split() if s]
for p in stores:
    errs = [i for i in lstar.validate(lstar.read(p)) if i.startswith("ERROR")]
    assert not errs, (p, errs)
print("  [py] all %d Seurat version-variant stores validate clean" % len(stores))
PY
