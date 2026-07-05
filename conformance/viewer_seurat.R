# Seurat -> extend_for_viewer conformance (the previously-untested seam). Builds a realistic Seurat
# object with the metadata that stresses the viewer's grouping/embedding/basis detection -- factor
# clusterings, a UMAP DimReduc, active Idents == the clustering, and a LOGICAL QC flag -- then runs
# read_seurat -> extend_for_viewer and asserts the viewer store is sane. In CI: a synthetic Seurat
# (no downloads). Locally: also runs on a real SeuratData object when available.
#   Rscript viewer_seurat.R <out_store>   [<seurat_dataset_name>]
suppressWarnings({
  args <- commandArgs(trailingOnly = TRUE); out <- args[1]; real <- if (length(args) >= 2) args[2] else ""
  rlib <- Sys.getenv("LSTAR_RLIB", ""); if (nzchar(rlib)) .libPaths(c(normalizePath(rlib), .libPaths()))
  ok <- suppressMessages(requireNamespace("SeuratObject", quietly = TRUE) &&
                         requireNamespace("Seurat", quietly = TRUE) &&
                         requireNamespace("lstar", quietly = TRUE))
})
if (!ok) { cat("  [skip] SeuratObject/Seurat/lstar not installed — skipping Seurat->viewer\n"); quit(status = 0) }
suppressMessages({ library(SeuratObject); library(Seurat); library(lstar) })

fail <- function(msg) { cat("  FAIL:", msg, "\n"); quit(status = 1) }

build_synthetic <- function() {
  set.seed(1); nc <- 200L; ng <- 60L
  cnt <- matrix(rpois(nc * ng, 2), ng, nc, dimnames = list(paste0("g", 1:ng), paste0("c", 1:nc)))
  so <- CreateSeuratObject(counts = cnt)
  so$seurat_clusters <- factor(paste0("k", (0:(nc - 1)) %% 5))     # clustering (factor)
  so$cell_type       <- factor(c("T", "B", "NK", "Mono", "DC")[(0:(nc - 1)) %% 5 + 1])
  so$qc_kept         <- (0:(nc - 1)) %% 7 != 0                     # a LOGICAL QC flag (must NOT become a grouping)
  Idents(so) <- "seurat_clusters"                                 # active idents == the clustering (common case)
  emb <- matrix(rnorm(nc * 2), nc, 2, dimnames = list(colnames(so), c("UMAP_1", "UMAP_2")))
  so[["umap"]] <- CreateDimReducObject(embeddings = emb, key = "UMAP_", assay = "RNA")
  so
}

check <- function(so, label) {
  cat("  --", label, "--\n")
  ds <- read_seurat(so)
  qc <- ds$fields[["qc_kept"]]
  if (!is.null(qc) && !is.logical(qc$values))
    fail(sprintf("a logical QC column came back as %s (should stay logical)", class(qc$values)[1]))
  dv <- extend_for_viewer(ds)
  gp <- sub("^stats_", "", sub("_sum$", "", grep("^stats_.*_sum$", names(dv$fields), value = TRUE)))
  cat("    groupings:", paste(gp, collapse = ", "), "\n")
  prim <- dv$fields[["counts_cellmajor_order"]]$provenance$group
  cat("    primary (reorder key):", prim, "\n")
  # (a) a boolean QC flag must NOT be a grouping (fix: read_seurat keeps it logical -> excluded)
  if ("qc_kept" %in% gp) fail("qc_kept (a boolean QC flag) was detected as a viewer grouping")
  # (a2) the active-idents mirror must NOT duplicate the clustering as a separate grouping
  if ("ident" %in% gp) fail("active idents ('ident') duplicated the clustering as a viewer grouping")
  # (b) the viewer must open on a real clustering, not on a QC/near-boolean field
  if (is.null(prim) || !nzchar(prim)) fail("no primary grouping keyed the reorder")
  if (identical(prim, "qc_kept")) fail("the viewer opened on the boolean qc_kept")
  # (c) the counts_cellmajor navigators are produced
  if (is.null(dv$fields[["counts_cellmajor"]])) fail("counts_cellmajor was not produced")
  if (is.null(dv$fields[["counts_cellmajor_order"]])) fail("counts_cellmajor_order (hybrid reorder) missing")
  cat("    OK: boolean excluded, primary =", prim, ", counts_cellmajor present\n")
  dv
}

dv <- check(build_synthetic(), "synthetic Seurat")
if (nzchar(real) && requireNamespace("SeuratData", quietly = TRUE)) {
  so <- tryCatch(suppressMessages(SeuratObject::UpdateSeuratObject(SeuratData::LoadData(real))),
                 error = function(e) { cat("    [skip] real", real, "unavailable:", conditionMessage(e), "\n"); NULL })
  if (!is.null(so)) {
    so <- subset(so, cells = colnames(so)[stats::complete.cases(so[[]])])
    so$qc_kept <- so$nFeature_RNA > stats::median(so$nFeature_RNA)   # inject a logical QC flag
    check(so, paste("real Seurat:", real))
  }
}

lstar_write(dv, out)   # hand the synthetic viewer store to the Python validator
cat("  wrote", out, "\n")
