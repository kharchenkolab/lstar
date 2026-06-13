# Sweep real 10x CITE-seq (RNA+ADT) .h5 files -> a Seurat object with an RNA Assay + an ADT Assay ->
# read_seurat -> write_seurat -> lstar_write. Public 10x feature-barcode output where "Antibody Capture"
# features are the ADT panel (surface proteins) -> the multimodal RNA+ADT case beyond cbmc/bmcite, on
# fresh tissue/panel variety (5k PBMC, 1k PBMC, 10k MALT lymphoma).
#
# Run: Rscript conformance/sweep/sweep_citeseq_10x.R   (writes /tmp/sweep_citeseq_10x.tsv)
.libPaths(c("/home/pkharchenko/p21/lstar/.Rlib", .libPaths()))
suppressMessages({library(Seurat); library(SeuratObject); library(lstar)})

files <- list(
  pbmc_5k_protein = "testdata/citeseq_10x/5k_pbmc_protein_v3.h5",
  pbmc_1k_protein = "testdata/citeseq_10x/pbmc_1k_protein_v3.h5",
  malt_10k_protein = "testdata/citeseq_10x/malt_10k_protein_v3.h5"
)
rep <- file("/tmp/sweep_citeseq_10x.tsv", "w")
cat("dataset\tstatus\tassays\tassay_classes\tnfeat\taxes\tfields\tnote\n", file = rep)
ok <- 0; fail <- 0

for (nm in names(files)) {
  fp <- files[[nm]]
  cat(sprintf("[cite10x] %-18s ... ", nm)); flush.console()
  r <- tryCatch({
    if (!file.exists(fp)) stop("missing: ", fp)
    d <- suppressWarnings(Read10X_h5(fp))
    # d is a named list with "Gene Expression" + "Antibody Capture"
    rna <- d[["Gene Expression"]]; adt <- d[["Antibody Capture"]]
    so <- CreateSeuratObject(counts = rna, assay = "RNA")
    so[["ADT"]] <- CreateAssayObject(counts = adt)
    cls <- paste(sapply(Assays(so), function(a) class(so[[a]])[1]), collapse = "+")
    nf  <- paste(sapply(Assays(so), function(a) nrow(so[[a]])), collapse = "+")
    ds <- read_seurat(so); so2 <- write_seurat(ds)              # write-back is the read_seurat round-trip
    errs <- character(0)
    p <- tempfile(fileext = ".zarr"); lstar_write(ds, p); unlink(p, recursive = TRUE)  # write IS validation
    list(s = if (length(errs)) "VALIDATE-ERR" else "PASS",
         as = paste(Assays(so), collapse = "+"), cl = cls, nf = nf,
         a = length(ds$axes), fld = length(ds$fields),
         n = if (length(errs)) substr(errs[1], 1, 90) else paste(head(ds$dropped, 3), collapse = ";"))
  }, error = function(e) list(s = "FAIL", as = "", cl = "", nf = "", a = "", fld = "",
                              n = substr(conditionMessage(e), 1, 100)))
  if (r$s == "PASS") ok <- ok + 1 else fail <- fail + 1
  cat(r$s, "\n"); flush.console()
  cat(sprintf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", nm, r$s, r$as, r$cl, r$nf, r$a, r$fld, r$n),
      file = rep); flush(rep)
}
close(rep)
cat(sprintf("10x CITE-seq sweep: %d PASS / %d FAIL\n", ok, fail))
