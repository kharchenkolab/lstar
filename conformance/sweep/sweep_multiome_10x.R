# Sweep real 10x multiome (RNA+ATAC) .h5 files -> a Seurat object with an RNA Assay + a Signac
# ChromatinAssay -> read_seurat -> write_seurat -> lstar_write. This is the multi-omics ATAC case
# beyond SeuratData's pbmcMultiome: fresh public 10x Cell Ranger ARC output (PBMC granulocyte 3k,
# human brain 3k), so the ChromatinAssay carries REAL peak coordinates (chr-start-end in the feature
# names) that read_seurat types as feature fields (seqnames/start/end) over the peaks axis.
#
# Run: Rscript conformance/sweep/sweep_multiome_10x.R   (writes /tmp/sweep_multiome_10x.tsv)
.libPaths(c("/home/pkharchenko/p21/lstar/.Rlib", .libPaths()))
suppressMessages({library(Seurat); library(SeuratObject); library(Signac); library(lstar)})

files <- list(
  pbmc_granulocyte_3k = "testdata/multiome_10x/pbmc_granulocyte_sorted_3k.h5",
  human_brain_3k      = "testdata/multiome_10x/human_brain_3k.h5"
)
rep <- file("/tmp/sweep_multiome_10x.tsv", "w")
cat("dataset\tstatus\tassays\tassay_classes\tnfeat\taxes\tfields\tnote\n", file = rep)
ok <- 0; fail <- 0

# Build a peaks GRanges from "chr-start-end" feature names, robust to chr names with a '-' (e.g. none here).
make_chromatin <- function(counts) {
  rn <- rownames(counts)
  # Cell Ranger ARC peaks are named "chrN:start-end" OR "chrN-start-end" depending on version.
  norm <- gsub(":", "-", rn)
  parts <- regmatches(norm, regexec("^(.*)-([0-9]+)-([0-9]+)$", norm))
  ok_rows <- vapply(parts, function(p) length(p) == 4L, logical(1))
  counts <- counts[ok_rows, , drop = FALSE]; parts <- parts[ok_rows]
  gr <- GenomicRanges::GRanges(
    seqnames = vapply(parts, `[`, "", 2L),
    ranges   = IRanges::IRanges(as.integer(vapply(parts, `[`, "", 3L)),
                                as.integer(vapply(parts, `[`, "", 4L))))
  Signac::CreateChromatinAssay(counts = counts, ranges = gr, sep = c(":", "-"))
}

for (nm in names(files)) {
  fp <- files[[nm]]
  cat(sprintf("[multiome10x] %-22s ... ", nm)); flush.console()
  r <- tryCatch({
    if (!file.exists(fp)) stop("missing: ", fp)
    d <- suppressWarnings(Read10X_h5(fp))
    rna <- d[["Gene Expression"]]; pk <- d[["Peaks"]]
    so <- CreateSeuratObject(counts = rna, assay = "RNA")
    so[["ATAC"]] <- make_chromatin(pk)
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
cat(sprintf("10x multiome sweep: %d PASS / %d FAIL\n", ok, fail))
