# Sweep the 10x RNA+ATAC multiome (Signac ChromatinAssay) -- the multi-omics ATAC case. Needs Signac.
.libPaths(c("/home/pkharchenko/p21/lstar/.Rlib", .libPaths()))
suppressMessages({library(Seurat); library(SeuratObject); library(Signac); library(SeuratData); library(lstar)})
if (!("pbmcMultiome" %in% tryCatch(InstalledData()$Dataset, error=function(e) character(0)))) {
  cat("[install] pbmcMultiome ...\n"); flush.console()
  tryCatch(InstallData("pbmcMultiome"), error = function(e) cat("install err:", conditionMessage(e), "\n"))
}
pkg <- "pbmcMultiome.SeuratData"
objs <- tryCatch(sub("\\s.*$", "", data(package = pkg)$results[, "Item"]), error = function(e) character(0))
cat("[multiome] objects:", paste(objs, collapse = ", "), "\n"); flush.console()
rep <- file("/tmp/sweep_multiome.tsv", "w")
cat("object\tstatus\tassays\tassay_classes\tnfeat\taxes\tnote\n", file = rep)
for (on in objs) {
  cat(sprintf("[multiome] %-14s ... ", on)); flush.console()
  r <- tryCatch({
    suppressWarnings(suppressMessages(data(list = on, package = pkg)))
    so <- UpdateSeuratObject(get(on))
    cls <- paste(sapply(Assays(so), function(a) class(so[[a]])[1]), collapse = "+")
    nf  <- paste(sapply(Assays(so), function(a) nrow(so[[a]])), collapse = "+")
    ds <- read_seurat(so); so2 <- write_seurat(ds)
    p <- tempfile(fileext = ".zarr"); lstar_write(ds, p)
    unlink(p, recursive = TRUE)
    list(s = "PASS", as = paste(Assays(so), collapse = "+"), cl = cls, nf = nf,
         a = length(ds$axes), n = paste(head(ds$dropped, 3), collapse = ";"))
  }, error = function(e) list(s = "FAIL", as = "", cl = "", nf = "", a = "", n = substr(conditionMessage(e), 1, 110)))
  cat(r$s, "\n"); flush.console()
  cat(sprintf("%s\t%s\t%s\t%s\t%s\t%s\t%s\n", on, r$s, r$as, r$cl, r$nf, r$a, r$n), file = rep); flush(rep)
}
close(rep); cat("MULTIOME SWEEP DONE\n")
