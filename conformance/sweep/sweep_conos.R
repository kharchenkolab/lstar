# Sweep real local Conos objects (the user's own research data) through write_conos -> L*.
.libPaths(c("/home/pkharchenko/p21/lstar/.Rlib", .libPaths()))
suppressMessages({library(conos); library(lstar)})
cands <- commandArgs(trailingOnly = TRUE)
if (!length(cands)) cands <- c("/home/pkharchenko/igor/adrenal/conI.rds", "/home/pkharchenko/cathy/katya/con.rds")
rep <- file("/tmp/sweep_conos.tsv", "w"); cat("object\tstatus\tsamples\taxes\tfields\tnote\n", file=rep)
ok <- 0; fail <- 0
for (f in cands) {
  if (!file.exists(f) || file.size(f) > 1.5e9) { cat(sprintf("%s\tSKIP\t\t\t\t(absent/too-big)\n", basename(f)), file=rep); next }
  r <- tryCatch({ co <- readRDS(f); ds <- suppressWarnings(write_conos(co))
    p <- tempfile(fileext=".zarr"); lstar_write(ds, p)
    list(s="PASS", n=length(co$samples), a=length(ds$axes), fl=length(ds$fields), note="") },
    error = function(e) list(s="FAIL", n="", a="", fl="", note=substr(conditionMessage(e), 1, 80)))
  if (r$s == "PASS") ok <- ok + 1 else fail <- fail + 1
  cat(sprintf("%s\t%s\t%s\t%s\t%s\t%s\n", basename(f), r$s, r$n, r$a, r$fl, r$note), file=rep); flush(rep)
}
close(rep); cat(sprintf("conos: %d PASS / %d FAIL\n", ok, fail))
