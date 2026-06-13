# Sweep installed SeuratData objects through read_seurat -> write_seurat.
.libPaths(c("/home/pkharchenko/p21/lstar/.Rlib", .libPaths()))
suppressMessages({library(Seurat); library(SeuratObject); library(SeuratData); library(lstar)})
# Each installed SeuratData dataset is its OWN package `<d>.SeuratData`, and its objects (e.g. `pbmc3k`
# and the lazily-loaded `pbmc3k.final`) are that package's data items -- NOT items of the SeuratData
# meta-package. Enumerate per-package, else we find 0 objects (the original bug).
inst <- tryCatch(InstalledData()$Dataset, error=function(e) character(0))
rep <- file("/tmp/sweep_seurat.tsv","w"); cat("dataset\tobject\tstatus\tassays\taxes\tnote\n", file=rep); ok<-0; fail<-0
for (d in inst) {
  pkg <- paste0(d, ".SeuratData")
  objs <- tryCatch(sub("\\s.*$", "", data(package=pkg)$results[,"Item"]), error=function(e) character(0))
  for (on in objs) {
    r <- tryCatch({ suppressWarnings(suppressMessages(data(list=on, package=pkg))); so <- UpdateSeuratObject(get(on))
      ds <- read_seurat(so); so2 <- write_seurat(ds); p<-tempfile(fileext=".zarr"); lstar_write(ds,p)
      list(s="PASS", as=paste(Assays(so),collapse="+"), a=length(ds$axes), n=paste(head(ds$dropped,2),collapse=";")) },
      error=function(e) list(s="FAIL", as="", a="", n=substr(conditionMessage(e),1,80)))
    if(r$s=="PASS")ok<-ok+1 else fail<-fail+1
    cat(sprintf("%s\t%s\t%s\t%s\t%s\t%s\n", d, on, r$s, r$as, r$a, r$n), file=rep); flush(rep)
  }
}
close(rep); cat(sprintf("SeuratData: %d PASS / %d FAIL (of %d installed)\n", ok, fail, length(inst)))
