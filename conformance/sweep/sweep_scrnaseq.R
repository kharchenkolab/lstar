.libPaths(c("/home/pkharchenko/p21/lstar/.Rlib", .libPaths()))
suppressMessages({library(scRNAseq); library(SingleCellExperiment); library(lstar)})
t <- listDatasets(); rep <- file("/tmp/sweep_scrnaseq.tsv","w")
cat("dataset\tstatus\tfields\taxes\tnote\n", file=rep)
ok<-0; fail<-0; loaderr<-0
for (i in seq_len(min(nrow(t), 30))) {
  nm <- t$Call[i]
  r <- tryCatch({ sce <- suppressMessages(eval(parse(text=nm)))
    if (is.list(sce) && !is(sce,"SummarizedExperiment")) sce <- sce[[1]]
    if (ncol(sce) > 80000) stop("too-big-skip")
    ds <- read_sce(sce); sce2 <- write_sce(ds); p<-tempfile(fileext=".zarr"); lstar_write(ds,p)
    list(s="PASS", f=length(ds$fields), a=length(ds$axes), n="") },
    error=function(e) list(s=if(grepl("there is no package|could not find|too-big|namespace",conditionMessage(e)))"LOADERR" else "FAIL",
                           f="",a="",n=substr(conditionMessage(e),1,90)))
  if(r$s=="PASS")ok<-ok+1 else if(r$s=="FAIL")fail<-fail+1 else loaderr<-loaderr+1
  cat(sprintf("%s\t%s\t%s\t%s\t%s\n", sub("\\(\\)","",nm), r$s, r$f, r$a, r$n), file=rep); flush(rep)
}
close(rep); cat(sprintf("scRNAseq: %d PASS / %d profile-FAIL / %d load-skip\n", ok, fail, loaderr))
