# Install the SeuratData breadth (all available datasets EXCEPT spatial -- stxBrain/stxKidney/ssHippo
# are deferred to the spatial tier per project decision) and sweep them through read_seurat/write_seurat.
# Heavy + local-only (GBs); resilient per-dataset so one bad install/load doesn't abort the batch.
.libPaths(c("/home/pkharchenko/p21/lstar/.Rlib", .libPaths()))
suppressMessages({library(Seurat); library(SeuratObject); library(SeuratData); library(lstar)})

SPATIAL <- c("stxBrain", "stxKidney", "ssHippo")          # deferred (Visium / Slide-seq) -> spatial tier
avail <- AvailableData()
want <- setdiff(avail$Dataset, SPATIAL)
inst0 <- tryCatch(InstalledData()$Dataset, error = function(e) character(0))
todo <- setdiff(want, inst0)
cat(sprintf("[install] %d datasets target; %d already installed; installing %d\n",
            length(want), length(inst0), length(todo)))
for (d in todo) {
  cat(sprintf("[install] %-18s ... ", d)); flush.console()
  ok <- tryCatch({ suppressWarnings(suppressMessages(InstallData(d))); "ok" },
                 error = function(e) substr(conditionMessage(e), 1, 70))
  cat(ok, "\n"); flush.console()
}

# --- sweep every installed (non-spatial) dataset's objects ---
inst <- setdiff(tryCatch(InstalledData()$Dataset, error = function(e) character(0)), SPATIAL)
rep <- file("/tmp/sweep_seurat.tsv", "w")
cat("dataset\tobject\tstatus\tassays\tassay_classes\taxes\tnote\n", file = rep)
ok <- 0; fail <- 0; nobj <- 0
for (d in inst) {
  pkg <- paste0(d, ".SeuratData")
  objs <- tryCatch(sub("\\s.*$", "", data(package = pkg)$results[, "Item"]), error = function(e) character(0))
  for (on in objs) {
    nobj <- nobj + 1
    r <- tryCatch({
      suppressWarnings(suppressMessages(data(list = on, package = pkg)))
      so <- UpdateSeuratObject(get(on))
      ds <- read_seurat(so); so2 <- write_seurat(ds)
      p <- tempfile(fileext = ".zarr"); lstar_write(ds, p); unlink(p, recursive = TRUE)
      cls <- paste(sapply(Assays(so), function(a) class(so[[a]])[1]), collapse = "+")
      list(s = "PASS", as = paste(Assays(so), collapse = "+"), cl = cls,
           a = length(ds$axes), n = paste(head(ds$dropped, 2), collapse = ";"))
    }, error = function(e) list(s = "FAIL", as = "", cl = "", a = "", n = substr(conditionMessage(e), 1, 80)))
    if (r$s == "PASS") ok <- ok + 1 else fail <- fail + 1
    cat(sprintf("%s\t%s\t%s\t%s\t%s\t%s\t%s\n", d, on, r$s, r$as, r$cl, r$a, r$n), file = rep); flush(rep)
    rm(list = intersect(on, ls())); gc(verbose = FALSE)         # free each object (some are large)
  }
}
close(rep)
cat(sprintf("\nSeuratData sweep: %d PASS / %d FAIL across %d objects (of %d installed datasets)\n",
            ok, fail, nobj, length(inst)))
