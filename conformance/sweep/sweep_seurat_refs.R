# Sweep the SeuratData *disk* datasets (Azimuth references) -- they have no lazy-load `data()` items, so
# the main sweep misses them; load via LoadData(). These are the SCTAssay + multi-reduction reference
# atlases (some CITE-seq, with ADT + spca), so they add the SCT/reference coverage the lazy datasets lack.
# Resilient + gc per ref (some are large). Appends to /tmp/sweep_seurat.tsv.
.libPaths(c("/home/pkharchenko/p21/lstar/.Rlib", .libPaths()))
suppressMessages({library(Seurat); library(SeuratObject); library(SeuratData); library(lstar)})

refs <- grep("ref$", tryCatch(InstalledData()$Dataset, error = function(e) character(0)), value = TRUE)
rep <- file("/tmp/sweep_seurat_refs.tsv", "w")
cat("dataset\tobject\tstatus\tassays\tassay_classes\taxes\tnote\n", file = rep)
ok <- 0; fail <- 0; skip <- 0
for (d in refs) {
  cat(sprintf("[ref] %-16s ... ", d)); flush.console()
  # Load FIRST and on its own: Azimuth references are disk datasets needing the Azimuth loader, so a
  # load failure is a SKIP (load-dep), NOT a profile FAIL. Only a read_seurat/write_seurat error is a bug.
  so <- tryCatch(suppressWarnings(suppressMessages(UpdateSeuratObject(LoadData(d)))),
                 error = function(e) structure(list(), loaderr = substr(conditionMessage(e), 1, 90)))
  if (!inherits(so, "Seurat")) {
    r <- list(s = "SKIP", as = "", cl = "", a = "", n = paste0("load-dep: ", attr(so, "loaderr")))
  } else r <- tryCatch({
    ds <- read_seurat(so); so2 <- write_seurat(ds)
    p <- tempfile(fileext = ".zarr"); lstar_write(ds, p); unlink(p, recursive = TRUE)
    cls <- paste(sapply(Assays(so), function(a) class(so[[a]])[1]), collapse = "+")
    list(s = "PASS", as = paste(Assays(so), collapse = "+"), cl = cls,
         a = length(ds$axes), n = paste(head(ds$dropped, 2), collapse = ";"))
  }, error = function(e) list(s = "FAIL", as = "", cl = "", a = "", n = substr(conditionMessage(e), 1, 90)))
  if (r$s == "PASS") ok <- ok + 1 else if (r$s == "SKIP") skip <- skip + 1 else fail <- fail + 1
  cat(r$s, "\n"); flush.console()
  cat(sprintf("%s\t%s\t%s\t%s\t%s\t%s\n", d, d, r$s, r$as, r$cl, r$a, r$n), file = rep); flush(rep)
  rm(list = intersect(c("so", "so2", "ds"), ls())); gc(verbose = FALSE)
}
close(rep)
cat(sprintf("\nSeuratData refs sweep: %d PASS / %d FAIL / %d SKIP(load-dep) across %d Azimuth references\n",
            ok, fail, skip, length(refs)))
