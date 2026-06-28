# Native R viewer prep: read a base store, extend_for_viewer, write it out. For conformance/viewer.sh
# leg (d) -- assert R extend_for_viewer == lstar Python extend_for_viewer on the same base store.
#   Rscript viewer_extend_r.R <base_store> <out_store>
args <- commandArgs(trailingOnly = TRUE); base <- args[1]; out <- args[2]
rlib <- Sys.getenv("LSTAR_RLIB", ""); if (nzchar(rlib)) .libPaths(c(normalizePath(rlib), .libPaths()))
suppressMessages(library(lstar))
ds <- lstar_read(base)
ds <- extend_for_viewer(ds)
lstar_write(ds, out)
cat("R extend_for_viewer ->", out, "\n")
