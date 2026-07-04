# Native R viewer prep: read a base store, extend_for_viewer, write it out. For conformance/viewer.sh
# leg (d) -- assert R extend_for_viewer == lstar Python extend_for_viewer on the same base store.
#   Rscript viewer_extend_r.R <base_store> <out_store>
args <- commandArgs(trailingOnly = TRUE); base <- args[1]; out <- args[2]
basis <- if (length(args) >= 3 && nzchar(args[3])) args[3] else NULL   # "lognorm" for corpus data w/o raw counts
primary <- if (length(args) >= 4 && nzchar(args[4])) args[4] else NULL # the grouping the viewer opens on (parity leg)
rlib <- Sys.getenv("LSTAR_RLIB", ""); if (nzchar(rlib)) .libPaths(c(normalizePath(rlib), .libPaths()))
suppressMessages(library(lstar))
ds <- lstar_read(base)
ds <- extend_for_viewer(ds, basis = basis, primary = primary)
lstar_write(ds, out)
cat("R extend_for_viewer ->", out, "\n")
