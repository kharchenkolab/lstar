# Write a viewer@0.1 store from a mock Pagoda2-shaped object via lstar's profile_pagoda2 path.
# Used by conformance/viewer.sh to check what the R side currently emits against the spec.
#   Rscript viewer_pagoda2.R <out_store> [grouping]
args  <- commandArgs(trailingOnly = TRUE)
out   <- args[1]; grouping <- if (length(args) >= 2) args[2] else "leiden"
rlib  <- Sys.getenv("LSTAR_RLIB", "")
if (nzchar(rlib)) .libPaths(c(normalizePath(rlib), .libPaths()))
suppressMessages(library(lstar))

set.seed(3); nc <- 90L; ng <- 18L
cnt <- as(Matrix::Matrix(rpois(nc * ng, 0.6), nc, ng, sparse = TRUE), "CsparseMatrix")
rownames(cnt) <- paste0("cell", 1:nc); colnames(cnt) <- paste0("g", 1:ng)
emb <- matrix(rnorm(nc * 2), nc, 2); rownames(emb) <- rownames(cnt)
meta <- data.frame(leiden = paste0("c", (0:(nc - 1)) %% 4),
                   cell_type = paste0("t", (0:(nc - 1)) %% 4),
                   row.names = rownames(cnt), stringsAsFactors = FALSE)
p2 <- list(getRawCounts = function(...) cnt, embeddings = list(PCA = list(UMAP = emb)),
           cellMeta = meta, misc = list())

if (dir.exists(out)) unlink(out, recursive = TRUE)
write_pagoda2(p2, out, grouping = grouping)
cat("wrote", out, "\n")
