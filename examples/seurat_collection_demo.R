#!/usr/bin/env Rscript
# Seurat v5 collection -> L*. A Seurat v5 integration workflow holds the samples *unintegrated*
# as a split Assay5: split(assay, f = sample) produces per-sample layers (counts.<sample>) that
# each cover only their sample's cells. That is a collection, not one aligned matrix, and L*
# ingests it as such: a samples axis, per-sample cells.<s> axes + counts.<s> measures, a union
# cells axis with a `sample` label, and any integrated reduction over the union cells.

suppressMessages({ library(lstar); library(SeuratObject); library(Matrix) })

OUT <- "/tmp/seurat_collection.lstar.zarr"
set.seed(1)
ng <- 200; nc <- 600
m <- matrix(rpois(ng * nc, 0.5), ng, nc, dimnames = list(paste0("g", 1:ng), paste0("c", 1:nc)))
obj <- CreateSeuratObject(counts = as(m, "dgCMatrix"))
obj$sample <- rep(c("donorA", "donorB", "donorC"), length.out = nc)

cat("joined assay layers: ", paste(Layers(obj[["RNA"]]), collapse = ", "), "\n")
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$sample)         # -> a collection of per-sample layers
cat("split  assay layers: ", paste(Layers(obj[["RNA"]]), collapse = ", "), "\n\n")

ds <- read_seurat(obj)
cat("read_seurat -> L*:\n"); print(ds)
stopifnot(ds$kind == "collection")

lstar_write(ds, OUT)
ds2 <- lstar_read(OUT)
ok <- TRUE
for (sn in c("donorA", "donorB", "donorC")) {
  nm <- paste0("counts.", sn)
  a <- ds$fields[[nm]]$values; b <- ds2$fields[[nm]]$values
  same <- nnzero(a) == nnzero(b) && abs(sum(a@x) - sum(b@x)) < 1e-9 * abs(sum(a@x) + 1)
  ok <- ok && same
  cat(sprintf("  %-16s %dx%d  nnz %d->%d  %s\n", nm, nrow(b), ncol(b),
              nnzero(a), nnzero(b), if (same) "OK" else "MISMATCH"))
}
cat(sprintf("\n%s\n", if (ok) "PASS: Seurat v5 collection round-trips faithfully" else "FAIL"))
quit(status = if (ok) 0 else 1)
