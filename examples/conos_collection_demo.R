#!/usr/bin/env Rscript
# Conos collection -> L* -> round-trip, on the real two-sample integration object.
#
# This is the headline test of the L* *collection* model: a Conos object holds two samples
# with DIFFERENT gene sets, their own PCA, plus a joint embedding / clustering / graph. An
# aligned cells x genes matrix cannot represent that without dropping the per-sample structure;
# L* keeps it (samples axis + per-sample axes + a union cells axis for the joint layer).

suppressMessages({ library(Matrix); library(lstar) })

RDS <- "/home/pkharchenko/p21/pagoda2/misc/conos_two_sample_integration/conos_two_sample.rds"
OUT <- "/tmp/conos_collection.lstar.zarr"

human <- function(n) { u <- c("B","KB","MB","GB"); i <- 1
  while (n >= 1024 && i < length(u)) { n <- n/1024; i <- i+1 }; sprintf("%.1f %s", n, u[i]) }
dir_size <- function(p) sum(file.info(list.files(p, recursive=TRUE, full.names=TRUE))$size, na.rm=TRUE)

t0 <- Sys.time(); co <- readRDS(RDS)
cat(sprintf("loaded Conos: %d samples (%s)  [%.1fs]\n",
            length(co$samples), paste(names(co$samples), collapse=", "),
            as.numeric(Sys.time()-t0, units="secs")))

t0 <- Sys.time(); ds <- write_conos(co)
cat(sprintf("\nwrite_conos -> L* collection  [%.1fs]\n", as.numeric(Sys.time()-t0, units="secs")))
print(ds)

# Show the heterogeneity the collection preserves: per-sample gene sets differ.
gene_axes <- grep("^genes\\.", names(ds$axes), value=TRUE)
gs <- lapply(gene_axes, function(a) ds$axes[[a]]$labels)
if (length(gs) >= 2) {
  inter <- length(Reduce(intersect, gs)); uni <- length(Reduce(union, gs))
  cat(sprintf("\nper-sample gene sets: sizes %s; shared %d / union %d (%.0f%% overlap)\n",
              paste(lengths(gs), collapse=", "), inter, uni, 100*inter/uni))
}

t0 <- Sys.time(); lstar_write(ds, OUT)
cat(sprintf("\nlstar_write -> zarr  [%.1fs]  store=%s\n",
            as.numeric(Sys.time()-t0, units="secs"), human(dir_size(OUT))))

t0 <- Sys.time(); ds2 <- lstar_read(OUT)
cat(sprintf("lstar_read  <- zarr  [%.1fs]\n", as.numeric(Sys.time()-t0, units="secs")))

# Fidelity: per-sample counts nnz + sums, joint graph edges, joint embedding shape.
cat("\nfidelity (original -> round-trip):\n")
ok <- TRUE
for (sn in names(co$samples)) {
  nm <- paste0("counts.", sn)
  a <- ds$fields[[nm]]$values; b <- ds2$fields[[nm]]$values
  same <- nnzero(a) == nnzero(b) && abs(sum(a@x) - sum(b@x)) < 1e-6 * abs(sum(a@x))
  ok <- ok && same
  cat(sprintf("  %-16s nnz %d->%d  sum %.6g->%.6g  %s\n", nm,
              nnzero(a), nnzero(b), sum(a@x), sum(b@x), if (same) "OK" else "MISMATCH"))
}
G <- ds2$fields$graph
if (!is.null(G)) {
  gg <- G$values
  cat(sprintf("  %-16s %dx%d  nnz(edges) %d  (relation over cells x cells)  %s\n",
              "graph", nrow(gg), ncol(gg), nnzero(gg),
              if (nnzero(gg) == nnzero(ds$fields$graph$values)) "OK" else "MISMATCH"))
}
E <- ds2$fields$embedding
if (!is.null(E)) cat(sprintf("  %-16s %s  (joint embedding over union cells)\n",
                             "embedding", paste(dim(as.matrix(E$values)), collapse="x")))
cat(sprintf("  %-16s %s\n", "sample label",
            paste(names(table(ds2$fields$sample$values)),
                  table(ds2$fields$sample$values), sep="=", collapse="  ")))

cat(sprintf("\n%s  store=%s\n", if (ok) "PASS: collection round-trips faithfully" else "FAIL", human(dir_size(OUT))))
quit(status = if (ok) 0 else 1)
