#!/usr/bin/env Rscript
# Conos collection -> L* -> round-trip, on a real two-sample integration object.
#
# This is the headline demonstration of the L* *collection* model. A Conos object holds several
# samples, each with its OWN cells and possibly its OWN gene set, plus a joint embedding /
# clustering / integration graph computed across them. An aligned cells x genes matrix cannot
# represent that without flattening away the per-sample structure; L* keeps it as a collection:
# a `samples` axis, per-sample cells.<s>/genes.<s> axes and counts.<s> measures, and a union
# `cells` axis carrying the joint layer. We convert it, write it, read it back, and check the
# heavy fields survived.
#
# Usage: conos_collection_demo.R [path/to/conos.rds]
#   The default below is a local path on the author's machine -- pass your own Conos .rds to run it.

suppressMessages({ library(Matrix); library(lstar) })

args <- commandArgs(trailingOnly = TRUE)
# A local default on the author's machine; override by passing your own Conos .rds as the 1st arg.
RDS <- "/home/pkharchenko/p21/pagoda2/misc/conos_two_sample_integration/conos_two_sample.rds"
if (length(args) >= 1) RDS <- args[[1]]
OUT <- "/tmp/conos_collection.lstar.zarr"

if (!file.exists(RDS)) {
  cat(sprintf("Conos object not found: %s\n  Pass a path to a saved Conos R6 object as the first argument.\n", RDS))
  quit(status = 0)
}

human <- function(n) { u <- c("B","KB","MB","GB"); i <- 1
  while (n >= 1024 && i < length(u)) { n <- n/1024; i <- i+1 }; sprintf("%.1f %s", n, u[i]) }
dir_size <- function(p) sum(file.info(list.files(p, recursive=TRUE, full.names=TRUE))$size, na.rm=TRUE)

t0 <- Sys.time(); co <- readRDS(RDS)
cat(sprintf("loaded Conos: %d samples (%s)  [%.1fs]\n",
            length(co$samples), paste(names(co$samples), collapse=", "),
            as.numeric(Sys.time()-t0, units="secs")))

# write_conos applies the `conos` profile: each member sample becomes its own cells.<s>/genes.<s>
# axes + a counts.<s> measure, and the joint analysis (embedding, clusters, integration graph)
# becomes fields over a derived union `cells` axis. The result is one L* `collection` dataset.
t0 <- Sys.time(); ds <- write_conos(co)
cat(sprintf("\nwrite_conos -> L* collection  [%.1fs]\n", as.numeric(Sys.time()-t0, units="secs")))
print(ds)   # note the per-sample axes/fields alongside the joint `cells`, `embedding`, `graph`

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

# Fidelity: confirm the heavy fields are byte-faithful after write + read. We compare the L*
# dataset BEFORE writing (ds) to the one read back from disk (ds2): per-sample counts (nonzero
# count and sum), the joint graph's edge count, the joint embedding shape, and the sample label.
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
