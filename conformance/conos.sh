#!/usr/bin/env bash
# Conos (graph-only integration) collection conversions. A Conos object integrates samples in GRAPH space
# -- there is NO corrected/batch-aligned expression matrix (computing one is expensive) -- so the joint
# layer is a graph + embedding + clustering over the union cells, atop per-sample raw counts with gene
# sets that overlap, differ, or are disjoint. This test proves that such a collection (a) round-trips
# Conos <-> L*, and (b) converts to the NATIVE collection formats WITHOUT fabricating a corrected matrix:
#   - Seurat v5: a split assay (per-sample raw layers) + Graphs(joint graph) + DimReduc(embedding) +
#                metadata(clusters, sample). read_seurat() reads it BACK as an L* collection.
#   - AnnData : a single object, X = raw joint counts (== conos getJointCountMatrix), obsp = the graph
#                (also aliased to `connectivities` + uns['neighbors'] so scanpy graph ops run), obsm =
#                embedding, obs = sample + clusters. (AnnData is one matrix, so this is a flattening.)
#
# The headline REAL leg builds an actual Conos object from conos::small_panel.preprocessed; it SKIPS if
# conos/pagoda2 are absent (Suggests-only -- not in CI). A SYNTHETIC divergent-genes leg always runs so CI
# covers the conversion logic (needs Seurat + anndata, both present in the r-cross-format job).
# Origin coverage: Conos-authored (real) | R collection_from-authored (synthetic) -- see conformance/README.md
set -eo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
RLIB="$ROOT/.Rlib"; export LSTAR_RLIB="$RLIB"
CS=/tmp/conos_real.lstar.zarr; SS=/tmp/conos_synth.lstar.zarr

# ── REAL Conos leg (skips if conos/pagoda2 not installed) ────────────────────────────────────────────
Rscript -e '.libPaths(c("'"$RLIB"'", .libPaths()))
have <- function(p) requireNamespace(p, quietly = TRUE)
if (!have("conos") || !have("pagoda2")) { cat("  [skip] real Conos leg: conos/pagoda2 absent (local-only)\n"); quit(status = 0) }
suppressMessages({library(conos); library(pagoda2); library(Matrix); library(SeuratObject); library(lstar)})
log <- file("/tmp/conos_chatter.log", "w")
con <- Conos$new(small_panel.preprocessed, n.cores = 1)
sink(log, type = "output")
suppressWarnings(con$runGraph(ncomps = 25))
suppressWarnings(con$runClustering(method = leiden.community))
suppressWarnings(con$runEmbedding(alpha = 0.001, sgd_batched = 1e8))
sink(type = "output"); close(log)

ds <- write_conos(con)                                            # Conos -> L* collection
stopifnot(identical(ds$kind, "collection"))
sn <- names(con$samples)
for (s in sn) stopifnot(!is.null(ds$axes[[paste0("cells.", s)]]), !is.null(ds$fields[[paste0("counts.", s)]]))
stopifnot(!is.null(ds$axes$cells), !is.null(ds$fields$graph),     # the joint layer is graph/embedding/cluster
          !is.null(ds$fields$embedding), identical(ds$fields$sample$subtype, "design"))
cat(sprintf("  [R ] write_conos: %d samples, joint graph(%d edges)+embedding+clusters, NO corrected matrix\n",
            length(sn), Matrix::nnzero(ds$fields$graph$values)))

co2 <- read_conos(ds)                                             # L* -> Conos (round-trip)
stopifnot(inherits(co2, "Conos"), length(co2$samples) == length(sn), !is.null(co2$graph))
cat("  [R ] read_conos round-trip: live Conos object (samples + joint graph restored)\n")

so <- write_seurat(ds)                                            # collection -> Seurat v5 split
stopifnot(inherits(so, "Seurat"), isTRUE(validObject(so)), inherits(so[["RNA"]], "Assay5"))
lyr <- SeuratObject::Layers(so, assay = "RNA")
stopifnot(length(lyr) == length(sn),                             # per-sample raw layers (split assay)
          "graph" %in% SeuratObject::Graphs(so),                # joint graph carried natively
          "embedding" %in% SeuratObject::Reductions(so),        # joint embedding as a DimReduc
          all(c("sample", "leiden") %in% colnames(so[[]])))     # clusters + sample in metadata
back <- read_seurat(so)                                          # native acceptance: reads BACK as collection
stopifnot(identical(back$kind, "collection"), length(grep("^cells\\.", names(back$axes))) == length(sn))
cat(sprintf("  [R ] -> Seurat v5: %d split layers + Graphs + DimReduc + meta; read_seurat -> collection\n", length(lyr)))
saveRDS(sn, "/tmp/conos_sn.rds")
lstar_write(ds, "'"$CS"'")' 2>&1 | grep -E "^  \[(R|skip)"

# Python leg of the real test: only if the store was produced (i.e. the R leg did not skip)
if [ -d "$CS" ]; then
PYTHONPATH="$ROOT/python/src" python3 - "$CS" <<'PY'
import sys, warnings; warnings.filterwarnings("ignore")
import numpy as np, lstar
from lstar import write_anndata, read_anndata
import scanpy as sc
c = lstar.read(sys.argv[1])
a = write_anndata(c)                                             # collection -> single AnnData (flattening)
assert a.X is not None and float(a.X.min()) >= 0                 # X = RAW joint counts (no corrected matrix)
assert {"sample", "leiden"} <= set(a.obs.columns)
assert any(k.startswith("X_") for k in a.obsm) and "connectivities" in a.obsp and "neighbors" in a.uns
drop = set(a.uns.get("lstar/dropped", []))
assert any(d.startswith("pca.") for d in drop), drop          # per-sample latent spaces honestly dropped
sc.tl.umap(a)                                                   # scanpy CONSUMES the conos graph, no extra prep
assert a.obsm["X_umap"].shape[0] == a.n_obs
print(f"  [py] -> AnnData: X={a.shape} raw joint counts, obsp graph+connectivities, obsm embedding, "
      f"obs sample/leiden; scanpy umap ran on the conos graph; dropped per-sample pca")
PY
fi

# ── SYNTHETIC divergent-genes leg (ALWAYS runs; covers the conversion logic in CI) ────────────────────
# A Conos-shaped collection assembled by R collection_from: 3 samples, divergent gene sets, a joint kNN
# graph + embedding + clustering over the union -- exactly write_conos's output shape, but with no Conos
# dependency. Convert to Seurat v5 (R) and AnnData (Python) and assert the structure survives.
Rscript -e '.libPaths(c("'"$RLIB"'", .libPaths())); suppressMessages({library(lstar); library(Matrix); library(SeuratObject)})
mk <- function(s, nc, genes) { d <- list(kind="sample", axes=list(), fields=list())
  d$axes$cells <- list(labels=paste0("c",seq_len(nc)), origin="observed", role="observation")
  d$axes$genes <- list(labels=genes, origin="observed", role="feature")
  d$fields$counts <- list(role="measure", span=c("cells","genes"), state="raw",
    values=as(Matrix::Matrix(matrix(rpois(nc*length(genes),2),nc),sparse=TRUE),"CsparseMatrix"))
  class(d)<-"lstar_dataset"; d }
shared <- paste0("g",1:200)                                      # underscore-free names (Seurat mangles "_")
samp <- list(A=mk("A",80, c(shared,paste0("Aonly",1:30))),
             B=mk("B",100,c(shared,paste0("Bonly",1:20))),
             C=mk("C",70, c(shared,paste0("Conly",1:25))))
ntot <- 80+100+70
g <- as(Matrix::Matrix(Matrix::rsparsematrix(ntot,ntot,0.05,symmetric=TRUE),sparse=TRUE),"CsparseMatrix")
g <- abs(g); diag(g) <- 0
col <- collection_from(samp, joint=list(umap=matrix(rnorm(ntot*2),ntot,2),
                                        graph=g,
                                        clusters=factor(sample(paste0("k",1:4),ntot,replace=TRUE))))
stopifnot(identical(col$kind,"collection"))
so <- write_seurat(col)                                          # divergent-genes collection -> Seurat v5
stopifnot(inherits(so,"Seurat"), isTRUE(validObject(so)), inherits(so[["RNA"]],"Assay5"))
lyr <- SeuratObject::Layers(so, assay="RNA")
ug <- length(union(union(c(shared,paste0("Aonly",1:30)),paste0("Bonly",1:20)),paste0("Conly",1:25)))
stopifnot(length(lyr)==3, nrow(so)==ug,                         # 3 split layers over the UNION gene set
          "graph" %in% SeuratObject::Graphs(so), "umap" %in% SeuratObject::Reductions(so),
          "clusters" %in% colnames(so[[]]))
# divergent genes unioned (absent-in-sample genes are 0 in that layer, present in the union features)
stopifnot("Aonly1" %in% rownames(so), "Conly1" %in% rownames(so))
cat(sprintf("  [R ] synthetic divergent collection -> Seurat v5: 3 layers over %d union genes + graph + umap\n", ug))
lstar_write(col, "'"$SS"'")' 2>&1 | grep -E "^  \[R"

PYTHONPATH="$ROOT/python/src" python3 - "$SS" <<'PY'
import sys, warnings; warnings.filterwarnings("ignore")
import numpy as np, lstar
from lstar import write_anndata
c = lstar.read(sys.argv[1])
assert c.kind == "collection" and "genes" not in c.axes        # per-sample genes.<s>, no shared genes axis
a = write_anndata(c)                                            # -> single union-genes AnnData
ng_union = len({g for s in "ABC" for g in np.asarray(c.axis(f"genes.{s}").labels)})
assert a.shape == (250, ng_union), (a.shape, ng_union)
assert a.X is not None and {"sample", "clusters"} <= set(a.obs.columns)
assert any(k.startswith("X_") for k in a.obsm) and "connectivities" in a.obsp
# the per-sample heterogeneity is unioned into one matrix (an AnnData is one matrix, by construction)
assert "Aonly1" in set(a.var_names) and "Conly1" in set(a.var_names)
print(f"  [py] synthetic divergent collection -> AnnData: X=(250 x {ng_union}) union, obsp graph, obsm embedding, obs sample/clusters")
PY
echo "conos conversion conformance PASSED."
