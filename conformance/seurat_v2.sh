#!/usr/bin/env bash
# Seurat v2 (pre-Assay) round-trip: a "very old" serialized Seurat object is the lowercase S4 class
# `seurat` -- it predates the Assay/Assay5 classes and the SeuratObject package entirely. SUPPORT.md
# long flagged this as untestable ("can't be synthetically tested without the ancient Seurat 2.x
# package, which won't co-install"). This fixture closes that gap WITHOUT ancient Seurat: it builds a
# structurally-FAITHFUL object from Seurat 2.3.4's *authoritative* S4 class definitions (the exact slot
# layout, lifted from the 2.3.4 source), then -- crucially -- removeClass()es them and readRDS()es the
# object back, so read_seurat() sees exactly what a modern R sees when handed an ancient .rds: an S4
# object whose `seurat` class is UNDEFINED, slots reachable only via attr(). (A genuine Seurat 2.3.4
# object built on the old R-3.6.2 install validates the same path locally -- see sweep notes.)
set -eo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
RLIB="$ROOT/.Rlib"

# Large-ish R program -> run from a FILE, not `Rscript -e` (a big -e expression overflows R's ~8KB
# command-line buffer and is silently ignored). Quoted heredoc keeps every `$` literal; </dev/null so R
# never blocks on the harness socket stdin.
RSRC="$(mktemp --suffix=.R)"; trap 'rm -f "$RSRC"' EXIT
cat > "$RSRC" <<'RSCRIPT'
args <- commandArgs(trailingOnly = TRUE); RLIB <- args[1]
.libPaths(c(RLIB, .libPaths())); suppressMessages({library(Matrix); library(lstar)})
set.seed(1)

# --- authoritative Seurat 2.3.4 S4 class definitions (pre-Assay), verbatim slot layout ---------------
setClass("dim.reduction", slots = list(cell.embeddings = "matrix", gene.loadings = "matrix",
  gene.loadings.full = "matrix", sdev = "numeric", key = "character", jackstraw = "ANY", misc = "ANY"))
setClass("assay", slots = list(raw.data = "ANY", data = "ANY", scale.data = "ANY", key = "character",
  misc = "ANY", var.genes = "vector", mean.var = "data.frame"))
setClass("seurat", slots = c(raw.data = "ANY", data = "ANY", scale.data = "ANY", var.genes = "vector",
  is.expr = "numeric", ident = "factor", meta.data = "data.frame", project.name = "character",
  dr = "list", assay = "list", hvg.info = "data.frame", imputed = "data.frame", cell.names = "vector",
  cluster.tree = "list", snn = "dgCMatrix", calc.params = "list", kmeans = "ANY", spatial = "ANY",
  misc = "ANY", version = "ANY"))

ng <- 40; nc <- 16; g <- paste0("Gene", 1:ng); ce <- paste0("Cell", 1:nc)
raw <- as(matrix(rpois(ng * nc, 3), ng, nc, dimnames = list(g, ce)), "CsparseMatrix")   # genes x cells
dat <- raw; dat@x <- log1p(dat@x)
vg  <- g[1:12]                                                              # variable genes (a subset)
scl <- matrix(rnorm(length(vg) * nc), length(vg), nc, dimnames = list(vg, ce))  # scaled over var.genes only
emb <- matrix(rnorm(nc * 5), nc, 5, dimnames = list(ce, paste0("PC", 1:5)))
ld  <- matrix(rnorm(length(vg) * 5), length(vg), 5, dimnames = list(vg, paste0("PC", 1:5)))  # loadings over var.genes
pca <- new("dim.reduction", cell.embeddings = emb, gene.loadings = ld,
           gene.loadings.full = matrix(nrow = 0, ncol = 0), sdev = sqrt(5:1), key = "PC")
np <- 8; pr <- paste0("ADT", 1:np)                                          # a v2 multimodal `assay` (CITE-seq)
adt <- as(matrix(rpois(np * nc, 5), np, nc, dimnames = list(pr, ce)), "CsparseMatrix")
cite <- new("assay", raw.data = adt, data = adt, scale.data = matrix(nrow = 0, ncol = 0),
            key = "CITE_", misc = list(), var.genes = character(0), mean.var = data.frame())
so <- new("seurat", raw.data = raw, data = dat, scale.data = scl, var.genes = vg, is.expr = 0,
          ident = factor(setNames(rep(c("Tcell", "Bcell"), length.out = nc), ce)),
          meta.data = data.frame(nGene = colSums(raw > 0),
                                 batch = factor(rep(c("b1", "b2"), length.out = nc)),  # v2-era meta.data: factors
                                 row.names = ce),
          project.name = "legacy", dr = list(pca = pca), assay = list(CITE = cite),
          hvg.info = data.frame(), imputed = data.frame(), cell.names = ce, cluster.tree = list(),
          snn = Matrix::sparseMatrix(i = 1:nc, j = 1:nc, x = 1, dims = c(nc, nc)), calc.params = list(),
          spatial = list(), misc = list(note = "legacy v2"), version = package_version("2.3.4"))
p0 <- "/tmp/sv2_legacy.rds"; saveRDS(so, p0)

# mimic the REAL scenario: the ancient classes are UNDEFINED when a modern R loads the old .rds
removeClass("seurat"); removeClass("dim.reduction"); removeClass("assay")
so2 <- readRDS(p0)
ds  <- read_seurat(so2)
stopifnot(
  identical(ds$kind, "sample"),
  any(grepl("object@seurat", ds$profiles)),                                 # v2 class recognized + recorded
  identical(as.character(ds$fields[["counts"]]$span), c("cells", "genes")), ds$fields[["counts"]]$state == "raw",
  ds$fields[["X"]]$state == "lognorm",
  identical(ds$fields[["scale.data"]]$coverage, "partial"),                 # scaled over var.genes -> partial
  length(ds$fields[["scale.data"]]$index) == length(vg),
  "pca" %in% names(ds$axes), identical(ds$fields[["pca"]]$role, "embedding"),
  length(ds$fields[["pca_stdev"]]$values) == 5,
  identical(as.character(ds$fields[["pca_loadings"]]$span), c("pca_features", "pca")),  # loadings over HVG subset
  length(ds$axes[["pca_features"]]$labels) == length(vg),
  identical(ds$fields[["ident"]]$subtype, "active_ident"), "ident" %in% names(ds$axes),
  "batch" %in% names(ds$axes),                                              # factor meta.data -> factor axis
  identical(ds$fields[["snn"]]$role, "relation"),
  identical(ds$fields[["variable_features"]]$role, "measure"),
  "CITE" %in% names(ds$axes), identical(ds$axes$CITE$role, "feature"),      # v2 multimodal assay -> 2nd feature axis
  identical(as.character(ds$fields[["CITE.counts"]]$span), c("cells", "CITE")))
cat("  [R] v2 seurat (pre-Assay) read    OK\n")

# old -> new conversion is the payoff: write_seurat emits a MODERN object (Assay/Assay5 + the reduction)
so3 <- write_seurat(ds)
stopifnot(inherits(so3[["RNA"]], "Assay") || inherits(so3[["RNA"]], "Assay5"),
          "pca" %in% SeuratObject::Reductions(so3), "CITE" %in% SeuratObject::Assays(so3))
cat("  [R] v2 -> modern Seurat (Assay/Assay5 + pca + CITE) OK\n")

p <- "/tmp/sv2_legacy.lstar.zarr"; if (dir.exists(p)) unlink(p, recursive = TRUE)
lstar_write(ds, p); cat(p, "\n", file = "/tmp/sv2_legacy_store.txt")
RSCRIPT
Rscript "$RSRC" "$RLIB" </dev/null 2>&1 | grep -E "^  \[R\]|Error|Execution halted|cannot|unable|no method"

# the v2-derived store validates in Python too (cross-language)
PYTHONPATH="$ROOT/python/src" python3 - <<'PY'
import warnings; warnings.filterwarnings("ignore")
import lstar
p = open("/tmp/sv2_legacy_store.txt").read().split()[0]
ds = lstar.read(p)
errs = [i for i in lstar.validate(ds) if i.startswith("ERROR")]
assert not errs, (p, errs)
assert "CITE" in ds.axes and ds.field("counts").state == "raw"
print("  [py] v2 (pre-Assay) store validates clean (%d axes, %d fields)" % (len(ds.axes), len(ds.fields)))
PY
