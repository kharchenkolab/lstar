#!/usr/bin/env bash
# "True collection" conformance (non-Conos): `collection_from()` assembles a collection of HETEROGENEOUS
# samples -- overlapping-but-divergent gene sets, and even fully DISJOINT (cross-species) ones -- which
# must round-trip AS a collection (never flattened to one cells x genes matrix), in BOTH Python and R, and
# support collection operations (streamed pseudobulk over the union). The differentiator a fixed tensor
# can't express: a linked set of heterogeneous samples + a joint union layer.
# Origin coverage: Py-authored ✓ | R-authored ✓ (each cross-read by the other) — see conformance/README.md
set -eo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
RLIB="$ROOT/.Rlib"; export LSTAR_RLIB="$RLIB"
PS=/tmp/coltrue_py.lstar.zarr; RS=/tmp/coltrue_r.lstar.zarr; DS=/tmp/coltrue_disjoint.lstar.zarr

# ── G1 (Python origin) + G3 (collection op): a 4-sample divergent-gene collection via collection_from,
#    streamed collection_pseudobulk over the joint clustering, written to a store. Realistic-ish size.
PYTHONPATH="$ROOT/python/src" python3 - "$PS" <<'PY'
import sys, warnings; warnings.filterwarnings("ignore")
import numpy as np, scipy.sparse as sp, lstar
from lstar import collection_from, collection_pseudobulk, Categorical
rng = np.random.default_rng(1)
shared = [f"g{i}" for i in range(350)]
samples = {}
for k, s in enumerate(["S1", "S2", "S3", "S4"]):
    nc = 300 + 50 * k
    genes = shared + [f"{s}_only{i}" for i in range(30 + 10 * k)]            # sample-specific genes
    d = lstar.Dataset(kind="sample")
    d.add_axis("cells", [f"c{i}" for i in range(nc)]); d.add_axis("genes", genes)
    d.add_field("counts", sp.random(nc, len(genes), density=0.1, format="csc").astype("float32"),
                role="measure", span=["cells", "genes"], state="raw")
    d.add_field("qc", rng.random(nc).astype("float32"), role="measure", span=["cells"])
    samples[s] = d
ntot = 300 + 350 + 400 + 450
joint = {"umap": rng.standard_normal((ntot, 2)).astype("float32"),
         "clusters": Categorical(rng.integers(0, 5, ntot), np.array([f"k{i}" for i in range(5)]))}
col = collection_from(samples, joint=joint)
assert col.kind == "collection" and not [e for e in lstar.validate(col) if e.startswith("ERROR")]
# never flattened: per-sample axes with DISTINCT lengths + sample-specific genes isolated
for s, nc in zip(["S1","S2","S3","S4"], [300,350,400,450]):
    assert col.axis(f"cells.{s}").labels.shape[0] == nc and f"counts.{s}" in col.fields
assert "S1_only0" in set(np.asarray(col.axis("genes.S1").labels))
assert "S1_only0" not in set(np.asarray(col.axis("genes.S2").labels))
assert col.axis("cells").labels.shape[0] == ntot and col.field("sample").subtype == "design"
# G3: streamed pseudobulk over the joint clustering, genes aligned by label across samples
pb = collection_pseudobulk(col, "clusters", field="counts")
ng_union = len({g for s in ["S1","S2","S3","S4"] for g in np.asarray(col.axis(f"genes.{s}").labels)})
assert col.field("pb.clusters.mean").values.shape == (5, ng_union), col.field("pb.clusters.mean").values.shape
lstar.write(col, sys.argv[1])
print(f"  [py] G1 built a 4-sample divergent collection (union {ntot} cells, {ng_union} genes); "
      f"G3 streamed pseudobulk -> (5 x {ng_union}); never-flattened; wrote store")
PY
Rscript -e '.libPaths(c("'"$RLIB"'", .libPaths())); suppressMessages(library(lstar))
c <- lstar_read("'"$PS"'")
stopifnot(identical(c$kind, "collection"))
for (s in c("S1","S2","S3","S4")) stopifnot(!is.null(c$axes[[paste0("cells.",s)]]), !is.null(c$fields[[paste0("counts.",s)]]))
gS1 <- as.character(c$axes[["genes.S1"]]$labels); gS2 <- as.character(c$axes[["genes.S2"]]$labels)
stopifnot("S1_only0" %in% gS1, !("S1_only0" %in% gS2))                       # never flattened (R cross-read)
stopifnot(!is.null(c$axes$cells), identical(c$fields$sample$subtype, "design"), !is.null(c$fields[["pb.clusters.mean"]]))
cat("  [R ] cross-read the Python collection: per-sample heterogeneity + union + pseudobulk intact\n")' \
  2>&1 | grep -E "^  \[R"

# ── G1b (R origin): R builds a divergent collection via collection_from -> Python cross-reads.
Rscript -e '.libPaths(c("'"$RLIB"'", .libPaths())); suppressMessages({library(lstar); library(Matrix)})
mk <- function(s, nc, genes) { d <- list(kind="sample", axes=list(), fields=list())
  d$axes$cells <- list(labels=paste0("c",seq_len(nc)), origin="observed", role="observation")
  d$axes$genes <- list(labels=genes, origin="observed", role="feature")
  d$fields$counts <- list(role="measure", span=c("cells","genes"), state="raw",
    values=as(Matrix::Matrix(matrix(rpois(nc*length(genes),2),nc),sparse=TRUE),"CsparseMatrix"))
  class(d)<-"lstar_dataset"; d }
shared <- paste0("g",1:300)
col <- collection_from(list(A=mk("A",200,c(shared,paste0("A_only",1:40))),
                            B=mk("B",250,c(shared,paste0("B_only",1:60))),
                            C=mk("C",300,c(shared,paste0("C_only",1:25)))),
                       joint=list(umap=matrix(rnorm(750*2),750,2)))
stopifnot(identical(col$kind,"collection"))
lstar_write(col, "'"$RS"'")
cat("  [R ] G1b authored a 3-sample divergent collection -> store\n")' 2>&1 | grep -E "^  \[R"
PYTHONPATH="$ROOT/python/src" python3 - "$RS" <<'PY'
import sys, warnings; warnings.filterwarnings("ignore")
import numpy as np, lstar
c = lstar.read(sys.argv[1])
assert c.kind == "collection" and not [e for e in lstar.validate(c) if e.startswith("ERROR")]
for s, nc in zip("ABC", [200,250,300]): assert c.axis(f"cells.{s}").labels.shape[0] == nc
gA = set(np.asarray(c.axis("genes.A").labels)); gB = set(np.asarray(c.axis("genes.B").labels))
assert "A_only1" in gA and "A_only1" not in gB and c.axis("cells").labels.shape[0] == 750
print("  [py] cross-read the R-authored collection: 3 divergent samples, never-flattened")
PY

# ── G2 (disjoint / cross-species): two samples with NON-overlapping gene namespaces. No union genes axis
#    is even possible; the joint analysis lives only on the union `cells`. Python origin -> R cross-read.
PYTHONPATH="$ROOT/python/src" python3 - "$DS" <<'PY'
import sys, warnings; warnings.filterwarnings("ignore")
import numpy as np, scipy.sparse as sp, lstar
from lstar import collection_from
rng = np.random.default_rng(2)
def sample(prefix, nc, ng):
    d = lstar.Dataset(kind="sample")
    d.add_axis("cells", [f"c{i}" for i in range(nc)])
    d.add_axis("genes", [f"{prefix}_{i}" for i in range(ng)])                # disjoint namespaces
    d.add_field("counts", sp.random(nc, ng, density=0.1, format="csc").astype("float32"),
                role="measure", span=["cells", "genes"], state="raw")
    return d
col = collection_from({"human": sample("HGNC", 300, 400), "mouse": sample("MGI", 250, 350)},
                      joint={"latent": rng.standard_normal((550, 10)).astype("float32")})   # integration latent
assert col.kind == "collection" and not [e for e in lstar.validate(col) if e.startswith("ERROR")]
gh = set(np.asarray(col.axis("genes.human").labels)); gm = set(np.asarray(col.axis("genes.mouse").labels))
assert gh.isdisjoint(gm), "cross-species gene sets must be disjoint"
assert "genes" not in col.axes, "a disjoint collection must NOT have a single shared genes axis"
assert col.field("latent").role == "embedding" and col.field("latent").values.shape == (550, 10)
lstar.write(col, sys.argv[1])
print("  [py] G2 built a cross-species collection (disjoint HGNC/MGI genes, joint latent over union)")
PY
Rscript -e '.libPaths(c("'"$RLIB"'", .libPaths())); suppressMessages(library(lstar))
c <- lstar_read("'"$DS"'")
gh <- as.character(c$axes[["genes.human"]]$labels); gm <- as.character(c$axes[["genes.mouse"]]$labels)
stopifnot(identical(c$kind,"collection"), length(intersect(gh, gm)) == 0, is.null(c$axes$genes),
          !is.null(c$fields$latent), ncol(c$fields$latent$values) == 10)
cat("  [R ] cross-read the disjoint cross-species collection: no shared gene axis; joint latent over the union\n")' \
  2>&1 | grep -E "^  \[R"
echo "true-collection conformance PASSED."
