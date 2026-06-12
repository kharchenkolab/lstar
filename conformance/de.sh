#!/usr/bin/env bash
# DE-bundle conformance: a differential-expression result typed from scanpy's `rank_genes_groups` is an
# ordinary bundle of measures over the induced (factor, genes) axes, so it round-trips across languages
# like any dense field. Python writes the bundle; R reads it back in the canonical (factor x genes)
# orientation and builds the tidy marker table (`lstar_markers`); Python regenerates rank_genes_groups.
set -e
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
RLIB="$ROOT/.Rlib"
PY=/tmp/de_py.lstar.zarr

PYTHONPATH="$ROOT/python/src" python3 - "$PY" <<'PY'
import sys, warnings; warnings.filterwarnings("ignore")
import numpy as np, anndata as ad, pandas as pd, scanpy as sc, scipy.sparse as sp, lstar
rng = np.random.default_rng(0); n, g = 120, 30
a = ad.AnnData(X=sp.csr_matrix(rng.poisson(0.4, size=(n, g)).astype("float32")),
               obs=pd.DataFrame({"leiden": pd.Categorical([str(i % 3) for i in range(n)])},
                                index=[f"cell{i}" for i in range(n)]),
               var=pd.DataFrame(index=[f"g{j}" for j in range(g)]))
sc.pp.normalize_total(a, target_sum=1e4); sc.pp.log1p(a)
sc.tl.rank_genes_groups(a, "leiden", method="t-test")
ds = lstar.read_anndata(a)
assert ds.field("de.leiden.lfc").span == ["leiden", "genes"]
lstar.write(ds, sys.argv[1])
# stash a couple of reference values for R to check
rgg = a.uns["rank_genes_groups"]
g0 = rgg["names"].dtype.names[0]
print("REF", g0, str(rgg["names"][g0][0]), float(rgg["scores"][g0][0]))
PY

REF=$(PYTHONPATH="$ROOT/python/src" python3 - "$PY" <<'PY'
import sys, warnings; warnings.filterwarnings("ignore")
import numpy as np, lstar
ds = lstar.read(sys.argv[1]); m = lstar.markers(ds, "leiden", top=3)
g0 = str(np.asarray(ds.axis("leiden").labels)[0])
sub = m[m["group"] == g0].sort_values("score", ascending=False).iloc[0]
print(g0, sub["gene"])
PY
)
echo "  [py] wrote DE bundle; markers() top gene for group $REF"

Rscript -e '.libPaths(c("'"$RLIB"'", .libPaths())); suppressMessages(library(lstar))
ds <- lstar_read("'"$PY"'")
f <- ds$fields[["de.leiden.score"]]
stopifnot(identical(as.character(f$span), c("leiden","genes")),
          nrow(as.matrix(f$values)) == length(ds$axes$leiden$labels),
          ncol(as.matrix(f$values)) == length(ds$axes$genes$labels))
mk <- lstar_markers(ds, "leiden", top = 3, sort_by = "score")
stopifnot(all(c("group","gene","score","lfc") %in% names(mk)), nrow(mk) == 3 * length(ds$axes$leiden$labels))
g0 <- as.character(ds$axes$leiden$labels)[1]
top <- mk[mk$group == g0, ][which.max(mk[mk$group == g0, "score"]), "gene"]
cat("  [R ] read (factor x genes) bundle; lstar_markers top gene for group", g0, "=", top, "\n")' \
  2>&1 | grep -vE "^Warning|deprecat|masked|following object|Attaching|^$"
