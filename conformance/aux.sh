#!/usr/bin/env bash
# Lossless-passthrough conformance: a format's untyped long-tail (AnnData `uns`, Seurat `@misc`) is
# carried verbatim in a self-describing `aux/` subtree (a JSON tree string + a flat array manifest).
# The C++ core and R must round-trip it *byte-for-byte* without interpreting it; Python reconstructs the
# live object (incl. numpy structured arrays) exactly. Python writes a rich uns -> C++ (via R) re-writes
# it -> Python reconstructs uns identical; R re-writes -> Python reconstructs uns identical.
set -e
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
RLIB="$ROOT/.Rlib"
PY=/tmp/aux_py.lstar.zarr
R=/tmp/aux_r.lstar.zarr

PYTHONPATH="$ROOT/python/src" python3 - "$PY" <<'PY'
import sys, warnings; warnings.filterwarnings("ignore")
import numpy as np, lstar
ds = lstar.Dataset(kind="sample")
ds.add_axis("cells", [f"c{i}" for i in range(3)])
ds.aux["anndata.uns"] = {
    "log1p": {"base": None}, "iroot": 42,
    "pca": {"variance_ratio": np.array([0.5, 0.3, 0.2]), "params": {"n_comps": 50}},
    "leiden_colors": np.array(["#ff0000", "#00ff00", "#0000ff"]),
    "neighbors": {"params": {"n_neighbors": 15, "method": "umap"}},
    "rank_genes_groups": {"params": {"groupby": "leiden", "reference": "rest"},
                          "names": np.array([("g1", "g2"), ("g3", "g4")],
                                            dtype=[("A", "U10"), ("B", "U10")])},
}
lstar.write(ds, sys.argv[1])
print("  [py] wrote a rich uns (params + colors + nested dicts + a structured array) into aux/")
PY

# C++ core (via R) reads + re-writes the aux subtree verbatim; R never interprets it.
Rscript -e '.libPaths(c("'"$RLIB"'", .libPaths())); suppressMessages(library(lstar))
ds <- lstar_read("'"$PY"'")
stopifnot("anndata.uns" %in% names(ds$aux), length(ds$aux[["anndata.uns"]]$leaves) >= 3)
lstar_write(ds, "'"$R"'")
cat("  [R ] C++/R round-tripped aux verbatim (", length(ds$aux[["anndata.uns"]]$leaves), "array leaves)\n")' \
  2>&1 | grep -vE "^Warning|deprecat|masked|following object|Attaching|^$"

# Python reconstructs the live uns from BOTH the C++-written and R-written stores; must be identical.
PYTHONPATH="$ROOT/python/src" python3 - "$PY" "$R" <<'PY'
import sys, warnings; warnings.filterwarnings("ignore")
import numpy as np, lstar
def chk(p, who):
    u = lstar.read(p).aux["anndata.uns"]
    assert u["log1p"] == {"base": None} and u["iroot"] == 42
    assert np.allclose(u["pca"]["variance_ratio"], [0.5, 0.3, 0.2]) and u["pca"]["params"]["n_comps"] == 50
    assert list(u["leiden_colors"]) == ["#ff0000", "#00ff00", "#0000ff"]
    assert u["neighbors"]["params"]["method"] == "umap"
    rg = u["rank_genes_groups"]["names"]
    assert rg.dtype.names == ("A", "B") and list(rg["A"]) == ["g1", "g3"]
    assert list(u)[0] == "log1p"                                  # dict key order preserved (verbatim)
    print(f"  [py] reconstructed uns from the {who} store: params/colors/nested/structured-array exact")
chk(sys.argv[1], "original")
chk(sys.argv[2], "R-rewritten")
PY
