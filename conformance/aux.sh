#!/usr/bin/env bash
# Lossless-passthrough conformance: a format's untyped long-tail (AnnData `uns`, Seurat `@misc`) is
# carried verbatim in a self-describing `aux/` subtree (a JSON tree string + a flat array manifest).
# The C++ core and R must round-trip it *byte-for-byte* without interpreting it; Python reconstructs the
# live object (incl. numpy structured arrays) exactly. Python writes a rich uns -> C++ (via R) re-writes
# it -> Python reconstructs uns identical; R re-writes -> Python reconstructs uns identical.
# Origin coverage: Py-authored ✓ | R-authored ✓ (R overwrites a leaf + Python cross-reads) — see conformance/README.md
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
lv <- ds$aux[["anndata.uns"]]$leaves                                       # check leaf KINDS, not just count
kinds <- vapply(lv, function(l) l$kind, character(1))
stopifnot("anndata.uns" %in% names(ds$aux), "utf8" %in% kinds, any(kinds != "utf8"), all(nzchar(names(lv))))
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

# R-AUTHORED aux content (origin coverage): R OVERWRITES a utf8 leaf's strings with R-native values +
# renames the namespace, then writes -> Python reconstructs the R-authored strings. Proves R-side data
# flows INTO the aux payload (not merely shuttled verbatim from what C++ handed it), distinct from the
# rewrite leg above. (Hand-authoring raw dense leaf bytes isn't a realistic R origin; mutating one is.)
RA=/tmp/aux_r_authored.lstar.zarr
Rscript -e '.libPaths(c("'"$RLIB"'", .libPaths())); suppressMessages(library(lstar))
ds <- lstar_read("'"$PY"'")
lv <- ds$aux[["anndata.uns"]]$leaves
isc <- vapply(lv, function(l) identical(l$kind,"utf8") && length(l$strings)==3, logical(1))  # the leiden_colors leaf
stopifnot(any(isc))
lv[[which(isc)[1]]]$strings <- c("R_red","R_green","R_blue")               # R authors the content
ds$aux <- list("r.authored" = list(attrs=ds$aux[["anndata.uns"]]$attrs, leaves=lv))  # same tree/manifest, renamed ns
lstar_write(ds, "'"$RA"'")
cat("  [R ] authored aux content (overwrote a utf8 leaf, renamed ns) -> wrote r.authored\n")' \
  2>&1 | grep -vE "^Warning|deprecat|masked|following object|Attaching|^$"

PYTHONPATH="$ROOT/python/src" python3 - "$RA" <<'PY'
import sys, warnings; warnings.filterwarnings("ignore")
import numpy as np, lstar
u = lstar.read(sys.argv[1]).aux["r.authored"]
assert list(np.asarray(u["leiden_colors"])) == ["R_red", "R_green", "R_blue"], u["leiden_colors"]
print("  [py] reconstructed R-authored aux: the utf8 leaf carries R-native values (R-origin cross-read)")
PY
