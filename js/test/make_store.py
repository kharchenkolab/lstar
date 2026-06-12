"""Generate a small, deterministic L* store + an expected-values JSON for the TS reader tests."""
import json
import os
import sys

import numpy as np
import scipy.sparse as sp

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "python", "src"))
import lstar
from lstar.lazy import stream_col_stats

HERE = os.path.dirname(__file__)
STORE = os.path.join(HERE, "data", "sample.lstar.zarr")
EXPECT = os.path.join(HERE, "data", "expected.json")


def main():
    os.makedirs(os.path.join(HERE, "data"), exist_ok=True)
    ncells, ngenes = 50, 20
    rng = np.random.default_rng(7)
    X = sp.csc_matrix(sp.random(ncells, ngenes, density=0.25, format="csc", random_state=rng))
    X.data = np.round(X.data * 9 + 1)                 # small integer-ish counts
    umap = rng.standard_normal((ncells, 2)).astype("f4")
    leiden = np.array([["A", "B", "C"][i % 3] for i in range(ncells)])
    n_umi = np.asarray(X.sum(axis=1)).ravel().astype("f4")
    qc = np.floor(n_umi).astype(np.int64)             # a nullable integer measure (a couple missing)
    qc_mask = np.zeros(ncells, dtype=np.uint8); qc_mask[[1, 3]] = 1

    ds = lstar.Dataset(kind="sample")
    ds.add_axis("cells", [f"cell{i}" for i in range(ncells)], role="observation")
    ds.add_axis("genes", [f"g{i}" for i in range(ngenes)], role="feature")
    ds.add_axis("umap", ["umap0", "umap1"], origin="derived", role="coordinate")
    ds.add_field("counts", X, role="measure", span=["cells", "genes"], state="raw")
    ds.add_field("umap", umap, role="embedding", span=["cells", "umap"])
    ds.add_field("leiden", leiden, role="label", span=["cells"])
    ds.add_field("n_umi", n_umi, role="measure", span=["cells"])
    ds.add_field("qc", qc, role="measure", span=["cells"], mask=qc_mask)   # nullable: 1 == missing
    ds.aux["anndata.uns"] = {                         # lossless passthrough: params + colors + nested
        "log1p": {"base": None},
        "pca": {"variance_ratio": np.array([0.6, 0.3, 0.1])},
        "leiden_colors": np.array(["#aa0000", "#00bb00", "#0000cc"]),
    }
    lstar.write(ds, STORE)                            # Zarr v2, single-chunk, consolidated metadata

    # Expected values for the TS tests.
    gcol = 3
    col = X[:, gcol].tocoo()
    mean, var, nnz = stream_col_stats(X, lognorm=True, engine="python")

    # A DE reference: per-gene log1p group means for leiden A vs B.
    cellsA = [i for i in range(ncells) if leiden[i] == "A"]
    cellsB = [i for i in range(ncells) if leiden[i] == "B"]
    Xl = X.copy(); Xl.data = np.log1p(Xl.data)
    meanA = np.asarray(Xl[cellsA].mean(axis=0)).ravel()
    meanB = np.asarray(Xl[cellsB].mean(axis=0)).ravel()

    expected = {
        "de_ref": {"cellsA": cellsA, "cellsB": cellsB,
                   "meanA": meanA.tolist(), "meanB": meanB.tolist()},
        "kind": ds.kind,
        "axes": list(ds.axes),
        "fields": list(ds.fields),
        "n_cells": ncells,
        "n_genes": ngenes,
        "umap": umap.ravel(order="C").tolist(),       # C-order, ncells x 2
        "leiden": leiden.tolist(),
        "n_umi": n_umi.tolist(),
        "qc": {"values": qc.tolist(), "mask": qc_mask.tolist()},
        "aux_ns": list(ds.aux),
        "uns": {"variance_ratio": [0.6, 0.3, 0.1], "colors": ["#aa0000", "#00bb00", "#0000cc"]},
        "gene_col": {"index": gcol, "rows": col.row.tolist(), "vals": col.data.tolist()},
        "colstats_lognorm": {"mean": mean.tolist(), "var": var.tolist(), "nnz": nnz.tolist()},
        "counts_sum": float(X.sum()),
        "counts_dense": X.toarray().ravel(order="C").tolist(),
    }
    json.dump(expected, open(EXPECT, "w"))
    print(f"wrote {STORE} ({ncells}x{ngenes}, nnz={X.nnz}) and {EXPECT}")


if __name__ == "__main__":
    main()
