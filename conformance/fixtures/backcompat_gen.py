#!/usr/bin/env python3
# (Re)generate the committed backwards-compatibility fixtures: a small MAXIMAL L* store (every encoding)
# written both as Zarr v2 and v3, plus expected.json — the value manifest derived from the INPUT arrays
# (NOT from reading the store back, so it's ground truth independent of the reader/writer). conformance/
# backcompat.sh reads the COMMITTED stores on every surface and asserts they still decode to expected.json.
#
# The v2 store is the backwards-compat guard: its bytes are what a PRE-flip lstar produced (the v2 write
# path is unchanged by the v2->v3 default flip), so "current readers still read old stores" is pinned.
# The v3 store freezes today's v3 layout so a future v3 change that breaks old-v3 reads is caught too.
#
# Run:  PYTHONPATH=python/src python3 conformance/fixtures/backcompat_gen.py   (from the repo root)
# then commit conformance/fixtures/{golden_v2,golden_v3}.lstar.zarr + expected.json.
import os, sys, json, shutil
import numpy as np, scipy.sparse as sp
import lstar

HERE = os.path.dirname(os.path.abspath(__file__))
NC, NG = 12, 8

def build():
    rng = np.random.default_rng(7)
    ds = lstar.Dataset(kind="sample")
    ds.profiles.append("backcompat@0.1")
    ds.add_axis("cells", [f"c{i}" for i in range(NC)], role="observation")
    ds.add_axis("genes", [f"g{j}" for j in range(NG)], role="feature")
    ds.add_axis("pca", [f"PC{k}" for k in range(3)], origin="derived", role="coordinate")

    counts = sp.csc_matrix((np.ceil(rng.random((NC, NG)) * 5) * (rng.random((NC, NG)) < 0.5)).astype("f4"))
    ds.add_field("counts", counts, role="measure", span=["cells", "genes"], state="raw")            # csc
    logn = sp.csr_matrix(counts.toarray().astype("f4"))
    ds.add_field("lognorm", logn, role="measure", span=["cells", "genes"], state="lognorm")          # csr
    pca = rng.standard_normal((NC, 3)).astype("f4")
    ds.add_field("pca", pca, role="embedding", span=["cells", "pca"])                                # dense
    leiden = np.array([f"cl{i % 3}" for i in range(NC)])
    ds.add_field("leiden", list(leiden), role="label", span=["cells"])                               # utf8
    import pandas as pd
    ct = pd.Categorical.from_codes(np.array([0, 1, 2, -1, 0, 1, 2, 0, 1, 2, 0, 1], dtype="int8"),
                                   categories=["T", "B", "NK"], ordered=True)
    ds.add_field("celltype", ct, role="label", span=["cells"])                                       # categorical (+ -1)
    qcv = (rng.random(NC) * 10).astype("f4"); qcmask = np.zeros(NC, "uint8"); qcmask[[2, 7]] = 1      # nullable
    ds.add_field("qc", qcv, role="measure", span=["cells"], mask=qcmask)
    cov = np.array([0, 2, 4, 6, 8, 10]); adt = (rng.random(len(cov)) * 3).astype("f4")               # partial
    ds.add_field("adt", adt, role="measure", span=["cells"], index=cov, index_axis="cells")
    G = sp.csr_matrix((rng.random((NC, NC)) < 0.25).astype("f4"))                                     # graph
    gf = ds.add_field("knn", G, role="graph", span=["cells", "cells"])
    ds.fields["knn"].directed = True; ds.fields["knn"].weighted = True
    ds.aux["anndata.uns"] = {"params": {"n_pca": 3, "method": "leiden"},                              # aux tree
                             "scores": np.arange(4.0), "names": ["a", "b", "c"]}
    return ds, counts, logn, pca, leiden, ct, qcv, qcmask, cov, adt, G

def manifest(counts, logn, pca, leiden, ct, qcv, qcmask, cov, adt, G):
    # ground truth from the INPUT arrays (independent of read/write)
    return {
        "counts": {"nnz": int(counts.nnz), "sum": round(float(counts.sum()), 4), "shape": [NC, NG]},
        "lognorm": {"nnz": int(logn.nnz), "sum": round(float(logn.sum()), 4), "shape": [NC, NG]},
        "pca": {"shape": [NC, 3], "sum": round(float(pca.sum()), 4)},
        "leiden": list(map(str, leiden)),
        "celltype": {"categories": list(ct.categories), "codes": [int(c) for c in ct.codes], "ordered": bool(ct.ordered)},
        "qc": {"missing_positions": [int(i) for i in np.nonzero(qcmask)[0]], "present_sum": round(float(qcv[qcmask == 0].sum()), 4)},
        "adt": {"index": [int(i) for i in cov], "sum": round(float(adt.sum()), 4)},
        "knn": {"nnz": int(G.nnz), "directed": True, "weighted": True},
        "aux": {"names": ["a", "b", "c"], "n_pca": 3},
    }

if __name__ == "__main__":
    built = build()
    ds = built[0]
    man = manifest(*built[1:])
    json.dump(man, open(os.path.join(HERE, "expected.json"), "w"), indent=2, sort_keys=True)
    for fmt, name in (("v2", "golden_v2.lstar.zarr"), ("v3", "golden_v3.lstar.zarr")):
        out = os.path.join(HERE, name)
        shutil.rmtree(out, ignore_errors=True)
        lstar.write(ds, out, format=fmt)                 # single-chunk (default) — keeps the committed tree small
    print("wrote golden_v2 + golden_v3 + expected.json under", HERE)
