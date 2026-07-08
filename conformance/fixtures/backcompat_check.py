#!/usr/bin/env python3
# Read a committed golden store (v2 or v3) with the CURRENT lstar and assert it still decodes to the frozen
# expected.json manifest (ground truth from the fixture's INPUT arrays). This is the backwards-compat guard:
# a store written by a pre-flip lstar (v2) — or today's v3 — must keep reading correctly forever.
# Usage: backcompat_check.py <store> <expected.json>
import sys, json
import numpy as np, scipy.sparse as sp, lstar

store, expected_path = sys.argv[1], sys.argv[2]
exp = json.load(open(expected_path))
ds = lstar.read(store)


def sp_of(f):
    v = ds.field(f).values
    return v if sp.issparse(v) else sp.csr_matrix(np.asarray(v))


def r4(x):
    return round(float(x), 4)


got = {}
for f in ("counts", "lognorm"):
    m = sp_of(f)
    got[f] = {"nnz": int(m.nnz), "sum": r4(m.sum()), "shape": list(m.shape)}
pca = np.asarray(ds.field("pca").values)
got["pca"] = {"shape": list(pca.shape), "sum": r4(pca.sum())}
got["leiden"] = [str(x) for x in np.asarray(ds.field("leiden").values)]
ct = ds.field("celltype")
cats = list(ct.categories) if getattr(ct, "categories", None) is not None else list(getattr(ct.values, "categories", []))
codes = ct.codes if getattr(ct, "codes", None) is not None else getattr(ct.values, "codes", None)
got["celltype"] = {"categories": [str(c) for c in cats], "codes": [int(c) for c in np.asarray(codes)],
                   "ordered": bool(getattr(ct, "ordered", getattr(ct.values, "ordered", False)))}
qc = ds.field("qc"); qv = np.asarray(qc.values); qm = np.asarray(qc.mask).astype(bool)
got["qc"] = {"missing_positions": [int(i) for i in np.nonzero(qm)[0]], "present_sum": r4(qv[~qm].sum())}
adt = ds.field("adt")
got["adt"] = {"index": [int(i) for i in np.asarray(adt.index)], "sum": r4(np.asarray(adt.values).sum())}
knn = ds.field("knn"); km = sp_of("knn")
got["knn"] = {"nnz": int(km.nnz), "directed": bool(knn.directed), "weighted": bool(knn.weighted)}
aux = ds.aux["anndata.uns"] if isinstance(ds.aux, dict) else dict(ds.aux)["anndata.uns"]
got["aux"] = {"names": [str(x) for x in aux["names"]], "n_pca": int(aux["params"]["n_pca"])}

bad = [k for k in exp if got.get(k) != exp[k]]
if bad:
    for k in bad:
        print(f"  MISMATCH {k}:\n    got  {got.get(k)}\n    want {exp[k]}", file=sys.stderr)
    sys.exit(1)
print(f"  [py] {store.split('/')[-1]}: reads to the frozen manifest (all encodings)")
