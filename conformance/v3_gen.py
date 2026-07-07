#!/usr/bin/env python3
# Build a maximal v2 lstar store (all encodings -- csc/csr/dense/categorical/utf8 + nullable/graph/partial
# + aux passthrough + viewer-derived fields and factor axes) as the seed for the v3-format conformance
# leg. Written v2 with gzip -- the format migration re-emits it as v3 and every surface reads it back.
import sys, numpy as np, scipy.sparse as sp, pandas as pd, lstar, numcodecs
from lstar.viewer import extend_for_viewer

out = sys.argv[1]
rng = np.random.default_rng(0)
n, g = 200, 50
ds = lstar.Dataset(kind="sample")
ds.add_axis("cells", ["c%d" % i for i in range(n)])
ds.add_axis("genes", ["g%d" % j for j in range(g)])
ds.add_axis("pca", ["PC%d" % k for k in range(10)], origin="derived")
C = sp.random(n, g, density=0.3, format="csc", random_state=1); C.data = np.ceil(C.data * 10).astype("float32")
ds.add_field("counts", C, role="measure", span=["cells", "genes"], state="raw")
ds.add_field("data", sp.random(n, g, density=0.3, format="csr", random_state=2).astype("float32"),
             role="measure", span=["cells", "genes"], state="lognorm")
ds.add_field("pca", rng.standard_normal((n, 10)).astype("float32"), role="embedding", span=["cells", "pca"])
ds.add_field("leiden", pd.Categorical(["cl%d" % (i % 5) for i in range(n)]), role="label", span=["cells"])
ds.add_field("barcode", np.array(["BC%04d" % i for i in range(n)]), role="label", span=["cells"])
ds.add_field("n_counts", rng.integers(0, 100, n).astype("int64"), role="measure", span=["cells"],
             mask=(rng.random(n) < 0.1).astype("uint8"))                     # nullable
gf = ds.add_field("knn", sp.random(n, n, density=0.05, format="csr", random_state=3).astype("float32"),
                  role="graph", span=["cells", "cells"]); gf.directed = True; gf.weighted = True  # graph
cov = np.sort(rng.choice(n, n // 2, replace=False))
ds.add_field("adt", rng.standard_normal(n // 2).astype("float32"), role="measure", span=["cells"],
             index=cov, index_axis="cells")                                  # partial
ds.aux["anndata.uns"] = {"params": {"n_pca": 10, "method": "leiden"},
                         "scores": np.arange(5.0), "names": ["a", "b", "c"]}  # aux tree + leaves
ds = extend_for_viewer(ds)                                                   # od_score, stats_*, markers_*, factor axes
lstar.write(ds, out, compressor=numcodecs.GZip(5), format="v2")             # explicit v2 seed (re-emitted to v3 downstream)
man = __import__("json").load(open(out + "/.zattrs"))["lstar"]
print("seed:", len(man["fields"]), "fields,", len(man["axes"]), "axes ->", out)
