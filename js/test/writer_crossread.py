"""Cross-language gate for the JS/WASM writer: read the store `writer_make.ts` wrote (chunked +
gzip-compressed, every encoding) with the **Python** lstar reader, assert it `validate()`s clean, and
check the values round-trip exactly. This is the JS-write -> Python-read leg the conformance matrix was
missing. Run after `node --experimental-strip-types js/test/writer_make.ts <dir>`.

Usage: PYTHONPATH=python/src python3 js/test/writer_crossread.py <store-dir>
"""
import sys

import numpy as np

import lstar
from lstar import Categorical

p = sys.argv[1] if len(sys.argv) > 1 else "/tmp/lstar-writer-cross.lstar.zarr"
ds = lstar.read(p)

errs = [e for e in lstar.validate(ds) if e.startswith("ERROR")]
assert not errs, ("validate", errs)

# CSC measure (chunked + gzip): reconstruct dense and check
import scipy.sparse as sp
counts = ds.field("counts")
assert counts.encoding == "csc" and counts.state == "raw"
m = sp.csc_matrix((counts.values.data, counts.values.indices, counts.values.indptr), shape=(10, 6)).toarray()
assert m[0, 0] == 1 and m[3, 0] == 2 and m[9, 4] == 8 and m[5, 5] == 0, m

# dense embedding (chunked + gzip)
umap = np.asarray(ds.field("umap").values)
assert umap.shape == (10, 2) and umap[1, 0] == 1.0 and umap[9, 1] == 9.5, umap

# categorical + induced factor axis
ct = ds.field("celltype").values
assert isinstance(ct, Categorical)
assert list(ct.codes) == [0, 1, 0, 2, 1, 0, -1, 2, 1, 0]
assert list(ct.categories) == ["T", "B", "NK"] and ct.ordered is False
ax = ds.axis("celltype")
assert ax.role == "factor" and ax.induced_by == "celltype" and list(np.asarray(ax.labels)) == ["T", "B", "NK"]

# nullable mask (1 == missing)
qc = ds.field("qc")
assert qc.mask is not None and list(np.asarray(qc.mask)) == [0, 0, 1, 0, 0, 0, 0, 1, 0, 0]
assert abs(float(np.asarray(qc.values)[0]) - 0.25) < 1e-6

# partial coverage
adt = ds.field("adt")
assert adt.coverage == "partial" and adt.index_axis == "cells"
assert list(np.asarray(adt.index)) == [0, 2, 4, 6, 8]
assert list(np.asarray(adt.values)) == [10, 11, 12, 13, 14]

# aux passthrough
aux = ds.aux["test.uns"]
assert aux["n_pca"] == 50 and aux["method"] == "leiden"
assert list(np.asarray(aux["scores"])) == [1.1, 2.2, 3.3]
assert list(aux["names"]) == ["foo", "bar"]

print("  [py] JS-written chunked+gzip store: every encoding round-trips + validate clean "
      "(csc, dense, categorical+factor, mask, partial, aux)")
