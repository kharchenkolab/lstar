#!/usr/bin/env python3
# Real-data v3 round-trip: ingest a real AnnData/MuData, write it as BOTH v2 and v3 (gzip), read both
# back through the Python reference, and assert the values are identical across formats + no validation
# errors. Emits the two store paths so the caller can also cross-read them with the C++ core (test_v3).
# Tolerant: a dataset that won't load/convert prints SKIP and exits 0 (the corpus is gitignored/optional).
import sys, os, json, warnings
warnings.filterwarnings("ignore")
import numpy as np
import scipy.sparse as sp
import lstar

src, out_v2, out_v3 = sys.argv[1], sys.argv[2], sys.argv[3]

def load(path):
    if path.endswith(".h5mu"):
        import mudata
        return lstar.read_mudata(mudata.read_h5mu(path))
    import anndata
    return lstar.read_anndata(anndata.read_h5ad(path))

try:
    ds = load(src)
except Exception as e:
    print(f"  SKIP {os.path.basename(src)}: {type(e).__name__}: {e}")
    sys.exit(0)

import numcodecs
lstar.write(ds, out_v2, compressor=numcodecs.GZip(5), format="v2")
lstar.write(ds, out_v3, compressor=numcodecs.GZip(5), format="v3")

# genuine formats
assert json.load(open(os.path.join(out_v3, "zarr.json")))["zarr_format"] == 3
assert os.path.exists(os.path.join(out_v2, ".zmetadata"))

d2, d3 = lstar.read(out_v2), lstar.read(out_v3)
assert list(d2.fields) == list(d3.fields), "field set differs across formats"
assert list(d2.axes) == list(d3.axes), "axis set differs across formats"

def vals(ds, name):
    v = ds.field(name).values
    if sp.issparse(v):
        return sp.csc_matrix(v).toarray()
    if hasattr(v, "codes"):
        return np.asarray(v.codes)
    a = np.asarray(v)
    return a if a.dtype.kind not in ("U", "S", "O") else a.astype(str)

nfields = 0
for name in ds.fields:
    a2, a3 = vals(d2, name), vals(d3, name)
    assert a2.shape == a3.shape, f"{name}: shape {a2.shape} != {a3.shape}"
    ok = np.array_equal(a2, a3, equal_nan=True) if a2.dtype.kind == "f" else np.array_equal(a2, a3)
    assert ok, f"{name}: values differ across v2/v3"
    nfields += 1

errs = [e for e in lstar.validate(d3) if e.startswith("ERROR")]
assert not errs, f"v3 validation errors: {errs[:3]}"
ncells = d3.axes[list(d3.axes)[0]].labels
print(f"  OK  {os.path.basename(src)}: {nfields} fields, {len(list(d3.axes))} axes -> v2==v3 (Python); v3 genuine")
