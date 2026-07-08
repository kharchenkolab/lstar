#!/usr/bin/env python3
# Real-data v3 round-trip: ingest a real AnnData/MuData, write it as BOTH v2 and v3 (uncompressed — the
# format-equality check is orthogonal to the codec), read both back through the Python reference, and assert
# the values are identical across formats + no validation errors. Then two asymmetric checks on real data: a
# v3->v2 downgrade, and a v3+zstd+SHARDED copy (the codec/sharding coverage). Emits the two store paths so
# the caller can also cross-read them with the C++ core (test_v3).
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

# Total DECOMPRESSED field bytes — the size that predicts the JS whole-array read's memory (io_dump
# base64s the decompressed arrays into one JSON; wasm_corpus materializes each in the WASM heap). The
# on-disk store size is compressed and a poor proxy; the caller gates the JS leg on this.
def _field_bytes(f):
    v = f.values
    if sp.issparse(v):
        return int(v.data.nbytes) + int(v.indices.nbytes) + int(v.indptr.nbytes)
    return int(np.asarray(v).nbytes) if not hasattr(v, "codes") else int(np.asarray(v.codes).nbytes)
print("DECOMP_MB=%d" % (sum(_field_bytes(f) for f in ds.fields.values()) // 1_000_000))

import numcodecs
# The v2==v3 equality check is about FORMAT, not the codec — compression is an orthogonal axis, covered by
# the zstd+shard leg below (real data) and synthetic chunked.sh. Write these two UNCOMPRESSED: on the large
# perturbation stores gzip-5 encode (~23 MB/s, single-threaded) dominated the sweep's wall-clock; none is ~7x.
lstar.write(ds, out_v2, format="v2")
lstar.write(ds, out_v3, format="v3")

# genuine formats
assert json.load(open(os.path.join(out_v3, "zarr.json")))["zarr_format"] == 3
assert os.path.exists(os.path.join(out_v2, ".zmetadata"))

d2, d3 = lstar.read(out_v2), lstar.read(out_v3)
assert list(d2.fields) == list(d3.fields), "field set differs across formats"
assert list(d2.axes) == list(d3.axes), "axis set differs across formats"

def _canon(v):
    # Comparable, MEMORY-CHEAP view of a field. A sparse matrix becomes a canonical CSC triplet (dedup +
    # drop explicit zeros + sorted indices) -- it is NEVER densified: a 41786^2 slide-seq cell graph is
    # 14 GB dense, and densify-to-compare made the sweep thrash (34 GB RSS, ~30 min) on one dataset. The
    # triplet compare is value-equivalent to a dense compare. Categorical -> codes; object/unicode -> str.
    if sp.issparse(v):
        m = sp.csc_matrix(v); m.sum_duplicates(); m.eliminate_zeros(); m.sort_indices()
        return ("csc", m.shape, m.indptr, m.indices, m.data)
    if hasattr(v, "codes"):
        return ("arr", np.asarray(v.codes))
    a = np.asarray(v)
    return ("arr", a if a.dtype.kind not in ("U", "S", "O") else a.astype(str))

def _eq(x, y):
    if x[0] != y[0]:
        return False                                        # sparse-vs-dense encoding mismatch is a real diff
    if x[0] == "csc":
        _, sh, ip, ix, da = x; _, sh2, ip2, ix2, da2 = y
        if sh != sh2 or not np.array_equal(ip, ip2) or not np.array_equal(ix, ix2):
            return False
    else:
        da, da2 = x[1], y[1]
        if da.shape != da2.shape:
            return False
    return np.array_equal(da, da2, equal_nan=True) if da.dtype.kind == "f" else np.array_equal(da, da2)

def cmp_fields(dsA, dsB, label):
    n = 0
    for name in ds.fields:
        assert _eq(_canon(dsA.field(name).values), _canon(dsB.field(name).values)), f"{label} {name}: values differ"
        n += 1
    return n

nfields = cmp_fields(d2, d3, "v2==v3")

# asymmetric on real data: (1) v3 -> v2 DOWNGRADE preserves values; (2) a v3 + zstd + SHARDED copy reads
# == the v2 store — real-data exercise of the downgrade direction + zstd + sharding (Python round-trips, so
# nullable masked values are preserved; no cross-surface 0-fill here).
down = out_v3 + ".down_v2"
lstar.write(d3, down, format="v2")                       # uncompressed: this leg tests the downgrade, not the codec
cmp_fields(lstar.read(down), d2, "v3->v2 downgrade")

zshard = out_v3 + ".zstd_shard"
lstar.write(ds, zshard, compressor=numcodecs.Zstd(3), format="v3", chunk_elems=50_000, shard_elems=400_000)
cmp_fields(lstar.read(zshard), d2, "v3+zstd+shard")

errs = [e for e in lstar.validate(d3) if e.startswith("ERROR")]
assert not errs, f"v3 validation errors: {errs[:3]}"
print(f"  OK  {os.path.basename(src)}: {nfields} fields, {len(list(d3.axes))} axes -> v2==v3, v3->v2 downgrade, v3+zstd+shard all equal")
