"""Property/fuzz tests for the encodings: randomized fields (dense / CSR / CSC / categorical / nullable)
+ deliberately nasty edge cases (empty axis, 0-nnz, all-missing mask, single element, int64 values past
2^31 -> JS BigInt boundary) round-trip through the store byte-faithfully and validate. Curated fixtures
test the typical shape; this sweeps the corners they miss.

Run: PYTHONPATH=python/src python3 python/tests/test_fuzz.py
"""
import os
import sys
import tempfile

import numpy as np
import scipy.sparse as sp

sys.path.insert(0, os.path.dirname(__file__))
import lstar  # noqa: E402


def _store():
    return os.path.join(tempfile.mkdtemp(), "fuzz.lstar.zarr")


def _roundtrip(ds):
    p = _store(); lstar.write(ds, p)
    assert not [i for i in lstar.validate(ds) if i.startswith("ERROR")], lstar.validate(ds)
    ds2 = lstar.read(p)
    assert not [i for i in lstar.validate(ds2) if i.startswith("ERROR")]
    return ds2


def test_random_encodings_roundtrip():
    rng = np.random.default_rng(12345)
    n_ok = 0
    for trial in range(40):
        nc = int(rng.integers(0, 40))          # includes 0 -> empty axis
        ng = int(rng.integers(0, 25))
        ds = lstar.Dataset(kind="sample")
        ds.add_axis("cells", [f"c{i}" for i in range(nc)])
        ds.add_axis("genes", [f"g{i}" for i in range(ng)])
        kind = rng.integers(0, 4)
        if kind == 0 and nc and ng:            # sparse CSC/CSR (density includes 0-nnz)
            fmt = "csc" if rng.random() < 0.5 else "csr"
            dens = float(rng.choice([0.0, 0.1, 0.6]))
            M = sp.random(nc, ng, density=dens, format=fmt, random_state=int(rng.integers(1e6))).astype("float32")
            ds.add_field("m", M, role="measure", span=["cells", "genes"], state="raw")
            ds2 = _roundtrip(ds)
            assert np.allclose(ds2.field("m").values.toarray(), M.toarray())
        elif kind == 1 and nc:                 # dense vector, with NaN
            v = rng.normal(size=nc).astype("float32")
            if nc:
                v[rng.integers(0, nc)] = np.nan
            ds.add_field("v", v, role="measure", span=["cells"])
            ds2 = _roundtrip(ds)
            assert np.allclose(np.asarray(ds2.field("v").values), v, equal_nan=True)
        elif kind == 2 and nc:                 # nullable int64 with a possibly all-missing mask + big values
            vals = (rng.integers(0, 5_000_000_000, nc)).astype(np.int64)   # > 2^31 (JS BigInt boundary)
            mask = (rng.random(nc) < rng.choice([0.0, 0.5, 1.0])).astype(np.uint8)  # incl. none/all missing
            ds.add_field("n", vals, role="measure", span=["cells"], mask=mask)
            ds2 = _roundtrip(ds)
            f = ds2.field("n")
            assert f.values.dtype.kind == "i" and (np.asarray(f.values) == vals).all()
            assert f.mask is not None and (np.asarray(f.mask) == mask).all()
        elif nc:                               # categorical incl. -1 missing + a single-category case
            k = int(rng.integers(1, 4))
            codes = rng.integers(-1, k, nc).astype(np.int64)
            cats = np.array([f"L{j}" for j in range(k)])
            ds.add_field("ct", lstar.Categorical(codes, cats), span=["cells"])
            ds2 = _roundtrip(ds)
            c2 = ds2.field("ct").values
            assert (np.asarray(c2.codes) == codes).all() and list(c2.categories) == list(cats)
        else:
            _roundtrip(ds)                     # empty axes / no field: must still validate + round-trip
        n_ok += 1
    print("fuzz: %d randomized encoding trials (empty axes, 0-nnz, all-missing masks, int64>2^31, "
          "-1 categoricals) round-trip + validate" % n_ok)


if __name__ == "__main__":
    test_random_encodings_roundtrip()
