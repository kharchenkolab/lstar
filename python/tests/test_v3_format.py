"""Zarr v3 on-disk format from the Python surface. lstar.write(..., format="v3") emits genuine v3
(zarr.json + inline-consolidated metadata) via the zarr-python 3 library, and lstar.read auto-detects
v2 vs v3. The two formats hold identical values, so the Python reference reads and writes both -- and
gives lstar a second, independent v3 writer (zarr-python) alongside the libzarr C++/R/JS cores.

Run standalone:  PYTHONPATH=python/src python3 python/tests/test_v3_format.py
"""
import json
import os
import tempfile

import numpy as np
import scipy.sparse as sp
import pandas as pd

import lstar
from lstar.zarr_io import write, read


def make_ds():
    rng = np.random.default_rng(0)
    n, g = 80, 30
    ds = lstar.Dataset(kind="sample")
    ds.add_axis("cells", ["c%d" % i for i in range(n)])
    ds.add_axis("genes", ["g%d" % i for i in range(g)])
    ds.add_axis("pc", ["PC%d" % i for i in range(8)], origin="derived")
    C = sp.random(n, g, density=0.3, format="csc", random_state=1)
    C.data = np.ceil(C.data * 10).astype("float32")
    ds.add_field("counts", C, role="measure", span=["cells", "genes"], state="raw")
    ds.add_field("pca", rng.standard_normal((n, 8)).astype("float32"), role="embedding", span=["cells", "pc"])
    ds.add_field("leiden", pd.Categorical(["cl%d" % (i % 4) for i in range(n)]), role="label", span=["cells"])
    ds.add_field("barcode", np.array(["BC%03d" % i for i in range(n)]), role="label", span=["cells"])
    ds.add_field("n_counts", rng.integers(0, 50, n).astype("int64"), role="measure", span=["cells"],
                 mask=(rng.random(n) < 0.1).astype("uint8"))
    return ds


def _field_vals(ds, name):
    v = ds.field(name).values
    if sp.issparse(v):
        return sp.csc_matrix(v).toarray()
    if hasattr(v, "codes"):
        return np.asarray(v.codes)
    return np.asarray(v)


def main():
    import numcodecs
    ds = make_ds()
    with tempfile.TemporaryDirectory() as d:
        v2, v3 = os.path.join(d, "s2.lstar.zarr"), os.path.join(d, "s3.lstar.zarr")
        write(ds, v2, compressor=numcodecs.GZip(5), format="v2")
        write(ds, v3, compressor=numcodecs.GZip(5), format="v3")

        # v3 is genuine v3; v2 is genuine v2
        rj = json.load(open(os.path.join(v3, "zarr.json")))
        assert rj["zarr_format"] == 3 and rj["node_type"] == "group", rj
        assert "consolidated_metadata" in rj, "v3 store must carry inline consolidated metadata"
        assert os.path.exists(os.path.join(v2, ".zmetadata")), "v2 store must carry .zmetadata"
        assert not os.path.exists(os.path.join(v2, "zarr.json")), "v2 store must not have zarr.json"

        # read both back (auto-detected) and compare values across every encoding
        d2, d3 = read(v2), read(v3)
        assert list(d2.fields) == list(d3.fields) == list(ds.fields)
        assert list(d2.axes) == list(d3.axes)
        def arreq(a, b):
            return a.shape == b.shape and (np.array_equal(a, b, equal_nan=True)
                                           if a.dtype.kind == "f" else np.array_equal(a, b))
        for name in ds.fields:
            assert arreq(_field_vals(d2, name), _field_vals(d3, name)), name
        # strings + mask
        assert list(d3.field("barcode").values) == list(d2.field("barcode").values)
        assert d3.field("leiden").values.categories.tolist() == d2.field("leiden").values.categories.tolist()
        # no validation errors on the v3 read
        assert not [e for e in lstar.validate(d3) if e.startswith("ERROR")]

        # invalid format rejected
        try:
            write(ds, os.path.join(d, "bad"), format="v4")
            assert False, "format='v4' should raise"
        except ValueError:
            pass

    print("v3 format OK: Python writes + reads v2 and v3; values identical across formats")


if __name__ == "__main__":
    main()
