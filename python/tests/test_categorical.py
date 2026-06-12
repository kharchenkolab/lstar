"""The `categorical` encoding: a label stored as integer codes + an ordered category set + a `-1`
missing sentinel, round-tripping with the order and missingness preserved. This is the Tier-1 dtype
fidelity fix and the substrate for factor axes (induction). Accepts both an L* Categorical and a
duck-typed pandas.Categorical.
"""
import os
import tempfile

import numpy as np

import lstar
from lstar import Categorical


def _store():
    return os.path.join(tempfile.mkdtemp(), "cat.lstar.zarr")


def test_categorical_roundtrip_ordered_and_missing():
    cats = np.array(["c0", "c1", "c2"])
    codes = np.array([0, 2, -1, 1, 1, -1, 0], dtype=np.int64)   # -1 = missing
    ds = lstar.Dataset(kind="sample")
    ds.add_axis("cells", [f"cell{i}" for i in range(len(codes))])
    ds.add_field("leiden", Categorical(codes, cats, ordered=True), span=["cells"])

    f = ds.field("leiden")
    assert f.role == "label" and f.encoding == "categorical"
    assert not lstar.validate(ds)

    p = _store()
    lstar.write(ds, p)
    rv = lstar.read(p).field("leiden").values
    assert isinstance(rv, Categorical)
    assert (rv.codes == codes).all()                 # codes preserved (incl. -1 missing)
    assert list(rv.categories) == list(cats)         # category set + order preserved
    assert rv.ordered is True
    # decoded view: missing -> ""
    dec = np.asarray(rv)
    assert list(dec) == ["c0", "c2", "", "c1", "c1", "", "c0"]
    print("categorical: codes / categories / ordered / -1 missing round-trip exact")


def test_categorical_from_pandas():
    import pandas as pd
    pc = pd.Categorical(["b", "a", "b", None, "c"], categories=["a", "b", "c"], ordered=False)
    ds = lstar.Dataset(kind="sample")
    ds.add_axis("cells", [f"cell{i}" for i in range(5)])
    ds.add_field("ct", pc, span=["cells"])           # pandas.Categorical accepted directly
    assert ds.field("ct").encoding == "categorical"
    p = _store()
    lstar.write(ds, p)
    rv = lstar.read(p).field("ct").values
    assert list(rv.categories) == ["a", "b", "c"]
    assert (rv.codes == np.asarray(pc.codes)).all()  # pandas -1 missing preserved
    print("pandas.Categorical -> L* categorical round-trips (codes/categories/missing)")


def test_categorical_validate_catches_bad_codes():
    ds = lstar.Dataset(kind="sample")
    ds.add_axis("cells", [f"c{i}" for i in range(3)])
    ds.add_field("g", Categorical(np.array([0, 1, 5]), np.array(["a", "b"])), span=["cells"])  # 5 >= 2
    issues = lstar.validate(ds)
    assert any("codes out of range" in i for i in issues), issues
    print("validate flags out-of-range categorical codes")


if __name__ == "__main__":
    test_categorical_roundtrip_ordered_and_missing()
    test_categorical_from_pandas()
    test_categorical_validate_catches_bad_codes()
