"""Nullable / extension dtypes: a field may carry an explicit `uint8` validity mask (`1 == missing`)
beside its values, so pandas nullable `Int64`/`boolean`/`string` columns round-trip with their
integer-ness and value-vs-missing distinction intact (coercing them to float-NaN is a silent
corruption, the same class P1 fixed for categoricals). Float keeps NaN (no mask).

Run: PYTHONPATH=python/src python3 python/tests/test_nullable.py
"""
import os
import tempfile

import numpy as np

import lstar


def _store():
    return os.path.join(tempfile.mkdtemp(), "nullable.lstar.zarr")


def test_mask_roundtrip_and_validate():
    ds = lstar.Dataset(kind="sample")
    ds.add_axis("cells", [f"c{i}" for i in range(5)])
    vals = np.array([10, 0, 7, 0, 3], dtype=np.int64)
    mask = np.array([0, 1, 0, 1, 0], dtype=np.uint8)        # positions 1,3 are missing (not zero)
    ds.add_field("n_counts", vals, role="measure", span=["cells"], mask=mask)
    assert not lstar.validate(ds)

    p = _store(); lstar.write(ds, p)
    f = lstar.read(p).field("n_counts")
    assert f.mask is not None and (f.mask == mask).all()
    assert (np.asarray(f.values) == vals).all()
    assert f.values.dtype.kind == "i"                       # integer-ness preserved (not float-NaN)
    print("nullable: integer values + validity mask round-trip exact (missing != 0)")


def test_validate_catches_bad_mask_length():
    ds = lstar.Dataset(kind="sample")
    ds.add_axis("cells", [f"c{i}" for i in range(4)])
    ds.add_field("x", np.arange(4), role="measure", span=["cells"],
                 mask=np.array([0, 1, 0], dtype=np.uint8))  # length 3 != 4
    assert any("mask length" in i for i in lstar.validate(ds))
    print("validate catches a mask whose length != axis length")


def test_anndata_nullable_columns():
    import anndata as ad
    import pandas as pd
    import scipy.sparse as sp
    from lstar import read_anndata, write_anndata, write, read

    n = 6
    obs = pd.DataFrame({
        "n_umi":   pd.array([100, None, 300, 400, None, 600], dtype="Int64"),
        "flagged": pd.array([True, None, False, True, None, False], dtype="boolean"),
        "donor":   pd.array(["d1", "d2", None, "d1", None, "d3"], dtype="string"),
        "frac":    pd.array([0.1, 0.2, None, 0.4, 0.5, 0.6], dtype="Float64"),
    }, index=[f"cell{i}" for i in range(n)])
    a = ad.AnnData(X=sp.random(n, 4, density=0.3, format="csr", random_state=0), obs=obs)

    ds = read_anndata(a)
    # the three nullable columns carry a mask; float keeps NaN (no mask)
    assert ds.field("n_umi").mask is not None and ds.field("n_umi").values.dtype.kind == "i"
    assert ds.field("flagged").mask is not None
    assert ds.field("donor").mask is not None
    assert ds.field("frac").mask is None
    assert not lstar.validate(ds)

    p = _store(); write(ds, p)
    a2 = write_anndata(read(p))                              # through the L* store + back to AnnData
    # nullable dtypes + exact NA positions reconstructed
    assert str(a2.obs["n_umi"].dtype) == "Int64"
    assert list(a2.obs["n_umi"].isna()) == list(a.obs["n_umi"].isna())
    assert list(a2.obs["n_umi"].dropna()) == list(a.obs["n_umi"].dropna())
    assert str(a2.obs["flagged"].dtype) == "boolean"
    assert list(a2.obs["flagged"].isna()) == list(a.obs["flagged"].isna())
    assert str(a2.obs["donor"].dtype) == "string"
    assert list(a2.obs["donor"].isna()) == list(a.obs["donor"].isna())
    assert list(a2.obs["donor"].dropna()) == list(a.obs["donor"].dropna())
    print("anndata nullable Int64/boolean/string round-trip type-faithful (NA positions + values exact)")


if __name__ == "__main__":
    test_mask_roundtrip_and_validate()
    test_validate_catches_bad_mask_length()
    test_anndata_nullable_columns()
