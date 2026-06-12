"""Nullable / extension dtypes: a field may carry an explicit `uint8` validity mask (`1 == missing`)
beside its values, so pandas nullable `Int64`/`boolean`/`string` columns round-trip with their
integer-ness and value-vs-missing distinction intact (coercing them to float-NaN is a silent
corruption, the same class P1 fixed for categoricals). Float keeps NaN (no mask).

Run: PYTHONPATH=python/src python3 python/tests/test_nullable.py
"""
import os
import sys
import tempfile

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
import corpus  # noqa: E402

import lstar  # noqa: E402


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
    """Grounded in a REAL dataset (pbmc3k): no public h5ad reliably ships pandas nullable extension
    dtypes, so -- the documented path -- we derive nullable obs columns from real obs *values* (real
    n_genes counts, a real QC threshold for missingness) and check the nullable Int64/boolean/string
    survive read->L*->write type-faithfully."""
    import warnings; warnings.filterwarnings("ignore")
    import pandas as pd
    from lstar import read_anndata, write_anndata, write, read

    a = corpus.pbmc3k_processed()
    if a is None:
        print("  SKIP test_anndata_nullable_columns (corpus unavailable)"); return
    a = a.copy()
    ng = np.asarray(a.obs["n_genes"], dtype=float)              # real per-cell gene counts
    miss = ng < np.percentile(ng, 5)                           # the lowest-QC cells: treat as missing
    a.obs["n_genes_nullable"] = pd.array(np.where(miss, pd.NA, ng.astype("int64")), dtype="Int64")
    a.obs["passes_qc"] = pd.array(np.where(miss, pd.NA, ng > np.median(ng)), dtype="boolean")
    donor = np.asarray(a.obs["louvain"].astype(str))           # real labels, some set missing
    a.obs["donor"] = pd.array(np.where(miss, pd.NA, donor), dtype="string")

    ds = read_anndata(a)
    assert ds.field("n_genes_nullable").mask is not None and ds.field("n_genes_nullable").values.dtype.kind == "i"
    assert ds.field("passes_qc").mask is not None and ds.field("donor").mask is not None
    assert ds.field("percent_mito").mask is None              # a real float column: NaN, no mask
    assert not lstar.validate(ds)

    a2 = write_anndata(read(_w(ds)))                           # through the store + back to AnnData
    assert str(a2.obs["n_genes_nullable"].dtype) == "Int64"
    assert list(a2.obs["n_genes_nullable"].isna()) == list(a.obs["n_genes_nullable"].isna())
    assert list(a2.obs["n_genes_nullable"].dropna()) == list(a.obs["n_genes_nullable"].dropna())
    assert str(a2.obs["passes_qc"].dtype) == "boolean"
    assert list(a2.obs["passes_qc"].isna()) == list(a.obs["passes_qc"].isna())
    assert str(a2.obs["donor"].dtype) == "string"
    assert list(a2.obs["donor"].isna()) == list(a.obs["donor"].isna())
    assert list(a2.obs["donor"].dropna()) == list(a.obs["donor"].dropna())
    print("anndata nullable Int64/boolean/string (real pbmc3k values): round-trip type-faithful")


def _w(ds):
    p = _store(); lstar.write(ds, p); return p


if __name__ == "__main__":
    test_mask_roundtrip_and_validate()
    test_validate_catches_bad_mask_length()
    test_anndata_nullable_columns()
