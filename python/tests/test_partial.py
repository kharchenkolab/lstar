"""Partial-coverage index arrays: a field may cover only a *subset* of a span axis, keyed by an `index`
of integer positions into that axis (`coverage="partial"`). This is the faithful representation of
partial overlap -- a modality measured on only some cells (10x multiome's per-modality barcode
whitelists, CITE-seq cells that fail one assay) -- without a separate `cells.<mod>` axis and without
zero/NA-padding the matrix to the full axis.

Run: PYTHONPATH=python/src python3 python/tests/test_partial.py
"""
import os
import sys
import tempfile

import numpy as np
import scipy.sparse as sp

sys.path.insert(0, os.path.dirname(__file__))
import lstar  # noqa: E402


def _store():
    return os.path.join(tempfile.mkdtemp(), "p.lstar.zarr")


def test_partial_measure_roundtrip():
    ds = lstar.Dataset(kind="sample")
    ds.add_axis("cells", [f"c{i}" for i in range(100)])      # 100 cells in the shared observation axis
    ds.add_axis("proteins", [f"p{i}" for i in range(8)])
    # the protein modality was measured on only 60 of the 100 cells -> a partial measure over (cells, proteins)
    covered = np.sort(np.random.default_rng(0).choice(100, 60, replace=False))
    vals = sp.random(60, 8, density=0.3, format="csr").astype(np.float32)
    ds.add_field("adt", vals, role="measure", span=["cells", "proteins"], state="raw",
                 index=covered, index_axis="cells")
    f = ds.field("adt")
    assert f.coverage == "partial" and f.index_axis == "cells"
    assert np.asarray(f.values).shape[0] == 60 if not sp.issparse(f.values) else f.values.shape[0] == 60
    assert not lstar.validate(ds)                            # partial shape (60) is valid, not 100

    p = _store(); lstar.write(ds, p)
    f2 = lstar.read(p).field("adt")
    assert f2.coverage == "partial" and f2.index_axis == "cells"
    assert (np.asarray(f2.index) == covered).all()           # the index round-trips
    assert f2.values.shape == (60, 8)
    assert np.allclose(f2.values.toarray(), vals.toarray())
    assert not lstar.validate(lstar.read(p))
    print("partial coverage: a (60-of-100 cells x 8 proteins) measure round-trips with its index; validates")


def test_validate_catches_out_of_range_index():
    ds = lstar.Dataset(kind="sample")
    ds.add_axis("cells", [f"c{i}" for i in range(10)])
    ds.add_axis("genes", [f"g{i}" for i in range(5)])
    bad = np.array([0, 3, 99])                                # 99 >= 10 -> out of range
    ds.add_field("x", sp.random(3, 5, density=0.5, format="csr"),
                 role="measure", span=["cells", "genes"], index=bad, index_axis="cells")
    issues = lstar.validate(ds)
    assert any("partial index out of range" in i for i in issues), issues
    print("validate catches a partial index that points outside its axis")


def test_validate_catches_bad_index_axis():
    ds = lstar.Dataset(kind="sample")
    ds.add_axis("cells", [f"c{i}" for i in range(10)])
    ds.add_axis("genes", [f"g{i}" for i in range(5)])
    ds.add_field("x", sp.random(3, 5, density=0.5, format="csr"),
                 role="measure", span=["cells", "genes"], index=np.array([0, 1, 2]), index_axis="samples")
    assert any("index_axis 'samples' not in span" in i for i in lstar.validate(ds))
    print("validate catches an index_axis that isn't in the field's span")


if __name__ == "__main__":
    test_partial_measure_roundtrip()
    test_validate_catches_out_of_range_index()
    test_validate_catches_bad_index_axis()
