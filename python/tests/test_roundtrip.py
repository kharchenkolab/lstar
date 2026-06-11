"""M1 round-trip test: build a small Sample dataset, write to Zarr, read back, compare.

Run standalone:  PYTHONPATH=python/src python3 python/tests/test_roundtrip.py
Or via pytest:   PYTHONPATH=python/src pytest python/tests/
"""
import os
import tempfile

import numpy as np
import scipy.sparse as sp

import lstar
from lstar.zarr_io import write, read


def make_ds():
    ds = lstar.Dataset(kind="sample")
    n_cells, n_genes, n_pc = 100, 50, 10
    ds.add_axis("cells", ["cell%d" % i for i in range(n_cells)])
    ds.add_axis("genes", ["g%d" % i for i in range(n_genes)])
    ds.add_axis("pca", ["PC%d" % i for i in range(n_pc)], origin="derived")

    rng = np.random.default_rng(0)
    counts = sp.random(n_cells, n_genes, density=0.1, format="csc", random_state=0)
    ds.add_field("counts", counts, role="measure", span=["cells", "genes"], state="raw")

    pca = rng.standard_normal((n_cells, n_pc)).astype("float32")
    ds.add_field("pca", pca, role="embedding", span=["cells", "pca"])

    leiden = np.array(["c%d" % (j % 5) for j in range(n_cells)])
    ds.add_field("leiden", leiden, role="label", span=["cells"])
    return ds


def run(path):
    ds = make_ds()
    write(ds, path)
    ds2 = read(path)

    assert set(ds2.axes) == set(ds.axes), (set(ds2.axes), set(ds.axes))
    assert set(ds2.fields) == set(ds.fields)
    assert (np.asarray(ds2.axis("cells").labels) == np.asarray(ds.axis("cells").labels)).all()
    assert ds2.axis("pca").origin == "derived"

    assert ds2.field("counts").encoding == "csc"
    assert ds2.field("counts").state == "raw"
    assert ds2.field("counts").role == "measure"
    assert (ds2.field("counts").values.toarray() == ds.field("counts").values.toarray()).all()

    assert ds2.field("pca").role == "embedding"
    assert (ds2.field("pca").values == ds.field("pca").values).all()

    assert ds2.field("leiden").role == "label"
    assert (ds2.field("leiden").values == ds.field("leiden").values).all()

    # inference check: omit role/span and confirm they are resolved
    ds3 = lstar.Dataset()
    ds3.add_axis("cells", ds.axis("cells").labels)
    ds3.add_axis("genes", ds.axis("genes").labels)
    f = ds3.add_field("counts2", make_ds().field("counts").values)  # span inferred by shape
    assert f.span == ["cells", "genes"], f.span
    assert f.role == "measure" and f.encoding == "csc"

    print("roundtrip OK:", ds2)
    return ds2


def test_roundtrip(tmp_path):
    run(str(tmp_path / "s.lstar.zarr"))


if __name__ == "__main__":
    d = tempfile.mkdtemp()
    run(os.path.join(d, "s.lstar.zarr"))
