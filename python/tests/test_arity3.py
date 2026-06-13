"""Arity-3 fields — lstar's "beyond cells × genes" claim, exercised. A field can span **three** axes:
a cell–cell **communication** tensor (`senders × receivers × lr_pairs`, the LIANA/CellPhoneDB shape) and
an **eQTL** table (`celltypes × genes × variants`). These have no native slot in AnnData/Seurat/SCE (they
go to the sidecar there), but the L* store holds them directly — typed, validated, and round-tripping
(including across languages: the R reconstruction used to collapse the 3rd axis to 2-D).

Run: PYTHONPATH=python/src python3 python/tests/test_arity3.py
"""
import os
import sys
import tempfile

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
import lstar  # noqa: E402


def _w(ds):
    p = os.path.join(tempfile.mkdtemp(), "a3.lstar.zarr"); lstar.write(ds, p); return p


def test_ccc_communication_tensor():
    rng = np.random.default_rng(0)
    groups = ["T", "B", "Mono", "NK"]
    lr = ["CD40LG_CD40", "TNF_TNFRSF1A", "IL2_IL2RA"]
    ds = lstar.Dataset(kind="sample")
    ds.add_axis("senders", groups); ds.add_axis("receivers", groups); ds.add_axis("lr_pairs", lr)
    score = rng.random((4, 4, 3)).astype("float32")
    pval = rng.random((4, 4, 3)).astype("float32")
    ds.add_field("ccc_score", score, role="measure", span=["senders", "receivers", "lr_pairs"],
                 subtype="communication")
    ds.add_field("ccc_pval", pval, role="measure", span=["senders", "receivers", "lr_pairs"],
                 subtype="communication", uncertainty="pval")
    assert not lstar.validate(ds)
    f = lstar.read(_w(ds)).field("ccc_score")
    assert f.span == ["senders", "receivers", "lr_pairs"]
    assert np.asarray(f.values).shape == (4, 4, 3) and np.allclose(f.values, score)
    print("arity-3 CCC tensor (senders×receivers×lr_pairs): typed measure + uncertainty; round-trips exact")


def test_eqtl_tensor():
    rng = np.random.default_rng(1)
    nct, ng, nv = 5, 8, 6
    ds = lstar.Dataset(kind="sample")
    ds.add_axis("celltypes", [f"ct{i}" for i in range(nct)])
    ds.add_axis("genes", [f"g{i}" for i in range(ng)])
    ds.add_axis("variants", [f"rs{i}" for i in range(nv)])
    beta = rng.normal(size=(nct, ng, nv)).astype("float32")
    ds.add_field("eqtl_beta", beta, role="measure", span=["celltypes", "genes", "variants"], subtype="eqtl")
    assert not lstar.validate(ds)
    f = lstar.read(_w(ds)).field("eqtl_beta")
    assert f.span == ["celltypes", "genes", "variants"] and np.asarray(f.values).shape == (nct, ng, nv)
    assert np.allclose(f.values, beta)
    print("arity-3 eQTL tensor (celltypes×genes×variants): typed measure; round-trips exact")


def test_validate_catches_arity3_shape_mismatch():
    ds = lstar.Dataset(kind="sample")
    ds.add_axis("a", ["a0", "a1"]); ds.add_axis("b", ["b0", "b1", "b2"]); ds.add_axis("c", ["c0"])
    ds.add_field("t", np.zeros((2, 3, 2), dtype="float32"), role="measure", span=["a", "b", "c"])  # c is len 1
    assert any("axis 'c'" in i for i in lstar.validate(ds))
    print("validate catches an arity-3 field whose dim disagrees with its axis length")


if __name__ == "__main__":
    test_ccc_communication_tensor()
    test_eqtl_tensor()
    test_validate_catches_arity3_shape_mismatch()
