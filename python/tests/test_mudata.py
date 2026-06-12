"""MuData multimodal profile: modalities (RNA / ADT / ATAC) become canonical feature axes
(genes/proteins/peaks) over one shared `cells` axis -- the same L* shape as a Seurat multi-assay
object. Grounded in real CITE-seq (minipbcite.h5mu) + a constructed RNA+ADT MuData.

Run: PYTHONPATH=python/src python3 python/tests/test_mudata.py
"""
import os
import sys
import tempfile
import warnings

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
import corpus  # noqa: E402

import lstar  # noqa: E402

warnings.filterwarnings("ignore")


def _store():
    return os.path.join(tempfile.mkdtemp(), "mu.lstar.zarr")


def test_synthetic_mudata_roundtrip():
    md = corpus.synthetic_mudata()
    if md is None:
        print("  SKIP test_synthetic_mudata_roundtrip (mudata not installed)"); return
    import scipy.sparse as sp
    ds = lstar.read_mudata(md)
    # RNA -> genes, ADT/prot -> proteins, both over the shared cells axis
    assert ds.axis("genes").role == "feature" and ds.axis("proteins").role == "feature"
    rna_m = [n for n, f in ds.fields.items() if f.role == "measure" and f.span[1:] == ["genes"]]
    prot_m = [n for n, f in ds.fields.items() if f.role == "measure" and f.span[1:] == ["proteins"]]
    assert rna_m and prot_m and ds.field(rna_m[0]).span[0] == "cells"
    # MuData lifts a modality's obs into global obs with a `<mod>:` prefix -> the categorical induces it
    assert ds.axis("rna:leiden").role == "factor"
    assert not lstar.validate(ds)

    ds2 = lstar.read(_w(ds))                                    # round-trip through the store
    assert not lstar.validate(ds2)
    md2 = lstar.write_mudata(ds2)
    assert set(md2.mod) == {"rna", "prot"}
    assert md2.mod["rna"].shape == md.mod["rna"].shape and md2.mod["prot"].shape == md.mod["prot"].shape
    def dense(x): return x.toarray() if sp.issparse(x) else np.asarray(x)
    assert np.allclose(dense(md2.mod["prot"].X), dense(md.mod["prot"].X), rtol=1e-4)
    print("MuData (synthetic RNA+ADT): modalities -> genes/proteins feature axes; round-trips + reconstructs")


def test_real_minipbcite():
    md = corpus.minipbcite()
    if md is None:
        print("  SKIP test_real_minipbcite (mudata absent / download failed)"); return
    ds = lstar.read_mudata(md)
    assert {"genes", "proteins"} <= set(n for n, a in ds.axes.items() if a.role == "feature")
    assert len(ds.axis("cells")) == md.n_obs
    # real per-modality feature counts preserved
    assert len(ds.axis("genes")) == md.mod["rna"].n_vars
    assert len(ds.axis("proteins")) == md.mod["prot"].n_vars
    # real global categoricals (celltype, leiden, ...) induced factor axes
    assert ds.axis("celltype").role == "factor"
    # per-modality uns preserved in the passthrough
    assert any(k.startswith("mudata.") for k in ds.aux)
    assert not lstar.validate(ds)

    ds2 = lstar.read(_w(ds))                                    # full pipeline through the store
    md2 = lstar.write_mudata(ds2)
    assert set(md2.mod) == set(md.mod)
    assert md2.mod["rna"].n_vars == md.mod["rna"].n_vars and md2.mod["prot"].n_vars == md.mod["prot"].n_vars
    print("MuData (real minipbcite CITE-seq): %d cells, rna(%d)+prot(%d) -> shared cells + feature axes; "
          "round-trips" % (md.n_obs, md.mod["rna"].n_vars, md.mod["prot"].n_vars))


def _w(ds):
    p = _store(); lstar.write(ds, p); return p


if __name__ == "__main__":
    test_synthetic_mudata_roundtrip()
    test_real_minipbcite()
