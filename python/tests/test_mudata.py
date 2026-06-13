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


def test_citeseq_mudata_roundtrip():
    md = corpus.citeseq_mudata()                               # REAL CITE-seq (subsampled minipbcite)
    if md is None:
        print("  SKIP test_citeseq_mudata_roundtrip (mudata/fixture unavailable)"); return
    import scipy.sparse as sp
    ds = lstar.read_mudata(md)
    # RNA -> genes, ADT/prot -> proteins, both over the shared cells axis
    assert ds.axis("genes").role == "feature" and ds.axis("proteins").role == "feature"
    rna_m = [n for n, f in ds.fields.items() if f.role == "measure" and f.span[1:] == ["genes"]]
    prot_m = [n for n, f in ds.fields.items() if f.role == "measure" and f.span[1:] == ["proteins"]]
    assert rna_m and prot_m and ds.field(rna_m[0]).span[0] == "cells"
    assert len(ds.axis("genes")) == 27 and len(ds.axis("proteins")) == 29
    assert not lstar.validate(ds)

    ds2 = lstar.read(_w(ds))                                    # round-trip through the store
    assert not lstar.validate(ds2)
    md2 = lstar.write_mudata(ds2)
    assert set(md2.mod) == {"rna", "prot"}
    assert md2.mod["rna"].shape == md.mod["rna"].shape and md2.mod["prot"].shape == md.mod["prot"].shape
    def dense(x): return x.toarray() if sp.issparse(x) else np.asarray(x)
    assert np.allclose(dense(md2.mod["prot"].X), dense(md.mod["prot"].X), rtol=1e-4)  # real ADT counts
    print("MuData (real CITE-seq RNA+ADT): modalities -> genes/proteins feature axes; round-trips exact")


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


def test_partial_overlap_mudata():
    """A modality measured on only a SUBSET of cells (10x multiome barcode whitelists, CITE-seq dropout)
    -> the measure is partial coverage over the *shared* cells axis (an index of covered positions), NOT
    a separate `cells.<mod>` axis and NOT zero/NA-padded. Round-trips to the cell subset on write-back."""
    try:
        import anndata as ad
        import mudata
        import pandas as pd
    except Exception:
        print("  SKIP test_partial_overlap_mudata (mudata unavailable)"); return
    import scipy.sparse as sp
    rng = np.random.default_rng(0)
    n = 100; cells = [f"c{i}" for i in range(n)]
    arna = ad.AnnData(sp.csr_matrix(rng.poisson(1, (n, 20)).astype("float32")),
                      obs=pd.DataFrame(index=cells), var=pd.DataFrame(index=[f"g{i}" for i in range(20)]))
    sub = sorted(rng.choice(n, 60, replace=False)); subcells = [cells[i] for i in sub]
    aprot = ad.AnnData(sp.csr_matrix(rng.poisson(2, (60, 8)).astype("float32")),
                       obs=pd.DataFrame(index=subcells), var=pd.DataFrame(index=[f"p{i}" for i in range(8)]))
    md = mudata.MuData({"rna": arna, "prot": aprot})

    ds = lstar.read_mudata(md)
    pm = [f for nm, f in ds.fields.items() if f.role == "measure" and f.span[1:] == ["proteins"]][0]
    assert pm.coverage == "partial" and pm.index_axis == "cells" and len(np.asarray(pm.index)) == 60
    assert len(ds.axis("cells")) == md.n_obs and "cells.prot" not in ds.axes   # shared axis, not a 2nd one
    assert not lstar.validate(ds)

    ds2 = lstar.read(_w(ds))                                    # through the store
    pm2 = [f for nm, f in ds2.fields.items() if f.coverage == "partial"][0]
    assert pm2.index is not None and not lstar.validate(ds2)
    md2 = lstar.write_mudata(ds2)                               # back to MuData on the covered subset
    assert md2.mod["rna"].n_obs == 100 and md2.mod["prot"].n_obs == 60
    assert set(md2.mod["prot"].obs_names) == set(subcells)
    print("MuData partial-overlap: prot on 60/100 cells -> partial coverage over the shared cells axis; round-trips")


def _w(ds):
    p = _store(); lstar.write(ds, p); return p


if __name__ == "__main__":
    test_citeseq_mudata_roundtrip()
    test_real_minipbcite()
    test_partial_overlap_mudata()
