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


def test_joint_method_storage_shapes():
    """How joint integration results are *stored/shaped* (we don't run the methods, we type their output):
      - WNN: a joint embedding (global `obsm`), per-cell modality weights (global `obs`), a joint graph
             (global `obsp`).
      - MOFA+: shared factor scores (global `obsm`) + **per-modality loadings** (`mod.varm`) that share
             **one factor axis** with the scores -- the lstar induction shape (one coord axis; an
             embedding over (cells, factor) + loadings over (feature, factor) per modality)."""
    try:
        import anndata as ad
        import mudata
        import pandas as pd
    except Exception:
        print("  SKIP test_joint_method_storage_shapes (mudata unavailable)"); return
    import scipy.sparse as sp
    rng = np.random.default_rng(1)
    n, ng, npr, k = 80, 20, 6, 5
    cells = [f"c{i}" for i in range(n)]
    arna = ad.AnnData(sp.csr_matrix(rng.poisson(1, (n, ng)).astype("float32")),
                      obs=pd.DataFrame(index=cells), var=pd.DataFrame(index=[f"g{i}" for i in range(ng)]))
    aprot = ad.AnnData(sp.csr_matrix(rng.poisson(2, (n, npr)).astype("float32")),
                       obs=pd.DataFrame(index=cells), var=pd.DataFrame(index=[f"p{i}" for i in range(npr)]))
    arna.varm["factors"] = rng.normal(size=(ng, k)).astype("float32")    # MOFA per-mod loadings (RNA)
    aprot.varm["factors"] = rng.normal(size=(npr, k)).astype("float32")  # MOFA per-mod loadings (ADT)
    aprot.layers["denoised"] = sp.csr_matrix(rng.random((n, npr)).astype("float32"))  # totalVI-style denoised
    md = mudata.MuData({"rna": arna, "prot": aprot})
    md.obsm["X_factors"] = rng.normal(size=(n, k)).astype("float32")     # MOFA shared factor scores
    md.obsm["X_wnn"] = rng.normal(size=(n, 2)).astype("float32")         # WNN joint embedding
    md.obs["rna:mod_weight"] = rng.random(n).astype("float32")           # WNN per-cell modality weights
    md.obs["prot:mod_weight"] = rng.random(n).astype("float32")
    md.obsp["wnn_connectivities"] = sp.random(n, n, density=0.1, format="csr").astype("float32")  # joint graph

    ds = lstar.read_mudata(md)
    # WNN joint embedding + modality weights + joint graph
    assert ds.field("X_wnn").role == "embedding" and ds.field("X_wnn").span == ["cells", "wnn"]
    assert ds.field("rna:mod_weight").role == "measure" and ds.field("rna:mod_weight").span == ["cells"]
    assert ds.field("wnn_connectivities").role == "relation" and ds.field("wnn_connectivities").span == ["cells", "cells"]
    # MOFA: ONE shared `factors` coord axis carries both the scores embedding and the per-mod loadings
    assert ds.axis("factors").role == "coordinate" and len(ds.axis("factors")) == k
    assert ds.field("X_factors").span == ["cells", "factors"]
    rl = ds.field("rna_factors_loadings"); pl = ds.field("prot_factors_loadings")
    assert rl.role == "loading" and rl.span == ["genes", "factors"]
    assert pl.role == "loading" and pl.span == ["proteins", "factors"]
    # facet-set provenance (S5): the factor-scores embedding records the feature axes that fed it
    assert ds.field("X_factors").provenance.get("input_axes") == ["genes", "proteins"]
    assert not lstar.validate(ds)
    ds2 = lstar.read(_w(ds))                                              # survives the store
    assert not lstar.validate(ds2)
    assert ds2.field("rna_factors_loadings").span == ["genes", "factors"]
    md2 = lstar.write_mudata(ds2)                                         # and back to a MuData
    assert "X_wnn" in md2.obsm and "X_factors" in md2.obsm
    assert "wnn_connectivities" in md2.obsp
    assert "factors" in md2.mod["rna"].varm and md2.mod["rna"].varm["factors"].shape == (ng, k)
    print("joint-method shapes: WNN (embedding + cell-weight measures + joint-graph relation) + MOFA "
          "(one shared factor axis: scores embedding + per-mod loadings) typed; round-trips to MuData")


def test_per_modality_pca_different_dim():
    """Regression for the sweep-caught minipbcite bug: two modalities each carry `varm['PCs']` but with
    DIFFERENT widths (rna 50 comps, prot 31). Naming the loadings' coordinate axis by the bare key (`PCs`)
    made the second modality span a length-mismatched axis -> validate error + lost loadings. The fix: a
    same-named coordinate axis is reused ONLY when the length matches (the genuine MOFA shared-factor case,
    test_joint_method_storage_shapes); a different length namespaces the axis per modality. The synthetic
    MOFA case uses EQUAL widths, so it can't guard this -- hence this distinct case."""
    try:
        import anndata as ad
        import mudata
        import pandas as pd
    except Exception:
        print("  SKIP test_per_modality_pca_different_dim (mudata unavailable)"); return
    rng = np.random.default_rng(3)
    n, ng, npr, k_rna, k_prot = 60, 40, 12, 50, 31      # DIFFERENT PCA widths per modality
    arna = ad.AnnData(rng.random((n, ng)).astype("float32"),
                      obs=pd.DataFrame(index=["c%d" % i for i in range(n)]),
                      var=pd.DataFrame(index=["g%d" % i for i in range(ng)]))
    aprot = ad.AnnData(rng.random((n, npr)).astype("float32"),
                       obs=pd.DataFrame(index=["c%d" % i for i in range(n)]),
                       var=pd.DataFrame(index=["p%d" % i for i in range(npr)]))
    arna.varm["PCs"] = rng.random((ng, k_rna)).astype("float32")
    aprot.varm["PCs"] = rng.random((npr, k_prot)).astype("float32")
    md = mudata.MuData({"rna": arna, "prot": aprot})
    ds = lstar.read_mudata(md)
    assert not lstar.validate(ds), lstar.validate(ds)    # the bug: validate flagged the length mismatch
    rl = ds.field("rna_PCs_loadings"); pl = ds.field("prot_PCs_loadings")
    # the two PCAs must NOT share one coordinate axis (different widths)
    assert rl.span[0] == "genes" and pl.span[0] == "proteins"
    assert rl.span[1] != pl.span[1], (rl.span, pl.span)
    assert len(ds.axis(rl.span[1])) == k_rna and len(ds.axis(pl.span[1])) == k_prot
    ds2 = lstar.read(_w(ds))                             # survives the store
    assert not lstar.validate(ds2)
    md2 = lstar.write_mudata(ds2)                        # and back to a MuData with both PCs intact
    assert md2.mod["rna"].varm["PCs"].shape == (ng, k_rna)
    assert md2.mod["prot"].varm["PCs"].shape == (npr, k_prot)
    print("per-modality PCA (rna 50 / prot 31, same `varm['PCs']` key): namespaced to distinct coordinate "
          "axes (no collision); round-trips both loadings -- regression for the minipbcite sweep bug")


def _w(ds):
    p = _store(); lstar.write(ds, p); return p


if __name__ == "__main__":
    test_citeseq_mudata_roundtrip()
    test_real_minipbcite()
    test_partial_overlap_mudata()
    test_joint_method_storage_shapes()
    test_per_modality_pca_different_dim()
