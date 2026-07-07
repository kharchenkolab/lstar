"""Cross-seam regression: converter output must be viewer-preppable, and `state` must reflect the
DATA, not the source slot/name. These tests cross the read_anndata/read_mudata -> extend_for_viewer
boundary that unit tests (which hand-build a `counts` field) never exercised -- the gap that let a
scanpy .h5ad (scaled .X + lognorm .raw, no counts) convert to a store that couldn't be viewer-prepped.
"""
import numpy as np
import scipy.sparse as sp
import pytest

import lstar

ad = pytest.importorskip("anndata")


def _counts(nc, ng, seed):
    rng = np.random.default_rng(seed)
    m = sp.random(nc, ng, density=0.35, format="csr", random_state=seed)
    m.data = np.rint(m.data * 10) + 1                      # non-negative integers -> raw
    return m.tocsr()


def _pbmc_like(nc=160, ng=50, seed=0):
    """scanpy-tutorial shape: scaled .X (negatives), lognorm .raw, NO counts anywhere."""
    counts = _counts(nc, ng, seed)
    ln = counts.copy().astype(float); ln.data = np.log1p(ln.data)      # lognorm
    scaled = np.asarray(ln.todense())
    scaled = (scaled - scaled.mean(0)) / (scaled.std(0) + 1e-8)        # z-scored -> negatives, dense
    a = ad.AnnData(X=ln)
    a.obs["leiden"] = np.array([str(i % 3) for i in range(nc)])
    a.obsm["X_umap"] = np.random.default_rng(seed).normal(size=(nc, 2))
    a.raw = a                                                          # .raw snapshot = lognorm
    a.X = scaled                                                       # .X = scaled
    return a


def _measure_states(ds):
    return {n: ds.field(n).state for n in ds.fields
            if ds.field(n).role == "measure" and ds.field(n).span
            and len(ds.field(n).span) == 2 and str(ds.field(n).span[0]).startswith("cells")}


# ---- state is inferred from CONTENT, not slot/name ----

def test_state_inferred_from_content():
    a = _pbmc_like()
    a.layers["counts"] = _counts(a.n_obs, a.n_vars, 7)        # a genuine raw counts layer
    ds = lstar.read_anndata(a)
    st = _measure_states(ds)
    assert st.get("X") == "scaled", st            # scaled .X was mislabeled None before
    assert st.get("raw") == "lognorm", st         # lognorm .raw was mislabeled "raw" before
    assert st.get("counts") == "raw", st          # integer layer -> raw


def test_raw_X_becomes_counts():
    a = ad.AnnData(X=_counts(120, 40, 1))                    # integer .X, nothing else
    ds = lstar.read_anndata(a)
    assert "counts" in ds.fields and ds.field("counts").state == "raw"
    # round-trips faithfully back to .X (no data moved to a counts layer)
    b = lstar.write_anndata(ds)
    xb = b.X.todense() if sp.issparse(b.X) else b.X
    assert np.allclose(np.asarray(xb), np.asarray(a.X.todense()))
    assert not b.layers


# ---- the reported bug: no-counts file must fail CLEARLY, and be preppable via basis="lognorm" ----

def test_auto_falls_back_to_lognorm():
    # scaled .X + lognorm .raw, NO counts: the default (basis='auto') now falls back to the lognorm
    # measure — with a warning — instead of erroring, so `--viewer` works on a scanpy .h5ad that kept
    # only normalized values. (var-of-lognorm stats are approximate; the warning says so.)
    ds = lstar.read_anndata(_pbmc_like())
    with pytest.warns(UserWarning, match="log-normalized"):
        lstar.extend_for_viewer(ds, order="none")
    assert "od_score" in ds.fields
    assert ds.field("od_score").provenance.get("basis") == "lognorm-input"
    assert ds.field("counts_cellmajor").state == "lognorm"


def test_scaled_only_gives_clear_error():
    # a store whose ONLY measure is scaled/z-scored (no raw, no lognorm) cannot be a viewer basis
    # (log1p on negatives is meaningless) -> a clear error that lists what is present and says why.
    rng = np.random.default_rng(0)
    a = ad.AnnData(X=rng.normal(size=(120, 40)).astype("float32"))     # negatives -> scaled, named "X"
    a.obs["leiden"] = np.array([str(i % 3) for i in range(120)])
    ds = lstar.read_anndata(a)
    with pytest.raises(ValueError) as e:
        lstar.extend_for_viewer(ds)
    msg = str(e.value)
    assert "no raw or log-normalized measure" in msg
    assert "scaled" in msg                                             # explains why it cannot be used


def test_lognorm_basis_fallback():
    ds = lstar.read_anndata(_pbmc_like())
    lstar.extend_for_viewer(ds, basis="lognorm", order="none")
    assert "od_score" in ds.fields
    assert ds.field("od_score").provenance.get("basis") == "lognorm-input"
    assert ds.field("counts_cellmajor").state == "lognorm"


# ---- the seam works when real counts exist, regardless of the measure NAME ----

def test_convert_seam_raw_counts_prepped():
    a = ad.AnnData(X=_counts(200, 60, 2))                    # counts in .X (-> named "counts")
    a.obs["leiden"] = np.array([str(i % 4) for i in range(200)])
    a.obsm["X_umap"] = np.random.default_rng(0).normal(size=(200, 2))
    ds = lstar.read_anndata(a)
    lstar.extend_for_viewer(ds)                              # default path, no kwargs
    for nav in ("counts_cellmajor", "od_score", "stats_leiden_sum", "markers_leiden_lfc"):
        assert nav in ds.fields, nav
    assert not [i for i in lstar.validate(ds) if i.startswith("ERROR")]


def test_basis_selected_by_state_not_name():
    # raw counts under a non-"counts" name -> still auto-selected (by state), no kwargs needed
    ds = lstar.Dataset(kind="sample")
    ds.add_axis("cells", [f"c{i}" for i in range(150)])
    ds.add_axis("genes", [f"g{j}" for j in range(40)])
    ds.add_field("X", _counts(150, 40, 3).tocsc(), role="measure", span=["cells", "genes"], state="raw")
    ds.add_field("leiden", np.array([str(i % 3) for i in range(150)]), role="label", span=["cells"])
    lstar.extend_for_viewer(ds, order="none")
    assert "od_score" in ds.fields and "stats_leiden_sum" in ds.fields


def test_scaled_field_named_counts_not_picked_as_raw():
    # regression (the name-shortcut hole): a measure literally named "counts" that is SCALED must not be
    # picked as raw and log1p'd -- the literal-"counts" fast path excludes scaled, like the lognorm name
    # fallback. Exercises _select_counts_basis (the shared basis picker) directly.
    from lstar.viewer import _select_counts_basis

    def _ds(fields):
        d = lstar.Dataset(kind="sample")
        d.add_axis("cells", [f"c{i}" for i in range(60)])
        d.add_axis("genes", [f"g{j}" for j in range(20)])
        for name, state in fields:
            d.add_field(name, _counts(60, 20, 5).tocsc(), role="measure", span=["cells", "genes"], state=state)
        return d
    # scaled "counts" + a real raw measure -> skip the scaled "counts", pick the raw one (log1p)
    assert _select_counts_basis(_ds([("counts", "scaled"), ("X", "raw")])) == ("X", True)
    # only a scaled "counts" -> clear error (a scaled measure is never a basis), not a silent log1p
    with pytest.raises(ValueError, match="no raw or log-normalized measure"):
        _select_counts_basis(_ds([("counts", "scaled")]))


def test_mudata_seam():
    mudata = pytest.importorskip("mudata")
    rna = ad.AnnData(X=_counts(140, 45, 4))
    rna.obs["leiden"] = np.array([str(i % 3) for i in range(140)])
    ds = lstar.read_mudata(mudata.MuData({"rna": rna}))
    # counts live under the modality name "rna" (state raw) -> auto-selected by state
    lstar.extend_for_viewer(ds, order="none")
    assert "od_score" in ds.fields
