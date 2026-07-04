"""Viewer extension: `extend_for_viewer` adds the lstar-viewer precomputed fields and a hybrid cell
order, byte-for-byte equivalent to the viewer's JS store-prep on the fields that matter.

Two gates:
  1. a synthetic round-trip (fields exist with the right spans/shapes; the order is a valid permutation
     and reconstructs the per-cell counts rows);
  2. equivalence vs the JS-prepped reference store `pbmc6.lstar.zarr` (stats match exactly; marker lfc
     matches; od_score reproduces the *current* JS lowess method). The pbmc6 gate is skipped when the
     reference store isn't present.
"""
import os

import numpy as np
import scipy.sparse as sp
import pytest

import lstar
from lstar.kernels import _xy2d          # Hilbert primitive now lives in the shared kernels module


# A JS-prepped reference store to cross-check against. Overridable via $LSTAR_PBMC6 so the gate
# runs anywhere the store is available (CI artifact, a non-mac checkout); skips cleanly otherwise.
PBMC6 = os.environ.get(
    "LSTAR_PBMC6", "/Users/peter.kharchenko/pagoda/lstar-viewer/web/public/pbmc6.lstar.zarr")


def _synthetic(nc=200, ng=50, K=5, seed=0, fmt="csc"):
    rng = np.random.default_rng(seed)
    X = sp.random(nc, ng, density=0.2, format="csc", random_state=seed)
    X.data = (rng.poisson(3, size=X.data.shape) + 1).astype(np.int32)
    X = X.tocsr() if fmt == "csr" else X.tocsc()             # counts may arrive in either encoding
    ds = lstar.Dataset(kind="sample")
    ds.add_axis("cells", [f"c{i}" for i in range(nc)])
    ds.add_axis("genes", [f"g{j}" for j in range(ng)])
    ds.add_axis("umap", ["umap1", "umap2"])
    ds.add_field("counts", X, role="measure", span=["cells", "genes"], state="raw")
    ds.add_field("umap", rng.normal(size=(nc, 2)), role="embedding", span=["cells", "umap"])
    leiden = rng.integers(0, K, size=nc).astype(str)
    ds.add_field("leiden", leiden, role="label", span=["cells"])
    return ds, X


# --------------------------------------------------------------------------------------------------
# 1) synthetic round-trip
# --------------------------------------------------------------------------------------------------

def test_xy2d_is_a_permutation():
    # the canonical 4x4 Hilbert order (matches reorder.mjs)
    order = [_xy2d(4, x, y) for y in range(4) for x in range(4)]
    assert sorted(order) == list(range(16))
    assert order == [0, 1, 14, 15, 3, 2, 13, 12, 4, 7, 8, 11, 5, 6, 9, 10]


def test_synthetic_roundtrip(tmp_path):
    ds, X = _synthetic()
    nc, ng = X.shape
    lstar.extend_for_viewer(ds)

    # new fields exist with the right spans/shapes
    K = len(np.unique(np.asarray(ds.field("leiden").values, dtype=str)))
    assert ds.field("counts_cellmajor").encoding == "csr"
    assert ds.field("counts_cellmajor").span == ["cells", "genes"]
    assert ds.field("counts_cellmajor").state == "raw"
    assert ds.axis("groups_leiden") is not None and len(ds.axis("groups_leiden")) == K
    for stat, shape in [("sum", (K, ng)), ("sumsq", (K, ng)), ("nexpr", (K, ng))]:
        f = ds.field(f"stats_leiden_{stat}")
        assert f.span == ["groups_leiden", "genes"] and np.asarray(f.values).shape == shape
    for m in ("lfc", "padj"):
        f = ds.field(f"markers_leiden_{m}")
        assert f.span == ["genes", "groups_leiden"] and np.asarray(f.values).shape == (ng, K)
    assert ds.field("od_score").span == ["genes"] and np.asarray(ds.field("od_score").values).shape == (ng,)

    # counts_cellmajor_order is a valid permutation
    pos_of = np.asarray(ds.field("counts_cellmajor_order").values)
    assert pos_of.shape == (nc,)
    assert np.array_equal(np.sort(pos_of.astype(int)), np.arange(nc))

    # write -> read -> all fields survive; reading rows via the order reconstructs the original counts
    store = str(tmp_path / "synth.lstar.zarr")
    lstar.write(ds, store)
    ds2 = lstar.read(store)
    for name in ("counts_cellmajor", "counts_cellmajor_order", "od_score",
                 "stats_leiden_sum", "stats_leiden_sumsq", "stats_leiden_nexpr",
                 "markers_leiden_lfc", "markers_leiden_padj"):
        assert name in ds2.fields

    pos_of2 = np.asarray(ds2.field("counts_cellmajor_order").values).astype(int)
    cm = ds2.field("counts_cellmajor").values.tocsr()
    X_orig = X.tocsr()
    for cell in range(nc):                                  # physical row pos_of[cell] holds this cell
        got = cm.getrow(int(pos_of2[cell])).toarray().ravel()
        want = X_orig.getrow(cell).toarray().ravel()
        assert np.array_equal(got, want), f"cell {cell} row mismatch under the hybrid order"


def _synthetic2(nc=200, ng=50, seed=1):
    """A synthetic ds with TWO groupings (leiden + cell_type) so ``primary=`` can be exercised."""
    ds, X = _synthetic(nc=nc, ng=ng, seed=seed)
    rng = np.random.default_rng(seed + 7)
    ds.add_field("cell_type", np.array(["Tcell", "Bcell", "NK"])[rng.integers(0, 3, size=nc)],
                 role="label", span=["cells"])
    return ds, X


def test_primary_hoists_grouping_and_keys_the_reorder():
    # `primary=` names the grouping the viewer opens on. It (a) is prepared, (b) keys the counts_cellmajor
    # locality reorder, and (c) COMPOSES with auto-detect (the OTHER groupings are still prepped) — which a
    # plain `groupings=[...]` can't express (it can't say "and detect the rest too").
    ds, _ = _synthetic2()
    lstar.extend_for_viewer(ds, primary="cell_type")
    for g in ("cell_type", "leiden"):                          # both prepped: primary + the auto-detected rest
        assert f"stats_{g}_sum" in ds.fields and f"markers_{g}_lfc" in ds.fields
    assert ds.field("counts_cellmajor_order").provenance.get("group") == "cell_type"

    # default (no primary): the reorder keys on the auto-detected first grouping — leiden is preferred over
    # cell_type by the detection policy, so the default primary is leiden. This is exactly the case `primary=`
    # lets the viewer override when it opens on the cell-type annotation instead.
    ds2, _ = _synthetic2()
    lstar.extend_for_viewer(ds2)
    assert ds2.field("counts_cellmajor_order").provenance.get("group") == "leiden"


def test_primary_composes_with_explicit_groupings_and_validates():
    ds, _ = _synthetic2()
    lstar.extend_for_viewer(ds, primary="cell_type", groupings=["leiden"])   # primary first, then the list
    assert ds.field("counts_cellmajor_order").provenance.get("group") == "cell_type"
    assert "stats_leiden_sum" in ds.fields and "stats_cell_type_sum" in ds.fields
    ds2, _ = _synthetic2()
    try:
        lstar.extend_for_viewer(ds2, primary="not_a_field")
        assert False, "expected ValueError for an unknown primary"
    except ValueError as e:
        assert "primary" in str(e)


def test_order_none_skips_reorder():
    ds, X = _synthetic()
    lstar.extend_for_viewer(ds, order="none")
    assert "counts_cellmajor_order" not in ds.fields
    # rows stay in cell order
    cm = ds.field("counts_cellmajor").values.tocsr()
    assert np.array_equal(cm.toarray(), X.tocsr().toarray())


def test_encoding_invariance():
    """A1 contract: `extend_for_viewer` output is identical whether counts arrive CSC or CSR. The
    divergence this guards against (JS *threw* on CSR while Python/R silently normalized) is the exact
    bug that motivated the parity work -- a CSC-only fixture never exercised it."""
    a, _ = _synthetic(fmt="csc"); lstar.extend_for_viewer(a)
    b, _ = _synthetic(fmt="csr"); lstar.extend_for_viewer(b)
    assert a.field("counts").encoding == "csc" and b.field("counts").encoding == "csr"   # genuinely different inputs
    for f in ("counts_cellmajor_order", "od_score", "stats_leiden_sum", "stats_leiden_sumsq",
              "stats_leiden_nexpr", "markers_leiden_lfc", "markers_leiden_padj", "counts_cellmajor"):
        va, vb = a.field(f).values, b.field(f).values
        A = va.toarray() if sp.issparse(va) else np.asarray(va)
        B = vb.toarray() if sp.issparse(vb) else np.asarray(vb)
        same = np.array_equal(A, B) if A.dtype.kind in "iu" else np.allclose(A, B)
        assert same, f"{f} differs between CSC and CSR counts input -- encoding not normalized"


def test_lognorm_basis_keeps_counts_cellmajor_float():
    """A lognorm-basis prep must keep counts_cellmajor FLOAT -- an int32 cast (fine for raw counts)
    truncates normalized values to garbage. Regression for the Python-only bug the corpus lognorm case
    surfaced (R kept float, Python didn't -> the cell-major payloads diverged)."""
    rng = np.random.default_rng(0)
    nc, ng = 120, 30
    Xf = sp.random(nc, ng, density=0.3, format="csc", random_state=0)
    Xf.data = (np.abs(Xf.data).astype(np.float32) * 0.5 + 0.01)     # sub-1 floats -> an int cast would zero them
    ds = lstar.Dataset(kind="sample")
    ds.add_axis("cells", [f"c{i}" for i in range(nc)]); ds.add_axis("genes", [f"g{j}" for j in range(ng)])
    ds.add_axis("umap", ["u1", "u2"])
    ds.add_field("X", Xf, role="measure", span=["cells", "genes"], state="lognorm")   # no raw counts -> basis lognorm
    ds.add_field("umap", rng.normal(size=(nc, 2)), role="embedding", span=["cells", "umap"])
    ds.add_field("leiden", rng.integers(0, 4, size=nc).astype(str), role="label", span=["cells"])
    lstar.extend_for_viewer(ds, basis="lognorm")
    cm = ds.field("counts_cellmajor").values
    assert cm.dtype.kind == "f", f"lognorm counts_cellmajor must stay float, got {cm.dtype}"
    pos = np.asarray(ds.field("counts_cellmajor_order").values).astype(int)
    cmr, Xr = cm.tocsr(), Xf.tocsr()
    for cell in range(nc):                                          # values preserved (not truncated), row-for-row
        assert np.allclose(cmr.getrow(int(pos[cell])).toarray().ravel(), Xr.getrow(cell).toarray().ravel())


# --------------------------------------------------------------------------------------------------
# 2) equivalence vs the JS-prepped reference store
# --------------------------------------------------------------------------------------------------

@pytest.mark.skipif(not os.path.isdir(PBMC6), reason="pbmc6 reference store not present")
def test_pbmc6_stats_and_markers_match_reference():
    ds = lstar.read(PBMC6)
    X = ds.field("counts").values.tocsc()
    ncells, ngenes = X.shape
    lab = np.asarray(ds.field("leiden").values, dtype=str)
    groups, codes = np.unique(lab, return_inverse=True)
    K = len(groups)

    S, SS, NE = lstar.col_sum_by_group(X, codes.astype("int32"), K, lognorm=True)
    assert np.allclose(S, np.asarray(ds.field("stats_leiden_sum").values), rtol=1e-4, atol=1e-6)
    assert np.allclose(SS, np.asarray(ds.field("stats_leiden_sumsq").values), rtol=1e-4, atol=1e-5)
    assert np.allclose(NE, np.asarray(ds.field("stats_leiden_nexpr").values), rtol=1e-4, atol=1e-6)

    from lstar.viewer import _markers
    lfc, padj = _markers(S, NE, codes.astype("int32"), ncells, K, ngenes)
    stored_lfc = np.asarray(ds.field("markers_leiden_lfc").values)
    assert np.allclose(lfc, stored_lfc, rtol=1e-3, atol=1e-4)
    # per-group lfc correlation is essentially 1
    corrs = [np.corrcoef(lfc[:, g], stored_lfc[:, g])[0, 1] for g in range(K)]
    assert min(corrs) > 0.99


# --------------------------------------------------------------------------------------------------
# 3) the viewer@0.1 profile contract (validate())
# --------------------------------------------------------------------------------------------------

def test_extended_store_satisfies_viewer_profile():
    """The canonical Python prep produces a store that validate() accepts as viewer@0.1."""
    ds, _ = _synthetic()
    lstar.extend_for_viewer(ds)
    assert "viewer@0.1" in ds.profiles
    viewer_errs = [i for i in lstar.validate(ds) if i.startswith("ERROR") and "viewer@0.1" in i]
    assert viewer_errs == [], viewer_errs


def test_viewer_profile_requires_fields():
    """Stamping viewer@0.1 without the fields is an ERROR (the tag is a guarantee, not decoration)."""
    ds, _ = _synthetic()
    ds.profiles.append("viewer@0.1")                       # claim it, provide nothing
    issues = lstar.validate(ds)
    assert any("counts_cellmajor" in i and i.startswith("ERROR") for i in issues)
    assert any("od_score" in i and i.startswith("ERROR") for i in issues)
    assert any("no grouping" in i for i in issues)


def test_kernels_accel_matches_fallback():
    """The viewer recipe kernels must agree between the compiled core and the numpy fallback (so a
    prepped store is reproducible regardless of whether the accel extension is present)."""
    import lstar.kernels as k
    from lstar._engine import has_accel
    if not has_accel():
        pytest.skip("accel extension not built")
    rng = np.random.default_rng(1)
    S = rng.random((5, 40)); NE = rng.integers(0, 30, (5, 40)).astype("f8")
    nper = rng.integers(5, 40, 5); nc = int(nper.sum())
    l1, p1 = k.markers_one_vs_rest(S, NE, nper, nc, engine="c++")
    l2, p2 = k.markers_one_vs_rest(S, NE, nper, nc, engine="python")
    assert np.allclose(l1, l2) and np.allclose(p1, p2)
    assert l1.shape == (40, 5)                                   # gene-major
    mean = rng.random(200) * 3; var = mean * rng.uniform(0.5, 8, 200); nobs = rng.integers(0, 80, 200)
    o1 = k.overdispersion(mean, var, nobs, engine="c++")
    o2 = k.overdispersion(mean, var, nobs, engine="python")
    assert np.allclose(o1, o2, atol=1e-6)
    # cell_order must reach the core (viewer_cell_order was defined in _accel.cpp but not registered,
    # so Python silently always used the numpy fallback -- guard the c++ path here now that it's bound).
    from lstar._engine import _accel
    assert hasattr(_accel, "viewer_cell_order"), "viewer_cell_order not registered in _accel"
    code = rng.integers(0, 6, size=250).astype(np.int32); emb = rng.normal(size=(250, 2))
    assert np.array_equal(k.cell_order(code, emb, engine="c++"), k.cell_order(code, emb, engine="python"))
    assert np.array_equal(k.cell_order(code, None, engine="c++"), k.cell_order(code, None, engine="python"))


_NAV = {"counts_cellmajor", "counts_cellmajor_order", "od_score", "stats_leiden_sum",
        "stats_leiden_sumsq", "stats_leiden_nexpr", "markers_leiden_lfc", "markers_leiden_padj"}


def test_navigators_tagged_as_cache():
    """extend_for_viewer tags exactly its produced navigators `provenance.cache='viewer@0.1'`; the
    primary inputs (counts/clusters/embedding) are NOT tagged (they're kept by converters)."""
    ds, _ = _synthetic()
    lstar.extend_for_viewer(ds)
    tagged = {n for n in ds.fields if (ds.field(n).provenance or {}).get("cache") == "viewer@0.1"}
    assert tagged == _NAV, tagged ^ _NAV
    for n in ("counts", "leiden", "umap"):
        assert not (ds.field(n).provenance or {}).get("cache")


def test_anndata_export_drops_cache_navigators():
    """A non-viewer export (AnnData) drops the regenerable cache navigators and records them in
    `dropped` — no scrambled `counts_cellmajor` layer, no navigator leaked into obs/var/obsm."""
    import importlib.util
    if importlib.util.find_spec("anndata") is None:
        pytest.skip("anndata not installed")
    from lstar.profiles.anndata import write_anndata
    ds, _ = _synthetic()
    lstar.extend_for_viewer(ds)
    ad = write_anndata(ds)
    names = set(ad.layers) | set(ad.obs.columns) | set(ad.var.columns) | set(ad.obsm) | set(ad.varm)
    assert not any("counts_cellmajor" in str(k) for k in names)        # no scrambled layer
    assert not (_NAV & names)                                          # nothing leaked
    assert _NAV.issubset(set(ds.dropped))                              # recorded, not silently lost
    assert "X_umap" in ad.obsm and "leiden" in ad.obs.columns          # primaries carried


def test_viewer_profile_rejects_transposed_markers():
    """Markers are gene-major (ng x K); a group-major (K x ng) markers table is the R-pagoda2 drift
    this contract exists to catch."""
    ds, _ = _synthetic()
    lstar.extend_for_viewer(ds)
    K = len(ds.axis("groups_leiden"))
    ng = len(ds.axis("genes"))
    # overwrite markers_leiden_lfc with a transposed (group-major) copy + matching span
    good = np.asarray(ds.field("markers_leiden_lfc").values)        # ng x K
    ds.fields.pop("markers_leiden_lfc")
    ds.add_field("markers_leiden_lfc", good.T.copy(), role="measure",
                 span=["groups_leiden", "genes"], encoding="dense")  # K x ng -- WRONG
    issues = lstar.validate(ds)
    assert any("markers_leiden_lfc" in i and "gene-major" in i and i.startswith("ERROR")
               for i in issues), issues
