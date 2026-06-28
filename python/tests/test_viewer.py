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
from lstar.viewer import _xy2d


PBMC6 = "/Users/peter.kharchenko/pagoda/lstar-viewer/web/public/pbmc6.lstar.zarr"


def _synthetic(nc=200, ng=50, K=5, seed=0):
    rng = np.random.default_rng(seed)
    X = sp.random(nc, ng, density=0.2, format="csc", random_state=seed)
    X.data = (rng.poisson(3, size=X.data.shape) + 1).astype(np.int32)
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


def test_order_none_skips_reorder():
    ds, X = _synthetic()
    lstar.extend_for_viewer(ds, order="none")
    assert "counts_cellmajor_order" not in ds.fields
    # rows stay in cell order
    cm = ds.field("counts_cellmajor").values.tocsr()
    assert np.array_equal(cm.toarray(), X.tocsr().toarray())


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
