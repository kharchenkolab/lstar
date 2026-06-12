"""Tier-1 promotions out of the lossless passthrough into typed, axis-bound fields -- grounded in
**real** datasets (pbmc68k_reduced's real `*_colors` + `uns['pca']`, and a real scVelo pancreas output
for velocity), not fabricated structures:
  - color palettes  uns['<key>_colors']     -> a label field over the factor axis (category order)
  - PCA variance     uns['pca']['variance*'] -> measures over the pca coordinate axis
  - RNA-velocity     spliced/unspliced layers -> measures; uns['velocity_graph'] -> a cell-cell relation
All round-trip through the store (ordinary fields) and regenerate their native uns/layers on write-back.

Run: PYTHONPATH=python/src python3 python/tests/test_tier1_promote.py
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
    return os.path.join(tempfile.mkdtemp(), "t1.lstar.zarr")


def test_color_palette_promotion():
    a = corpus.pbmc68k_reduced()                            # REAL bulk_labels_colors + louvain_colors
    if a is None:
        print("  SKIP test_color_palette_promotion (corpus unavailable)"); return
    ds = lstar.read_anndata(a)
    cf = ds.field("bulk_labels_colors")
    assert cf.subtype == "color" and cf.span == ["bulk_labels"]            # bound to the factor axis
    assert len(np.asarray(cf.values)) == len(ds.axis("bulk_labels"))       # one color per category
    assert "bulk_labels_colors" not in ds.aux.get("anndata.uns", {})       # promoted out of the tail
    assert not lstar.validate(ds)
    a2 = lstar.write_anndata(lstar.read(_w(ds)))
    assert list(a2.uns["bulk_labels_colors"]) == list(a.uns["bulk_labels_colors"])  # regenerated
    print("colors (real pbmc68k): uns['bulk_labels_colors'] -> field over the factor axis; round-trips")


def test_pca_variance_promotion():
    a = corpus.pbmc68k_reduced()                            # REAL uns['pca']['variance_ratio']
    if a is None:
        print("  SKIP test_pca_variance_promotion (corpus unavailable)"); return
    ds = lstar.read_anndata(a)
    vr = ds.field("pca_variance_ratio")
    assert vr.subtype == "pca_var" and vr.span == ["pca"]                  # measure over the pca axis
    assert np.allclose(np.asarray(vr.values), a.uns["pca"]["variance_ratio"], rtol=1e-5)
    assert "variance_ratio" not in ds.aux.get("anndata.uns", {}).get("pca", {})
    a2 = lstar.write_anndata(lstar.read(_w(ds)))
    assert np.allclose(a2.uns["pca"]["variance_ratio"], a.uns["pca"]["variance_ratio"], rtol=1e-5)
    print("pca (real pbmc68k): uns['pca']['variance_ratio'] -> measure over the pca axis; round-trips")


def test_velocity_real_scvelo():
    a = corpus.pancreas_velocity()                          # REAL scVelo output (committed fixture)
    if a is None:
        print("  SKIP test_velocity_real_scvelo (fixture unavailable)"); return
    ds = lstar.read_anndata(a)
    # spliced/unspliced layers come free as measures
    assert ds.field("spliced").role == "measure" and ds.field("spliced").span == ["cells", "genes"]
    assert ds.field("unspliced").role == "measure"
    # the velocity graph lives in uns (NOT obsp) -> typed as a cell-cell relation, not dropped
    vg = ds.field("velocity_graph")
    assert vg.role == "relation" and vg.span == ["cells", "cells"] and vg.subtype == "uns_graph"
    assert "velocity_graph" in ds.fields and "velocity_graph_neg" in ds.fields
    assert "velocity_graph" not in ds.aux.get("anndata.uns", {})           # not left in passthrough
    # real clusters categorical + colors also promoted
    assert ds.axis("clusters").role == "factor" and "clusters_colors" in ds.fields
    assert not lstar.validate(ds)

    import scipy.sparse as sp
    nnz0 = a.uns["velocity_graph"].nnz
    a2 = lstar.write_anndata(lstar.read(_w(ds)))            # round-trip through the store
    assert "spliced" in a2.layers and "unspliced" in a2.layers
    assert sp.issparse(a2.uns["velocity_graph"]) and a2.uns["velocity_graph"].shape == a.uns["velocity_graph"].shape
    assert a2.uns["velocity_graph"].nnz == nnz0            # the velocity graph survives byte-faithfully
    print("velocity (real scVelo): spliced/unspliced measures + uns velocity_graph -> cell-cell relation; round-trips")


def _w(ds):
    p = _store(); lstar.write(ds, p); return p


if __name__ == "__main__":
    test_color_palette_promotion()
    test_pca_variance_promotion()
    test_velocity_real_scvelo()
