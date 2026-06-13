"""Spatial -- **conceptual** support only (per project decision): spatial coordinates become a *named
observed coordinate axis* (`spatial`), the lstar concept the multimodal design recommended. Serious
spatial support (images, vendor coordinate frames, molecule tables) is deferred to its own tier, so the
`uns['spatial']` image blob is left in the lossless passthrough rather than typed.

This uses a small synthetic Visium-like object (spots x 2-D coords); no real spatial data.

Run: PYTHONPATH=python/src python3 python/tests/test_spatial.py
"""
import os
import sys
import tempfile
import warnings

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
import lstar  # noqa: E402

warnings.filterwarnings("ignore")


def _synth_visium(n=64):
    import anndata as ad
    import pandas as pd
    import scipy.sparse as sp
    rng = np.random.default_rng(0)
    a = ad.AnnData(sp.csr_matrix(rng.poisson(1.0, (n, 30)).astype("float32")),
                   obs=pd.DataFrame(index=[f"spot{i}" for i in range(n)]),
                   var=pd.DataFrame(index=[f"g{i}" for i in range(30)]))
    side = int(np.sqrt(n))
    xy = np.array([[i % side, i // side] for i in range(n)], dtype="float32")   # a tissue grid of spots
    a.obsm["spatial"] = xy
    a.uns["spatial"] = {"libA": {"images": {"hires": np.zeros((8, 8, 3), "float32")},
                                 "scalefactors": {"spot_diameter_fullres": 12.3}}}   # vendor blob (deferred)
    return a, xy


def test_spatial_named_axis_and_roundtrip():
    a, xy = _synth_visium()
    ds = lstar.read_anndata(a)
    # spatial coordinates -> a named, *observed* coordinate axis (not a derived embedding like pca/umap)
    assert "spatial" in ds.axes
    ax = ds.axis("spatial")
    assert ax.role == "coordinate" and ax.origin == "observed"
    sf = ds.field("spatial")
    assert sf.role == "embedding" and sf.span == ["cells", "spatial"] and sf.subtype == "spatial"
    # the vendor image blob is deferred -> kept verbatim in the passthrough, not silently dropped
    assert "spatial" in ds.aux.get("anndata.uns", {})
    assert not lstar.validate(ds)

    p = os.path.join(tempfile.mkdtemp(), "sp.lstar.zarr"); lstar.write(ds, p)
    a2 = lstar.write_anndata(lstar.read(p))
    assert "spatial" in a2.obsm and "X_spatial" not in a2.obsm   # back to obsm['spatial'], not X_spatial
    assert np.allclose(np.asarray(a2.obsm["spatial"]), xy)
    assert "spatial" in a2.uns                                   # image blob regenerated from passthrough
    print("spatial (conceptual): obsm['spatial'] -> a named observed `spatial` coordinate axis; "
          "round-trips to obsm['spatial']; uns image blob preserved in passthrough (images deferred)")


if __name__ == "__main__":
    test_spatial_named_axis_and_roundtrip()
