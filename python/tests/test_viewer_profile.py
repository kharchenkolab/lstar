"""viewer@0.1 profile: precomputed fields are correct and round-trip cleanly."""
import os, tempfile
import numpy as np
import scipy.sparse as sp
import lstar


def _toy(nc=120, ng=20, seed=0):
    rng = np.random.default_rng(seed)
    X = sp.csc_matrix(rng.poisson(0.6, (nc, ng)).astype("f4"))
    leiden = np.array(["c%d" % (i % 4) for i in range(nc)])
    ds = lstar.Dataset(kind="sample")
    ds.add_axis("cells", ["cell%d" % i for i in range(nc)], role="observation")
    ds.add_axis("genes", ["g%d" % j for j in range(ng)], role="feature")
    ds.add_field("counts", X, role="measure", span=["cells", "genes"], state="raw")
    ds.add_field("leiden", list(leiden), role="label", span=["cells"])
    return ds, X, leiden


def test_write_viewer_fields_and_stats():
    ds, X, leiden = _toy()
    lstar.write_viewer(ds, "leiden")
    assert "viewer@0.1" in ds.profiles
    for f in ["stats_leiden_sum", "stats_leiden_sumsq", "stats_leiden_nexpr",
              "markers_leiden_lfc", "markers_leiden_padj", "cell_order", "counts_cellmajor"]:
        assert f in ds.fields, f
    assert "od_genes" not in ds.axes                  # gene scope is on-the-fly, never precomputed
    # cluster sufficient stats == numpy per-group colSums(log1p)
    Xl = X.copy().astype("f8"); Xl.data = np.log1p(Xl.data); Xlr = Xl.tocsr()
    groups = sorted(set(leiden.tolist())); code = np.array([groups.index(l) for l in leiden])
    S = np.array([np.asarray(Xlr[code == g].sum(0)).ravel() for g in range(len(groups))])
    assert np.max(np.abs(np.asarray(ds.field("stats_leiden_sum").values) - S)) < 1e-4
    # counts_cellmajor is the full counts (cells, genes) in CSR (cell-major), raw == counts re-oriented
    dp = ds.field("counts_cellmajor")
    assert tuple(dp.values.shape) == (X.shape[0], X.shape[1]) and dp.encoding == "csr"
    assert np.max(np.abs(np.asarray(dp.values.tocsc().todense()) - np.asarray(X.todense()))) == 0


def test_write_viewer_roundtrip():
    ds, _, _ = _toy(seed=1)
    lstar.write_viewer(ds, "leiden", n_od=10)
    assert not [e for e in lstar.validate(ds) if e.startswith("ERROR")]
    p = os.path.join(tempfile.mkdtemp(), "v.lstar.zarr")
    lstar.write(ds, p)
    ds2 = lstar.read(p)
    assert "viewer@0.1" in ds2.profiles
    assert "counts_cellmajor" in ds2.fields and "stats_leiden_sum" in ds2.fields


def test_write_viewer_recomputes_no_stale():
    # the profile is always recomputed/overwritten from current inputs — a same-named field that
    # was stale (e.g. computed from a different clustering) gets corrected, not preserved.
    ds, X, leiden = _toy(seed=2)
    lstar.write_viewer(ds, "leiden")
    n_fields = len(ds.fields)
    # corrupt a derived field, then re-run: it must be restored to the correct values
    ds.field("stats_leiden_sum").values[:] = -999.0
    lstar.write_viewer(ds, "leiden")
    assert len(ds.fields) == n_fields                # idempotent (overwrites, no duplicates)
    assert ds.profiles.count("viewer@0.1") == 1
    Xl = X.copy().astype("f8"); Xl.data = np.log1p(Xl.data); Xlr = Xl.tocsr()
    groups = sorted(set(leiden.tolist())); code = np.array([groups.index(l) for l in leiden])
    S = np.array([np.asarray(Xlr[code == g].sum(0)).ravel() for g in range(len(groups))])
    assert np.max(np.abs(np.asarray(ds.field("stats_leiden_sum").values) - S)) < 1e-4
