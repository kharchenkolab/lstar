"""Conformance: the compiled lstar._accel kernels (libstar core, the same C++ bound to WASM and R)
match the numpy reference to machine precision. Skipped when the extension isn't built."""
import numpy as np
import scipy.sparse as sp
import pytest
import lstar

pytestmark = pytest.mark.skipif(not lstar.has_accel(), reason="lstar._accel not built")


def _toy(nc=200, ng=15, K=4, seed=0):
    rng = np.random.default_rng(seed)
    X = sp.csc_matrix(rng.poisson(0.6, (nc, ng)).astype("f4"))
    code = (np.arange(nc) % K).astype("i4")
    return X, code, K


def test_col_sum_by_group_matches_numpy():
    from lstar import _accel
    X, code, K = _toy()
    nc, ng = X.shape
    S, SS, NE = _accel.col_sum_by_group(X.data, X.indptr, X.indices, nc, ng, code, K, True, 0)
    Xlr = X.astype("f8").copy(); Xlr.data = np.log1p(Xlr.data); Xlr = Xlr.tocsr()
    for g in range(K):
        sub = Xlr[code == g]
        assert np.abs(S[g] - np.asarray(sub.sum(0)).ravel()).max() < 1e-9
        assert np.abs(SS[g] - np.asarray(sub.multiply(sub).sum(0)).ravel()).max() < 1e-9
        assert np.abs(NE[g] - np.asarray((sub > 0).sum(0)).ravel()).max() < 1e-9


def test_subsample_de_rank_matches_numpy():
    from lstar import _accel
    X, code, K = _toy()
    nc, ng = X.shape
    Xcsr = X.tocsr()                                  # the kernel is CSR cell-major (the de_panel)
    mem = np.full(nc, -1, "i4"); mem[code == 0] = 0; mem[code == 1] = 1
    mA, mB, lfc, nA, nB = _accel.subsample_de_rank(Xcsr.data, Xcsr.indptr, Xcsr.indices, nc, ng, mem, True)
    Xl = Xcsr.astype("f8").copy(); Xl.data = np.log1p(Xl.data)
    mA0 = np.asarray(Xl[code == 0].mean(0)).ravel(); mB0 = np.asarray(Xl[code == 1].mean(0)).ravel()
    assert (nA, nB) == (int((code == 0).sum()), int((code == 1).sum()))
    assert np.abs(mA - mA0).max() < 1e-9
    assert np.abs(mB - mB0).max() < 1e-9
    assert np.abs(lfc - (mA0 - mB0)).max() < 1e-9


def test_write_viewer_cpp_equals_python():
    # the exporter's cluster stats are identical on either engine (so stores are engine-agnostic).
    from lstar.profiles.viewer import write_viewer

    def build():
        rng = np.random.default_rng(3)
        X = sp.csc_matrix(rng.poisson(0.7, (180, 24)).astype("f4"))
        ds = lstar.Dataset(kind="sample")
        ds.add_axis("cells", ["c%d" % i for i in range(180)], role="observation")
        ds.add_axis("genes", ["g%d" % j for j in range(24)], role="feature")
        ds.add_field("counts", X, role="measure", span=["cells", "genes"], state="raw")
        ds.add_field("leiden", ["k%d" % (i % 5) for i in range(180)], role="label", span=["cells"])
        return ds

    a = write_viewer(build(), "leiden", n_od=12, engine="c++")
    b = write_viewer(build(), "leiden", n_od=12, engine="python")
    for f in ["stats_leiden_sum", "stats_leiden_sumsq", "stats_leiden_nexpr"]:
        assert np.abs(np.asarray(a.field(f).values) - np.asarray(b.field(f).values)).max() < 1e-5
