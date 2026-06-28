"""Public compute kernels — the shared libstar primitives, exposed so downstream tools (e.g. the
pagoda3 viewer's store prep) build on lstar's fast path instead of reimplementing reductions. Each
uses the compiled C++ accelerator when present and an identical numpy fallback otherwise (results
match; see tests/test_accel.py)."""
import numpy as np
import scipy.sparse as sp

from ._engine import resolve_engine, _accel


def col_sum_by_group(X, code, ngroups, lognorm=True, engine="auto"):
    """Per-(group, gene) sufficient stats over a CSC measure: returns (sum, sumsq, n_expr), each a
    dense (ngroups, ngenes) array, computed over log1p(X) when ``lognorm`` (else raw X).

    X is a (cells, genes) matrix (densified/copied to CSC if needed); ``code`` is a length-ncells
    int array mapping each cell to a group in [0, ngroups). This is the reduction cluster stats and
    marker tables are built from."""
    X = sp.csc_matrix(X) if not sp.issparse(X) else X.tocsc()
    ncells, ngenes = X.shape
    code = np.asarray(code)
    if resolve_engine(engine) == "c++" and hasattr(_accel, "col_sum_by_group"):
        return _accel.col_sum_by_group(X.data, X.indptr, X.indices, ncells, ngenes,
                                       code.astype("int32"), int(ngroups), bool(lognorm), 0)
    Xl = X.astype("f8").copy()
    if lognorm:
        Xl.data = np.log1p(Xl.data)
    Xlr = Xl.tocsr()
    S = np.zeros((ngroups, ngenes)); SS = np.zeros((ngroups, ngenes)); NE = np.zeros((ngroups, ngenes))
    for g in range(ngroups):
        sub = Xlr[code == g]
        S[g] = np.asarray(sub.sum(0)).ravel()
        SS[g] = np.asarray(sub.multiply(sub).sum(0)).ravel()
        NE[g] = np.asarray((sub > 0).sum(0)).ravel()
    return S, SS, NE


def markers_one_vs_rest(S, NE, nper, ncells, engine="auto"):
    """1-vs-rest marker table from per-(group,gene) sufficient stats (`viewer.markers/1-vs-rest`).

    ``S``, ``NE`` are (ngroups, ngenes) dense arrays (group-major, over log1p — the output of
    :func:`col_sum_by_group`); ``nper`` = group sizes (length ngroups); ``ncells`` = total cells.
    Returns ``(lfc, padj)``, each **(ngenes, ngroups)** — the viewer profile's gene-major orientation.
    """
    S = np.ascontiguousarray(S, dtype="f8"); NE = np.ascontiguousarray(NE, dtype="f8")
    nper = np.ascontiguousarray(nper, dtype="i8"); ngroups, ngenes = S.shape
    if resolve_engine(engine) == "c++" and hasattr(_accel, "markers_one_vs_rest"):
        return _accel.markers_one_vs_rest(S, NE, nper, int(ncells))
    grand = S.sum(0)
    lfc = np.empty((ngenes, ngroups)); padj = np.empty((ngenes, ngroups))
    for g in range(ngroups):
        d = S[g] / max(nper[g], 1) - (grand - S[g]) / max(ncells - nper[g], 1)
        lfc[:, g] = d
        padj[:, g] = np.clip(np.exp(-np.abs(d * np.sqrt(NE[g] + 1.0))), 1e-12, 1.0)
    return lfc, padj


def _lowess_predict(xs, ys, span=0.3, n_anchor=200):
    """tricube local-linear LOWESS evaluated at ``n_anchor`` anchors and linearly interpolated back to
    each ``xs`` (constant edge extrapolation). The numpy reference the core's ``lowess_predict`` matches."""
    xs = np.asarray(xs, dtype="f8"); ys = np.asarray(ys, dtype="f8"); n = xs.shape[0]
    if n < 3:
        return np.full(n, ys.mean() if n else 0.0)
    ordr = np.argsort(xs); sx = xs[ordr]; sy = ys[ordr]
    win = max(2, int(span * n))
    ax = np.empty(n_anchor); ay = np.empty(n_anchor)
    for a in range(n_anchor):
        x0 = sx[0] + (sx[n - 1] - sx[0]) * a / (n_anchor - 1)
        l = max(0, int(np.searchsorted(sx, x0, side="left")) - (win >> 1))
        r = min(n, l + win); l = max(0, r - win)
        seg_x = sx[l:r]; seg_y = sy[l:r]
        maxd = max(1e-9, float(np.abs(seg_x - x0).max()))
        w = (1.0 - (np.abs(seg_x - x0) / maxd) ** 3) ** 3
        sw = w.sum(); swx = (w * seg_x).sum(); swy = (w * seg_y).sum()
        swxx = (w * seg_x * seg_x).sum(); swxy = (w * seg_x * seg_y).sum()
        den = sw * swxx - swx * swx
        if abs(den) < 1e-12:
            ay[a] = swy / sw
        else:
            b1 = (sw * swxy - swx * swy) / den
            ay[a] = (swy - b1 * swx) / sw + b1 * x0
        ax[a] = x0
    return np.interp(xs, ax, ay)


def overdispersion(mean, var, nobs, engine="auto"):
    """Per-gene overdispersion score (pagoda2 ``adjustVariance``): residual of log(var) about a lowess
    fit of log(var)~log(mean), scored by the upper-tail variance-ratio F-test ``-log P(F>exp(res);
    nobs,nobs)``. Genes with ``nobs<3`` or mean/var ≤0 score 0. ``mean``/``var``/``nobs`` are per-gene."""
    mean = np.ascontiguousarray(mean, dtype="f8"); var = np.ascontiguousarray(var, dtype="f8")
    nobs = np.ascontiguousarray(nobs, dtype="i8")
    if resolve_engine(engine) == "c++" and hasattr(_accel, "overdispersion"):
        return _accel.overdispersion(mean, var, nobs)
    from scipy.stats import f as _f                       # numpy/scipy fallback (matches the core)
    ng = mean.shape[0]; od = np.zeros(ng)
    ok = (nobs >= 3) & (mean > 0) & (var > 0)
    if ok.sum() > 10:
        xs = np.log(mean[ok]); ys = np.log(var[ok])
        res = ys - _lowess_predict(xs, ys)
        od[ok] = -_f.logsf(np.exp(res), nobs[ok], nobs[ok])
    return od
