"""Public compute kernels — the shared libstar primitives, exposed so downstream tools (e.g. the
pagoda3 viewer's store prep) build on lstar's fast path instead of reimplementing reductions. Each
uses the compiled C++ accelerator when present and an identical numpy fallback otherwise (results
match; see tests/test_accel.py)."""
import numpy as np

from ._engine import resolve_engine, _accel
from ._sparse import as_csc


def col_sum_by_group(X, code, ngroups, lognorm=True, engine="auto"):
    """Per-(group, gene) sufficient stats over a CSC measure: returns (sum, sumsq, n_expr), each a
    dense (ngroups, ngenes) array, computed over log1p(X) when ``lognorm`` (else raw X).

    X is a (cells, genes) matrix (densified/copied to CSC if needed); ``code`` is a length-ncells
    int array mapping each cell to a group in [0, ngroups). This is the reduction cluster stats and
    marker tables are built from."""
    X = as_csc(X)
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


# --- viewer@0.1 canonical cell order (the physical row reorder key) --------------------------------
_N_GRID = 1024                                            # Hilbert grid resolution (power of two)


def _xy2d(N, x, y):
    """Canonical Hilbert index of grid cell (x, y) on an N x N curve (reflection uses N-1) -- the
    numpy-fallback twin of ``lstar::hilbert_xy2d``."""
    d = 0; s = N >> 1
    while s > 0:
        rx = 1 if (x & s) else 0; ry = 1 if (y & s) else 0
        d += s * s * ((3 * rx) ^ ry)
        if ry == 0:
            if rx == 1:
                x = N - 1 - x; y = N - 1 - y
            x, y = y, x
        s >>= 1
    return d


def _hilbert_index(emb, grid=_N_GRID):
    """Per-cell Hilbert index over a ``grid`` x ``grid`` grid of the min-max-scaled first 2 embedding
    dims -- the numpy-fallback twin of ``lstar::hilbert_index``."""
    xy = np.asarray(emb, dtype=np.float64)[:, :2]
    x, y = xy[:, 0], xy[:, 1]; n = xy.shape[0]
    xr = (x.max() - x.min()) or 1.0; yr = (y.max() - y.min()) or 1.0
    gx = np.minimum(grid - 1, np.floor((x - x.min()) / xr * (grid - 1))).astype(np.int64)
    gy = np.minimum(grid - 1, np.floor((y - y.min()) / yr * (grid - 1))).astype(np.int64)
    out = np.empty(n, dtype=np.int64)
    for i in range(n):
        out[i] = _xy2d(grid, int(gx[i]), int(gy[i]))
    return out


def cell_order(primary_code, embedding=None, grid=_N_GRID, engine="auto"):
    """viewer@0.1 canonical physical cell order shared by Python/R/JS (C++ core ``viewer_cell_order``).

    Returns ``pos_of`` (int64, length ncells): each cell's physical row after a stable sort by
    (cluster ``primary_code``, then Hilbert index of ``embedding`` when given, else cell index). The
    numpy fallback matches the core exactly, so a store prepped without the accelerator is byte-identical
    to one prepped with it (and to the R/JS surfaces)."""
    primary_code = np.ascontiguousarray(primary_code, dtype=np.int32)
    n = primary_code.size
    emb = None if embedding is None else np.ascontiguousarray(np.asarray(embedding, dtype=np.float64))
    if resolve_engine(engine) == "c++" and hasattr(_accel, "viewer_cell_order"):
        return np.asarray(_accel.viewer_cell_order(primary_code, emb, int(grid)), dtype=np.int64)
    pc64 = primary_code.astype(np.int64)
    if emb is None:                                      # sort by (cluster, cell index) -- lexsort last key = primary
        perm = np.lexsort((np.arange(n), pc64))
    else:
        perm = np.lexsort((np.arange(n), _hilbert_index(emb, grid), pc64))
    pos = np.empty(n, dtype=np.int64); pos[perm] = np.arange(n)
    return pos
