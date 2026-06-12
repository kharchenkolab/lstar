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
