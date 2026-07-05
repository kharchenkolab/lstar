"""Shared sparse-matrix normalization.

`as_csc` / `as_csr` turn a measure — a dense ndarray or any scipy sparse matrix — into CSC / CSR. The
"convert if already sparse, else densify into sparse" idiom lived inline in several modules (viewer, de,
kernels); centralizing it means a measure read can't quietly become dense-unsafe in one place while its
siblings handle both. This is the Python twin of the JS reader's `fieldAsCsc` and R's
`as(., "CsparseMatrix")` coercion — keeping the surfaces' measure normalization aligned is the point of
the cross-language parity contract (docs/parity.md); the JS side once diverged here and threw on a dense
primary measure.
"""
import scipy.sparse as sp


def as_csc(x):
    """`x` (a dense array or any scipy sparse matrix) as a CSC matrix."""
    return x.tocsc() if sp.issparse(x) else sp.csc_matrix(x)


def as_csr(x):
    """`x` (a dense array or any scipy sparse matrix) as a CSR matrix."""
    return x.tocsr() if sp.issparse(x) else sp.csr_matrix(x)
