"""Old serialized-format recognition, distilled to a synthetic test (no real old data, no old anndata
install needed -- we synthesize the *on-disk layout* with h5py). lstar's backed/streaming h5ad reader
must recognize the **legacy** sparse attributes (`h5sparse_format`/`h5sparse_shape`, anndata < 0.7) as
well as the modern ones (`encoding-type`/`shape`, anndata >= 0.7) -- the graceful-version-recognition
rule applied to the on-disk format, so an old h5ad streams into an L* store the same as a new one.

Run: PYTHONPATH=python/src python3 python/tests/test_legacy_format.py
"""
import os
import sys
import tempfile

import numpy as np
import scipy.sparse as sp

sys.path.insert(0, os.path.dirname(__file__))
from lstar.profiles.anndata import _BackedH5Sparse  # noqa: E402


def _write_sparse_group(path, M, legacy):
    import h5py
    with h5py.File(path, "w") as f:
        g = f.create_group("X")
        g.create_dataset("data", data=M.data)
        g.create_dataset("indices", data=M.indices)
        g.create_dataset("indptr", data=M.indptr)
        if legacy:                                  # anndata < 0.7: h5sparse_* attributes
            g.attrs["h5sparse_format"] = "csr" if sp.isspmatrix_csr(M) else "csc"
            g.attrs["h5sparse_shape"] = np.array(M.shape)
        else:                                       # anndata >= 0.7: encoding-type / shape
            g.attrs["encoding-type"] = "csr_matrix" if sp.isspmatrix_csr(M) else "csc_matrix"
            g.attrs["encoding-version"] = "0.1.0"
            g.attrs["shape"] = np.array(M.shape)


def _check(legacy):
    M = sp.random(30, 12, density=0.3, format="csr", random_state=7).astype("float32")
    p = os.path.join(tempfile.mkdtemp(), "x.h5")
    _write_sparse_group(p, M, legacy=legacy)
    b = _BackedH5Sparse(p, "X")
    assert b.fmt == "csr" and b.shape == (30, 12), (b.fmt, b.shape)
    # blocked read (the streaming path lstar.write(stream=True) drives) matches the full matrix
    full = sp.vstack([blk for _a, _b, blk in b.blocks(8)]).tocsr()
    assert full.shape == M.shape and np.allclose(full.toarray(), M.toarray())
    b.close()


def test_legacy_and_modern_h5sparse():
    _check(legacy=True)                             # anndata < 0.7 on-disk layout
    _check(legacy=False)                            # anndata >= 0.7 on-disk layout
    print("legacy-format recognition: both `h5sparse_format` (<0.7) and `encoding-type` (>=0.7) on-disk "
          "sparse layouts are recognized + stream block-by-block identically")


if __name__ == "__main__":
    test_legacy_and_modern_h5sparse()
