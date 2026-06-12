"""Lazy / streaming field values.

`read(path, lazy=True)` returns field values as these proxies instead of materialized
numpy/scipy objects. A proxy holds the open zarr arrays and reads nothing until asked, so a
large store opens in milliseconds. The sparse proxy additionally supports *streaming* by
column block, so per-gene reductions (variance modeling, HVG selection) run over a measure
that never fully lands in memory -- the pagoda2 lazy-normalized-view pattern, in Python.

All proxies satisfy `np.asarray(proxy)` (materialization) and expose `.shape`/`.dtype`, so
they drop into existing code paths; `.materialize()` is the explicit form.
"""
import numpy as np
import scipy.sparse as sp


class LazyDense:
    """A dense field value backed by an unread zarr array."""

    def __init__(self, z):
        self._z = z
        self.shape = tuple(z.shape)
        self.dtype = z.dtype

    def __getitem__(self, idx):
        return self._z[idx]                 # zarr streams only the requested selection

    def materialize(self):
        return self._z[...]

    def __array__(self, dtype=None):
        a = self._z[...]
        return a.astype(dtype) if dtype is not None else a

    def __len__(self):
        return self.shape[0]

    def sum(self):
        # stream row blocks so a big dense array is never fully resident
        total = 0.0
        n = self.shape[0]
        step = max(1, 1_000_000 // max(1, int(np.prod(self.shape[1:])) or 1))
        for i in range(0, n, step):
            total += float(np.asarray(self._z[i:i + step]).sum())
        return total

    def __repr__(self):
        return "LazyDense(shape=%s, dtype=%s)" % (self.shape, self.dtype)


class LazyCSX:
    """A CSR/CSC field value backed by unread zarr `data`/`indices`/`indptr` arrays.

    `indptr` (length n+1) is small and read eagerly so any slice's extent is known without
    touching `data`/`indices`. For CSC, "outer" = columns (genes for a cells x genes measure);
    for CSR, "outer" = rows.
    """

    def __init__(self, fmt, data_z, indices_z, indptr_z, shape):
        self.fmt = fmt                      # "csc" | "csr"
        self._data = data_z
        self._indices = indices_z
        self.indptr = np.asarray(indptr_z[...])
        self.shape = tuple(int(s) for s in shape)
        self.dtype = data_z.dtype
        self.idtype = indices_z.dtype       # so it can also be a streaming-WRITE source (zarr_io)
        self.nnz = int(self.indptr[-1])
        self.n_outer = self.shape[1] if fmt == "csc" else self.shape[0]

    # ---- streaming access over the compressed (outer) axis ----
    def outer_block(self, start, stop):
        """The slice [start:stop] along the compressed axis, as a small scipy matrix."""
        a, b = int(self.indptr[start]), int(self.indptr[stop])
        data = np.asarray(self._data[a:b])
        indices = np.asarray(self._indices[a:b])
        indptr = self.indptr[start:stop + 1] - self.indptr[start]
        n = stop - start
        inner = self.shape[0] if self.fmt == "csc" else self.shape[1]
        cls = sp.csc_matrix if self.fmt == "csc" else sp.csr_matrix
        shape = (inner, n) if self.fmt == "csc" else (n, inner)
        return cls((data, indices, indptr), shape=shape)

    def blocks(self, block=2048):
        """Iterate (start, stop, submatrix) over the compressed axis in blocks."""
        for start in range(0, self.n_outer, block):
            stop = min(start + block, self.n_outer)
            yield start, stop, self.outer_block(start, stop)

    def materialize(self):
        cls = sp.csc_matrix if self.fmt == "csc" else sp.csr_matrix
        return cls((np.asarray(self._data[...]), np.asarray(self._indices[...]),
                    self.indptr), shape=self.shape)

    def __array__(self, dtype=None):
        m = self.materialize().toarray()
        return m.astype(dtype) if dtype is not None else m

    def sum(self):
        # stream the data array in blocks; never hold all nonzeros at once
        total = 0.0
        n = self.nnz
        step = 4_000_000
        for a in range(0, n, step):
            total += float(np.asarray(self._data[a:min(a + step, n)]).sum())
        return total

    def __repr__(self):
        return "LazyCSX(%s, shape=%s, nnz=%d)" % (self.fmt, self.shape, self.nnz)


def _block_col_stats(sub, nrows, lognorm):
    """Per-column mean/var/nnz of one CSC block, zero-aware (the same form as libstar's
    `csc_col_mean_var`).

    Memory-lean: the nonzero values stay in their stored dtype (float32 measures stay float32 --
    no widening copy), while the moments accumulate in float64 because `np.bincount` always sums
    weights in double. That is the low-precision-storage / high-precision-accumulation pattern, so
    a float32 measure costs no extra memory yet the per-gene mean/var are float64-accurate (the
    `Σx^2 − (Σx)^2/n` identity is well-conditioned here: relative error ~ eps·(1 + mean^2/var)).
    """
    sub = sub.tocsc()
    nb = sub.shape[1]
    counts = np.diff(sub.indptr).astype(np.int64)        # nonzeros per column
    data = sub.data                                      # native dtype, no upcast
    if lognorm:
        data = np.log1p(data)
    col = np.repeat(np.arange(nb, dtype=np.int32), counts)   # column index per nonzero
    s = np.bincount(col, weights=data, minlength=nb)         # float64 accumulation (free)
    sq = np.bincount(col, weights=data * data, minlength=nb)
    mu = s / nrows
    ss = sq - s * s / nrows                              # zero-aware: implicit zeros contribute 0
    var = ss / (nrows - 1) if nrows > 1 else np.zeros(nb)
    return mu, var, counts


def stream_col_stats(field_value, lognorm=False, block=2048, n_threads=1, engine="auto"):
    """Zero-aware per-column mean/variance of a CSC measure, streamed by column block.

    Mirrors libstar's `csc_col_mean_var`: implicit zeros are accounted for, and `lognorm` applies
    log1p per nonzero on the fly (the normalized matrix is never built). Accepts a `LazyCSX`
    (streams from disk, constant memory) or a materialized scipy CSC.

    engine selects the compute path: 'auto' (default) uses the compiled C++/OpenMP accelerator when
    present and the pure-Python kernel otherwise; 'c++'/'python' force one (see `lstar.has_accel`,
    `lstar.show_config`). Results are identical across engines and thread counts.

    n_threads is the threading policy from this call: 1 = serial (default), N = N threads, 0/None =
    all cores. With the C++ engine the OpenMP kernel parallelizes columns within each block; with
    the Python engine a thread pool runs blocks concurrently (the numpy kernels release the GIL).
    """
    import os
    from ._engine import resolve_engine, _accel
    eng = resolve_engine(engine)
    nt = (os.cpu_count() or 1) if n_threads in (0, None) else n_threads

    if isinstance(field_value, LazyCSX):
        assert field_value.fmt == "csc", "stream_col_stats needs a CSC field"
        ncols, nrows = field_value.n_outer, field_value.shape[0]
        def get_block(a, b):
            return field_value.outer_block(a, b)
    else:
        mat = field_value.tocsc()
        nrows, ncols = mat.shape
        def get_block(a, b):
            return mat[:, a:b]

    mean = np.empty(ncols)
    var = np.empty(ncols)
    nnz = np.empty(ncols, dtype=np.int64)
    ranges = [(a, min(a + block, ncols)) for a in range(0, ncols, block)]

    if eng == "c++":
        # Stream blocks (bounded memory); the OpenMP kernel reduces each block's columns in
        # parallel. Pass raw values + the lognorm flag -- log1p happens inside the kernel.
        for a, b in ranges:
            sub = get_block(a, b).tocsc()
            m, v, c = _accel.col_mean_var(sub.data, sub.indptr.astype(np.int64, copy=False),
                                          int(nrows), int(nt), bool(lognorm))
            mean[a:b] = m
            var[a:b] = v
            nnz[a:b] = c
        return mean, var, nnz

    def work(rng):
        a, b = rng
        mu, v, counts = _block_col_stats(get_block(a, b), nrows, lognorm)
        return a, b, mu, v, counts

    if nt > 1 and len(ranges) > 1:
        from concurrent.futures import ThreadPoolExecutor
        with ThreadPoolExecutor(max_workers=nt) as ex:
            results = list(ex.map(work, ranges))
    else:
        results = (work(r) for r in ranges)

    for a, b, mu, v, counts in results:
        mean[a:b] = mu
        var[a:b] = v
        nnz[a:b] = counts
    return mean, var, nnz
