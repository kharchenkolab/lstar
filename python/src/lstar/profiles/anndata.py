"""AnnData profile (scverse): read_anndata / write_anndata.

Maps an AnnData object to/from an L* Dataset. The native location of each field is
recorded in field.provenance["anndata"] so write-back is exact; field *names* follow the
shared core vocabulary (X_pca -> 'pca', etc.) so cross-format conversion is meaningful.

anndata is imported lazily, so `import lstar` never requires it.
"""
import numpy as np

from ..model import Categorical, Dataset, OBSERVED, DERIVED, _is_categorical

PROFILE = "anndata@0.1"


def _anndata_version():
    """Detected anndata version, recorded so a reader knows which schema produced the store.

    The in-memory object normalizes most cross-version differences on read (e.g. anndata >=0.7
    migrates uns['neighbors']['distances'/'connectivities'] into .obsp), but .raw, the uns
    layout, and dtype conventions still vary, so we recognize rather than assume.
    """
    try:
        import anndata
        return "anndata@%s" % anndata.__version__
    except Exception:
        return "anndata@?"

# known varm -> coordinate-axis pairings (shared vocabulary)
_VARM_PAIR = {"PCs": "pca"}


def _strip_x(key):
    return key[2:] if key.startswith("X_") else key


def _guess_state(layer_name):
    n = layer_name.lower()
    if "count" in n or n == "raw":
        return "raw"
    if n in ("data", "lognorm", "log1p", "lognorm", "normalized", "norm"):
        return "lognorm"
    if "scale" in n:
        return "scaled"
    return None


def _guess_subtype(key):
    n = key.lower()
    if "dist" in n:
        return "distance"
    if "conn" in n or "knn" in n or "neighbor" in n or "snn" in n:
        return "similarity"
    return None


def _by_dtype_series(s):
    import pandas as pd
    if pd.api.types.is_bool_dtype(s):
        return np.asarray(s.values), "label"
    if isinstance(s.dtype, pd.CategoricalDtype):       # preserve category order + NaN (-> -1 missing)
        c = s.values
        return Categorical(np.asarray(c.codes), np.asarray(c.categories, dtype=str),
                           ordered=bool(c.ordered)), "label"
    if pd.api.types.is_numeric_dtype(s):
        return np.asarray(s.values), "measure"
    return np.asarray(s.astype(str).values, dtype=str), "label"


def _to_pandas_col(f):
    """An L* obs/var field as a column value for write-back: a `Categorical` field is rebuilt as a
    `pd.Categorical` (categories, order, and `-1`->NaN missing preserved), else a plain array."""
    if _is_categorical(f.values):
        import pandas as pd
        c = f.values
        return pd.Categorical.from_codes(np.asarray(c.codes), categories=list(c.categories),
                                         ordered=bool(c.ordered))
    return np.asarray(f.values)


def _coord_axis(ds, name, ncol, observed=False):
    if name not in ds.axes:
        ds.add_axis(name, ["%s%d" % (name, i) for i in range(ncol)],
                    origin=(OBSERVED if observed else DERIVED), role="coordinate")
    return name


def _pair_coord_axis(ds, varm_key, ncol):
    target = _VARM_PAIR.get(varm_key, _strip_x(varm_key))
    if target in ds.axes and len(ds.axes[target]) == ncol:
        return target
    if target not in ds.axes:
        ds.add_axis(target, ["%s%d" % (target, i) for i in range(ncol)],
                    origin=DERIVED, role="coordinate")
    return target


def _sparse_attrs(g):
    """(fmt, shape) for an h5ad sparse group across format versions, or None if it isn't sparse.

    Recognizes the modern `encoding-type`/`shape` attributes (anndata >= 0.7) AND the legacy
    `h5sparse_format`/`h5sparse_shape` attributes (older h5ad) -- the same graceful version handling
    the in-memory profile uses, applied to the on-disk layout."""
    a = getattr(g, "attrs", {})
    et = str(a.get("encoding-type", ""))
    if "csr" in et or "csc" in et:
        return ("csc" if "csc" in et else "csr"), tuple(int(x) for x in a["shape"])
    if "h5sparse_format" in a:
        hf = str(a["h5sparse_format"])
        return ("csc" if "csc" in hf else "csr"), tuple(int(x) for x in a["h5sparse_shape"])
    return None


class _BackedH5Sparse:
    """A streaming sparse source over an on-disk h5ad sparse group (`data`/`indices`/`indptr`).
    Yields scipy blocks straight from disk so `lstar.write(..., stream=True)` can copy a large
    matrix into an L* store without ever materializing it. Holds an open h5py handle; call `.close()`
    when done (or let it be garbage-collected)."""

    def __init__(self, filename, key):
        import h5py
        self._f = h5py.File(filename, "r")
        g = self._f[key]
        self.fmt, self.shape = _sparse_attrs(g)               # handles modern + legacy attrs
        self._data, self._indices = g["data"], g["indices"]
        self.indptr = g["indptr"][:]                          # the compressed-axis pointer; small
        self.nnz = int(self._data.shape[0])
        self.dtype = self._data.dtype
        self.idtype = self._indices.dtype
        self.n_outer = self.shape[0] if self.fmt == "csr" else self.shape[1]
        self._inner = self.shape[1] if self.fmt == "csr" else self.shape[0]

    def outer_block(self, a, b):
        """Outer slice [a:b) as a small scipy CSR/CSC matrix (reads only that block from disk)."""
        import scipy.sparse as sp
        lo, hi = int(self.indptr[a]), int(self.indptr[b])
        iptr = (self.indptr[a:b + 1] - self.indptr[a]).astype(np.int64)
        cls = sp.csr_matrix if self.fmt == "csr" else sp.csc_matrix
        shape = (b - a, self._inner) if self.fmt == "csr" else (self._inner, b - a)
        return cls((self._data[lo:hi], self._indices[lo:hi], iptr), shape=shape)

    def blocks(self, block=4096):
        for a in range(0, self.n_outer, block):
            b = min(a + block, self.n_outer)
            yield a, b, self.outer_block(a, b)

    def close(self):
        try:
            self._f.close()
        except Exception:
            pass


def _backed_sparse(filename, *keys):
    """Return a `_BackedH5Sparse` for the first of `keys` that is a sparse group, else None.

    `keys` lets callers try location variants across h5ad versions (e.g. raw at `raw/X` vs `raw.X`)."""
    if not filename:
        return None
    try:
        import h5py
        with h5py.File(filename, "r") as f:
            hit = next((k for k in keys if k in f and _sparse_attrs(f[k]) is not None), None)
        return _BackedH5Sparse(filename, hit) if hit else None
    except Exception:
        return None


def read_anndata(adata, kind="sample"):
    """Read a live AnnData object into an L* Dataset.

    If `adata` was opened in backed mode (`anndata.read_h5ad(path, backed="r")`), the large sparse
    matrices (`X`, `layers`, `.raw`) are held as on-disk streaming sources rather than read into
    memory, so a subsequent `lstar.write(..., stream=True)` performs a bounded-memory conversion.
    """
    ds = Dataset(kind=kind)
    ds.profiles = [PROFILE, _anndata_version()]
    cells = np.asarray(adata.obs_names.to_numpy(), dtype=str)
    genes = np.asarray(adata.var_names.to_numpy(), dtype=str)
    ds.add_axis("cells", cells, origin=OBSERVED, role="observation")
    ds.add_axis("genes", genes, origin=OBSERVED, role="feature")

    # In backed mode the heavy matrices stay on disk: wrap them as streaming sources keyed by their
    # h5ad location, so `write(stream=True)` copies them block-by-block. `fn` is None when not backed.
    fn = getattr(adata, "filename", None) if getattr(adata, "isbacked", False) else None

    if adata.X is not None:
        try:
            state = adata.uns.get("lstar/state")
        except Exception:
            state = None
        x = _backed_sparse(fn, "X") or adata.X      # backed -> streaming source; else the in-memory X
        ds.add_field("X", x, role="measure", span=["cells", "genes"],
                     state=state, provenance={"anndata": "X"})

    # .raw: older pipelines stash pre-HVG raw counts here, frequently over a *larger* gene set.
    # Recognize it and keep its own gene axis when the vocabulary differs (a within-object
    # collection of two feature spaces), rather than forcing it onto `genes` or dropping it.
    raw = getattr(adata, "raw", None)
    if raw is not None and getattr(raw, "X", None) is not None:
        raw_genes = np.asarray(np.asarray(raw.var_names), dtype=str)
        if raw_genes.shape[0] == genes.shape[0] and np.array_equal(raw_genes, genes):
            gax = "genes"
        else:
            gax = "genes_raw"
            ds.add_axis(gax, raw_genes, origin=OBSERVED, role="feature")
        rawx = _backed_sparse(fn, "raw/X", "raw.X") or raw.X    # modern vs legacy raw location
        ds.add_field("raw", rawx, role="measure", span=["cells", gax],
                     state="raw", provenance={"anndata": "raw/X"})

    for k in list(adata.layers.keys()):
        lk = _backed_sparse(fn, "layers/%s" % k) or adata.layers[k]
        ds.add_field(k, lk, role="measure", span=["cells", "genes"],
                     state=_guess_state(k), provenance={"anndata": "layers/%s" % k})

    for col in adata.obs.columns:
        vals, role = _by_dtype_series(adata.obs[col])
        ds.add_field(str(col), vals, role=role, span=["cells"],
                     provenance={"anndata": "obs/%s" % col})
    for col in adata.var.columns:
        vals, role = _by_dtype_series(adata.var[col])
        ds.add_field(str(col), vals, role=role, span=["genes"],
                     provenance={"anndata": "var/%s" % col})

    for k in list(adata.obsm.keys()):
        v = np.asarray(adata.obsm[k])
        name = _strip_x(k)
        cax = _coord_axis(ds, name, v.shape[1], observed=(name == "spatial"))
        ds.add_field(name, v, role="embedding", span=["cells", cax],
                     provenance={"anndata": "obsm/%s" % k})

    for k in list(adata.varm.keys()):
        v = np.asarray(adata.varm[k])
        cax = _pair_coord_axis(ds, k, v.shape[1])
        ds.add_field("%s_loadings" % cax, v, role="loading", span=["genes", cax],
                     provenance={"anndata": "varm/%s" % k})

    for k in list(adata.obsp.keys()):
        ds.add_field(k, adata.obsp[k], role="relation", span=["cells", "cells"],
                     subtype=_guess_subtype(k), weighted=True,
                     provenance={"anndata": "obsp/%s" % k})
    for k in list(adata.varp.keys()):
        ds.add_field("%s_varp" % k, adata.varp[k], role="relation", span=["genes", "genes"],
                     subtype=_guess_subtype(k), weighted=True,
                     provenance={"anndata": "varp/%s" % k})

    # uns is not imported in M2; record the loss (never silent)
    ds.dropped = ["uns/%s" % k for k in adata.uns.keys() if not str(k).startswith("lstar/")]
    return ds


def _vocab_location(ds, name, f):
    """Fallback native location for an L* field with no `anndata` provenance."""
    sp1 = f.span
    if f.role == "measure" and sp1 == ["cells", "genes"]:
        if name in ("counts",) or f.state == "raw":
            return "layers/counts" if name != "X" else "X"
        return "X"
    if f.role == "embedding" and len(sp1) == 2 and sp1[0] == "cells":
        return "obsm/X_%s" % name
    if f.role == "loading" and len(sp1) == 2 and sp1[0] == "genes":
        base = name[:-9] if name.endswith("_loadings") else name
        return "varm/%s" % ("PCs" if base == "pca" else base)
    if f.role in ("label", "measure") and sp1 == ["cells"]:
        return "obs/%s" % name
    if f.role in ("label", "measure") and sp1 == ["genes"]:
        return "var/%s" % name
    if f.role == "relation" and sp1 == ["cells", "cells"]:
        return "obsp/%s" % name
    if f.role == "relation" and sp1 == ["genes", "genes"]:
        return "varp/%s" % name
    return None


def _route_fields(ds):
    """Group an L* Dataset's fields by their AnnData destination -- the routing shared by the eager
    (`write_anndata`) and streamed (`write_anndata_streamed`) writers, so the two never disagree on
    where a field lands. Values are the L* `Field` objects (not yet materialized), keyed by their
    native sub-name; `X`/`raw` are single fields, the rest are dicts; unplaceable fields go to
    `dropped`. The destination is the field's recorded `anndata` provenance, else the shared-vocabulary
    fallback (`_vocab_location`)."""
    r = {"X": None, "raw": None, "layers": {}, "obs": {}, "var": {},
         "obsm": {}, "varm": {}, "obsp": {}, "varp": {}, "dropped": []}
    for name, f in ds.fields.items():
        loc = (f.provenance or {}).get("anndata") or _vocab_location(ds, name, f)
        if loc == "X":
            r["X"] = f
        elif loc == "raw/X":
            r["raw"] = f
        elif loc and "/" in loc and loc.split("/", 1)[0] in r and isinstance(r[loc.split("/", 1)[0]], dict):
            grp, key = loc.split("/", 1)
            r[grp][key] = f
        else:
            r["dropped"].append(name)
    return r


def write_anndata(ds):
    """Write an L* Dataset back to an AnnData object (lossy where no slot fits)."""
    import anndata as ad
    import pandas as pd

    cells = np.asarray(ds.axis("cells").labels, dtype=str)
    genes = np.asarray(ds.axis("genes").labels, dtype=str)
    r = _route_fields(ds)
    X = r["X"].values if r["X"] is not None else None
    raw_field = r["raw"]
    layers = {k: f.values for k, f in r["layers"].items()}
    obs = {k: _to_pandas_col(f) for k, f in r["obs"].items()}
    var = {k: _to_pandas_col(f) for k, f in r["var"].items()}
    obsm = {k: np.asarray(f.values) for k, f in r["obsm"].items()}
    varm = {k: np.asarray(f.values) for k, f in r["varm"].items()}
    obsp = {k: f.values for k, f in r["obsp"].items()}
    varp = {k: f.values for k, f in r["varp"].items()}
    dropped = list(r["dropped"])

    obs_df = pd.DataFrame(obs, index=cells) if obs else pd.DataFrame(index=cells)
    var_df = pd.DataFrame(var, index=genes) if var else pd.DataFrame(index=genes)
    adata = ad.AnnData(X=X, obs=obs_df, var=var_df,
                       layers=layers or None, obsm=obsm or None, varm=varm or None,
                       obsp=obsp or None, varp=varp or None)
    if raw_field is not None:
        rg = np.asarray(ds.axis(raw_field.span[1]).labels, dtype=str)
        raw_var = pd.DataFrame(index=rg)
        adata.raw = ad.AnnData(X=raw_field.values, obs=pd.DataFrame(index=cells), var=raw_var)
    if dropped:
        adata.uns["lstar/dropped"] = list(dropped)
    return adata


def convert_anndata(h5ad_path, store_path, **write_kwargs):
    """Convert an `.h5ad` to an L* store with **bounded memory**.

    Reads the source in backed mode (its `X`/layers/`.raw` stay on disk) and streams them block-by-
    block into the store, so even a multi-million-cell matrix never lands in RAM. The small parts
    (`obs`/`var`/`obsm`/graphs) are read normally. `write_kwargs` are forwarded to `lstar.write`
    (e.g. `compressor=numcodecs.GZip(5)`, `chunk_elems=...`); `stream=True` is set by default.

    Returns `store_path`. For the in-memory (fast, unbounded) path, use `read_anndata` +
    `lstar.write` on a non-backed AnnData instead.
    """
    import anndata as ad

    from ..zarr_io import write as _write

    adata = ad.read_h5ad(h5ad_path, backed="r")
    ds = read_anndata(adata)
    write_kwargs.setdefault("stream", True)
    try:
        _write(ds, store_path, **write_kwargs)
    finally:
        for f in ds.fields.values():            # close the on-disk handles the streaming sources hold
            if hasattr(f.values, "close"):
                f.values.close()
        try:
            adata.file.close()
        except Exception:
            pass
    return store_path


def _stream_sparse_to_h5(g, source, chunk_elems):
    """Stream a CSR/CSC source block-by-block into an open h5py group as a *modern* h5ad sparse
    group (`data`/`indices`/`indptr` + `encoding-type`/`encoding-version`/`shape`). `data`/`indices`
    are resizable and grown per block; `indptr` (small) is filled incrementally. The matrix is never
    fully resident, and its native orientation is preserved (no transpose). Mirrors the zarr sparse
    sink (`zarr_io._write_sparse_streaming`) -- same `iter_sized_blocks` policy, different store."""
    from ..lazy import iter_sized_blocks
    fmt = source.fmt
    n_outer = int(source.n_outer)
    nnz = int(getattr(source, "nnz", np.asarray(source.indptr)[-1]))
    data_dtype = np.dtype(getattr(source, "dtype", np.float32))
    idx_dtype = np.dtype(getattr(source, "idtype", np.int32))
    ce = int(chunk_elems) if chunk_elems else 1_000_000
    hchunk = (max(1, min(ce, nnz or 1)),)               # 1-D HDF5 chunk along the nonzero axis
    d = g.create_dataset("data", shape=(0,), maxshape=(None,), dtype=data_dtype, chunks=hchunk)
    ix = g.create_dataset("indices", shape=(0,), maxshape=(None,), dtype=idx_dtype, chunks=hchunk)
    indptr_dtype = np.int32 if 0 <= nnz < 2 ** 31 else np.int64   # anndata's convention (int32)
    out_indptr = np.empty(n_outer + 1, dtype=indptr_dtype)
    out_indptr[0] = 0
    pos = 0
    for a, b, sub in iter_sized_blocks(source, ce):
        sub = sub.tocsr() if fmt == "csr" else sub.tocsc()
        k = int(sub.nnz)
        d.resize((pos + k,)); d[pos:pos + k] = np.asarray(sub.data, dtype=data_dtype)
        ix.resize((pos + k,)); ix[pos:pos + k] = np.asarray(sub.indices, dtype=idx_dtype)
        out_indptr[a + 1:b + 1] = pos + sub.indptr[1:]
        pos += k
    g.create_dataset("indptr", data=out_indptr)
    g.attrs["encoding-type"] = "%s_matrix" % fmt
    g.attrs["encoding-version"] = "0.1.0"
    g.attrs["shape"] = np.asarray([int(x) for x in source.shape], dtype="int64")


def write_anndata_streamed(ds, path, chunk_elems=None):
    """Write an L* Dataset to an `.h5ad` with **bounded memory** -- the L*->native counterpart of
    `convert_anndata`.

    The small parts (obs/var/obsm/varm/obsp/varp/uns) are written through anndata, so their fiddly
    on-disk encoding (categoricals, dataframes) is reused rather than re-implemented. Then each large
    sparse measure (`X`, `raw/X`, every `layers/*`) is streamed straight into its h5ad sparse group
    block-by-block via h5py, so the whole matrix never lands in RAM. A measure that is dense (or
    otherwise not sparse) falls back to the in-memory skeleton. Peak memory is bounded by one block
    (~`chunk_elems` nonzeros, default 1e6).

    Sparse cell-cell / gene-gene graphs (`obsp`/`varp`) are streamed too, since on a large atlas a
    k-NN graph is itself big; only a (rare) dense graph falls back to the in-memory skeleton.

    Typical use is a fully bounded L* store -> h5ad conversion:
        write_anndata_streamed(lstar.read("atlas.lstar.zarr", lazy=True), "atlas.h5ad")
    where the lazy read leaves measures (and graphs) as on-disk `LazyCSX` sources that stream
    through untouched. Returns `path`.
    """
    import anndata as ad
    import h5py
    import pandas as pd

    from ..lazy import as_source

    cells = np.asarray(ds.axis("cells").labels, dtype=str)
    genes = np.asarray(ds.axis("genes").labels, dtype=str)
    r = _route_fields(ds)

    # Partition the big measures into streamable (sparse) and eager (dense/other). `stream_jobs` is a
    # list of (h5ad group path, source); eager pieces are built into the skeleton AnnData below.
    stream_jobs = []
    X_eager = None
    layers_eager = {}
    raw_eager = None

    if r["X"] is not None:
        s = as_source(r["X"].values)
        if s is not None:
            stream_jobs.append(("X", s))
        else:
            X_eager = np.asarray(r["X"].values)
    for k, f in r["layers"].items():
        s = as_source(f.values)
        if s is not None:
            stream_jobs.append(("layers/%s" % k, s))
        else:
            layers_eager[k] = np.asarray(f.values)

    obs = {k: _to_pandas_col(f) for k, f in r["obs"].items()}
    var = {k: _to_pandas_col(f) for k, f in r["var"].items()}
    obsm = {k: np.asarray(f.values) for k, f in r["obsm"].items()}
    varm = {k: np.asarray(f.values) for k, f in r["varm"].items()}

    # Graphs (obsp/varp) are sparse and, on a large atlas, big -- so stream the sparse ones too;
    # only a (rare) dense graph falls back to the in-memory skeleton.
    obsp, varp = {}, {}
    for grp, dest in (("obsp", obsp), ("varp", varp)):
        for k, f in r[grp].items():
            s = as_source(f.values)
            if s is not None:
                stream_jobs.append(("%s/%s" % (grp, k), s))
            else:
                dest[k] = np.asarray(f.values)

    obs_df = pd.DataFrame(obs, index=cells) if obs else pd.DataFrame(index=cells)
    var_df = pd.DataFrame(var, index=genes) if var else pd.DataFrame(index=genes)
    adata = ad.AnnData(X=X_eager, obs=obs_df, var=var_df,
                       layers=layers_eager or None, obsm=obsm or None, varm=varm or None,
                       obsp=obsp or None, varp=varp or None)

    # .raw: stream it too. anndata's writer needs a raw AnnData to lay down `raw/var`; give it a
    # zero-nonzero placeholder X of the right shape, then overwrite `raw/X` with the streamed data.
    raw_field = r["raw"]
    raw_src = as_source(raw_field.values) if raw_field is not None else None
    if raw_field is not None:
        rg = np.asarray(ds.axis(raw_field.span[1]).labels, dtype=str)
        if raw_src is not None:
            import scipy.sparse as sp
            placeholder = sp.csr_matrix((len(cells), len(rg)), dtype=raw_src.dtype)
            adata.raw = ad.AnnData(X=placeholder, obs=pd.DataFrame(index=cells),
                                   var=pd.DataFrame(index=rg))
        else:
            raw_eager = raw_field
            adata.raw = ad.AnnData(X=np.asarray(raw_field.values),
                                   obs=pd.DataFrame(index=cells), var=pd.DataFrame(index=rg))
    if r["dropped"]:
        adata.uns["lstar/dropped"] = list(r["dropped"])

    adata.write_h5ad(path)                  # lay down the skeleton (small parts + any eager measures)

    if raw_src is not None:
        stream_jobs.append(("raw/X", raw_src))   # replace the placeholder raw/X below

    # Re-open and stream each big sparse measure straight into its h5ad sparse group.
    try:
        with h5py.File(path, "a") as f:
            for gpath, src in stream_jobs:
                if gpath in f:
                    del f[gpath]                  # drop the placeholder/empty group anndata wrote
                _stream_sparse_to_h5(f.create_group(gpath), src, chunk_elems)
    finally:
        for _, src in stream_jobs:                # release any on-disk handles the sources held
            if hasattr(src, "close"):
                src.close()
    return path


def convert_to_h5ad(store_path, h5ad_path, chunk_elems=None):
    """Convert an L* store to an `.h5ad` with **bounded memory** -- the reverse of `convert_anndata`.

    Reads the store lazily (its measures stay on disk as `LazyCSX` sources) and streams them into
    the h5ad block-by-block, so a multi-gigabyte store converts in the memory of one block. Returns
    `h5ad_path`.
    """
    from ..zarr_io import read as _read
    return write_anndata_streamed(_read(store_path, lazy=True), h5ad_path, chunk_elems=chunk_elems)
