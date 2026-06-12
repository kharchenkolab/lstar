"""AnnData profile (scverse): read_anndata / write_anndata.

Maps an AnnData object to/from an L* Dataset. The native location of each field is
recorded in field.provenance["anndata"] so write-back is exact; field *names* follow the
shared core vocabulary (X_pca -> 'pca', etc.) so cross-format conversion is meaningful.

anndata is imported lazily, so `import lstar` never requires it.
"""
import numpy as np

from ..model import Dataset, OBSERVED, DERIVED

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
    if pd.api.types.is_numeric_dtype(s):
        return np.asarray(s.values), "measure"
    return np.asarray(s.astype(str).values, dtype=str), "label"


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


def write_anndata(ds):
    """Write an L* Dataset back to an AnnData object (lossy where no slot fits)."""
    import anndata as ad
    import pandas as pd

    cells = np.asarray(ds.axis("cells").labels, dtype=str)
    genes = np.asarray(ds.axis("genes").labels, dtype=str)
    X = None
    raw_field = None
    layers, obs, var, obsm, varm, obsp, varp = {}, {}, {}, {}, {}, {}, {}
    dropped = []

    for name, f in ds.fields.items():
        loc = (f.provenance or {}).get("anndata") or _vocab_location(ds, name, f)
        if loc is None:
            dropped.append(name)
            continue
        if loc == "X":
            X = f.values
        elif loc == "raw/X":
            raw_field = f
        elif loc.startswith("layers/"):
            layers[loc.split("/", 1)[1]] = f.values
        elif loc.startswith("obs/"):
            obs[loc.split("/", 1)[1]] = np.asarray(f.values)
        elif loc.startswith("var/"):
            var[loc.split("/", 1)[1]] = np.asarray(f.values)
        elif loc.startswith("obsm/"):
            obsm[loc.split("/", 1)[1]] = np.asarray(f.values)
        elif loc.startswith("varm/"):
            varm[loc.split("/", 1)[1]] = np.asarray(f.values)
        elif loc.startswith("obsp/"):
            obsp[loc.split("/", 1)[1]] = f.values
        elif loc.startswith("varp/"):
            varp[loc.split("/", 1)[1]] = f.values
        else:
            dropped.append(name)

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
