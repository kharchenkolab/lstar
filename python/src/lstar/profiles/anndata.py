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
    """An obs/var column -> (values, role, mask). pandas **nullable extension** dtypes (Int64/boolean/
    string) carry a validity mask so integer-ness and the value-vs-missing distinction survive -- the
    same silent-corruption class P1 fixed for categoricals. Float keeps NaN (no mask)."""
    import pandas as pd
    if isinstance(s.dtype, pd.CategoricalDtype):       # preserve category order + NaN (-> -1 missing)
        c = s.values
        return Categorical(np.asarray(c.codes), np.asarray(c.categories, dtype=str),
                           ordered=bool(c.ordered)), "label", None
    if pd.api.types.is_extension_array_dtype(s):       # nullable Int/boolean/string: values + null mask
        mask = np.asarray(s.isna(), dtype=np.uint8)
        if pd.api.types.is_bool_dtype(s):
            return np.asarray(s.fillna(False)).astype(bool), "label", mask
        if pd.api.types.is_integer_dtype(s):
            return np.asarray(s.fillna(0)).astype(np.int64), "measure", mask
        if pd.api.types.is_float_dtype(s):
            return np.asarray(s.astype("float64")), "measure", None     # NaN already encodes missing
        return np.asarray(s.fillna("").astype(str).values, dtype=str), "label", mask
    if pd.api.types.is_bool_dtype(s):
        return np.asarray(s.values), "label", None
    if pd.api.types.is_numeric_dtype(s):
        return np.asarray(s.values), "measure", None
    return np.asarray(s.astype(str).values, dtype=str), "label", None


def _to_pandas_col(f):
    """An L* obs/var field as a column value for write-back: a `Categorical` -> `pd.Categorical`
    (categories/order/`-1`->NaN preserved); a masked field -> the matching pandas **nullable** dtype
    (Int64/boolean/string), so a round-trip is type-faithful; else a plain array."""
    if _is_categorical(f.values):
        import pandas as pd
        c = f.values
        return pd.Categorical.from_codes(np.asarray(c.codes), categories=list(c.categories),
                                         ordered=bool(c.ordered))
    arr = np.asarray(f.values)
    if getattr(f, "mask", None) is not None:
        import pandas as pd
        m = np.asarray(f.mask).astype(bool)
        if arr.dtype.kind == "b":
            return pd.arrays.BooleanArray(arr.astype(bool), m)
        if arr.dtype.kind in ("i", "u"):
            return pd.arrays.IntegerArray(arr.astype("int64"), m)
        if arr.dtype.kind in ("U", "S", "O"):
            obj = arr.astype(object); obj[m] = pd.NA
            return pd.array(obj, dtype="string")
    return arr


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
        vals, role, mask = _by_dtype_series(adata.obs[col])
        ds.add_field(str(col), vals, role=role, span=["cells"], mask=mask,
                     provenance={"anndata": "obs/%s" % col})
    for col in adata.var.columns:
        vals, role, mask = _by_dtype_series(adata.var[col])
        ds.add_field(str(col), vals, role=role, span=["genes"], mask=mask,
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

    # uns: preserved *verbatim* via lossless passthrough (params, colors, dendrograms, DE tables, ...)
    # rather than recorded name-only -- the `aux/` subtree round-trips it and a reader can later promote
    # recognized structures out of the tail. lstar-internal markers (e.g. lstar/state) are excluded.
    # Capture uns from its *raw backing dict*, not the OverloadedDict view: `adata.uns['neighbors']`
    # is overloaded to inject the obsp graphs (connectivities/distances) into the returned dict -- we
    # don't want those (they're already captured as `relation`s, and they'd become None through the
    # passthrough and break anndata's neighbors-setter on restore). `_safe_uns_copy` then spine-copies
    # so promotions can pop safely without touching the caller (and never follows OverloadedDict cycles).
    raw_uns = getattr(adata.uns, "data", None)
    src = raw_uns if isinstance(raw_uns, dict) else adata.uns
    uns = {k: _safe_uns_copy(v) for k, v in src.items() if not str(k).startswith("lstar/")}
    if uns:
        ds.aux["anndata.uns"] = uns
    _read_rank_genes_groups(adata, ds)              # type one-vs-rest DE; leave pairwise in passthrough
    _read_uns_promotions(adata, ds)                 # promote color palettes + PCA variance out of the tail
    _read_uns_graphs(adata, ds)                     # type uns cell/gene graphs (scVelo velocity_graph)
    return ds


# Subtypes of fields that are regenerated into `uns` on write (so they aren't routed as ordinary fields).
_TYPED_UNS_SUBTYPES = {"de", "color", "pca_var", "uns_graph"}


def _read_uns_graphs(adata, ds):
    """Type top-level `uns` sparse **square** matrices as relations over the matching axis. scVelo writes
    `velocity_graph` / `velocity_graph_neg` (cell x cell transition graphs) to `uns`, not `obsp` -- so
    without this they'd land in the passthrough and be dropped (the aux tree can't hold a sparse leaf).
    A (n_cells x n_cells) matrix becomes a `relation` over (cells, cells); (n_genes x n_genes) over
    (genes, genes). Removed from the passthrough; regenerated on write. (A PAGA cluster x cluster graph
    is nested under `uns['paga']` and sized to the cluster axis, so it is not matched here.)"""
    import scipy.sparse as sp
    uns = (getattr(ds, "aux", None) or {}).get("anndata.uns")
    if not uns:
        return
    ncell = len(ds.axes["cells"]) if "cells" in ds.axes else -1
    ngene = len(ds.axes["genes"]) if "genes" in ds.axes else -1
    for key in [k for k in list(uns) if sp.issparse(uns[k])]:
        m = uns[key]
        ax = "cells" if m.shape == (ncell, ncell) else "genes" if m.shape == (ngene, ngene) else None
        if ax is None:
            continue
        ds.add_field(str(key), m, role="relation", span=[ax, ax], subtype="uns_graph",
                     directed=True, weighted=True, provenance={"anndata": "uns/%s" % key})
        uns.pop(key, None)


def _write_uns_graphs(adata, ds):
    """Regenerate the typed `uns` graphs (e.g. velocity_graph) from their relation fields."""
    for name, f in ds.fields.items():
        if f.subtype == "uns_graph":
            adata.uns[str(name)] = f.values


def _safe_uns_copy(v, _depth=0):
    """A mutation-safe copy of an `uns` value that copies the dict/list *spine* and references
    array/scalar leaves. Unlike `copy.deepcopy` it never follows anndata's `OverloadedDict` internal
    back-references (to the AnnData / obsp), so it doesn't explode on real `uns['neighbors']`."""
    if _depth < 32 and isinstance(v, dict):
        return {k: _safe_uns_copy(val, _depth + 1) for k, val in v.items()}
    if _depth < 32 and isinstance(v, (list, tuple)):
        return [_safe_uns_copy(x, _depth + 1) for x in v]
    return v


def _read_uns_promotions(adata, ds):
    """Promote recognized structures out of the lossless passthrough into typed fields **bound to their
    axes** -- the payoff of capturing `uns` verbatim first, then typing incrementally:

    - **color palettes** `uns['<key>_colors']` -> a `<key>_colors` label field over the factor axis
      `<key>` (one color per category, in category order; reordering the factor re-permutes the palette);
    - **PCA variance** `uns['pca']['variance'|'variance_ratio']` -> a measure over the `pca` coordinate
      axis (the same axis as the `X_pca` scores / `PCs` loadings).

    Promoted entries are removed from the passthrough so they aren't double-stored; everything else stays
    in the tail."""
    uns = (getattr(ds, "aux", None) or {}).get("anndata.uns")
    if not uns:
        return
    for key in [k for k in list(uns) if str(k).endswith("_colors")]:
        base = str(key)[:-len("_colors")]
        ax = ds.axes.get(base)
        if ax is not None and ax.role == "factor":
            colors = np.asarray(uns[key], dtype=str)
            if colors.ndim == 1 and colors.shape[0] == len(ax):
                ds.add_field(str(key), colors, role="label", span=[base], subtype="color",
                             provenance={"anndata": "uns/%s" % key})
                uns.pop(key, None)
    pca = uns.get("pca")
    if isinstance(pca, dict) and "pca" in ds.axes:
        n = len(ds.axes["pca"])
        for src, fld in (("variance", "pca_variance"), ("variance_ratio", "pca_variance_ratio")):
            v = pca.get(src)
            if v is not None and np.asarray(v).ndim == 1 and np.asarray(v).shape[0] == n:
                ds.add_field(fld, np.asarray(v, dtype=float), role="measure", span=["pca"],
                             subtype="pca_var", provenance={"anndata": "uns/pca/%s" % src})
                pca.pop(src, None)
        if not pca:
            uns.pop("pca", None)


def _write_uns_promotions(adata, ds):
    """Regenerate the promoted `uns` structures (color palettes, PCA variance) from their typed fields."""
    for name, f in ds.fields.items():
        if f.subtype == "color":
            adata.uns[str(name)] = np.asarray(f.values, dtype=str)
        elif f.subtype == "pca_var":
            src = str((f.provenance or {}).get("anndata", "")).rsplit("/", 1)[-1] or str(name)
            pca = adata.uns.get("pca")
            if not isinstance(pca, dict):
                pca = {}
                adata.uns["pca"] = pca
            pca[src] = np.asarray(f.values, dtype=float)


def _read_rank_genes_groups(adata, ds):
    """Type a **one-vs-rest** `uns['rank_genes_groups']` into a DE bundle over `(factor, gene-axis)` and
    remove it from the passthrough (so it isn't stored twice). Pairwise / reference-group results are
    left in the passthrough verbatim -- the documented asymmetry (only the representable shape is typed).
    The per-group ranking `names` is *not* stored: it is recoverable by argsort of the scores on write."""
    from ..de import de_field_name
    rgg = adata.uns.get("rank_genes_groups")
    if not isinstance(rgg, dict) or "names" not in rgg:
        return
    params = dict(rgg.get("params", {}))
    if str(params.get("reference", "rest")) != "rest":      # pairwise / ref-group -> passthrough only
        return
    factor = str(params.get("groupby", ""))
    if not factor:
        return
    names = rgg["names"]
    groups = list(names.dtype.names)
    use_raw = bool(params.get("use_raw", False))
    gene_axis = "genes_raw" if (use_raw and "genes_raw" in ds.axes) else "genes"
    glabels = np.asarray(ds.axis(gene_axis).labels, dtype=str)
    gidx = {g: i for i, g in enumerate(glabels)}
    if factor in ds.axes and ds.axes[factor].role == "factor":   # reuse the induced clustering axis
        forder = list(np.asarray(ds.axis(factor).labels, dtype=str))
    else:                                                        # or derive one from the DE groups
        forder = [str(g) for g in groups]
        ds.add_axis(factor, forder, origin=DERIVED, role="factor")
    grow = {g: i for i, g in enumerate(forder)}
    n_f, n_g = len(forder), len(glabels)
    stat_keys = {"scores": "score", "logfoldchanges": "lfc", "pvals": "pval", "pvals_adj": "padj"}
    made = False
    for key, stat in stat_keys.items():
        if key not in rgg:
            continue
        rec = rgg[key]
        arr = np.full((n_f, n_g), np.nan, dtype=np.float64)     # NaN where a group didn't rank a gene
        for grp in groups:
            r = grow.get(str(grp))
            if r is None:
                continue
            cols = np.array([gidx.get(nm, -1) for nm in np.asarray(names[grp]).astype(str)])
            ok = cols >= 0
            arr[r, cols[ok]] = np.asarray(rec[grp], dtype=float)[ok]
        ds.add_field(de_field_name(factor, stat), arr, role="measure", span=[factor, gene_axis],
                     subtype="de", provenance={"anndata": "uns/rank_genes_groups", "de_factor": factor,
                                               "de_stat": stat, "method": params.get("method"),
                                               "reference": "rest", "use_raw": use_raw})
        made = True
    if made and "anndata.uns" in ds.aux:
        ds.aux["anndata.uns"].pop("rank_genes_groups", None)   # typed now -> drop from the passthrough


def _write_rank_genes_groups(adata, ds):
    """Regenerate `uns['rank_genes_groups']` from the typed DE bundle: rank each group's genes by score
    (argsort) to rebuild the per-group `names`, and emit the structured arrays + `params`."""
    from ..de import de_bundle, de_factors
    for factor in de_factors(ds):
        bundle = de_bundle(ds, factor)
        any_f = next(iter(bundle.values()))
        gene_axis = any_f.span[1]
        groups = [str(g) for g in np.asarray(ds.axis(factor).labels, dtype=str)]
        glabels = np.asarray(ds.axis(gene_axis).labels, dtype=str)
        prov = any_f.provenance or {}
        basis = bundle.get("score") or bundle.get("lfc") or any_f
        bvals = np.asarray(basis.values, dtype=float)
        orders = [np.argsort(-np.nan_to_num(bvals[gi], nan=-np.inf), kind="stable")
                  for gi in range(len(groups))]
        names_rec = np.empty(len(glabels), dtype=[(g, glabels.dtype) for g in groups])
        for gi, g in enumerate(groups):
            names_rec[g] = glabels[orders[gi]]
        rgg = {"params": {"groupby": factor, "reference": "rest",
                          "method": prov.get("method") or "lstar", "use_raw": bool(prov.get("use_raw", False))},
               "names": names_rec}
        for st, key in {"score": "scores", "lfc": "logfoldchanges",
                        "pval": "pvals", "padj": "pvals_adj"}.items():
            if st not in bundle:
                continue
            vals = np.asarray(bundle[st].values, dtype=float)
            rec = np.empty(len(glabels), dtype=[(g, "f4") for g in groups])
            for gi, g in enumerate(groups):
                rec[g] = vals[gi][orders[gi]].astype("f4")
            rgg[key] = rec
        adata.uns["rank_genes_groups"] = rgg


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
        if f.subtype in _TYPED_UNS_SUBTYPES or (f.provenance or {}).get("anndata") == "uns/rank_genes_groups":
            continue                            # DE / colors / pca-var -> regenerated into uns, not routed
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


def _restore_uns(adata, ds):
    """Reproduce the passthrough `uns` captured by `read_anndata` (lossless round-trip of params, color
    palettes, dendrograms, DE tables, ...). Typed fields already placed elsewhere win; uns only fills
    keys not otherwise set."""
    uns = (getattr(ds, "aux", None) or {}).get("anndata.uns")
    if uns:
        for k, v in uns.items():
            adata.uns[k] = v


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
    _restore_uns(adata, ds)                         # reproduce the passthrough uns (the untyped tail)
    _write_rank_genes_groups(adata, ds)             # regenerate the typed DE bundle into uns
    _write_uns_promotions(adata, ds)                # regenerate promoted colors + PCA variance into uns
    _write_uns_graphs(adata, ds)                    # regenerate uns cell/gene graphs (velocity_graph)
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
    _restore_uns(adata, ds)                 # reproduce the passthrough uns
    _write_rank_genes_groups(adata, ds)     # regenerate the typed DE bundle into uns
    _write_uns_promotions(adata, ds)        # regenerate promoted colors + PCA variance into uns
    _write_uns_graphs(adata, ds)            # regenerate uns cell/gene graphs (velocity_graph)
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
