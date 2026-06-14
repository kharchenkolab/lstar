"""Package-free AnnData reader — read an ``.h5ad``'s on-disk HDF5 layout directly via ``h5py`` (no
``anndata`` package), producing the SAME core L* dataset ``read_anndata`` builds. This is the fallback
``lstar convert`` uses when ``anndata`` isn't installed (``--backend direct``/``auto``).

It targets the **modern group-based encoding** (anndata ≥ 0.7: ``encoding-type`` attrs, dataframe groups,
sparse ``data``/``indices``/``indptr`` groups). A layout it doesn't recognize — a legacy compound
``obs``/``var`` dataset, an unknown column encoding — raises :class:`NeedsPackage` naming ``anndata``, so
the user knows exactly what to install. Coverage is the shared **core** (X, raw, layers, obs/var fields
with categorical→factor induction, obsm/varm/obsp/varp, and ``uns`` preserved verbatim in ``aux``); the
``uns``→typed-field promotions (rank_genes_groups, ``*_colors``, velocity graphs) remain a native-only
enhancement — the data is still carried in ``aux``, just not promoted.
"""
from __future__ import annotations

import numpy as np

from ..model import Categorical, Dataset, OBSERVED
from .anndata import (_coord_axis, _guess_state, _guess_subtype, _pair_coord_axis,
                      _sparse_attrs, _strip_x, _uniq)


def _needs(thing, package="anndata", install="pip install anndata"):
    from ..cli import NeedsPackage
    return NeedsPackage(thing, package, install)


def _str_arr(a) -> np.ndarray:
    """An HDF5 string dataset/array → a numpy ``str`` array (handles fixed-S, var-len object, unicode)."""
    a = np.asarray(a)
    if a.dtype.kind == "S":
        return np.char.decode(a, "utf-8")
    if a.dtype.kind == "O":
        flat = [x.decode("utf-8") if isinstance(x, bytes) else str(x) for x in a.ravel()]
        return np.array(flat, dtype=str).reshape(a.shape)
    return a.astype(str)


def _read_matrix(node):
    """A dense ndarray or a scipy sparse matrix from an h5ad ``X``/layer node (dataset or sparse group)."""
    import h5py
    if isinstance(node, h5py.Dataset):
        return node[...]
    sp_attrs = _sparse_attrs(node)
    if sp_attrs:
        import scipy.sparse as sp
        fmt, shape = sp_attrs
        cls = sp.csr_matrix if fmt == "csr" else sp.csc_matrix
        return cls((node["data"][...], node["indices"][...], node["indptr"][...]), shape=shape)
    raise _needs(f"an unrecognized matrix encoding ({dict(node.attrs)!r})")


def _read_column(group, col):
    """One obs/var column → ``(values, role, mask)``, mirroring the native ``_by_dtype_series`` mapping."""
    import h5py
    node = group[col]
    if isinstance(node, h5py.Group):
        et = str(node.attrs.get("encoding-type", ""))
        if et == "categorical":
            codes = np.asarray(node["codes"][...]).astype(np.int64)
            cats = _str_arr(node["categories"][...])
            return (Categorical(codes, cats, ordered=bool(node.attrs.get("ordered", False))),
                    "label", None)
        if et in ("nullable-integer", "nullable-boolean"):
            vals = np.asarray(node["values"][...])
            mask = np.asarray(node["mask"][...]).astype(np.uint8)        # anndata mask: True == missing
            if et == "nullable-boolean":
                return vals.astype(bool), "label", mask
            return vals.astype(np.int64), "measure", mask
        raise _needs(f"obs/var column {col!r} has an unsupported encoding {et!r}")
    a = node[...]
    if a.dtype.kind == "b":
        return np.asarray(a), "label", None
    if a.dtype.kind in ("i", "u"):
        return np.asarray(a).astype(np.int64), "measure", None
    if a.dtype.kind == "f":
        return np.asarray(a).astype(np.float64), "measure", None
    return _str_arr(a), "label", None


def _df_columns(group):
    """Column names of an h5ad dataframe group, in stored order, excluding the index."""
    idx = group.attrs.get("_index", "_index")
    order = group.attrs.get("column-order")
    if order is not None:
        cols = [c.decode() if isinstance(c, bytes) else str(c) for c in np.asarray(order).ravel()]
    else:
        cols = [k for k in group.keys() if k != idx]
    return [c for c in cols if c != idx and c in group]


def _read_uns(node):
    """Read a uns group/dataset recursively into a plain nested dict/array for the ``aux`` passthrough."""
    import h5py
    if isinstance(node, h5py.Group):
        if _sparse_attrs(node):
            return _read_matrix(node)
        if str(node.attrs.get("encoding-type", "")) == "lstar-record":   # reconstruct a structured array
            cols = {k: np.asarray(_read_uns(node[k])) for k in node.keys()}
            n = int(node.attrs.get("length", len(next(iter(cols.values()))) if cols else 0))
            out = np.empty(n, dtype=[(k, v.dtype) for k, v in cols.items()])
            for k, v in cols.items():
                out[k] = v
            return out
        return {k: _read_uns(node[k]) for k in node.keys() if not str(k).startswith("lstar")}
    a = node[...]
    if np.ndim(a) == 0:
        v = a.item() if hasattr(a, "item") else a
        return v.decode("utf-8") if isinstance(v, bytes) else v
    if np.asarray(a).dtype.kind in ("S", "O"):
        return _str_arr(a)
    return np.asarray(a)


def _require_modern(f, group):
    import h5py
    if group not in f:
        raise _needs(f"this .h5ad has no '{group}'")
    if not isinstance(f[group], h5py.Group):
        raise _needs(f"this .h5ad uses a legacy compound '{group}' (pre-0.7 layout)")


def read_h5ad_direct(path: str) -> Dataset:
    """Read *path* into an L* :class:`Dataset` using only ``h5py`` (no ``anndata``)."""
    try:
        import h5py
    except ImportError:
        raise _needs("reading .h5ad without the anndata package", "h5py",
                     "pip install h5py   (the package-free path; or: pip install anndata)")

    with h5py.File(path, "r") as f:
        _require_modern(f, "obs")
        _require_modern(f, "var")
        ds = Dataset(kind="sample")
        cells = _str_arr(f["obs"][f["obs"].attrs.get("_index", "_index")][...])
        genes = _str_arr(f["var"][f["var"].attrs.get("_index", "_index")][...])
        ds.add_axis("cells", cells, origin=OBSERVED, role="observation")
        ds.add_axis("genes", genes, origin=OBSERVED, role="feature")

        state = None
        if "uns/lstar/state" in f:                                       # lstar's own state marker
            sv = f["uns/lstar/state"][()]
            state = sv.decode("utf-8") if isinstance(sv, bytes) else sv

        if "X" in f:
            ds.add_field("X", _read_matrix(f["X"]), role="measure", span=["cells", "genes"],
                         state=state, provenance={"anndata": "X"})

        if "raw" in f and "raw/X" in f:                                  # pre-HVG raw, possibly larger gene set
            raw_genes = _str_arr(f["raw/var"][f["raw/var"].attrs.get("_index", "_index")][...]) \
                if "raw/var" in f else genes
            if raw_genes.shape[0] == genes.shape[0] and np.array_equal(raw_genes, genes):
                gax = "genes"
            else:
                gax = "genes_raw"
                ds.add_axis(gax, raw_genes, origin=OBSERVED, role="feature")
            ds.add_field("raw", _read_matrix(f["raw/X"]), role="measure", span=["cells", gax],
                         state="raw", provenance={"anndata": "raw/X"})

        if "layers" in f:
            for k in f["layers"].keys():
                ds.add_field(_uniq(ds, str(k), "layer"), _read_matrix(f["layers"][k]), role="measure",
                             span=["cells", "genes"], state=_guess_state(k),
                             provenance={"anndata": "layers/%s" % k})

        for axis, group in (("cells", "obs"), ("genes", "var")):
            g = f[group]
            for col in _df_columns(g):
                vals, role, mask = _read_column(g, col)
                ds.add_field(_uniq(ds, str(col), axis), vals, role=role, span=[axis], mask=mask,
                             provenance={"anndata": "%s/%s" % (group, col)})

        if "obsm" in f:
            for k in f["obsm"].keys():
                v = np.asarray(f["obsm"][k][...])
                name = _strip_x(k)
                cax = _coord_axis(ds, name, v.shape[1], observed=(name == "spatial"))
                ds.add_field(name, v, role="embedding", span=["cells", cax],
                             subtype=("spatial" if name == "spatial" else None),
                             provenance={"anndata": "obsm/%s" % k})
        if "varm" in f:
            for k in f["varm"].keys():
                v = np.asarray(f["varm"][k][...])
                cax = _pair_coord_axis(ds, k, v.shape[1])
                ds.add_field("%s_loadings" % cax, v, role="loading", span=["genes", cax],
                             provenance={"anndata": "varm/%s" % k})
        if "obsp" in f:
            for k in f["obsp"].keys():
                ds.add_field(k, _read_matrix(f["obsp"][k]), role="relation", span=["cells", "cells"],
                             subtype=_guess_subtype(k), weighted=True,
                             provenance={"anndata": "obsp/%s" % k})
        if "varp" in f:
            for k in f["varp"].keys():
                ds.add_field("%s_varp" % k, _read_matrix(f["varp"][k]), role="relation",
                             span=["genes", "genes"], subtype=_guess_subtype(k), weighted=True,
                             provenance={"anndata": "varp/%s" % k})

        if "uns" in f:
            uns = {k: _read_uns(f["uns"][k]) for k in f["uns"].keys() if not str(k).startswith("lstar")}
            if uns:
                ds.aux["anndata.uns"] = uns

        ds.profiles = ["anndata@direct"]
    return ds


# ── package-free writer (h5py): emit the pinned modern anndata encoding ───────────────────────────────
# Pinned versions: anndata 0.1.0 / array 0.2.0 / csr|csc_matrix 0.1.0 / dataframe 0.2.0 / categorical
# 0.2.0 / string-array 0.2.0 / nullable-* 0.1.0 / dict 0.1.0. "latest/recent so everything fits"; native
# anndata >= 0.8 reads it (verified by the native-acceptance check + the direct-vs-native cross-validation).
_EV = {"anndata": "0.1.0", "array": "0.2.0", "dataframe": "0.2.0", "categorical": "0.2.0",
       "string-array": "0.2.0", "nullable-integer": "0.1.0", "nullable-boolean": "0.1.0",
       "dict": "0.1.0", "string": "0.2.0", "numeric-scalar": "0.2.0"}


def _vlen():
    import h5py
    return h5py.special_dtype(vlen=str)


def _w_strings(parent, name, arr):
    a = np.asarray(arr)
    data = np.array([("" if x is None else str(x)) for x in a.ravel()], dtype=object).reshape(a.shape)
    d = parent.create_dataset(name, data=data, dtype=_vlen())
    d.attrs["encoding-type"], d.attrs["encoding-version"] = "string-array", _EV["string-array"]
    return d


def _w_array(parent, name, arr):
    arr = np.asarray(arr)
    if arr.dtype.kind in ("U", "S", "O"):
        return _w_strings(parent, name, arr)
    d = parent.create_dataset(name, data=arr)
    d.attrs["encoding-type"], d.attrs["encoding-version"] = "array", _EV["array"]
    return d


def _w_matrix(parent, name, m):
    import scipy.sparse as sp
    if sp.issparse(m):
        fmt = "csr" if sp.isspmatrix_csr(m) else ("csc" if sp.isspmatrix_csc(m) else None)
        if fmt is None:
            m, fmt = m.tocsr(), "csr"
        g = parent.create_group(name)
        g.attrs["encoding-type"], g.attrs["encoding-version"] = f"{fmt}_matrix", "0.1.0"
        g.attrs["shape"] = np.asarray(m.shape, dtype="int64")
        g.create_dataset("data", data=m.data)
        g.create_dataset("indices", data=m.indices)
        g.create_dataset("indptr", data=m.indptr)
        return g
    return _w_array(parent, name, np.asarray(m))


def _w_column(group, name, fld):
    v = fld.values
    if isinstance(v, Categorical):
        g = group.create_group(name)
        g.attrs["encoding-type"], g.attrs["encoding-version"] = "categorical", _EV["categorical"]
        g.attrs["ordered"] = bool(v.ordered)
        cd = g.create_dataset("codes", data=np.asarray(v.codes).astype("int32"))
        cd.attrs["encoding-type"], cd.attrs["encoding-version"] = "array", _EV["array"]
        _w_strings(g, "categories", np.asarray(v.categories, dtype=str))
        return
    arr = np.asarray(v)
    mask = getattr(fld, "mask", None)
    if mask is not None and arr.dtype.kind in ("b", "i", "u"):       # nullable Int/boolean (True == missing)
        et = "nullable-boolean" if arr.dtype.kind == "b" else "nullable-integer"
        g = group.create_group(name)
        g.attrs["encoding-type"], g.attrs["encoding-version"] = et, _EV[et]
        _w_array(g, "values", arr.astype(bool) if arr.dtype.kind == "b" else arr.astype("int64"))
        g.create_dataset("mask", data=np.asarray(mask).astype(bool))
        return
    _w_array(group, name, arr)                                       # plain numeric / bool / string


def _w_dataframe(root, name, index, items):
    import h5py
    g = root.create_group(name)
    g.attrs["encoding-type"], g.attrs["encoding-version"] = "dataframe", _EV["dataframe"]
    g.attrs["_index"] = "_index"
    cols = [c for c, _ in items]
    g.attrs.create("column-order", data=np.array(cols, dtype=object),
                   shape=(len(cols),), dtype=h5py.special_dtype(vlen=str))
    _w_strings(g, "_index", np.asarray(index, dtype=str))
    for cname, fld in items:
        _w_column(g, cname, fld)


def _w_uns(parent, name, obj):
    if isinstance(obj, dict):
        g = parent.create_group(name)
        g.attrs["encoding-type"], g.attrs["encoding-version"] = "dict", _EV["dict"]
        for k, v in obj.items():
            _w_uns(g, str(k), v)
    elif isinstance(obj, str):
        d = parent.create_dataset(name, data=obj)
        d.attrs["encoding-type"], d.attrs["encoding-version"] = "string", _EV["string"]
    elif isinstance(obj, (bool, np.bool_, int, float, np.integer, np.floating)):
        d = parent.create_dataset(name, data=obj)
        d.attrs["encoding-type"], d.attrs["encoding-version"] = "numeric-scalar", _EV["numeric-scalar"]
    elif isinstance(obj, np.ndarray) and obj.dtype.names:           # structured/record array (e.g.
        g = parent.create_group(name)                              # rank_genes_groups['names']): a group of
        g.attrs["encoding-type"], g.attrs["encoding-version"] = "lstar-record", "0.1.0"   # per-field arrays
        g.attrs["length"] = int(obj.shape[0]) if obj.ndim else 0   # (h5py can't write a numpy U-dtype compound)
        for fld in obj.dtype.names:
            _w_uns(g, str(fld), np.asarray(obj[fld]))
    elif isinstance(obj, np.ndarray):
        _w_array(parent, name, obj)
    else:
        _w_array(parent, name, np.asarray(obj))                     # best effort for the long tail


def write_h5ad_direct(ds, path: str) -> str:
    """Write an L* :class:`Dataset` to an ``.h5ad`` using only ``h5py`` (no ``anndata``), in the pinned
    modern encoding. Reuses ``write_anndata``'s field routing so a field lands in the same slot as native."""
    try:
        import h5py
    except ImportError:
        raise _needs("writing .h5ad without the anndata package", "h5py",
                     "pip install h5py   (the package-free path; or: pip install anndata)")
    from .anndata import _route_fields
    r = _route_fields(ds)
    cells = np.asarray(ds.axis("cells").labels, dtype=str)
    genes = np.asarray(ds.axis("genes").labels, dtype=str)

    with h5py.File(path, "w") as f:
        f.attrs["encoding-type"], f.attrs["encoding-version"] = "anndata", _EV["anndata"]
        if r["X"] is not None:
            _w_matrix(f, "X", r["X"].values)
        _w_dataframe(f, "obs", cells, list(r["obs"].items()))
        _w_dataframe(f, "var", genes, list(r["var"].items()))

        for grp, is_matrix in (("layers", True), ("obsm", False), ("varm", False),
                               ("obsp", True), ("varp", True)):
            items = r[grp]
            if not items:
                continue
            g = f.create_group(grp)
            g.attrs["encoding-type"], g.attrs["encoding-version"] = "dict", _EV["dict"]
            for k, fld in items.items():
                _w_matrix(g, k, fld.values) if is_matrix else _w_array(g, k, np.asarray(fld.values))

        if r["raw"] is not None:                                    # pre-HVG raw (its own gene axis)
            rg = f.create_group("raw")
            _w_matrix(rg, "X", r["raw"].values)
            _w_dataframe(rg, "var", np.asarray(ds.axis(r["raw"].span[1]).labels, dtype=str), [])

        uns = dict((getattr(ds, "aux", None) or {}).get("anndata.uns", {}))
        xstate = r["X"].state if (r["X"] is not None and r["X"].state) else None
        if r["dropped"]:
            uns["lstar/dropped"] = list(r["dropped"])
        if uns or xstate:
            ug = f.create_group("uns")
            ug.attrs["encoding-type"], ug.attrs["encoding-version"] = "dict", _EV["dict"]
            for k, v in uns.items():
                _w_uns(ug, str(k), v)
            if xstate:                                              # lstar's own X-state marker
                lg = ug.create_group("lstar")
                lg.attrs["encoding-type"], lg.attrs["encoding-version"] = "dict", _EV["dict"]
                _w_uns(lg, "state", xstate)
    return path
