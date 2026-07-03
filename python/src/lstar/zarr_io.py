"""L* Zarr serialization (Python reference implementation).

Writes and reads the axes/ fields/ models/ registry of the proposal (Appendix A). All L*
metadata lives under an "lstar" key in each group's attributes.

Cross-language notes (so the C++ core can read these stores):
  - strings (axis labels, string-valued fields) use a UTF-8 `data` (uint8) + `offsets`
    (int64, length n+1) encoding rather than fixed-width unicode arrays;
  - `compressor=None` by default writes uncompressed chunks (no codec dependency in C++);
  - targets Zarr v2 here (Python 3.8); written to be v3-ready.
"""
import json
import os

import numpy as np
import scipy.sparse as sp
import zarr

from .lazy import is_stream_source as _is_stream_source
from .lazy import iter_sized_blocks
from .model import Dataset, Field, Categorical, _is_categorical, as_categorical

LSTAR = "lstar"
SPEC_VERSION = "0.1"


def write(ds, path, compressor=None, chunk_elems=None, stream=False, viewer=False):
    """Serialize an L* Dataset to a Zarr store at `path`.

    compressor=None (default) writes uncompressed chunks; pass a numcodecs codec (e.g.
    numcodecs.GZip(5)) to compress. chunk_elems=None writes each array as a single chunk
    (the portable default); set it (e.g. 1_000_000) to chunk large arrays along their first
    axis, which is what makes lazy streaming read only the touched blocks. Both chunked and
    gzip-compressed stores are readable by the C++ core.

    stream=True writes with bounded memory: any sparse field whose values are a streaming source
    (a backed AnnData, or a `LazyCSX` from `read(..., lazy=True)`) is copied block-by-block instead
    of being materialized, so a large `h5ad`->L* or L*->L* conversion never holds the whole matrix.
    (Such sources are streamed even without `stream=True`; the flag also chunks the rest of the
    store, since streaming output is inherently chunked.)

    viewer=True first calls `lstar.extend_for_viewer(ds)` to add the lstar-viewer precomputed fields
    (counts_cellmajor, per-group stats + marker tables, od_score, and a hybrid cell order) so the
    resulting store is ready for fast differential-expression / variable-gene / dotplot browsing.
    """
    if viewer:
        from .viewer import extend_for_viewer
        extend_for_viewer(ds)
    if stream and chunk_elems is None:
        chunk_elems = 1_000_000
    is_zip = str(path).endswith(".zip")
    _zip_target = None
    if is_zip:                                         # single-file .lstar.zarr.zip: write a dir, pack STORED
        import tempfile
        _zip_target = str(path)
        path = tempfile.mkdtemp(suffix=".lstar.zarr")
    root = zarr.open_group(path, mode="w")
    axg = root.create_group("axes")
    flg = root.create_group("fields")
    root.create_group("models")

    for name, ax in ds.axes.items():
        g = axg.create_group(name)
        _write_strings(g, "labels", ax.labels, compressor, chunk_elems)
        g.attrs[LSTAR] = {"kind": "axis", "origin": ax.origin, "role": ax.role,
                          "induced_by": ax.induced_by, "provenance": ax.provenance}

    for name, fl in ds.fields.items():
        g = flg.create_group(name)
        meta = {"kind": "field", "role": fl.role, "span": fl.span, "state": fl.state,
                "encoding": fl.encoding, "coverage": fl.coverage, "directed": fl.directed,
                "weighted": fl.weighted, "subtype": fl.subtype, "uncertainty": fl.uncertainty,
                "provenance": fl.provenance}
        _write_values(g, fl, meta, compressor, chunk_elems)
        if fl.mask is not None:                        # nullable Int/bool/string: an explicit validity mask
            _ds(g, "mask", np.asarray(fl.mask, dtype=np.uint8), compressor, chunk_elems)
            meta["nullable"] = True
        if getattr(fl, "index", None) is not None:     # partial coverage: int positions into index_axis
            _ds(g, "index", np.asarray(fl.index, dtype=np.int64), compressor, chunk_elems)
            meta["coverage"] = "partial"
            meta["index_axis"] = fl.index_axis
        g.attrs[LSTAR] = meta

    aux = getattr(ds, "aux", None) or {}
    if aux:                                            # verbatim passthrough (uns/@misc) -> passthrough/<ns>
        from .passthrough import to_store as _aux_to_store
        auxg = root.create_group("passthrough")
        for ns, obj in aux.items():
            g = auxg.create_group(ns)
            tree, leaves = _aux_to_store(obj)
            manifest = []
            for a in leaves:
                if a["kind"] == "utf8":
                    _write_strings(g, a["id"], np.asarray(a["data"], dtype=str), compressor, chunk_elems)
                else:
                    _ds(g, a["id"], np.asarray(a["data"]), compressor, chunk_elems)
                manifest.append({"id": a["id"], "kind": a["kind"]})
            # tree is stored as an opaque JSON *string*: zarr sorts attribute object keys, which would
            # scramble the passthrough's dict order; a string is preserved verbatim (and lets the
            # C++/R readers round-trip it without parsing JSON).
            g.attrs[LSTAR] = {"kind": "passthrough", "tree": json.dumps(tree), "arrays": manifest}

    root.attrs[LSTAR] = {"spec_version": ds.spec_version or SPEC_VERSION, "kind": ds.kind,
                         "profiles": list(ds.profiles), "dropped": list(ds.dropped),
                         "axes": list(ds.axes), "fields": list(ds.fields), "passthrough": list(aux)}
    zarr.consolidate_metadata(path)
    if is_zip:                                         # pack the finished store into ONE file, every entry STORED
        import shutil
        _pack_stored_zip(path, _zip_target)
        shutil.rmtree(path, ignore_errors=True)
        return _zip_target
    return path


def _pack_stored_zip(src_dir, zippath):
    """Pack a directory store into `zippath` as ONE file with every entry STORED (no deflate).

    A `.lstar.zarr.zip` is always STORED: zarr chunks are already codec-compressed, so re-deflating
    them wastes CPU for ~no gain, and — the load-bearing reason — only a STORED entry stays
    byte-range-readable inside the archive, which is the whole point of a hosted single file (a
    reader issues one HTTP Range into the zip for a chunk; a deflated entry would force fetching +
    inflating the whole entry). Going via a directory (rather than zarr's ZipStore) keeps the write
    path identical to a normal store — the streaming writer's array resizes can't leave duplicate
    zip entries, and a `compressor=` argument only ever compresses the inner chunks, never the zip."""
    import zipfile
    src_dir = str(src_dir)
    entries = []
    for root, _dirs, files in os.walk(src_dir):
        for fn in files:
            fp = os.path.join(root, fn)
            arc = os.path.relpath(fp, src_dir).replace(os.sep, "/")
            entries.append((arc, fp))
    # deterministic archive; metadata (.z*) first so a reader hits the manifest early
    entries.sort(key=lambda e: (not os.path.basename(e[0]).startswith(".z"), e[0]))
    with zipfile.ZipFile(zippath, "w", compression=zipfile.ZIP_STORED, allowZip64=True) as zf:
        for arc, fp in entries:
            zf.write(fp, arcname=arc)


def _open_root(path):
    """Open the store root, preferring consolidated metadata (one read, no listing)."""
    if str(path).endswith(".zip"):
        return _open_zip_root(path)
    if os.path.exists(os.path.join(str(path), ".zmetadata")):
        try:
            return zarr.open_consolidated(path, mode="r")
        except Exception:
            pass
    return zarr.open_group(path, mode="r")


def _open_zip_root(path):
    """Open a single-file `.lstar.zarr.zip`. Reject a DEFLATE-packed archive with a clear message: a
    hosted single-file store must be STORED so its chunks stay byte-range-readable (a deflated entry
    silently defeats the range access that is the point of the single file). The open ZipStore is held
    alive by the returned group's arrays, so lazy reads stream chunks straight from the zip."""
    import zipfile
    with zipfile.ZipFile(str(path)) as z:
        bad = [i.filename for i in z.infolist() if i.compress_type != zipfile.ZIP_STORED]
    if bad:
        raise ValueError(
            f"{path}: this .lstar.zarr.zip is DEFLATE-compressed ({len(bad)} of its entries, "
            f"e.g. {bad[0]!r}) — a hosted single-file store must be written STORED so its chunks "
            f"stay byte-range-readable. Repack it STORED (via `lstar convert`, or `zip -0 -r`).")
    store = zarr.ZipStore(str(path), mode="r")
    try:
        return zarr.open_consolidated(store=store, mode="r")
    except Exception:
        return zarr.open_group(store=store, mode="r")


def read(path, lazy=False):
    """Read an L* Dataset from a Zarr store at `path`.

    lazy=False (default) materializes every field. lazy=True leaves field values as proxies
    (`lstar.lazy.LazyDense` / `LazyCSX`) backed by the open zarr arrays: the store opens without
    reading the heavy arrays, and a CSC measure can be reduced by streaming column blocks
    (`lstar.lazy.stream_col_stats`) without ever fully materializing.
    """
    root = _open_root(path)
    rmeta = dict(root.attrs[LSTAR])
    ds = Dataset(kind=rmeta.get("kind", "sample"),
                 spec_version=rmeta.get("spec_version", SPEC_VERSION))
    ds.profiles = list(rmeta.get("profiles", []))
    ds.dropped = list(rmeta.get("dropped", []))

    for name in rmeta["axes"]:
        g = root["axes"][name]
        m = dict(g.attrs[LSTAR])
        ds.add_axis(name, _read_strings(g, "labels"), origin=m.get("origin", "observed"),
                    role=m.get("role"), induced_by=m.get("induced_by"),
                    provenance=m.get("provenance", {}))

    for name in rmeta["fields"]:
        g = root["fields"][name]
        m = dict(g.attrs[LSTAR])
        vals = _lazy_values(g, m) if lazy else _read_values(g, m)
        mask = np.asarray(g["mask"], dtype=np.uint8) if m.get("nullable") and "mask" in g else None
        index = np.asarray(g["index"], dtype=np.int64) if "index" in g else None
        ds.fields[name] = Field(
            name, vals, role=m.get("role"), span=m.get("span"), state=m.get("state"),
            encoding=m.get("encoding"), coverage=m.get("coverage", "full"),
            directed=m.get("directed"), weighted=m.get("weighted"),
            subtype=m.get("subtype"), uncertainty=m.get("uncertainty"),
            mask=mask, index=index, index_axis=m.get("index_axis"), provenance=m.get("provenance", {}))

    for ns in rmeta.get("passthrough", []):            # verbatim passthrough -> reconstruct the object
        from .passthrough import from_store as _aux_from_store
        g = root["passthrough"][ns]
        am = dict(g.attrs[LSTAR])
        leaves = []
        for a in am.get("arrays", []):
            data = _read_strings(g, a["id"]) if a["kind"] == "utf8" else np.asarray(g[a["id"]])
            leaves.append({"id": a["id"], "kind": a["kind"], "data": data})
        tree = am.get("tree")
        ds.aux[ns] = _aux_from_store(json.loads(tree) if isinstance(tree, str) else tree, leaves)
    return ds


# ---- field value encodings ----

def _chunks_for(shape, chunk_elems):
    """A chunk shape splitting the first axis so each chunk holds ~chunk_elems elements.

    Returns None (single chunk) when chunk_elems is None or the array already fits.
    """
    if chunk_elems is None or len(shape) == 0 or shape[0] == 0:
        return None
    inner = 1
    for s in shape[1:]:
        inner *= int(s)
    rows = max(1, chunk_elems // max(1, inner))
    if rows >= shape[0]:
        return None
    return (rows,) + tuple(int(s) for s in shape[1:])


def _ds(g, name, arr, compressor, chunk_elems):
    arr = np.asarray(arr)
    g.create_dataset(name, data=arr, compressor=compressor,
                     chunks=_chunks_for(arr.shape, chunk_elems))


def _write_sparse_streaming(g, source, meta, compressor, chunk_elems):
    """Write a CSR/CSC field block-by-block from a streaming source (a backed h5ad sparse group, or
    a LazyCSX from another store), so the whole matrix is never resident. `data`/`indices` grow by
    `append`; `indptr` (small) is filled incrementally. Outer blocks are sized to ~chunk_elems
    nonzeros each (using the source's small full `indptr`), so peak memory is bounded regardless of
    how dense the data is. Streaming implies a chunked output (the final size isn't known up front)."""
    fmt = source.fmt                                   # "csr" | "csc"
    n_outer = int(source.n_outer)
    nnz = int(getattr(source, "nnz", -1))
    data_dtype = np.dtype(getattr(source, "dtype", np.float64))
    idx_dtype = np.dtype(getattr(source, "idtype", np.int32))
    indptr_dtype = np.int32 if 0 <= nnz < 2 ** 31 else np.int64
    ce = int(chunk_elems) if chunk_elems else 1_000_000
    data_arr = g.create_dataset("data", shape=(0,), chunks=(ce,), dtype=data_dtype, compressor=compressor)
    indices_arr = g.create_dataset("indices", shape=(0,), chunks=(ce,), dtype=idx_dtype, compressor=compressor)
    out_indptr = np.empty(n_outer + 1, dtype=indptr_dtype)
    out_indptr[0] = 0
    pos = 0
    for a, b, sub in iter_sized_blocks(source, ce):    # blocks sized to ~ce nonzeros (bounded memory)
        sub = sub.tocsr() if fmt == "csr" else sub.tocsc()
        data_arr.append(np.asarray(sub.data, dtype=data_dtype))
        indices_arr.append(np.asarray(sub.indices, dtype=idx_dtype))
        out_indptr[a + 1:b + 1] = pos + sub.indptr[1:]
        pos += int(sub.nnz)
    g.create_dataset("indptr", data=out_indptr, compressor=compressor,
                     chunks=_chunks_for(out_indptr.shape, chunk_elems))
    meta["encoding"] = fmt
    meta["shape"] = [int(x) for x in source.shape]


def _write_values(g, fl, meta, compressor, chunk_elems=None):
    enc = fl.encoding
    if _is_stream_source(fl.values):                   # backed/lazy sparse -> stream block-by-block
        _write_sparse_streaming(g, fl.values, meta, compressor, chunk_elems)
    elif enc in ("csr", "csc"):
        m = fl.values.tocsr() if enc == "csr" else fl.values.tocsc()
        _ds(g, "data", m.data, compressor, chunk_elems)
        _ds(g, "indices", m.indices, compressor, chunk_elems)
        _ds(g, "indptr", m.indptr, compressor, chunk_elems)
        meta["shape"] = [int(x) for x in m.shape]
    elif enc == "coo":                                 # coo is a Python-only on-disk form (C++/R/JS have no
        m = fl.values.tocsc()                          # coo reader) -> normalize to csc so every surface reads it
        _ds(g, "data", m.data, compressor, chunk_elems)
        _ds(g, "indices", m.indices, compressor, chunk_elems)
        _ds(g, "indptr", m.indptr, compressor, chunk_elems)
        meta["encoding"] = "csc"
        meta["shape"] = [int(x) for x in m.shape]
    elif enc == "categorical" or _is_categorical(fl.values):
        cat = as_categorical(fl.values)                # codes (-1 = missing) + categories + ordered
        _ds(g, "codes", cat.codes.astype(np.int32), compressor, chunk_elems)
        if meta.get("categories") is None:             # inline categories (P1); an axis ref (P2) skips this
            _write_strings(g, "categories", cat.categories, compressor, chunk_elems)
        meta["encoding"] = "categorical"
        meta["ordered"] = bool(cat.ordered)
        meta["shape"] = [int(len(cat))]
    else:
        arr = np.asarray(fl.values)
        if arr.dtype.kind in ("U", "S", "O"):          # string field -> utf8 + offsets
            _write_strings(g, "values", arr, compressor, chunk_elems)
            meta["encoding"] = "utf8"
            meta["shape"] = [int(arr.shape[0])]
        else:
            _ds(g, "values", arr, compressor, chunk_elems)
            meta["encoding"] = "dense"
            meta["shape"] = [int(x) for x in arr.shape]   # dense: shape in the manifest too (parity with sparse)


def _read_values(g, m):
    enc = m.get("encoding")
    if enc in ("csr", "csc"):
        cls = sp.csr_matrix if enc == "csr" else sp.csc_matrix
        return cls((np.asarray(g["data"]), np.asarray(g["indices"]),
                    np.asarray(g["indptr"])), shape=tuple(m["shape"]))
    if enc == "coo":
        return sp.coo_matrix((np.asarray(g["weight"]),
                              (np.asarray(g["row"]), np.asarray(g["col"]))),
                             shape=tuple(m["shape"]))
    if enc == "utf8":
        return _read_strings(g, "values")
    if enc == "categorical":
        cats = _read_strings(g, "categories") if "categories" in g else m.get("_categories")
        return Categorical(np.asarray(g["codes"]), cats, ordered=bool(m.get("ordered", False)))
    return np.asarray(g["values"])


def _lazy_values(g, m):
    """Field value as a lazy proxy (sparse/dense) or, for small string fields, eager."""
    from .lazy import LazyDense, LazyCSX
    enc = m.get("encoding")
    if enc in ("csr", "csc"):
        return LazyCSX(enc, g["data"], g["indices"], g["indptr"], tuple(m["shape"]))
    if enc == "utf8":                       # labels are small; materialize
        return _read_strings(g, "values")
    if enc in ("coo", "categorical"):       # rare / small; materialize
        return _read_values(g, m)
    return LazyDense(g["values"])


# ---- string encoding (utf8 bytes + offsets) ----

def _write_strings(g, name, values, compressor, chunk_elems=None):
    arr = np.asarray(values)
    bs = [str(x).encode("utf-8") for x in arr.tolist()]
    offs = np.zeros(len(bs) + 1, dtype=np.int64)
    for i, b in enumerate(bs):
        offs[i + 1] = offs[i] + len(b)
    data = (np.frombuffer(b"".join(bs), dtype=np.uint8).copy()
            if bs else np.zeros(0, dtype=np.uint8))
    _ds(g, name, data, compressor, chunk_elems)
    _ds(g, name + "_offsets", offs, compressor, chunk_elems)


def _read_strings(g, name):
    buf = np.asarray(g[name], dtype=np.uint8).tobytes()
    offs = np.asarray(g[name + "_offsets"])
    return np.array([buf[int(offs[i]):int(offs[i + 1])].decode("utf-8")
                     for i in range(len(offs) - 1)], dtype=str)
