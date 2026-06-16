# Python API (reference)

Package `lstar` (under `python/src/lstar`). Pure-Python on numpy/scipy/zarr, with an optional
compiled C++ accelerator that is used automatically when present.

## Construction

```python
import lstar
ds = lstar.Dataset(kind="sample")            # or kind="collection"
ds.add_axis(name, labels, origin="observed", role=None, induced_by=None, provenance=None)
ds.add_field(name, values, role=None, span=None, state=None, encoding=None,
             coverage="full", directed=None, weighted=None, subtype=None, provenance=None)
ds.axis(name)        # -> Axis (has .name, .labels, .origin, .role; len(axis) == #labels)
ds.field(name)       # -> Field (.values, .role, .span, .state, .encoding, .provenance, ...)
ds.fields_over(axis_name)   # fields whose span includes an axis
ds.profiles, ds.dropped     # provenance / recorded losses (lists)
```
`add_field` infers `span`/`role`/`encoding` from `values` if omitted (see `reference/model.md`).

## Serialization

```python
lstar.write(ds, path, compressor=None, chunk_elems=None)
#   compressor: a numcodecs codec, e.g. numcodecs.GZip(5); None = uncompressed
#   chunk_elems: chunk arrays along axis 0 (~that many elements/chunk); None = single chunk.
#                Chunking is what lets a lazy read touch only the blocks it needs.
ds = lstar.read(path, lazy=False)
#   lazy=True -> field .values are proxies (lstar.LazyDense / lstar.LazyCSX), nothing materialized;
#                prefers zarr.open_consolidated when .zmetadata exists.
errs = lstar.validate(ds)                    # [] if valid; checks spans/shapes/relation arity
```

## Lazy / streaming

```python
from lstar.lazy import LazyDense, LazyCSX
d = lstar.read("big.lstar.zarr", lazy=True)
v = d.field("counts").values                 # LazyCSX(csc, shape=(C,G), nnz=...)
v.shape, v.nnz, v.dtype
v.materialize()                              # full scipy matrix
np.asarray(v)                                # dense (careful: may be large)
for start, stop, sub in v.blocks(block=2048): ...   # stream column blocks (CSC)
v.outer_block(a, b)                          # one block as a small scipy matrix
```

### `stream_col_stats` — zero-aware per-column mean/var, streamed

```python
mean, var, nnz = lstar.stream_col_stats(field_value, lognorm=False, block=2048,
                                        n_threads=1, engine="auto")
```
- Accepts a `LazyCSX` (streams from disk, bounded memory) or a materialized scipy CSC.
- `lognorm=True` applies `log1p` per nonzero on the fly (the normalized matrix is never built).
- Zero-aware: implicit zeros are accounted for; values stay in their stored dtype, moments
  accumulate in float64 (lean + accurate).
- `n_threads`: `1` serial (default), `N` threads, `0`/`None` all cores.
- `engine`: `"auto"` (C++ if available else Python), `"c++"`, `"python"`. Results identical.

## The compute engine (fast by default)

```python
lstar.has_accel()     # True if the compiled C++ accelerator imported
lstar.show_config()   # prints which path is active + OpenMP threads
# env override (debug/bench only): LSTAR_ENGINE=python|c++
```
The accelerator (`lstar._accel`, a pybind11 binding over libstar) is autodetected at import. A wheel
install gives it automatically; absent it, everything still works on the pure-Python path. The C++
reduction scales to the core count (≈20× at 56 threads on a 77M-nonzero measure); the pure-Python
streamer peaks ≈2× around 4 threads (GIL on the orchestration) — for heavy reductions the C++ engine
is the path.

## AnnData profile

```python
from lstar.profiles.anndata import read_anndata, write_anndata
ds = read_anndata(adata)        # X, layers, obs/var, obsm(embeddings), varm(loadings),
                                # obsp/varp(relations), .raw (own genes_raw axis if divergent)
adata = write_anndata(ds)       # writes back to the recorded native locations; ds.dropped -> adata.uns["lstar/dropped"]
```
- Records `anndata@<version>` in `ds.profiles`; uns is dropped (recorded in `ds.dropped`).
- Field names follow the shared vocabulary (`X_pca` → `pca`, `PCs` → `pca` loadings) so cross-format
  conversion is meaningful.
- A `kind="collection"` dataset is **flattened** to one union-genes AnnData (X = raw joint counts;
  embedding → `obsm`, graph → `obsp` aliased to `connectivities` + `uns['neighbors']`, clustering/sample
  → `obs`; per-sample PCA → `dropped`). An AnnData is one matrix, so this loses the per-sample structure —
  keep the L\* store to retain it.

## Collections

```python
col = lstar.collection_from(samples, joint=None, sample_field="sample", prefix_cells=True)
#   samples: dict {name: Dataset|AnnData|MuData} or a list (auto-named s0, s1, ...). Each keeps its own
#            cells.<name>/genes.<name> axes + <field>.<name> measures; gene sets may overlap, differ, or
#            be entirely disjoint. Builds the `samples` + union `cells` axes and the `sample` design label.
#   joint:   dict over the UNION cells -- 2-D array -> embedding; (n×n) sparse -> graph relation;
#            Categorical/label -> clustering (induces a factor axis).
pb = lstar.collection_pseudobulk(col, "clusters", field="counts")   # streamed pseudobulk over the union
```
Assembles the same canonical shape `write_conos` / a split Seurat v5 assay produce, from any list of
per-sample objects. See `reference/model.md` and `conformance/collection_true.sh`.

## Packaging (PyPI)

- `pyproject.toml` ([project] metadata + build-system + `[tool.cibuildwheel]`), `setup.py` builds the
  **optional** extension (install never fails for lack of a compiler → pure-Python fallback),
  `MANIFEST.in` grafts `core/include` into the sdist.
- `.github/workflows/wheels.yml` runs `cibuildwheel` (manylinux/macOS/Windows) + sdist, trusted-publish
  to PyPI on release. macOS OpenMP via `brew install libomp`.
- Build locally: `pip install .` (build isolation). The extension needs Python dev headers + a C++17
  compiler with OpenMP; `lstar.show_config()` confirms it loaded.

## Files

`model.py` (Dataset/Axis/Field), `zarr_io.py` (read/write + chunk/codec), `lazy.py` (proxies +
`stream_col_stats`), `_engine.py` (autodetect/dispatch), `_accel.cpp` (pybind11 binding),
`validate.py`, `profiles/anndata.py`. Tests: `python/tests/test_*.py`.
