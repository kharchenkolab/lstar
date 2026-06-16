# L★ model & store layout (reference)

## Axes

An **Axis** is a labelled set you index by. Attributes:

| field | meaning |
|---|---|
| `name` | unique key (`cells`, `genes`, `pca`, `samples`, `cells.GSM123`) |
| `labels` | the ordered labels (strings) — the identity of each position |
| `origin` | `observed` (measured) or `derived` (computed, e.g. a PCA coordinate axis) |
| `role` | optional: `observation \| feature \| coordinate \| sample \| ...` |
| `induced_by` | optional: the field/operation that induced a derived axis |
| `provenance` | optional dict (native location, parameters) |

Per-sample axes are named with a suffix convention `cells.<sample>` / `genes.<sample>`.

## Fields

A **Field** is typed data over a tuple of axes (its `span`).

| field | meaning |
|---|---|
| `name` | unique key |
| `values` | the data: dense `ndarray`, sparse `Matrix`, or strings |
| `role` | `measure \| embedding \| loading \| relation \| label \| sequence \| design \| transform` |
| `span` | ordered list of axis names, e.g. `["cells","genes"]`, `["cells","cells"]` (relation) |
| `state` | optional measure state: `raw \| lognorm \| scaled` |
| `encoding` | `dense \| csr \| csc \| coo \| utf8` (inferred if omitted) |
| `coverage`/`index`/`index_axis` | `full` (default) or `partial` — covers a subset of one span axis (`index_axis`) via an `index` of int positions into it (implemented; round-trips Py/C++/R) |
| `subtype` | role-specific tag (`distance`, `similarity`, `knn`, `design`, …) |
| `directed`,`weighted` | relation flags |
| `provenance` | native location for exact write-back (e.g. `{"anndata":"obsm/X_pca"}`) |

### Roles — what each is for

- **measure** — values over `(cells, genes)`: counts, normalized expression. `state` distinguishes raw/lognorm/scaled.
- **embedding** — coordinates over `(cells, <coord-axis>)`: PCA/UMAP/tSNE.
- **loading** — feature weights over `(genes, <coord-axis>)`: PCA loadings.
- **relation** — a graph over `(axis, axis)`, usually `(cells, cells)`: kNN/SNN, distances, connectivities. `directed`/`weighted`/`subtype` qualify it.
- **label** — categorical annotation over one axis: clusters, cell types. Encoding `utf8`.
- **design** — a `label`/covariate used as experimental design (e.g. the per-cell `sample`).
- **sequence**, **transform** — ordered data / a model (a `transform` field *is* a model, e.g. a PCA rotation).

### Role/encoding inference

If you pass only `values`, lstar infers: `span` by matching each dimension to a unique axis length;
`role` (`label` for string 1-D, else `measure`); `encoding` from the sparse format or `dense`. Pass
explicit `span=`/`role=` when ambiguous (e.g. a square `(cells, cells)` relation).

## Kinds: sample vs collection

- `kind="sample"` — one uniform measurement (the pagoda2 unit): shared `cells`/`genes`, aligned facets.
- `kind="collection"` — a **collection of heterogeneous samples** (the conos unit). Canonical shape:
  - a `samples` axis (the members);
  - per-sample `cells.<s>` and `genes.<s>` axes + `counts.<s>` measures (samples may differ in cells *and* genes);
  - a union `cells` axis for the joint layer;
  - a `sample` **design** label over union cells;
  - the joint `embedding`, joint cluster `label`(s), and the integration **graph** as a `relation` over `(cells, cells)`.

This is the core differentiator: alignment is legitimate only *within* a sample; across samples you
keep a collection joined by a graph, not a concatenated matrix.

**Building one:** `collection_from(samples, joint=...)` (Python `lstar.collection_from`, R `collection_from`)
assembles this canonical shape from any **list/dict of per-sample objects** — `Dataset`/`AnnData`/`MuData`
(Py) or `lstar_dataset`/`Seurat`/`SingleCellExperiment` (R) — namespacing each sample's axes/fields as
`<x>.<s>`, building the `samples` + union `cells` axes and the `sample` design label. `joint=` fields land
over the union cells: a 2-D array → embedding, a `(cells×cells)` sparse matrix → graph `relation`, a
factor/categorical → clustering label. `write_conos` and a split Seurat v5 assay produce the *same* shape,
so every collection — however assembled — has one structure. (`conformance/collection_true.sh`,
`conformance/conos.sh`.)

## Store layout (Zarr v2)

```
store.lstar.zarr/
  .zattrs            -> {"lstar": {spec_version, kind, profiles, dropped, axes:[...], fields:[...]}}
  .zmetadata         -> consolidated metadata (one read)
  axes/<name>/       -> .zattrs {"lstar":{kind:"axis", origin, role, provenance}}, labels(+_offsets)
  fields/<name>/     -> .zattrs {"lstar":{kind:"field", role, span, state, encoding, shape, ...}}
                        dense:  values
                        csc/csr: data, indices, indptr
                        coo:    row, col, weight
                        utf8:   values(uint8) + values_offsets(int64)
  models/            -> (transform fields / models)
```

- **Strings** are encoded as UTF-8 bytes (`uint8`) + `int64` offsets (length n+1), not fixed-width
  unicode — so C++ reads them without an encoding dependency.
- `provenance` records the exact native location so a profile can write back losslessly.
- `dropped` lists native locations a profile could not represent (loss is recorded, never silent).

## Lossless round-trip & fixed point

A value can leave its native format, pass through L★ any number of hops (even across formats and
languages), and return unchanged: native → L★ → native is a **fixed point**. What L★ cannot hold for
a given target (e.g. `uns` for AnnData, graphs through Seurat) is recorded in `dropped`.

See `docs/model.md` and `docs/format.md` for the normative spec, `reference/profiles.md` for the
native↔L★ mappings.
