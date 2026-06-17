# The on-disk format

An L★ dataset is serialized as a [Zarr](https://zarr.dev) group tree, named `*.lstar.zarr` by
convention. This page is the practical reading of Appendix A of the proposal
([`../misc/Lstar_proposal.md`](../misc/Lstar_proposal.md)); consult that for the normative schema. The
guiding rule: **all L★ metadata lives under an `"lstar"` key** in each group's `.zattrs`, so the store
is *also* a plain Zarr group that scanpy/Vitessce/zarrita can open and read the common parts of.

> **Format stability.** The on-disk format is at `spec_version` **0.1** and **may change in
> backward-incompatible ways before 1.0**. Every store records its `spec_version` and readers check it,
> so the format is self-describing; reaching **1.0** is what turns that version line into a compatibility
> commitment. Until then, regenerate stores from source rather than relying on them for long-term archival.

> **Implemented today.** The Python and C++ writers emit the tree below in **Zarr v2** with a
> consolidated `.zmetadata`; the C++/R/Python readers all read it (chunked + gzip included). A nullable
> validity `mask` and **partial coverage** (`coverage="partial"` + an `index` into `index_axis`) are
> written and read by all three. The spec targets **Zarr v3 + sharding** (the request-economy path for
> million-cell remote stores); that and the `recipe`/`ragged`/`raster` encodings and `models/` are not
> yet written.

## The tree

```text
store.lstar.zarr/
├── .zgroup                     {"zarr_format": 2}
├── .zattrs                     {"lstar": <root metadata>}            # see below
├── .zmetadata                  consolidated metadata (one read; required for HTTP access)
├── axes/
│   ├── .zgroup
│   └── <axis>/
│       ├── .zattrs             {"lstar": {kind:"axis", origin, role, induced_by, provenance}}
│       ├── labels/             the element labels (utf8: a uint8 byte array …)
│       └── labels_offsets/     … + an int64 offsets array, length n+1
├── fields/
│   ├── .zgroup
│   └── <field>/
│       ├── .zattrs             {"lstar": {kind:"field", role, span, state, encoding, shape, …}}
│       └── <value arrays>      depend on the encoding (below)
├── aux/                        lossless passthrough of a format's untyped long-tail
│   ├── .zgroup
│   └── <namespace>/            e.g. "anndata.uns"
│       ├── .zattrs             {"lstar": {kind:"aux", tree:"<JSON string>", arrays:[{id, kind}]}}
│       └── <id>               the array leaves the tree references (dense, or utf8 bytes+offsets)
└── models/                     fitted transforms (apply contract + weights)   [spec; not yet emitted]
```

`axes` and `fields` in the root metadata give the read order; `span` references axes *by name*, so a
`cells × cells` relation and a `cells × genes` measure are unambiguous even when the two counts are
equal — placement and reading are deterministic lookups, never shape inference.

## Root metadata (`.zattrs` → `lstar`)

```json
{
  "spec_version": "0.1",
  "kind": "sample",                          // or "collection"
  "profiles": ["anndata@0.1", "anndata@0.8.0"],   // who wrote it, with versions
  "dropped": ["uns/velocity_graph"],         // native locations no profile rule could hold (recorded loss)
  "axes":   ["cells", "genes", "pca", "umap", "leiden"],
  "fields": ["counts", "data", "pca", "umap", "knn", "leiden"]
}
```

`profiles` lets a reader interpret idiosyncrasies and know which conformance suite applied. `dropped`
is the loss manifest — **what a writer could not represent is recorded, never silently lost** (see the
worked conversion in the proposal §4.3).

## Lossless passthrough (`aux/`)

A format's untyped long-tail — AnnData `uns`, Seurat `@misc`/`@commands` (params, color palettes,
dendrograms, DE tables) — is preserved **verbatim** under `aux/<namespace>/` rather than recorded
name-only in `dropped`. Each namespace is a *self-describing* subtree: a JSON `tree` (stored as an
opaque **string** so the store's key order survives zarr's key-sorting) whose array leaves are
references into a flat `arrays` manifest; the manifest's `dense`/`utf8` leaves are ordinary L* arrays.
Leaf grammar: a JSON scalar, `{"$obj":…}`, `{"$list":…}`, `{"$array":id}`, `{"$strings":id}`,
`{"$record":…}` (a numpy structured array, e.g. `rank_genes_groups`), or `{"$dropped":…}` for a
genuinely unrepresentable leaf (recorded, never silent). The C++ and R cores **round-trip the subtree
without interpreting it** (they carry the `tree` string + the array leaves); only the originating
profile walks the tree to rebuild the live object. This keeps the tail *inspectable*, so a recognized
structure can be promoted to a typed field later.

## Field value arrays, by encoding

The `lstar.encoding` attribute selects which value arrays a field group holds:

| encoding | arrays under `fields/<name>/` | for |
|---|---|---|
| `dense` | `values` (C-order, N-D matching `span`) | dense matrices/vectors/embeddings |
| `csc` | `data`, `indices`, `indptr` (+ `lstar.shape`) | sparse, **gene-compressed** — color any gene cheaply |
| `csr` | `data`, `indices`, `indptr` (+ `lstar.shape`) | sparse, **cell-compressed** — fetch a cell's genes cheaply |
| `coo` / `edge_list` | `row`/`source`, `col`/`target`, `weight` | unordered sparse triples; graphs |
| `utf8` | `values` (uint8), `values_offsets` (int64) | string / `label` fields with no category set |
| `categorical` | `codes` (int, `-1`=missing) + inline `categories` (utf8); `lstar.ordered` | factor/categorical `label` fields — preserves category order + missingness; the substrate for **factor axes** (a label's categories *are* a derived axis) and dictionary-encoded labels |
| `ragged` | `values`, `offsets` | variable-length `sequence` data *(spec; not yet written)* |
| `raster` | an OME-NGFF multiscale group | images over `(y, x, channel)` *(spec; not yet written)* |
| `recipe` | *no value arrays* — `lstar.recipe = {op, inputs:[…], params}` | a **virtual** field computed on read *(spec; not yet written)* |

**Strings (`utf8`).** Labels and string-valued fields are stored as concatenated UTF-8 *bytes*
(`uint8`) plus an `int64` *offsets* array of length n+1: item *i* is `bytes[offsets[i] :
offsets[i+1]]`. This avoids fixed-width unicode arrays and is trivially decodable from C++ and
JavaScript (`TextDecoder`). The same encoding holds axis `labels`.

**Nullability (validity mask).** Any field may carry an optional `mask` array (`uint8`, same length as
its values, `1 = missing`) plus `lstar.nullable = true`. This is how pandas **nullable** `Int64` /
`boolean` / `string` columns keep their integer-ness and their value-vs-missing distinction instead of
collapsing to float-NaN. It is distinct from the categorical `-1` sentinel (which is built into `codes`)
and from float `NaN` (which already encodes missing, so a float field needs no mask). An absent `mask`
means no nulls — existing stores read unchanged.

**Partial coverage** *(implemented).* A field that covers only some of one span axis sets
`coverage="partial"` and `index_axis=<axis>` in its `lstar` attrs and adds an `index` array under
`fields/<name>/index` — `int64` **positions** into `axes/<axis>/labels` (one per value row along that
axis), so the covered subset is `axes/<axis>/labels[index]`. The field's stored shape is `len(index)`
along that axis (not the full axis length); uncovered elements are *absent* (not zero/NA-padded). Used
e.g. for a modality measured on a subset of cells, or a Seurat `scale.data` over the variable features
only. Round-trips across the Python, C++, and R readers/writers.

**Models** *(spec; not yet emitted).* `models/<model>/` holds a fitted `transform`: an `apply` contract
(`{in:[axes], out:[axes]}`) plus learned `weights` (inline or a URI), with provenance.

## dtypes

- values: `<f4` (float32 — common, e.g. AnnData/Seurat) or `<f8`; integer measures at their stored
  width. Readers keep values in their stored dtype and accumulate in float64 (lean + accurate).
- sparse `indices`/`indptr`: `<i4` or `<i8` by size — scipy emits int32 until a dimension or nnz
  exceeds 2³¹. Readers normalize to int64 for computation (`as_i64` in the C++ core).
- `labels`/`utf8` data: `|u1`; offsets: `<i8`. Little-endian, C-order.

## Chunking & compression

- Arrays may be **single-chunk** (the portable default) or **chunked** along the first axis (set
  `chunk_elems` on write). Chunking is what lets a lazy reader fetch only the blocks a query touches —
  e.g. one gene's CSC column.
- Chunks may be **uncompressed** (default) or **gzip/zlib** compressed (`compressor=` on write). The
  C++ core reads any chunk grid (C-order, fill-padded edge chunks; a missing chunk reads as
  `fill_value` 0) and decodes gzip/zlib when built with zlib; `blosc` and Zarr v3 sharding are planned.

## Consolidated metadata

`.zmetadata` collects every `.zgroup`/`.zattrs`/`.zarray` into one JSON document
(`{"zarr_consolidated_format": 1, "metadata": {…}}`), so a reader makes **one request** instead of
many small ones — required for access over HTTP. Both the Python and C++ writers emit it; readers
prefer it (`zarr.open_consolidated`, zarrita's consolidated path) and fall back to walking the tree.

## Packaging — one model, four layouts

Because a collection references its members by identity (axes are namespaced), packaging is decoupled
from logical structure (proposal §3.2):

1. **Standalone sample** — one self-contained store, valid on its own: `human_1.lstar.zarr/`.
2. **Collection, bundled** — one store physically containing its members (portable; one `.zip`):
   `samples/<id>/` sub-stores + collection-level `axes/ fields/ models/`.
3. **Collection, by reference** — a thin overlay; members live elsewhere (local or `s3://`), pinned by
   `id` + content `hash` in `lstar.members`; loads lazily over HTTP.
4. **Mixed** — some members embedded, some referenced.

A sample store always opens on its own; the same sample may be referenced by several collections
without duplication; bundling and splitting are lossless. *(lstar today writes standalone stores;
collections are written as a single store with per-sample axes/fields — the `samples/<id>/` sub-store
and `members` reference layouts are specified but not yet emitted.)*

## Cross-implementation guarantee

The same store reads byte-faithfully from **Python, R, C++, and the browser (WASM/zarrita)**. A store
written by any one reads in the others with identical field values (counts sums, graph nnz,
embeddings); chunked + gzip stores round-trip across all of them. This is enforced by the conformance
suite ([`../conformance/`](../conformance)).
