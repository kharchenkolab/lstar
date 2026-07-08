# The on-disk format

An L★ dataset is serialized as a [Zarr](https://zarr.dev) group tree, named `*.lstar.zarr` by
convention. This page is the practical reading of Appendix A of the proposal
([`../misc/Lstar_proposal.md`](../misc/Lstar_proposal.md)); consult that for the normative schema. The
guiding rule: **all L★ metadata lives under an `"lstar"` key** in each node's attributes (the v3
`zarr.json` `attributes`, or the v2 `.zattrs`), so the store is *also* a plain Zarr group that any
standard Zarr reader (scanpy, Vitessce, a generic v2/v3 client) can open and read the common parts of.

> **Format stability.** The on-disk format is at `spec_version` **0.1** and **may change in
> backward-incompatible ways before 1.0**. Every store records its `spec_version` and readers check it,
> so the format is self-describing; reaching **1.0** is what turns that version line into a compatibility
> commitment. Until then, regenerate stores from source rather than relying on them for long-term archival.

> **Implemented today.** The Python, C++, R, and JS writers emit the tree below in **Zarr v3 by
> default** (`zarr.json` per node + inline consolidated metadata); legacy **Zarr v2** (`.zarray`/
> `.zgroup`/`.zattrs` + a consolidated `.zmetadata`) is still available as an opt-in (`format="v2"`).
> All four surfaces read *both* formats (chunked + gzip + zstd included), and **v3 sharding** (the
> `sharding_indexed` codec — the request-economy path for million-cell remote stores) is written and
> read across surfaces. A nullable validity `mask` and **partial coverage** (`coverage="partial"` + an
> `index` into `index_axis`) are written and read everywhere. The `recipe`/`ragged`/`raster` encodings
> and `models/` are not yet written.

## The tree

```text
store.lstar.zarr/               Zarr v3 by default (per-node zarr.json); v2 layout noted below
├── zarr.json                   {node_type:"group", attributes:{lstar:<root metadata>},
│                                consolidated_metadata:{…}}   # inline: one read, required for HTTP access
├── axes/
│   ├── zarr.json
│   └── <axis>/
│       ├── zarr.json           attributes.lstar = {kind:"axis", origin, role, induced_by, provenance}
│       ├── labels/             the element labels (utf8: a uint8 byte array …)
│       └── labels_offsets/     … + an int64 offsets array, length n+1
├── fields/
│   ├── zarr.json
│   └── <field>/
│       ├── zarr.json           attributes.lstar = {kind:"field", role, span, state, encoding, shape, …}
│       └── <value arrays>      depend on the encoding (below); chunk keys are c/<i>
├── passthrough/                lossless passthrough of a format's untyped long-tail
│   ├── zarr.json
│   └── <namespace>/            e.g. "anndata.uns"
│       ├── zarr.json           attributes.lstar = {kind:"passthrough", tree:"<JSON string>", arrays:[{id, kind}]}
│       └── <id>                the array leaves the tree references (dense, or utf8 bytes+offsets)
└── models/                     fitted transforms (apply contract + weights)   [spec; not yet emitted]
```

Legacy **Zarr v2** (`format="v2"` / `--zarr-format v2`) writes the *same* tree with v2's metadata files
instead: each node's `lstar` metadata lives in `.zattrs` (a group also has `.zgroup`, an array `.zarray`),
consolidated metadata is a separate root `.zmetadata` document, and chunk keys are `<i>` (dotted for N-D).
The L★ content is identical across formats — only the container differs — and all four surfaces read both.

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

## Lossless passthrough (`passthrough/`)

A format's untyped long-tail — AnnData `uns`, Seurat `@misc`/`@commands` (params, color palettes,
dendrograms, DE tables) — is preserved **verbatim** under `passthrough/<namespace>/` rather than recorded
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

## Chunking, compression & sharding

- Arrays may be **single-chunk** (the portable default) or **chunked** along the first axis (set
  `chunk_elems` on write). Chunking is what lets a lazy reader fetch only the blocks a query touches —
  e.g. one gene's CSC column.
- Chunks may be **uncompressed** or compressed with **gzip** or **zstd** (`compressor=` on write; zstd is
  a Zarr v3 codec and zarr-python 3's own default). All four surfaces **read** both codecs — C++/R via
  libzstd (the build degrades to gzip-only if it's absent), the JS/WASM reader via a decode-only zstd
  build — and Python/CLI/R/JS all **write** both. The C++ core reads any chunk grid (C-order, fill-padded
  edge chunks; a missing chunk reads as `fill_value` 0). (`blosc` is read by zarr-python but no lstar
  writer emits it.)
- A **compressed** array stays byte-range-readable: the reader decodes only the chunk(s) covering the
  requested range, not the whole array — so compression and the sub-chunk fast path (one gene / one cell)
  coexist. Smaller chunks trade more over-read per read against more objects/index entries.
- **v3 sharding** (`sharding_indexed`; `shard_elems` on write) packs many inner chunks into fewer shard
  *objects* — fewer files to host — while each inner chunk stays byte-range-readable via the shard index.
  Written by Python/CLI/R/JS, read by all surfaces.
- **Per-field write layout.** Chunking, compression, and sharding can be set **per field** (xarray-
  `encoding` style — a `write` override on a field), so hot and bulk arrays in one store get different
  layouts. The `viewer@0.1` prep uses this (see the viewer section).

## Consolidated metadata

So a reader makes **one request** instead of many small ones (required for access over HTTP), the store
carries consolidated metadata. In **v3** it is inline under the root `zarr.json`'s `consolidated_metadata`
key (every node's metadata, keyed by path); in **v2** it is a separate root `.zmetadata` document
(`{"zarr_consolidated_format": 1, "metadata": {…}}`). Every surface's writer emits it; readers prefer it
(`zarr.open_consolidated`, or the libstar/WASM reader's consolidated open) and fall back to walking the tree.

## The viewer profile (`viewer@0.1`)

A *profile* is a contract: a writer stamps `profiles` with a name@version, and a reader may then
**rely** on the fields that profile guarantees. `viewer@0.1` is the profile the `lstar-viewer` web
app consumes — precomputed summaries so the browser renders differential-expression / variable-gene /
dotplot views and does low-latency selection reads without recomputing client-side. Any surface
(Python `extend_for_viewer`, R `extend_for_viewer`, the JS prep, a profile writer like `pagoda2`) that
stamps `viewer@0.1` MUST satisfy this contract; `validate()` enforces it.

**Precompute is optional.** A plain (un-prepped) store opens and is fully usable — the viewer computes
markers, stats, selection DE and overdispersion **on the fly** (the same kernels). The profile only
*precomputes* these as global navigators so a large/remote store opens instantly. So this contract
binds a store **only when it declares `viewer@0.1`**; a bare store is valid and viewable without it.
Because the precompute and the on-the-fly path call the **same** core kernels, a prepped field equals
its live-computed value (within the conformance tolerance).

**Cross-surface consistency.** The Python, R, and JS/WASM preps drive the *same* C++ core for both the
heavy kernels **and** the cell reorder (`viewer_cell_order`: cluster-contiguous, then a Hilbert curve
over the embedding when present), and each normalizes the input `counts` to the layout the kernels
expect — so counts may arrive **CSR or CSC**, and a store prepped by any surface is field-for-field
identical to one prepped by another (byte-identical on `counts_cellmajor` / `counts_cellmajor_order`;
stats/markers/`od_score` to the conformance tolerance). Grouping auto-detection is single-sourced across
the surfaces. The `conformance/viewer*.sh` legs (incl. a corpus-driven cross-surface check) and
`conformance/policy_linter.py` enforce this.

A `viewer@0.1` store contains, for **at least one** categorical cell grouping `<g>` (e.g. `leiden`,
`cell_type`), with `cells`/`genes` standing for the observation/feature axes and `K` = number of
groups in `<g>` (an induced factor axis `groups_<g>` or the grouping's own factor axis):

| field | required | encoding | span (orientation) | meaning |
|---|---|---|---|---|
| `counts_cellmajor` | **yes** | `csr` (`state` follows the basis: raw counts → int; `lognorm` → float) | `[cells, genes]` | cell-major copy of the count basis; substrate for per-cell row reads + scope compute |
| `stats_<g>_sum` | **yes** | dense | `[<g>, genes]` = **K×ng** | Σ `log1p(counts)` per (group, gene) |
| `stats_<g>_sumsq` | **yes** | dense | `[<g>, genes]` = K×ng | Σ `log1p(counts)²` per (group, gene) |
| `stats_<g>_nexpr` | **yes** | dense | `[<g>, genes]` = K×ng | count of nonzero cells per (group, gene) |
| `markers_<g>_lfc` | **yes** | dense | `[genes, <g>]` = **ng×K** | 1-vs-rest log-fold-change (see below) |
| `markers_<g>_padj` | **yes** | dense | `[genes, <g>]` = ng×K | 1-vs-rest significance surrogate (see below) |
| `od_score` | **yes** | dense | `[genes]` | per-gene overdispersion score (pagoda2 F-test; variable-gene ranking) |
| `counts_cellmajor_order` | optional | dense (`state="permutation"`) | `[cells]` | cell → physical row in `counts_cellmajor` after a locality reorder |

**Orientation (load-bearing, do not transpose).** Stats are **group-major** (`K×ng`): a row is one
cluster's full gene profile (contiguous). Markers are **gene-major** (`ng×K`): a row is one gene's
per-cluster values — the path the viewer takes when coloring an embedding by a marker. The two
orientations are deliberately different and are part of the contract.

**On-disk compression (default per-field layout).** A `viewer@0.1` store is compressed by default, each
array's codec/chunking/sharding chosen for its access pattern (via the per-field write layout): the
gene-major count basis stays **raw, single-chunk** so gene-coloring reads exact bytes with no decode;
`counts_cellmajor` is **zstd, chunked + sharded** so a cell-subset read touches ~one chunk (per-chunk
decompress); the dense `stats_*` / `markers_*` / `od_score` (read whole) are **zstd, single-chunk** for
the best ratio at no read penalty. `compress_primary` also compresses the gene-major counts (smaller
store, a small decode on gene-coloring); `compress=false` writes it all uncompressed. This is a *default*,
not part of the contract — the field names/encodings/orientations above are the contract; codecs are free.

**Marker definition (`viewer.markers/1-vs-rest`, a fast surrogate — not a calibrated test).** Over
`log1p(counts)`, with group size `n_g` and rest size `n − n_g`:
`lfc[gene, g] = mean_log_in_g − mean_log_in_rest`;
`padj[gene, g] = clip(exp(−|lfc| · sqrt(nexpr_g + 1)), 1e-12, 1)`.

**Overdispersion (`od_score`, pagoda2-style).** Per gene, take `mean`/`var` of `log1p(counts)` over
all cells and `nobs` = number of expressing cells. Fit `log(var) ~ log(mean)` with a tricube lowess
(span 0.3, 200 anchors, linear interpolation, constant edge extrapolation), giving residual
`r = log(var) − trend`. The score is `−log P(F > exp(r); df1 = df2 = nobs)` — the upper-tail variance-
ratio F-test from pagoda2's `adjustVariance`, so a sparsely-expressed gene (small `nobs`) can't reach a
high score. Genes with `nobs < 3`, `mean ≤ 0`, or `var ≤ 0` get `0`. The F-test score (rather than the
raw residual `r`) is canonical because it is what the viewer computes **live**, so prepped == live.

**Hybrid order (`counts_cellmajor_order`).** When present, `counts_cellmajor`'s rows are physically
permuted by `lexsort(hilbert_index, primary_cluster_code)` (primary key = the first grouping's cluster
code, secondary = a Hilbert index over a 1024×1024 grid of the min-max-scaled 2-D embedding). The
field stores, for each cell, its physical row (`f8`, exact integers), and records the grouping it was
keyed on in `provenance.group`. The reorder-key grouping is the first *detected* grouping by default;
`extend_for_viewer(primary=…)` names it explicitly — set it to the grouping the viewer opens on so the
locality reorder matches the first view. A reader keys on the `_order` **sibling** of a field name so
that a cluster/lasso selection coalesces into a few byte-range reads.

These quantities have a single C++-core implementation bound to every surface; pure-language
fallbacks (Python without the accel extension, JS before WASM) must match the core within the
conformance tolerance (`stats` exact; `lfc`/`od` to ~1e-3).

**The navigators are caches.** Each of the eight fields above carries `provenance.cache = "viewer@0.1"`:
they are *regenerable* re-encodings/summaries of `counts` (+ the grouping), carrying no decision a user
made — unlike clusters, embeddings, or graphs (which are *primary* analysis results and are never
tagged). The tag is **advisory metadata, not an access gate**: read-by-name (the viewer, `lstar.read`,
byte-range reads) reads cache fields normally. What acts on the tag is **format-mapping conversion** —
exporting a store to a non-viewer object (`write_anndata`, `write_seurat`, `write_sce`,
`Pagoda2$fromLstar`) **drops `cache` fields and records them in `dropped`** rather than carrying a
redundant or mis-aligned copy (notably, `counts_cellmajor` is physically row-reordered, so a generic
consumer must not read it as an aligned matrix — the tag prevents that). Eager `lstar.read` still loads
everything (use `lazy=True` to avoid materializing unused caches).

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

### Single-file `.lstar.zarr.zip` (STORED)

A store can be packaged as **one file** — a `.lstar.zarr.zip`: the same directory store, zipped with
every entry **STORED (no deflate)**. This is the artifact to host or hand around; the directory form
stays the default *working* format (a zip is append-only, so it can't back an "open → add a field →
resave" or stream-in-place workflow).

**STORED is a requirement, not a preference — for two reasons:**

1. **No double compression.** Zarr chunks are already codec-compressed (gzip/…), so re-deflating them
   in the zip layer costs CPU for essentially no size gain.
2. **Byte-range readability.** Only a STORED entry is readable by byte range *inside* the archive. A
   hosted single file is read by issuing an HTTP `Range` into the zip for one chunk (lstar's JS
   `ZipStore` / `httpZipSource`); a DEFLATE-compressed entry would force
   fetching and inflating the whole entry, silently defeating range access. lstar therefore **forces
   STORED on every `.zip` write** (regardless of a chunk `compressor=`) and **rejects a DEFLATE-packed
   `.lstar.zarr.zip` on read** with an actionable message (repack it STORED). Large stores (>4 GB or
   >65535 entries) use **ZIP64** transparently.

**Produce one** (all forced STORED): `lstar convert dir.lstar.zarr out.lstar.zarr.zip` (repackage a
directory store), `lstar convert … --viewer out.lstar.zarr.zip` (a single-file *viewer* store — the
hosted-viewer artifact), or from a library — Python `lstar.write(ds, "x.lstar.zarr.zip")`, R
`lstar_write(ds, "x.lstar.zarr.zip")`, JS `writeStoreZip("x.lstar.zarr.zip", ds)`.

**Read one** on any surface: Python `lstar.read("x.lstar.zarr.zip")` (over zarr's `ZipStore`), R
`lstar_read(...)` and C++ `lstar::read(...)` (extract the STORED entries locally, then read — a copy
with no decompression), and JS `ZipStore.open(nodeFileSource(path))` for a local file or
`ZipStore.open(httpZipSource(url))` for a **hosted** zip (the seek-into-zip path that makes a
single URL a range-served store).

## Cross-implementation guarantee

The same store reads byte-faithfully from **Python, R, C++, and the browser** (the last via the libstar
core compiled to WebAssembly — the same reader, not a separate JS one). A store written by any one reads
in the others with identical field values (counts sums, graph nnz, embeddings); v2 and v3 stores, chunked
+ gzip/zstd, and sharded, all round-trip across every surface, as does the single-file `.lstar.zarr.zip`
packaging. This is enforced by the conformance suite
([`../conformance/`](../conformance) — e.g. `zip.sh` / `zip_r.sh` / `zip_js.sh`).
