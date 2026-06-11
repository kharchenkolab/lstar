# Extending lstar for viewers

*Audience: an agent (or developer) extending the lstar implementations — C++ core `libstar`, the
Python/R packages, or the `js/` browser layer — to serve an interactive viewer. The reference viewer
is the next-generation pagoda2 app (`~/p21/pagoda2/misc/app_prop2.md`), but the guidance generalizes.*

This document catalogs the kinds of **flexibility and extension** a viewer needs from lstar, maps each
to the L★ model, says **where in the codebase** it belongs and **what already exists**, and — most
importantly — states the **principles** that keep such extensions from fragmenting the model. Read the
principles first; they are the part you must not violate.

---

## 0. Orientation: what a viewer asks of lstar

A viewer never loads a dataset whole. It opens a manifest, then issues many small, specific reads and
a handful of compute reductions, interactively, over HTTP or local disk. Concretely (from
app_prop2.md), it must: navigate an embedding of 1M+ cells; color cells by **one gene** or a metadata
column; show **grouped heatmaps** and cluster markers; **select** cells (lasso/metadata) and run
**differential expression on arbitrary selections**; and handle **multi-sample collections** as a
first-class case. Each of these is either an *access pattern* (read a slice sparingly), a *precomputed
derived field* (so a query needs no matrix access), or a *compute kernel* (a reduction run in WASM).

lstar's job is to make those three things cheap, consistent across languages, and expressible without
new schema slots.

---

## 1. Principles for extending (do not violate)

1. **Extend with fields and roles, never with schema slots.** Everything a viewer needs to *store* —
   precomputed cluster statistics, a cell-major DE panel, a coherent cell order, a dendrogram,
   quantized embeddings — is an ordinary L★ field over ordinary (possibly induced) axes. A summary is
   a `measure` over `(cluster, gene)`; pseudobulk is a `measure` over `(sample, cluster, gene)`; a cell
   order is a `measure`/permutation over `(cells)`. **Do not invent a `viewer/` slot namespace with
   bespoke semantics.** This is the whole reason L★ exists (proposal §1; see `docs/principles.md` §1a).

2. **Generic readers must keep working.** A viewer-oriented store is still a valid L★ store. A reader
   that doesn't know about "summaries" must see them as plain `measure` fields and ignore them; nothing
   a viewer adds may break `lstar.read` / `lstar_read` / `openLstar`. Viewer-specific *meaning* lives in
   a **profile** (a `viewer`/`pagoda2.app` profile), not in the core schema.

3. **One kernel, every runtime.** Compute belongs in the C++ core (`core/include/lstar/lstar.hpp`,
   the translation-primitives section), exposed to R (cpp11), Python (pybind11, `python/src/lstar/_accel.cpp`),
   and the browser (Emscripten, `js/wasm/lstar_wasm.cpp`). Never reimplement a numeric kernel in
   TypeScript or R that already exists (or should exist) in libstar — drift between runtimes is a bug.
   Every new kernel gets a **cross-language conformance test** (`conformance/`, e.g. `js.sh` already
   checks WASM == Python).

4. **Lean and fast by default.** Keep stored dtypes (float32 stays float32; quantize *for the eye* with
   uint8/int16 but retain full precision in `counts`); accumulate in float64. The fast path (the
   compiled accelerator / WASM) is used automatically when present, with a pure-language fallback —
   see `feedback-memory-lean-precision` and the `_engine.py` autodetect pattern.

5. **Record loss; report approximation.** A profile that can't represent something writes it to
   `dropped`. A viewer that downsamples or returns ranking-grade (not exact) statistics must *say so* in
   the result — silent approximation misleads biology (app_prop2.md §4, §5).

6. **The store schema is normative; the in-memory API is not.** Conformance is the on-disk store + the
   profiles (proposal §3.2). A viewer extension is "done" when the store round-trips and the readers in
   all relevant languages agree, not when one binding works.

---

## 2. Access-pattern flexibility (reading sparingly)

What the viewer needs to *fetch*, and how lstar should serve it.

| Need | L★ mapping | Where / status |
|---|---|---|
| **Manifest only** (axes/fields, no data) | root `.zattrs` + per-axis/field `.zattrs`, via consolidated `.zmetadata` | **Done.** `openLstar` (js), `read` (py/cpp) read the manifest first. |
| **One gene's expression** (color by gene) | a slice of a **gene-major CSC** measure: `indptr[g..g+1]`, then `data`/`indices` ranges | **Done** in `js/core/reader.ts::cscColumn`. Needs the measure stored **CSC** and **chunked** (a single chunk forces reading the whole array). |
| **One cell's genes** (inspect a cell; cell-major DE) | a slice of a **cell-major CSR** measure | **Gap.** Add a `csrRow(field, row)` mirroring `cscColumn`. Requires a CSR copy or a transpose (`csc_to_csr`). The viewer's cell-major DE panel (§3) is the intended store form. |
| **An embedding** (positions) | a dense `embedding` field, often int16-quantized | **Done** (`fieldDense`/`view.embedding`). Quantization dequant is a §4 gap. |
| **A metadata column** | a `label` (utf8 codes+categories) or numeric `measure` over `(cells)` | **Done** (`view.metadata`). |
| **Few requests on remote stores** | consolidated metadata + **Zarr v3 sharding** (two range reads per shard) | **Partial.** Consolidated `.zmetadata` done; **v3 + sharding is a gap** — the Python writer emits v2. This is *the* scaling lever for million-cell remote stores (app_prop2.md §3); zarrita.js already supports v3 sharding on read. |
| **Store transports** | directory-Zarr over HTTP (`FetchStore`), local dir / drag-dropped zip (File System Access / `ZipFileStore`), single-file reference store | **Partial.** The js reader takes any zarrita store; wire up `FetchStore`/`ZipFileStore` adapters next to `NodeFSStore`. The kerchunk-style reference store is later. |

**Guidance.** The unit of "read sparingly" is the *chunk*. For any field a viewer slices, ensure the
exporter writes it **chunked along the sliced axis** (CSC `data`/`indices` chunked by nnz; embeddings
chunked by cells). Reading economy is a property of how the store was *written* as much as how it's
read — see the `chunk_elems` parameter and `docs/format.md`.

---

## 3. Precomputed derived fields (the viewer profile)

The biggest single win: precompute the reductions so common queries need **no matrix access**. Each is
an ordinary L★ field; together they constitute a **viewer profile** (call it `pagoda2.app` /
`viewer@0.1`) whose *exporter* is the main new piece of work. None of these require a model change.

| Derived field | L★ form | Enables |
|---|---|---|
| **Cluster sufficient stats** | `measure` over `(grouping, genes)` holding `sum`, `sumsq`, `n_expr` (three fields, or an arity-3 `(grouping, genes, stat)`) | instant grouped heatmaps + cluster markers, no matrix read. Computed by a per-group kernel (§4). |
| **Pseudobulk** | `measure` over `(samples, grouping, genes)` (arity-3) | condition DE + composition, statistically correct for multi-sample. |
| **Marker tables** | `markers.<grouping>.lfc` / `.padj` : `measure` over `(genes, group)` + `uncertainty` | default ranked gene tables. |
| **Cell order** | a permutation `measure` over `(cells)` (dendrogram/embedding-coherent) | coherent storage order so an embedding selection maps to near-contiguous chunks → cheap subsample DE. |
| **Cell-major DE panel** | a `measure` (csr or dense uint8) over `(cells, od_genes)`, OD genes only, in cell order | subsample DE on arbitrary selections at constant cost (read ~hundreds of rows). |
| **Dendrogram** | a tree: a derived node axis + a `relation` parent field + a `membership` leaf-map (proposal §2.3) | heatmap ordering + cluster navigation. |
| **Normalization recipe** | a `recipe`/virtual field, or provenance on `data` | reproduce the display transform on the fly (pagoda2 matrix-views; §4). |

**Guidance.** Implement the exporter as a **profile writer** in Python and/or R (`write_viewer(ds,
…)`), reusing existing kernels: the cluster stats are exactly `csc_col_mean_var` generalized to
per-group (§4); the cell-major panel is `csc_to_csr` of the OD-gene submatrix. Mark the profile in
`ds.profiles` (`viewer@0.1`) so a reader knows the precomputed fields are present and trustworthy.
Because these are plain fields, a non-viewer reader still loads them as measures — principle 2 holds
for free.

---

## 4. Compute kernels the viewer needs (add to libstar)

Kernels live in `core/include/lstar/lstar.hpp` (templated on the value dtype, n_threads parameter,
double accumulation) and are bound to all three runtimes. Status and gaps:

- **Per-gene mean/variance** (`csc_col_mean_var`) — **done.** HVG ranking, the basis of cluster stats.
- **CSC↔CSR transpose** (`csc_to_csr`) — **done.** Orientation flips; building the cell-major panel.
- **Per-group sufficient statistics** — **gap.** A `csc_col_sum_by_group(data, indptr, indices,
  group_of_cell, n_groups, …)` returning `sum`/`sumsq`/`n_expr` per `(group, gene)` (pagoda2's
  `colSumByFac`/`colMeanVarView`). This powers grouped heatmaps and cluster markers and is the
  exporter's workhorse. Parallelize over genes (columns); accumulate per-group in thread-local buffers.
- **Subsample-DE ranking** — **gap.** An AUC/Wilcoxon (pagoda2 specificity) or fold-change ranker over
  a small cell-major submatrix. The current `js view.subsampleDE` does a fold-change pass in JS over the
  whole measure; promote the statistic to a libstar kernel and feed it the panel rows (§3).
- **On-the-fly normalization** of a sparse vector from a recipe (library size + log1p) — **partial.**
  `view.geneExpression` does log1p in JS; once a stored `recipe` exists, apply it consistently in a
  kernel so coloring matches what the analysis used.
- **Quantize / dequantize** (float ↔ uint8/int16-to-bbox) — **gap.** Small kernels; keep the bbox/scale
  as field attributes (a recipe) so a generic reader can recover physical values.
- **Crossfilter bitops** (AND across facet bitsets) — **partial.** `js view.Crossfilter` is fine in
  TS for now; move to WASM only if profiling at very large N demands it.

**Guidance.** When you add a kernel: template it on the value dtype, take `n_threads`, accumulate in
float64, and add a conformance row (a fixture comparing WASM/Python/C++ — see
`js/test/fixture_python.json` and `conformance/js.sh`). The kernel is the source of truth; the bindings
are thin.

---

## 5. Encodings, quantization, and virtual fields

- **Quantized arrays for the eye.** Embeddings as int16-to-bbox, expression as uint8 — these are just
  dense fields of that dtype (already storable). What's missing is a *recorded dequantization*: store
  the bbox/scale (a `recipe` or field attributes) so any reader recovers physical coordinates. Keep
  full precision in `counts`.
- **`recipe` (virtual) fields** — specified but **not implemented** (`docs/format.md`). The viewer's
  "store raw counts + a normalization recipe, not a materialized normalized matrix" is exactly this:
  `data` is a `recipe` over `counts` + library sizes. Implementing `recipe` (compute-on-read, with the
  op/inputs/params in `.zattrs`) benefits viewers *and* the pagoda2 matrix-views story. Start by
  reading a recipe and materializing in the kernel; lazy/virtual evaluation later.
- **`ragged`/`raster`** — specified, not implemented; needed for spatial/imaging viewers and sequence
  data, not for the first pagoda2 app.

**Guidance.** A new encoding is a cross-cutting change: add it to the spec (`misc/Lstar_proposal.md`
Appendix A + `docs/format.md`), to every writer and reader (Python `zarr_io.py`, C++ `read_array`/IO,
js `reader.ts`), and to the conformance suite. Until all readers support it, gate it behind a profile
capability so older readers degrade rather than fail.

---

## 6. Scale: streaming, LOD, threads

- **Coarse-first / downsample** — a precomputed subsample index (a `measure` over `(cells)` marking a
  ~50–100k representative subset) lets the viewer render instantly, then stream the rest. Store it as a
  field; it's also pyramid level 0.
- **Embedding tile pyramids** (LOD, for >5M cells) — a multiscale, quantized representation
  (app_prop2.md §6). Model as a derived per-level axis + quantized embedding fields, or reuse the
  OME-NGFF multiscale convention. This is a large, optional, late-phase item.
- **WASM threads** — the WASM kernels are **single-threaded** today (no `-pthread`). Threaded reductions
  need Emscripten pthreads + `SharedArrayBuffer`, which requires the page to be **cross-origin
  isolated** (COOP/COEP headers). Worth it only once a reduction is the interaction bottleneck; until
  then, stream and keep blocks small. The kernels already take `n_threads`, so enabling threads is a
  build-flag + deploy-header change, not a rewrite.

---

## 7. Collections in viewers

The collection model (proposal §2, Appendix B.5; implemented for Conos + split Seurat) is what
multi-dataset viewing rests on. A viewer reads: per-sample `counts.<s>` measures; the `samples` axis
and a `sample` design label over the union `cells`; and joint embedding/clusters/graph. From these it
does small-multiples (filter by sample — `O(cells)` bitset), composition (per-(sample × cluster) counts
straight from labels), and condition DE (pseudobulk, §3). **Guidance:** don't add a separate "multi"
code path — a collection is the same axes/fields with per-sample namespacing; the viewer profile's
summaries simply gain a `samples` dimension (pseudobulk).

---

## 8. What stays in the viewer (not lstar)

So agents don't push the boundary the wrong way: **rendering** (deck.gl/regl), the **crossfilter/
coordination state**, **linked-view wiring**, lasso/point-in-polygon, color maps, and the **application
shell** (the framework decision in app_prop2.md §7) are the viewer's, not lstar's. lstar provides the
*data primitives* — fields, lazy reads, label codes, bitset-ready arrays, and the compute kernels — and
stops at the typed arrays. `js/core/view.ts` is the seam: it returns `Float32Array`/codes/stats; the
app turns those into pixels.

---

## 9. How to add an extension (recipes)

**Add a compute kernel (e.g. per-group stats).** Implement `template<class T> … (…, int n_threads)` in
`core/include/lstar/lstar.hpp` (double accumulation). Bind it: pybind11 in `python/src/lstar/_accel.cpp`,
Emscripten in `js/wasm/lstar_wasm.cpp`, cpp11 in R if R needs it. Add a Python reference + a fixture and
a `conformance/` check that WASM == Python == C++. Re-vendor the header into `R/inst/include/`.

**Add a precomputed field (viewer profile).** Write a profile exporter (`write_viewer` in Python/R)
that computes the field with an existing kernel and adds it via `add_field` with a clear name and role
(`measure`/`relation`/`label`), recording `viewer@0.1` in `ds.profiles`. No model or reader change —
it's a plain field. Add a reader-side convenience in `js/core/view.ts` that *uses* the precomputed
field when present and falls back to computing it otherwise.

**Add an encoding.** Spec (`misc/Lstar_proposal.md` Appendix A) → `docs/format.md` → every writer/reader
→ conformance. Gate behind a profile capability until all readers support it.

**Add a store transport.** Implement the zarrita `get(key)` store (or reuse `FetchStore`/`ZipFileStore`)
and pass it to `openLstar`; no reader change. Put node-only stores behind a separate import (as
`node-store.ts` is) so the browser bundle stays clean.

---

## 10. Pointers

- Viewer design: `~/p21/pagoda2/misc/app_prop2.md` (and `app_prop1.md`).
- The model + the normative profile catalog: `misc/Lstar_proposal.md`; the docs: `docs/`.
- Current implementation: `core/` (C++ kernels + IO), `python/src/lstar/`, `R/`, `js/`
  (`js/README.md`, `js/core/reader.ts`, `js/core/view.ts`, `js/wasm/lstar_wasm.cpp`).
- The phased WASM/viewer plan (A–E) and what's done: this repo's commit log and `js/README.md` —
  A (WASM kernels), B (reader), C (view API) are done; **D (viewer profile + exporter + Zarr v3/
  sharding) and E (a deck.gl slice) are the next substantive work**, and most of §3–§5 above lands in D.
