# lstar 0.2.0

This release brings the Zarr **v3** on-disk format to every surface (C++/Python/R/JS) and makes
compressed, range-readable viewer stores the default.

## Default on-disk format is now Zarr v3

* `lstar_write()` now defaults to `format = "v3"` (was `"v2"`): stores are written with a per-node
  `zarr.json` + inline consolidated metadata. The legacy Zarr v2 layout (`.zarray`/`.zgroup`/`.zattrs`
  + a consolidated `.zmetadata`) remains available via `format = "v2"`, and `lstar_read()` reads both
  formats transparently — so the change is invisible to readers; only newly *written* stores change
  layout. All four surfaces share the default.

## Zstd compression and sharded writes

* Stores can be written and read with **zstd** — Zarr v3's standard codec — in addition to gzip, across
  all four surfaces. `lstar_write()` gained a `shard_elems` argument for **sharded** v3 writes, which
  pack many inner chunks into fewer store objects while staying range-readable, so a many-chunk array
  can be hosted without a file-per-chunk explosion.

## Compressed, range-readable viewer stores by default

* `extend_for_viewer()` now compresses the viewer store **per field** by default (zstd): the gene-major
  count basis stays a single raw chunk for exact-byte gene-color reads, the cell-major counts are zstd
  chunked + sharded, and every other array is zstd single-chunk. The reader resolves compressed arrays
  at **chunk granularity**, so a hosted viewer fetches only the chunks it displays instead of whole
  arrays. Pass `compress = FALSE` for the previous all-raw layout, or `compress_primary = TRUE` to trade
  gene-color latency for a smaller store.

# lstar 0.1.7

## `extend_for_viewer()` auto-selects the count basis (no longer errors on normalized-only inputs)

* `extend_for_viewer()` now auto-selects the basis for the viewer's counts instead of erroring when an
  object kept only normalized values. It prefers **raw** counts (`log1p`-transformed); failing that it
  falls back to a **log-normalized** measure (used as-is, with a warning that HVG / marker rankings are
  then approximate); failing that it raises a clear error. A **scaled / z-scored** measure is never
  chosen — previously a name-based fallback could pick a scaled `X`, corrupting the ranking statistics.
  So a converted 'scanpy' object that dropped its raw layer now yields a working viewer store. The
  selection contract is identical across R (`.viewer_counts_basis`), Python (`_select_counts_basis`) and
  JS (`selectCountsBasis`), and is enforced by cross-surface parity tests.

# lstar 0.1.6

The R package version jumps 0.1.0 -> 0.1.6 to align with the companion Python package (`lstar-sc` on
PyPI) and the shared on-disk format; the entries below cover everything the R package gained since the
0.1.0 CRAN release.

## Seurat → viewer prep (the previously-untested seam)

* A Seurat object converted with `read_seurat()` then run through `extend_for_viewer()` now yields a
  clean viewer store: a **logical** `meta.data` column (a QC flag like `qc_kept`) stays boolean and is
  **not** detected as a viewer grouping (was coerced to a `"TRUE"/"FALSE"` string and became a noise
  grouping), and the **active identity** (`Idents()`, captured as the `ident` field) no longer duplicates
  the clustering it mirrors — the viewer's grouping detection skips the `active_ident` mirror on all
  surfaces (Python/R/JS). The active ident is still preserved for the Seurat round-trip. New
  `conformance/viewer_seurat.sh` covers the seam (synthetic Seurat in CI; a real `SeuratData` object
  locally): boolean QC excluded, opens on a real clustering, viewer@0.1-clean.

## `extend_for_viewer(primary=)` — align the prep with the viewer's default open

* `extend_for_viewer()` gains a `primary` argument: the grouping the viewer opens on. It is hoisted to the
  front of the prepared groupings, so it keys the `counts_cellmajor` locality reorder AND is summarized
  first — the eager-prepare a fast launch waits on. Unlike ordering groupings by hand, `primary` **composes
  with auto-detect** (`primary="cell_type"` with `grouping=NULL` preps *every* detected grouping but keys the
  reorder on `cell_type`) — which matters because the auto-detect policy prefers clusterings while the viewer
  may open on a cell-type annotation. `counts_cellmajor_order` now records `provenance$group` (the reorder
  key), matching Python/JS. Same `primary=` option added to the Python and JS/WASM `extend_for_viewer`. A
  `primary` that isn't a grouping over the cell axis is rejected with a clear error (was a cryptic reorder
  crash). Cross-surface parity (Py==R==JS reorder for a given `primary`) is enforced by
  `conformance/viewer_primary.sh`.

## Single-file `.lstar.zarr.zip` packaging

* `lstar_read()` / `lstar_write()` now accept a single-file `*.lstar.zarr.zip` (a store packed into ONE
  file with every entry **STORED**, so its already-compressed chunks stay byte-range-readable when
  hosted — the point of a single file). Writing forces STORED (never DEFLATE) and is ZIP64-aware;
  reading a DEFLATE-packed `.lstar.zarr.zip` is rejected with a clear message. R rides the C++ core's
  `.zip` dispatch, so it reads/writes the same artifact as Python, C++, and the browser (JS reads a
  *hosted* zip by HTTP range). See `docs/format.md` §Packaging; enforced by `conformance/zip_r.sh`.

## viewer@0.1 cross-language parity

* `extend_for_viewer` now yields a store field-for-field identical to the Python and JS/WASM preps.
  The cell reorder is the shared C++ core (`viewer_cell_order`: cluster-contiguous, then a Hilbert curve
  over the embedding) instead of a cluster-only sort; grouping auto-detection returns **all** groupings
  ranked by a single-sourced preferred-name policy (was: a single grouping, ranked differently); and a
  `basis = "lognorm"` prep keeps `counts_cellmajor` float.
* Enforced by cross-surface conformance legs (`conformance/viewer*.sh`, including a corpus-driven check
  over `corpus.py`/`synth.py`) and `conformance/policy_linter.py`. See `docs/parity.md`.

## Cross-surface fidelity (parity audit)

* R now preserves a float32 value dtype (was widened f4→f8) and a graph relation's `directed`/`weighted`
  flags (were dropped) across a `lstar_read`/`lstar_write` round-trip. Guarded by `conformance/r_fidelity.sh`.
* `extend_for_viewer` gains `order=` and `markers=` (parity with Python); grouping detection restricted to
  string-like labels over the cell axis; the lognorm measure-name fallback now picks in field order (was
  name-list order). All viewer policy constants are single-sourced (`viewer_policy.json` + `policy_linter.py`).
* Store fidelity: coo matrices normalize to csc on write (portable to every reader); dense fields carry
  their shape in the manifest; the C++ core preserves a field's `uncertainty`; the `.h5ad` direct backend
  infers `state` from content like the native backend.

# lstar 0.1.0

First release. lstar is a uniform data model (L\*) and a Zarr interchange format for single-cell /
spatial omics, with a shared C++ core (`libstar`) and bindings in R, Python and C++.

## Data model & store

* Datasets of **axes** (labelled sets you index by) and **fields** (typed data over a tuple of axes,
  with a role / span / encoding / state / coverage / provenance).
* Read / write `.lstar.zarr` stores (`lstar_read`, `lstar_write`), multi-chunk and gzip-compressed,
  byte-compatible with the Python and C++ readers.
* **Collections** of heterogeneous samples (`kind = "collection"`) — per-sample `cells.<s>`/`genes.<s>`
  axes over gene sets that may overlap, differ, or be disjoint, plus a union `cells` axis carrying the
  joint embedding / clustering / integration graph. Build one from any list of per-sample objects with
  `collection_from()`.

## Format converters (profiles)

* **Seurat** (legacy v2 through v5; a split v5 assay is read as a collection, and a collection is written
  back as a split assay) — `read_seurat` / `write_seurat`. Cell-cell graphs round-trip as `Graphs()`.
  **Multimodal** objects (CITE-seq RNA+ADT, multiome RNA+ATAC, ECCITE …) round-trip every assay, each on
  a **canonical feature axis** (`proteins`, `peaks`, …) shared with MuData/pagoda2 — so a modality is the
  same L\* feature space regardless of source format; the original assay name is kept in provenance.
* **SingleCellExperiment** — `read_sce` / `write_sce`.
* **Conos** — `write_conos` (Conos → L\*) and `read_conos` (L\* → a live Conos), preserving the
  per-sample data and the joint graph / embedding / clustering.
* The shared vocabulary keeps the common core (counts, data/X, pca + loadings, umap, labels, metadata)
  lossless across conversions; what a target cannot hold is recorded in `dropped`, never lost silently.

## Viewer profile (`viewer@0.1`)

* `extend_for_viewer(ds)` and the `lstar viewer <store>` CLI add the **viewer profile**: a cell-major
  `counts_cellmajor` (physically reordered cluster-contiguous, with a `counts_cellmajor_order`
  permutation for locality reads), per-grouping cluster stats (`stats_<g>_*`, group-major), 1-vs-rest
  marker tables (`markers_<g>_*`, gene-major), and a pagoda2-style `od_score` (lowess + F-test). The
  profile is specified in `docs/format.md` and enforced by `validate()`.
* The recipe math (`markers_one_vs_rest`, `overdispersion`) lives in the shared `libstar` core and is
  bound to R, Python and WebAssembly, so a store prepped from any surface — and the browser viewer's
  on-the-fly compute — agree (a cross-language conformance gate checks it). `write_pagoda2` now emits
  a fully conformant `viewer@0.1` store.

## Performance

* Streamed, multi-threaded (OpenMP) per-gene/per-cell reductions over the compiled core
  (`lstar_read_block`, `lstar_stream_col_sum_by_group`, `col_sum_by_group`), with bounded memory and
  thread-invariant results.
