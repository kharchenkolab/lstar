# lstar (development version)

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
