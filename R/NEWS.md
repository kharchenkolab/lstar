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
