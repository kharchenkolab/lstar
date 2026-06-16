# lstar 0.0.1

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

## Performance

* Streamed, multi-threaded (OpenMP) per-gene/per-cell reductions over the compiled core
  (`lstar_read_block`, `lstar_stream_col_sum_by_group`, `col_sum_by_group`), with bounded memory and
  thread-invariant results.
