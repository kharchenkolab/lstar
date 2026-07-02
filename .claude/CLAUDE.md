# lstar (L★)

## What lstar is

L★ is a general model for single-cell/spatial omics data, built from **axes** (the entities you
index by — cells, genes, samples, clusters) and **fields** (typed data over them — counts,
embeddings, loadings, graphs, labels, designs). Because everything is just axes and fields, one
small model spans cases a fixed `cells × genes` container strains on: multi-sample (even
cross-species) integrations kept as a *collection* of heterogeneous samples, CITE-seq with a second
feature axis, case-control cohorts carrying a statistical *design*.

Its most immediate use is **lossless, explicit conversion** between the formats people already use —
AnnData, Seurat, SingleCellExperiment, pagoda/conos — routing each through one shared-vocabulary L★
store, preserving the meaning of every piece and *reporting* (never silently dropping) whatever a
target can't hold. It reads/writes a portable Zarr-based format and streams heavy operations in
bounded memory, so multi-gigabyte datasets convert on a laptop.

## Components (one model, four surfaces)

| Surface | Path | Role |
|---|---|---|
| C++ core | `core/` | `libstar` — header-only: the model, chunked+gzip Zarr IO, the fast kernels |
| Python | `python/` | the `lstar` package on zarr-python, binding the optional compiled C++ accelerator |
| R | `R/` | the `lstar` package; the format profiles (Seurat, SCE, Conos) live here |
| Browser/Node | `js/` | a TypeScript reader (zarrita) + the kernels compiled to WebAssembly, for viewers |

`conformance/` holds the shared round-trip / cross-format / cross-language test suite.

## Core aim: cross-language consistency

**We aim to implement consistent logic and feature sets across C++, Python, R, and JS/WASM — ideally
by minimizing code duplication** (share the C++ core; bind it from Python and R; compile it to WASM
for the browser) rather than reimplementing the same logic four times and letting the versions drift.
When adding or changing a feature, treat all four surfaces as one system: keep the behavior, the API
shape, and the results identical across languages, and extend the `conformance/` suite so any drift
is caught.
