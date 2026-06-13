---
name: lstar
description: >-
  Use when working with the lstar package or the L-star (L*) data model and Zarr interchange format
  for single-cell / spatial omics, ESPECIALLY to convert single-cell data between formats (AnnData /
  h5ad, Seurat, SingleCellExperiment, Conos, pagoda2) and languages (Python, R, C++). Covers
  converting/exporting/importing between those formats via profiles (read_anndata/write_anndata,
  read_seurat/write_seurat, read_sce/write_sce, write_conos), building datasets of axes and fields,
  reading/writing .lstar.zarr stores, collections of heterogeneous samples, lazy/streaming reads, the
  C++ accelerator (libstar), per-gene reductions, and format/version recognition. Keywords: lstar, L*,
  L-star, convert, conversion, glue, interchange, h5ad, AnnData, Seurat, SingleCellExperiment, SCE,
  Conos, pagoda2, profile, export, import, axes, fields, measure, embedding, loading, relation, label,
  collection, zarr, csc, csr, lazy, streaming, stream_col_stats, libstar, accelerator, single-cell.
---

# lstar

lstar is a lightweight, fast "glue" library: a uniform data model (**L\***) and a **Zarr**
interchange format for single-cell/spatial omics, with bindings in **Python**, **R**, and **C++**
(shared core `libstar`), and bidirectional converters for **AnnData, Seurat, SingleCellExperiment,
Conos, and pagoda2**. Repo: `~/p21/lstar` (local git only — do **not** push to GitHub).

## The model in three sentences

- A dataset is a set of **Axes** (labelled sets you index by — `cells`, `genes`, `pca`, `samples`)
  and **Fields** (typed data over a tuple of axes — counts, embeddings, graphs, labels).
- A field has a **role** (`measure | embedding | loading | relation | label | sequence | design |
  transform`), a **span** (which axes it lives over), an **encoding** (`dense | csr | csc | coo |
  utf8`), and optional `state` (`raw | lognorm | scaled`), `coverage`, `provenance`.
- A **collection of heterogeneous samples is a collection, not one aligned `cells × genes` tensor**:
  a `samples` axis + per-sample `cells.<s>`/`genes.<s>` axes & measures + a union `cells` axis for
  the joint embedding/clusters/graph.

## When to use this skill

Reach for it when the task involves: converting between single-cell formats; reading/writing
`.lstar.zarr`; representing a multi-sample collection without flattening; computing per-gene/per-cell
statistics at scale (lazy, streamed, multithreaded); making format conversion version-robust; or
building/packaging the lstar Python/R libraries.

## Main usage patterns

**Convert between formats (the near-term selling point).** `convert(X → Y) = write_Y(read_X(obj))`,
with the L★ dataset (or an on-disk `.lstar.zarr` store) as the universal intermediate. Readers/writers:
`read_anndata`/`write_anndata` (Python), `read_seurat`/`write_seurat`, `read_sce`/`write_sce`,
`write_conos` (R). The shared vocabulary makes it lossless on the common core (counts, data/X, pca +
**pca_loadings**, umap, labels, metadata); what a target can't hold goes to `ds.dropped`, not silently.
```r
# R, in memory: Seurat -> SingleCellExperiment
sce <- write_sce(read_seurat(seurat_obj))
```
```bash
# Cross-language: AnnData (Python) -> Seurat (R), bridged by the on-disk store
python3 -c 'import anndata as ad, lstar; from lstar.profiles.anndata import read_anndata
lstar.write(read_anndata(ad.read_h5ad("x.h5ad")), "x.lstar.zarr")'
Rscript  -e 'library(lstar); saveRDS(write_seurat(lstar_read("x.lstar.zarr")), "x.rds")'
```
Full guide: `reference/conversions.md` and `docs/conversions.md`.

**Python — build, write, read, validate**
```python
import numpy as np, scipy.sparse as sp, lstar
ds = lstar.Dataset(kind="sample")
ds.add_axis("cells", [f"c{i}" for i in range(100)])
ds.add_axis("genes", [f"g{i}" for i in range(50)])
ds.add_field("counts", sp.random(100, 50, density=0.1, format="csc"),
             role="measure", span=["cells", "genes"], state="raw")
lstar.write(ds, "s.lstar.zarr")              # add chunk_elems=, compressor= for big/compressed
ds2 = lstar.read("s.lstar.zarr")             # add lazy=True to defer the heavy arrays
assert not lstar.validate(ds2)
```

**Python — AnnData round-trip (returns to the original format)**
```python
from lstar import read_anndata, write_anndata, write, read
write(read_anndata(adata), "a.lstar.zarr")
adata2 = write_anndata(read("a.lstar.zarr"))  # X/obs/obsm/obsp restored; uns -> ds.dropped
```

**Python — lazy + streamed + multithreaded per-gene stats (fast path automatic)**
```python
ds = lstar.read("big.lstar.zarr", lazy=True)             # +MBs, not the whole matrix
mean, var, nnz = lstar.stream_col_stats(                 # bounded memory, never densified
    ds.field("counts").values, lognorm=True, n_threads=0)  # n_threads: 1=serial, N, 0=all
lstar.show_config()                                       # is the C++ accelerator active?
```

**R — read/write and profiles**
```r
library(lstar)
ds  <- lstar_read("a.lstar.zarr")
so  <- write_seurat(ds);  ds1 <- read_seurat(so)     # Seurat legacy v2 -> v5; split assay -> collection
sce <- write_sce(ds);     ds2 <- read_sce(sce)
dsc <- write_conos(conos_obj)                        # a Conos collection -> L*
lstar_write(ds2, "out.lstar.zarr")
```

**C++ — read a store, reduce with OpenMP**
```cpp
#include "lstar/lstar.hpp"
auto ds = lstar::read("a.lstar.zarr");                 // multi-chunk + gzip/zlib
auto* f = ds.field("counts");
auto ip = lstar::as_i64(f->indptr);                    // index width-normalize
auto s  = lstar::csc_col_mean_var(f->data.as<float>(), ip.data(),
                                  f->shape[1], f->shape[0], /*n_threads*/0, /*lognorm*/true);
```

## Key principles (do not violate)

- **Collection ≠ tensor.** Never flatten a multi-sample dataset into one matrix; model it as a
  collection (see `reference/model.md`).
- **Memory-lean.** Don't widen stored dtypes; float32 measures stay float32, accumulate moments in
  float64. (`reference/python.md`, `reference/cpp.md`)
- **Recognize versions gracefully.** Detect Seurat legacy v2 → v5, pagoda2 accessor-vs-slot, AnnData
  `.raw`/uns layout; record `<format>@<version>`; route the unrepresentable to `dropped`, never
  silently lose it. (`reference/conversions.md`, `reference/r.md`)
- **Fast by default.** Python auto-uses the compiled C++ accelerator when present and falls back to
  pure Python; results are identical. Don't make users opt in. (`reference/python.md`)

## Reference files (read the one you need)

- `reference/conversions.md` — **format glue**: the readers/writers, the conversion matrix, what is
  preserved vs. dropped, version recognition. (The near-term selling point.)
- `reference/model.md` — axes, fields, roles, encodings, collections, the store layout.
- `reference/python.md` — full Python API, lazy/streaming, profiles, packaging.
- `reference/r.md` — full R API, profiles, CRAN packaging.
- `reference/cpp.md` — libstar core: structs, Zarr IO, translation primitives, the accelerator.
- `reference/recipes.md` — task-oriented recipes (conversions, collections, perf, round-trips).

Deeper background lives in `docs/` (principles, model & format specs, worked examples) and
`misc/Lstar_proposal.md` (the full design proposal).
