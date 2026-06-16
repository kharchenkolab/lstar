---
name: lstar
description: >-
  Use when working with the lstar package or the L-star (L*) data model and Zarr interchange format
  for single-cell / spatial omics, ESPECIALLY to convert single-cell data between formats (AnnData /
  h5ad, Seurat, SingleCellExperiment, Conos, pagoda2) and languages (Python, R, C++). Covers
  converting/exporting/importing between those formats via profiles (read_anndata/write_anndata,
  read_seurat/write_seurat, read_sce/write_sce, write_conos/read_conos) or the one-command `lstar convert`
  CLI (with a fidelity report, native-acceptance check, and a package-free `--backend` fallback that
  converts .h5ad via h5py / Seurat & SCE .rds via base R, without anndata/SeuratObject/SingleCellExperiment),
  building datasets of axes and fields, assembling collections of heterogeneous samples from per-sample
  objects (collection_from) and converting a collection to a Seurat v5 split assay or a single AnnData,
  reading/writing .lstar.zarr stores, lazy/streaming reads, the
  C++ accelerator (libstar), per-gene reductions, and format/version recognition. Keywords: lstar, L*,
  L-star, convert, conversion, glue, interchange, h5ad, AnnData, Seurat, SingleCellExperiment, SCE,
  Conos, pagoda2, profile, export, import, lstar convert, CLI, native-acceptance, mapping, backend,
  package-free, h5py, axes, fields, measure, embedding, loading, relation, label,
  collection, collection_from, read_conos, split assay, zarr, csc, csr, lazy, streaming,
  stream_col_stats, libstar, accelerator, single-cell.
---

# lstar

lstar is a lightweight, fast "glue" library: a uniform data model (**L\***) and a **Zarr**
interchange format for single-cell/spatial omics, with bindings in **Python**, **R**, and **C++**
(shared core `libstar`), and bidirectional converters for **AnnData, Seurat, SingleCellExperiment,
Conos, and pagoda2**. Repo: `~/p21/lstar` → public GitHub remote `kharchenkolab/lstar`.

## The model in three sentences

- A dataset is a set of **Axes** (labelled sets you index by — `cells`, `genes`, `pca`, `samples`)
  and **Fields** (typed data over a tuple of axes — counts, embeddings, graphs, labels).
- A field has a **role** (`measure | embedding | loading | relation | label | sequence | design |
  transform`), a **span** (which axes it lives over), an **encoding** (`dense | csr | csc | coo |
  utf8`), and optional `state` (`raw | lognorm | scaled`), `coverage`, `provenance`.
- A **collection of heterogeneous samples is a collection, not one aligned `cells × genes` tensor**:
  a `samples` axis + per-sample `cells.<s>`/`genes.<s>` axes & measures + a union `cells` axis for
  the joint embedding/clusters/graph. Build one from any list of per-sample objects with
  `collection_from(...)` (Python & R), or get one from `write_conos` / a split Seurat v5 assay.

## When to use this skill

Reach for it when the task involves: converting between single-cell formats; reading/writing
`.lstar.zarr`; representing a multi-sample collection without flattening; computing per-gene/per-cell
statistics at scale (lazy, streamed, multithreaded); making format conversion version-robust; or
building/packaging the lstar Python/R libraries.

## Main usage patterns

**Convert between formats (the near-term selling point).** `convert(X → Y) = write_Y(read_X(obj))`,
with the L★ dataset (or an on-disk `.lstar.zarr` store) as the universal intermediate. The shared
vocabulary makes it lossless on the common core (counts, data/X, pca + **pca_loadings**, umap, labels,
metadata); what a target can't hold goes to `ds.dropped`, not silently.

The quickest path is the **`lstar convert` CLI** — it detects formats by path, routes through the store
(in-process for Python formats, an `Rscript` bridge for Seurat/SCE), and reports what crossed:
```bash
lstar convert x.h5ad x.rds                 # AnnData -> Seurat (.rds); --to sce for SingleCellExperiment
lstar convert x.h5ad x.lstar.zarr --report # -> store + fidelity report (fields + provenance + `dropped`)
lstar inspect x.h5ad                        # read + report its L★ structure, no write
```
`--check` (default on; `--strict` to gate the exit code) opens the result in its native library and runs
a canonical-ops smoke (scanpy/Seurat/scran) — verifies native tools accept it, not just that bytes
round-tripped. The deterministic role→slot contract is `docs/mapping.md`.

`--backend auto|native|direct` picks the codec: `auto` (default) uses the format's native package when
present, else lstar's **package-free** fallback — so you don't *need* the domain packages. Without them:
`.h5ad` ↔ store (read+write) needs only **`h5py`**; **Seurat `.rds`** ↔ store (read+write) and **SCE
`.rds`** → store (**read**) need only **base R + the `lstar` R package** (no SeuratObject/SingleCellExperiment).
Native-only: **SCE write** and `.h5mu`. At a wall (unknown on-disk version, `BPCells`-backed matrix) it
names the package to install. The analysis packages (scanpy/Seurat/scran) are only for `--check`.

Under the hood it's the readers/writers — `read_anndata`/`write_anndata` (Python),
`read_seurat`/`write_seurat`, `read_sce`/`write_sce`, `write_conos` (R) — which you can also call directly:
```r
sce <- write_sce(read_seurat(seurat_obj))   # R, in memory: Seurat -> SingleCellExperiment
```
Full guide: `reference/conversions.md` and `docs/conversions.md`; the mapping contract: `docs/mapping.md`.

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
so  <- write_seurat(ds);  ds1 <- read_seurat(so)     # Seurat legacy v2 -> v5; split assay <-> collection
sce <- write_sce(ds);     ds2 <- read_sce(sce)
dsc <- write_conos(conos_obj); con <- read_conos(dsc) # Conos collection <-> L* (round-trips)
lstar_write(ds2, "out.lstar.zarr")
```

**Assemble or convert a collection of heterogeneous samples**
```r
col <- collection_from(list(s1 = ds_a, s2 = ds_b),   # per-sample objects (lstar_dataset / Seurat / SCE),
                       joint = list(umap = emb, graph = knn))   # divergent or disjoint gene sets OK
so  <- write_seurat(col)   # -> Seurat v5 SPLIT assay (per-sample layers + Graphs + DimReduc); reads back as a collection
```
```python
import lstar
col = lstar.collection_from({"s1": ad_a, "s2": ad_b}, joint={"umap": emb})  # AnnData/MuData/Dataset
a   = lstar.write_anndata(col)   # -> ONE AnnData (flattens to union genes; graph->obsp, embedding->obsm)
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
- **Native-valid, not just round-trip-faithful.** A conversion's target must be canonical enough that
  the destination's *own* tools accept it (Seurat `_`-terminated keys, scanpy categorical `groupby`, SCE
  `logcounts`). The role→slot mapping is deterministic (`docs/mapping.md`); verify with
  `lstar convert --check` (open in the native lib + a canonical-ops smoke).
- **Fast by default.** Python auto-uses the compiled C++ accelerator when present and falls back to
  pure Python; results are identical. Don't make users opt in. (`reference/python.md`)

## Reference files (read the one you need)

- `reference/conversions.md` — **format glue**: the `lstar convert` CLI, the readers/writers, the
  conversion matrix, what is preserved vs. dropped, version recognition, and the deterministic role→slot
  mapping + native-acceptance (`docs/mapping.md`). (The near-term selling point.)
- `reference/model.md` — axes, fields, roles, encodings, collections, the store layout.
- `reference/python.md` — full Python API, lazy/streaming, profiles, packaging.
- `reference/r.md` — full R API, profiles, CRAN packaging.
- `reference/cpp.md` — libstar core: structs, Zarr IO, translation primitives, the accelerator.
- `reference/recipes.md` — task-oriented recipes (conversions, collections, perf, round-trips).

Deeper background lives in `docs/` (principles, model & format specs, worked examples) and
`misc/Lstar_proposal.md` (the full design proposal).
