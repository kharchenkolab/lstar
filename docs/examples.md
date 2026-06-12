# Examples

These are worked, runnable examples chosen to teach the model by doing. Each one opens with **what it
shows** and why it matters, and the code carries running comments explaining each step — read the
comments, not just the calls. Full end-to-end scripts live in [`../examples/`](../examples); the
conformance suite that exercises the same paths is in [`../conformance/`](../conformance).

If you haven't read [the model](model.md) yet, the one-line version: a dataset is **axes** (the things
you index by) and **fields** (typed data over axes), and a field's `role` says what it means
(`measure`, `embedding`, `loading`, `relation`, `label`, …).

Converting between formats (AnnData ↔ Seurat ↔ SCE ↔ Conos) is the most common reason to reach for
lstar — examples 2 and 6–10 below cover it, and **[conversions.md](conversions.md)** is the dedicated
guide (the function table, what's preserved vs. recorded as dropped).

Contents:

- [1. Build a dataset by hand, write, read, validate (Python)](#1-build-write-read-validate-python)
- [2. AnnData round-trip — and why it's a *fixed point* (Python)](#2-anndata-round-trip-python)
- [3. Lazy open + streamed, threaded per-gene stats (Python)](#3-lazy-streaming-python)
- [4. Chunking + compression, and what they buy you (Python)](#4-chunking--compression-python)
- [5. Build a collection by hand — heterogeneous samples (Python)](#5-build-a-collection-python)
- [6. R: read/write and the Seurat profile](#6-r-readwrite-and-seurat)
- [7. R: a split Seurat v5 assay *is* a collection](#7-r-seurat-v5-split-collection)
- [8. R: SingleCellExperiment](#8-r-singlecellexperiment)
- [9. R: a Conos collection](#9-r-a-conos-collection)
- [10. A round-trip across formats *and* languages](#10-cross-language-round-trip)
- [11. C++: read a store and reduce with OpenMP](#11-c-read-and-reduce)
- [12. Browser/Node: read an L★ store and answer viewer queries (WASM)](#12-browser-wasm)

---

## 1. Build, write, read, validate (Python)

**What this shows.** The whole model in code: you create *axes*, add *fields* over them with a `role`,
and serialize. Note that a 2-D embedding *induces its own coordinate axis* (`umap`) — that's
[induction rule 1](model.md#three-induction-rules-why-the-model-stays-small) — and that a string
vector over `cells` is inferred to be a `label`. You never touch a fixed "slot"; you declare what each
field *is*.

```python
import numpy as np, scipy.sparse as sp, lstar

# A dataset is created empty; "sample" means one uniform measurement (vs. "collection", example 5).
ds = lstar.Dataset(kind="sample")

# Axes are labelled, ordered sets. The labels ARE the identity of each position (label-keyed).
ds.add_axis("cells", [f"cell{i}" for i in range(100)], role="observation")  # the observations
ds.add_axis("genes", [f"g{i}" for i in range(50)], role="feature")          # the features

# A measure over (cells, genes). state="raw" lets later steps refuse, e.g., clustering raw counts.
# A scipy CSC matrix maps to L*'s 'csc' encoding (gene-compressed) automatically.
counts = sp.random(100, 50, density=0.1, format="csc")
ds.add_field("counts", counts, role="measure", span=["cells", "genes"], state="raw")

# An embedding places `cells` into a 2-D coordinate space. Those 2 columns ARE the 'umap' axis
# (induction rule 1) — so we declare the axis, then the field over (cells, umap).
ds.add_axis("umap", ["umap0", "umap1"], origin="derived", role="coordinate")
ds.add_field("umap", np.random.randn(100, 2), role="embedding", span=["cells", "umap"])

# A categorical assignment of cells. A string vector is inferred to be a `label` if you omit role=.
ds.add_field("leiden", np.array([f"c{i % 5}" for i in range(100)]),
             role="label", span=["cells"])

lstar.write(ds, "sample.lstar.zarr")          # -> a Zarr v2 group with consolidated metadata
ds2 = lstar.read("sample.lstar.zarr")

# validate() returns [] when the store is well-formed: spans reference real axes, shapes match, etc.
assert lstar.validate(ds2) == []
print(ds2.field("counts").encoding,            # 'csc' — inferred from the scipy format
      ds2.field("counts").values.shape)        # (100, 50)
```

**What to notice.** No `obsm`/`obsp`/`uns` — just fields with roles. Adding a *new* kind of result
later (a gene–gene network, a fate-probability embedding) is the same three lines, never a schema
change.

---

## 2. AnnData round-trip (Python)

**What this shows.** A profile (`anndata`) is a *bidirectional* mapping. Reading an AnnData yields L★
fields named in the shared vocabulary (`X`→`X`, `obsm['X_pca']`→`pca`, `obsp['connectivities']`→a
relation); writing reverses it to the original slots. Anything L★ can't put back (`uns`) is **recorded
in `dropped`, not silently lost**. And because every field records its native location, the loop is a
*fixed point*: native → L★ → native returns the original, for a chain of any length.

```python
import anndata as ad, numpy as np, lstar
from lstar.profiles.anndata import read_anndata, write_anndata

adata = ad.read_h5ad("pbmc.h5ad")
ds = read_anndata(adata)                  # AnnData slots -> L* fields (X, layers, obs/var, obsm, obsp, .raw)
print(ds.profiles)                        # e.g. ['anndata@0.1', 'anndata@0.8.0'] — who wrote it + version
print(ds.dropped)                         # e.g. ['uns/neighbors'] — recorded losses, never silent

lstar.write(ds, "pbmc.lstar.zarr")
adata2 = write_anndata(lstar.read("pbmc.lstar.zarr"))
# X / obs / obsm (PCA, UMAP) / obsp (graphs) are restored to their original slots;
# anything dropped is re-surfaced under adata2.uns["lstar/dropped"] so it's visible.

# The fixed-point property: repeated native->L*->native conversions don't drift. A conversion CHAIN
# of any length therefore returns to the original format unchanged on the representable content.
cur = adata
for _ in range(4):
    cur = write_anndata(read_anndata(cur))
assert np.allclose(np.asarray(cur.obsm["X_pca"]), np.asarray(adata.obsm["X_pca"]))
```

**What to notice.** The shared vocabulary is what makes this *meaningful* across formats: `X_pca` and a
Seurat PCA both become `pca`, so example 10 can hand the same fields to a different ecosystem. Full
script with a real 40k-cell dataset: [`examples/roundtrip_chain.py`](../examples/roundtrip_chain.py).

---

## 3. Lazy streaming (Python)

**What this shows.** L★ stores are meant to be read *sparingly*. `read(lazy=True)` opens a store
without materializing the heavy arrays; a CSC measure is then reduced by streaming column blocks, so
per-gene statistics run in bounded memory and the matrix is *never densified*. The threading degree is
controllable from the call, and the fast path (a compiled C++ kernel) is used automatically when
present.

```python
import numpy as np, lstar

# Lazy: the store opens in megabytes, not the whole matrix. Field .values are proxies (LazyCSX),
# backed by the open zarr arrays — nothing is read until you touch it.
ds = lstar.read("big.lstar.zarr", lazy=True)
v = ds.field("counts").values                          # LazyCSX(csc, shape=(C, G), nnz=...)

# Per-gene mean/variance over log1p-normalized values, computed by STREAMING column blocks:
#   - lognorm=True applies log1p per nonzero on the fly (the normalized matrix is never built),
#   - n_threads: 1 = serial (default), N = N threads, 0 = all cores,
#   - engine="auto" uses the C++/OpenMP accelerator if installed, else pure Python (identical results).
mean, var, nnz = lstar.stream_col_stats(v, lognorm=True, n_threads=0)
hvg = np.argsort(var)[::-1][:2000]                     # the top highly-variable genes

lstar.show_config()                                    # prints whether the C++ accelerator is active
print("accelerated:", lstar.has_accel())               # True if a wheel/build provided it
```

**What to notice (measured on real data, 40,220 × 20,138, 77.6M nonzeros).** Lazy open costs ~9 MB vs.
~780 MB eager; the C++ reduction is ~7.6× faster at 16 threads with the thread count controllable; a
float32 measure stays float32 (no widening) while the moments accumulate in float64 — *memory-lean,
float64-accurate*. Scripts: [`examples/lazy_streaming_demo.py`](../examples/lazy_streaming_demo.py),
[`examples/real_perf.py`](../examples/real_perf.py).

---

## 4. Chunking + compression (Python)

**What this shows.** *Why* chunking matters: a single-chunk array can't be read sparingly (any slice
pulls the whole chunk), so to stream a gene's column you must chunk. Compression shrinks the store; the
ratio tracks the data's entropy.

```python
import numcodecs, lstar
ds = lstar.read("sample.lstar.zarr")

# chunk_elems chunks each array along its first axis (~that many elements/chunk), so a lazy read of
# one gene's CSC column touches only a few chunks; GZip compresses each chunk.
lstar.write(ds, "sample.gz.lstar.zarr",
            chunk_elems=1_000_000,                 # without this, arrays are a single chunk
            compressor=numcodecs.GZip(5))
```

**What to notice.** gzip shrinks a *raw-count* store ~4.8× (integers compress well) but a *normalized
float* matrix only ~1.6× (high entropy) — the ratio is a property of the data, not a bug. The C++ and
browser readers read these chunked, compressed stores directly.

---

## 5. Build a collection (Python)

**What this shows.** The model's signature feature: a *collection of heterogeneous samples* is **not**
one aligned `cells × genes` matrix. Each sample keeps its **own** namespaced `cells.<s>` (and
`genes.<s>`) axis and its own `counts.<s>` measure; the joint analysis lives over a *derived union*
`cells` axis carrying a `sample` label and the integration graph as a `relation`. (Compare to the
proposal's conos example, Appendix B.5.)

```python
import numpy as np, scipy.sparse as sp, lstar

ds = lstar.Dataset(kind="collection")               # the collection kind
A = sp.csc_matrix(sp.random(40, 25, density=0.2))   # sample S1: 40 cells x 25 genes
B = sp.csc_matrix(sp.random(60, 25, density=0.2))   # sample S2: 60 cells (DIFFERENT cell count)

# Per-sample, namespaced axes + measures. Samples may differ in cells AND genes — nothing forces a
# shared axis. This is the heterogeneity an aligned tensor erases.
for s, M in (("S1", A), ("S2", B)):
    ds.add_axis(f"cells.{s}", [f"{s}_c{i}" for i in range(M.shape[0])], role="observation")
    ds.add_axis(f"genes.{s}", [f"g{i}" for i in range(M.shape[1])], role="feature")
    ds.add_field(f"counts.{s}", M, role="measure",
                 span=[f"cells.{s}", f"genes.{s}"], state="raw")

# The `samples` axis is the collection itself.
ds.add_axis("samples", ["S1", "S2"], role="sample")

# The joint analysis layer lives over a DERIVED union `cells` axis (induction rule 3).
ucells = [f"S1_c{i}" for i in range(40)] + [f"S2_c{i}" for i in range(60)]
ds.add_axis("cells", ucells, origin="derived", role="observation")

# `sample` is a design label: which sample each union cell came from (the partition).
ds.add_field("sample", np.array(["S1"] * 40 + ["S2"] * 60),
             role="label", span=["cells"], subtype="design")

# The integration graph is a `relation` over (cells, cells) — not a special slot, just a field.
G = sp.csc_matrix(sp.random(100, 100, density=0.05))
ds.add_field("graph", G, role="relation", span=["cells", "cells"], subtype="knn")

lstar.write(ds, "collection.lstar.zarr")
assert lstar.read("collection.lstar.zarr").kind == "collection"
```

**What to notice.** There is no concatenated matrix. Alignment is legitimate *within* a sample;
*across* samples you keep a collection joined by the graph. A portable, dependency-free version is the
conformance test [`conformance/collection.sh`](../conformance/collection.sh).

---

## 6. R: read/write and Seurat

**What this shows.** The R package binds the same C++ core; an `lstar_dataset` round-trips to a Seurat
object via the `seurat` profile. Note the **orientation flip**: L★ measures are `cells × genes`, but
Seurat assays are `genes × cells`, so the profile transposes (a `dgCMatrix` CSC maps directly to L★
`csc`).

```r
library(lstar)
ds  <- lstar_read("pbmc.lstar.zarr")
so  <- write_seurat(ds)        # L* -> Seurat: measures -> assay layers (transposed to genes x cells),
                               #                embeddings -> DimReduc, arity-1 cell fields -> meta.data
ds1 <- read_seurat(so)         # Seurat -> L*; ds1$profiles records SeuratObject@<v> + assay@<v3|v5>
lstar_write(ds1, "pbmc.from_seurat.lstar.zarr")

field_value(ds, "leiden")      # accessor: a field's values (a vector / matrix / Matrix)
```

**What to notice.** `read_seurat` *detects the Seurat version* (v3/v4 `Assay` vs. v5 `Assay5`) and
adapts — it doesn't assume one object layout.

---

## 7. R: Seurat v5 split collection

**What this shows.** Seurat v5's integration workflow holds samples *unintegrated* as a split
`Assay5` — `split(assay, f = sample)` makes per-sample layers (`counts.<sample>`). That is a
collection, and `read_seurat` ingests it as one (`kind="collection"`), exactly like example 5.

```r
library(lstar); library(SeuratObject)
obj$sample <- rep(c("donorA", "donorB"), length.out = ncol(obj))
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$sample)   # per-sample layers: counts.donorA, counts.donorB
ds <- read_seurat(obj)
stopifnot(ds$kind == "collection")                    # samples axis + per-sample cells.<s>/counts.<s>
```

Full script: [`examples/seurat_collection_demo.R`](../examples/seurat_collection_demo.R).

---

## 8. R: SingleCellExperiment

**What this shows.** The same profile pattern for Bioconductor's SCE: assays ↔ measures, `reducedDims`
↔ embeddings, `colData` ↔ arity-1 cell fields.

```r
library(lstar); library(SingleCellExperiment)
ds  <- lstar_read("pbmc.lstar.zarr")
sce <- write_sce(ds)           # measures -> assays (genes x cells), embeddings -> reducedDims, fields -> colData
ds2 <- read_sce(sce)           # and back to L*
```

---

## 9. R: a Conos collection

**What this shows.** The reference collection: a Conos object is a set of Pagoda2 samples linked by a
joint graph. `write_conos` restates it as L★ axes and fields — exactly the structure of example 5,
read from a real object. The pagoda2 version is detected (the modern `getRawCounts()` accessor vs. the
legacy `$counts` slot).

```r
library(lstar)
co <- readRDS("conos_two_sample.rds")     # a Conos R6 object (a collection of Pagoda2 samples)
ds <- write_conos(co)                     # -> an L* collection:
#   samples axis;
#   per-sample cells.<s>/genes.<s> axes + counts.<s> measures + pca.<s> embeddings;
#   a union `cells` axis; a `sample` design label;
#   the joint embedding, joint leiden label, and integration graph as a (cells x cells) relation.
print(ds)
lstar_write(ds, "conos.lstar.zarr")       # reads back identically in Python / C++ / the browser
```

**What to notice.** On real two-sample data this is a ~16k-cell collection with the joint graph as a
16,022² relation — none of which an aligned matrix could hold. Full script:
[`examples/conos_collection_demo.R`](../examples/conos_collection_demo.R).

---

## 10. Cross-language round-trip

**What this shows.** Because the shared vocabulary aligns field names across profiles, a round-trip can
cross **formats *and* languages** and return to the origin. The chain length is arbitrary; the
shared-vocabulary core (`counts`/`X`, `pca`, `umap`, a label) survives, and what a given hop can't
carry (neighbor graphs through Seurat, `uns`) is recorded.

```text
AnnData (Python) → L* → [ Seurat (R) → L* → SCE (R) → L* ] × N → AnnData (Python)
```

```bash
# Drives a real AnnData through R's Seurat and SCE profiles and back, checking the core survives:
bash examples/roundtrip_xlang.sh 1 path/to/file.h5ad
```

Source: [`examples/roundtrip_xlang.sh`](../examples/roundtrip_xlang.sh). The synthetic version runs in
CI as [`conformance/cross_format.sh`](../conformance/cross_format.sh).

---

## 11. C++: read and reduce

**What this shows.** The header-only C++ core (`libstar`) reads any L★ store (multi-chunk + gzip) and
runs the translation primitives directly on the stored arrays — the same kernels the R and Python
packages use, so the numbers match.

```cpp
#include "lstar/lstar.hpp"
#include <cstdio>

int main() {
    auto ds = lstar::read("big.lstar.zarr");        // reads a chunk grid; decodes gzip/zlib
    auto* f = ds.field("counts");                   // a CSC measure (cells x genes)
    auto ip = lstar::as_i64(f->indptr);             // normalize the index dtype (int32 or int64 by size)

    // Per-gene log1p mean/variance, OpenMP over columns; the float32 data is read IN PLACE (no copy).
    auto s = lstar::csc_col_mean_var(f->data.as<float>(), ip.data(),
                                     f->shape[1], f->shape[0], /*n_threads*/0, /*lognorm*/true);
    printf("gene0: mean=%.4f var=%.4f\n", s.mean[0], s.var[0]);

    // Orientation flip (cells x genes -> genes x cells), preserving the value dtype:
    auto idx = lstar::as_i64(f->indices);
    auto csr = lstar::csc_to_csr(f->data.as<float>(), idx.data(), ip.data(),
                                 f->shape[0], f->shape[1]);
    return 0;
}
```

Build: `-Icore/include -std=c++17 -fopenmp` (add `-DLSTAR_HAVE_ZLIB -lz` for gzip). Working test:
[`core/test/test_chunked.cpp`](../core/test/test_chunked.cpp).

---

## 12. Browser WASM

**What this shows.** The data layer for a web viewer: TypeScript reads the L★ store over HTTP or local
disk via zarrita.js — fetching only the chunks a view needs — and the libstar **WASM** kernels do the
compute, so the browser shows the same numbers as R/Python. This is the [`js/`](../js) package.

```ts
import { openLstar, LstarView, scalarToRGBA } from "lstar-js";
import { NodeFSStore } from "lstar-js/node-store";          // browser: new zarrita FetchStore(url)

const ds = await openLstar(new NodeFSStore("sample.lstar.zarr"));
ds.kind; ds.axisNames(); ds.fieldNames();                   // the manifest (one consolidated read)

// The hot path: fetch ONE gene's CSC column as a slice (not the whole matrix), normalize on the fly.
const view = new LstarView(ds);
const expr = await view.geneExpression("g3", { lognorm: true });  // per-cell scalar, 0 where unexpressed
const colors = scalarToRGBA(expr.values, expr.max);               // -> RGBA for a deck.gl color attribute

const { data, n, dim } = await view.embedding("umap");      // positions for a point layer
const md = await view.metadata("leiden");                   // {kind:'categorical', codes, categories}
const { mean, var: variance } = await view.colStats({ lognorm: true });  // per-gene stats via WASM
```

**What to notice.** `colStats` runs the WASM kernel and matches Python's `stream_col_stats` exactly —
one C++ core, three runtimes. The gene-column fetch reads only that column's chunks, which is the whole
point of reading "sparingly." See [`js/README.md`](../js/README.md) and the tests in `js/test/`.
