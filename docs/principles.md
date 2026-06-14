# Principles

L★ exists to make single-cell/spatial data **move** — between formats (AnnData, Seurat,
SingleCellExperiment, Conos, pagoda2) and between languages (Python, R, C++) — without loss and
without forcing everything into one rigid container. This page explains the ideas behind it.

## 1. A dataset is axes + fields

Most formats hard-code a few special matrices (an expression matrix, an embedding, a graph) and a
grab-bag of metadata. L★ generalizes that to two primitives:

- an **Axis** is a labelled set you index by — `cells`, `genes`, a `pca` coordinate axis, a `samples`
  axis;
- a **Field** is typed data over a *tuple* of axes, tagged with a **role** that says what it means.

So a count matrix is a `measure` over `(cells, genes)`; a UMAP is an `embedding` over
`(cells, umap)`; PCA loadings are a `loading` over `(genes, pca)`; a kNN graph is a `relation` over
`(cells, cells)`; a clustering is a `label` over `(cells)`. One algebra, not a dozen special cases.
New kinds of data don't need new container slots — they're just fields over axes.

> **Example.** In AnnData, `X` → measure `(cells, genes)`; `obsm['X_pca']` → embedding;
> `varm['PCs']` → loading; `obsp['connectivities']` → relation; an `obs` column → a label or a
> 1-arity measure over `(cells)`. The mapping is mechanical because the model is general.

## 1a. The long tail has a home

This is the argument that motivates the whole construction (proposal §1.1). Every fixed schema gives
you a few named slots, and the *routine* workflow fits them. But analysis keeps producing structures
that don't:

- **RNA velocity** — the vectors fit in `layers`, but the velocity *graph* becomes an opaque
  `uns['velocity_graph']`, and fate probabilities get jammed into `obsm` as if they were an embedding.
- **Gene regulatory networks** — a TF→target graph lives on the *gene* axis, which no slot names, so
  SCENIC ships a separate `.loom`.
- **Cell–cell communication** — a `cluster × cluster × ligand-receptor` tensor is defined on the
  *cell-type* axis, which the schema doesn't have.
- **Fitted models** — a trained scVI or CellTypist isn't data; it's a function applied to data.

Each of these is a structure over a *different* pair (or triple) of axes — `cells×cells`,
`genes×genes`, `celltype×celltype×lr`. Adding an AnnData slot for each never converges. The
axis-and-field construction expresses all of them uniformly: a new method adds a **role** and maybe an
**axis**, not a schema slot, so the format does not change as methods do. That is why L★ has a small,
general core instead of a growing list of slots — and why the off-main-path work that today lands in
`uns`/`misc` becomes first-class, queryable, provenance-bearing fields.

## 2. A collection is not a tensor

This is the central design choice. A study is usually **many heterogeneous samples** — different
donors, ages, conditions, sometimes different species or feature sets. Two ways to represent that:

- **The aligned-tensor way** (scverse/AnnData/MuData): one global `cells × genes` matrix, samples
  reduced to a categorical column, modalities forced to share cells. Clean, but it *erases* the
  per-sample structure and assumes a common feature space.
- **The collection way** (Conos): a *linked set* of samples, each with its own cells (and possibly
  genes), integrated by a **joint graph**, not by concatenation.

L★ takes the collection view as first-class. A `kind="collection"` dataset has:

- a `samples` axis (the members),
- per-sample `cells.{s}` and `genes.{s}` axes and `counts.{s}` measures — samples may genuinely
  differ in cells *and* gene sets,
- a union `cells` axis for the joint analysis layer, with a `sample` **design** label,
- the joint `embedding`, joint cluster `label`s, and the integration **graph** as a `relation` over
  `(cells, cells)`.

Alignment is legitimate only *within* a sample (its facets share cells); *across* samples you keep a
collection joined by a graph. A single Pagoda2 object is the `sample` unit; a Conos object is the
`collection` unit.

> **Example.** A real two-sample Conos object (16,022 cells) ingests as a collection: two
> `counts.{s}` measures with their own cell axes, a `sample` label, a joint UMAP, a joint Leiden
> labelling, and the integration graph as a 16,022² relation — and it reads back identically in
> Python. Flattened to one matrix, all of that per-sample structure would be gone.

## 3. Lossless, recorded conversion

A converter you can't trust is useless. L★'s contract: **native → L★ → native is a fixed point** —
a value can leave its format, pass through L★ any number of hops (even across formats and languages),
and return unchanged. The shared *vocabulary* (`counts`, `X`, `pca`, `umap`, `leiden`, …) keeps field
names aligned across formats so the round-trip is meaningful, and each field records its **native
location** in `provenance` so write-back is exact.

When a target format genuinely can't hold something — AnnData `uns`, or a neighbor graph passing
through Seurat — it is recorded in `dropped`, **never silently changed**. Loss you can see is
recoverable; silent loss is a bug.

> **Example.** AnnData → L★ → AnnData over a real 40k-cell dataset is a fixed point across repeated
> cycles: `X`, `obs`, `obsm`, and `obsp` come back identical; `uns` shows up in `ds.dropped`.

A round-trip back to the *same* format, though, only proves L★ preserved its own representation — it
says nothing about whether the object handed to Seurat is *canonical Seurat that Seurat's own tools
accept*. So conversion is also **deterministic and native-valid by construction**: a field's typed
`(role, state, span)` fixes its canonical slot in each target (a `DimReduc` with a `_`-terminated key, a
categorical `groupby`, the `logcounts` assay scran looks for), and `lstar convert --check` verifies it by
opening the result in its *own* library and running a canonical-ops smoke. The full role→slot contract is
**[mapping.md](mapping.md)**.

This is the practical, near-term reason to use lstar — moving data between AnnData, Seurat, and SCE
without losing your analysis. The dedicated guide is **[conversions.md](conversions.md)**.

## 4. Recognize versions gracefully

External formats change shape across releases, and a glue library that assumes one layout breaks the
day a user upgrades. L★ **detects** the version and adapts:

- a legacy Seurat **v2** object (the pre-`Assay` lowercase `seurat` S4 class — slots read via `attr()`,
  so the ancient class need not be defined) through v3/v4 (`Assay`, fixed slots) vs v5 (`Assay5`,
  layers), with a `GetAssayData` fallback for older SeuratObject that lacks the `Layers()` API;
- pagoda2's removed `$counts` slot vs the `getRawCounts()` accessor;
- AnnData's `.raw` (kept on its own gene axis when it diverges) and the `uns['neighbors']` → `obsp`
  migration.

The detected `<format>@<version>` is recorded in the store, and anything unrepresentable goes to
`dropped`. Adapt, don't assume.

## 5. Fast and memory-lean by default

Glue sits on hot paths — large sparse matrices, remote stores. So:

- **Lazy & streaming.** `read(lazy=True)` opens a store without materializing the heavy arrays; a CSC
  measure is reduced by streaming column blocks in bounded memory (the matrix is never densified).
- **Memory-lean precision.** A float32 measure (the common case) **stays float32** end to end — no
  widening copy — while moments accumulate in float64. Low-precision storage, high-precision
  accumulation: lean *and* accurate. Don't upcast to gain precision; accumulate in float64 instead.
- **Fast by default.** Python uses a compiled C++/OpenMP accelerator automatically when present and
  the pure-Python path otherwise — identical results, no opt-in. Reductions scale to the core count
  in C++ (≈20× at 56 threads on a 77M-nonzero measure).
- **Determinism contract (stated and tested).** The streaming reducers are **thread-count-invariant and
  bit-identical**: each column is accumulated independently in float64, column-parallel, with **no
  cross-thread reduction**, so a column's result is a pure function of its data — not of how many threads
  ran. The same input gives byte-for-byte the same mean/variance/nnz at 1, 2, 4, or 8 threads (asserted
  exactly, `== 0`, in `conformance/stream_reduce.sh` and `python/tests/test_determinism.py`). This is
  both a reproducibility guarantee and the property that lets a *different* library reuse the same kernel
  and get identical summaries — the precondition for sharing one tuned core across lstar, a viewer, or a
  consumer like pagoda2.

> **Example.** Lazily opening a 77.6M-nonzero measure costs ~9 MB instead of ~780 MB; per-gene
> mean/variance streams in bounded memory and matches a dense computation; the C++ reduction is
> ~7.6× faster at 16 threads, with the thread count controllable from the call.

## 6. One core, three languages

The numerics live once, in the header-only C++ core `libstar`. R binds it via cpp11; Python binds it
via pybind11 (and falls back to a pure-Python reference when the extension isn't built). The same
`.lstar.zarr` store is read byte-faithfully from all three. This keeps the model honest (a single
source of truth) and lets each ecosystem use its native objects.

---

Next: the [model spec](model.md), the [format spec](format.md), and the [examples](examples.md).
