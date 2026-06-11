# lstar

A lightweight, fast library implementing **L★** — a uniform model and Zarr interchange format for
single-cell and spatial omics — with bindings in **Python**, **R**, and **C++**, and bidirectional
conversion to/from the common formats (AnnData, Seurat, SingleCellExperiment, the pagoda-verse /
Conos).

L★ represents a dataset as a set of **axes** (the entities one indexes by — cells, genes, clusters,
samples, …) and **fields** (typed data over axes — counts, embeddings, graphs, labels, …). Existing
formats are expressed as **profiles**: precise, bidirectional mappings into L★. See
[`misc/Lstar_proposal.md`](misc/Lstar_proposal.md) for the model and
[`misc/plan1.md`](misc/plan1.md) for the implementation plan.

> Status: early development. Tri-language (Python/C++/R) read+write of the same store; profiles for
> AnnData, Seurat (v3/v4/v5), SingleCellExperiment, and Conos; collection model. Not yet released.

## Collections, not just matrices

A central design point: a *collection of heterogeneous samples* is a collection, not one aligned
`cells × genes` tensor. L★ represents it with a `samples` axis, **per-sample** `cells.<s>`/`genes.<s>`
axes and measures (samples may differ in cells *and* gene sets), and a **union** `cells` axis for the
joint analysis layer — joint embedding, clusters, and the integration **graph** as a `relation`. The
R package ingests both a **Conos** object (`write_conos`) and a split **Seurat v5** assay
(`split(assay, f = sample)` → a collection) this way. See
[`examples/conos_collection_demo.R`](examples/conos_collection_demo.R).

## Version recognition

External formats ship many versions with different object layouts, so profiles **detect** the
version and degrade gracefully rather than assuming one shape: Seurat v3/v4 `Assay` vs v5 `Assay5`
(and a `GetAssayData` fallback for SeuratObject < 5); pagoda2's removed `$counts` slot vs the
`getRawCounts()` accessor; AnnData's `.raw` slot (kept on its own gene axis when it diverges). The
detected `<format>@<version>` is recorded in the store, and anything a profile can't represent is
recorded in `dropped` rather than silently lost.

## Documentation

- [`docs/principles.md`](docs/principles.md) — the idea and the design philosophy (start here)
- [`docs/model.md`](docs/model.md) — the model spec (axes, fields, roles, encodings, collections)
- [`docs/format.md`](docs/format.md) — the Zarr store layout (on-disk spec)
- [`docs/examples.md`](docs/examples.md) — worked, runnable examples (Python, R, C++, cross-language)

An agent skill (keywords, usage patterns, dense reference) lives in
[`.claude/skills/lstar/`](.claude/skills/lstar/SKILL.md). The full design proposal is in
[`misc/Lstar_proposal.md`](misc/Lstar_proposal.md).

## Layout

```
docs/            principles, model & format specs, worked examples
core/            libstar — the C++ core (model, chunked+gzip Zarr IO, fast translation)
python/          the `lstar` Python package (zarr-python + optional C++ accelerator)
R/               the `lstar` R package (cpp11 → libstar)
conformance/     shared round-trip / cross-profile / chunked test suite
examples/        runnable end-to-end demos
misc/            proposal + plans
```

## Quickstart (Python)

```python
import numpy as np, scipy.sparse as sp, lstar

ds = lstar.Dataset(kind="sample")
ds.add_axis("cells", [f"cell{i}" for i in range(100)])
ds.add_axis("genes", [f"g{i}" for i in range(50)])
ds.add_field("counts", sp.random(100, 50, density=0.1, format="csc"),
             role="measure", span=["cells", "genes"], state="raw")

lstar.write(ds, "sample.lstar.zarr")
ds2 = lstar.read("sample.lstar.zarr")
```

## Lazy & streaming

`read(path, lazy=True)` opens a store without materializing the heavy arrays, and a CSC measure
can be reduced by streaming column blocks — bounded memory, never densified:

```python
ds = lstar.read("big.lstar.zarr", lazy=True)        # opens without reading the heavy arrays
mean, var, nnz = lstar.stream_col_stats(             # per-gene stats, bounded memory
    ds.field("counts").values, lognorm=True,         # on-the-fly log1p, zero-aware variance
    n_threads=8)                                     # threading controllable from the call

# chunk + compress so lazy reads touch only the blocks they need:
import numcodecs
lstar.write(ds, "big.lstar.zarr", chunk_elems=1_000_000, compressor=numcodecs.GZip(5))
```

Measured on real data — Tabula Muris Marrow, 40,220 × 20,138, **77.6M nonzeros, float32** (see
[`misc/plan1.md`](misc/plan1.md) §12): lazy open **+9 MB vs +779 MB** eager (~87× less); per-gene
mean/var streamed in **~165 MB**, matching dense `np.var`; the C++ `log1p` reduction **7.6× at 16
threads** (thread count controllable, results thread-invariant); gzip-5 shrinks a raw-count store
**4.8×**. Float32 measures stay float32 end to end (no widening copy) while moments accumulate in
float64 — memory-lean, float64-accurate. AnnData → L★ → AnnData is a **fixed point** (any chain
length returns to the original; `uns` it can't hold is recorded in `dropped`). The C++ core reads
multi-chunk + gzip/zlib stores and writes a consolidated `.zmetadata`, so Python / C++ / R all read
the same chunked, compressed store.

## License

MIT.
