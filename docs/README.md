# lstar documentation

**L★** is a uniform data model and a [Zarr](https://zarr.dev) interchange format for single-cell and
spatial omics, with bindings in **Python**, **R**, and **C++** (a shared core, `libstar`) and
bidirectional converters for **AnnData, Seurat, SingleCellExperiment, Conos, and pagoda2**. It is
meant to be the lightweight, fast *glue* between formats and languages — not another monolithic
container.

## Contents

- [**Principles**](principles.md) — the idea, the design philosophy, the long tail, and why a
  collection is not a tensor. Start here.
- [**Model**](model.md) — a precise, worked description of axes, fields, roles, the induction rules,
  and collections, building on the proposal.
- [**Format**](format.md) — the Zarr store layout (the on-disk spec).
- [**Examples**](examples.md) — worked, runnable, **commented** examples in Python, R, C++, and the
  browser. *Read this to learn by doing.*

The full normative specification — the model, the Zarr schema, and the bidirectional **profile rule
catalog** for AnnData, Seurat, pagoda2, Conos, and cacoa — is the proposal,
[`../misc/Lstar_proposal.md`](../misc/Lstar_proposal.md). These docs explain and exemplify it; that
document is the source of truth. Each docs page flags where lstar's current implementation covers only
part of the spec.

## Thirty-second tour

A dataset is a set of **axes** (labelled sets you index by) and **fields** (typed data over a tuple
of axes). Existing formats are expressed as **profiles** — precise, bidirectional mappings into L★.

```python
import scipy.sparse as sp, lstar
ds = lstar.Dataset(kind="sample")
ds.add_axis("cells", [f"c{i}" for i in range(100)])
ds.add_axis("genes", [f"g{i}" for i in range(50)])
ds.add_field("counts", sp.random(100, 50, density=0.1, format="csc"),
             role="measure", span=["cells", "genes"], state="raw")
lstar.write(ds, "sample.lstar.zarr")
ds2 = lstar.read("sample.lstar.zarr")
```

## Design tenets (the short version)

1. **A collection of heterogeneous samples is a collection, not one aligned `cells × genes` tensor.**
2. **Lossless, recorded conversion** — round-trips return to the original format; what can't be
   carried is recorded in `dropped`, never silently lost.
3. **Recognize versions gracefully** — adapt to Seurat v3/v4/v5, pagoda2 accessor-vs-slot, AnnData
   `.raw`/uns layouts, rather than assuming one object shape.
4. **Fast by default, memory-lean** — the same store reads from three languages; a compiled C++ core
   accelerates Python automatically; float32 stays float32 with float64 accumulation.

See [principles.md](principles.md) for the reasoning, and `../misc/Lstar_proposal.md` for the full
design proposal.
