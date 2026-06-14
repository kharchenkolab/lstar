# The conversion mapping — a deterministic, native-valid contract

**Round-trip fidelity ≠ native validity.** A `X → store → X` round-trip proves lstar preserved *its own*
representation; it says nothing about whether the object handed back to Seurat is *canonical Seurat that
Seurat's own tools accept*. The destination ecosystem's tools (`FindClusters`, `sc.tl.rank_genes_groups`,
`scran::modelGeneVar`) are a **third reader** that round-trip tests never exercise. This document is the
explicit contract for producing output those tools won't choke on — and `lstar convert --check` is the
gate that verifies it (open the result in its native library + run a canonical-ops smoke).

## Why it is deterministic

The "best translation" is **not** derived from first principles — it is a finite, explicit contract in
three layers:

1. **Role → slot is a function of the L★ type.** Once a field is typed `(role, state, span-shape)`, its
   canonical destination slot in each target is fixed. The typed model *is* the disambiguator — this is
   the work L★ does that a direct `as.Seurat()` cannot.
2. **A per-format convention table** for the residue that actually makes native tools choke (key
   suffixes, name munging, dimname alignment, dtype). Concrete and enumerable, below.
3. **Recorded defaults for genuine ambiguity.** Where the source under-determines the target, pick the
   canonical default *and record the choice in `provenance`* — so it is reproducible and inspectable, not
   silent.

What cannot be derived is **enumerated and validated against the native API** (the `--check` smoke).

## Layer 1 — role + state + span → canonical slot

| L★ field | AnnData | Seurat | SingleCellExperiment |
|---|---|---|---|
| measure `state=raw` (cells×genes) | `layers["counts"]` (X if primary) | `Assay@counts` (genes×cells) | assay `counts` |
| measure `state=lognorm` (the "primary") | `X` | `Assay@data` | assay `logcounts` |
| measure `state=scaled` (HVG subset) | `layers["scaled"]` | `Assay@scale.data` (subset) | — (recorded if uncoercible) |
| second feature space (ADT/HTO) | a modality / `obsm` (MuData) | a second `Assay` | `altExp` |
| embedding `[cells, <red>]` | `obsm["X_<red>"]` | `DimReduc "<red>"`, `cell.embeddings` | `reducedDim("<RED>")` |
| loading `[features, <red>]` | `varm["<RED>s"]` | `DimReduc@feature.loadings` | `attr(reducedDim, "rotation")` |
| relation `[cells, cells]` | `obsp` | `Graphs` / `Neighbor` | `colPair` |
| relation `[genes, genes]` | `varp` | (assay-level graph) | `rowPair` |
| label categorical `[cells]` | `obs` (categorical) | `meta.data` factor (+ active `Idents`) | `colData` factor |
| label `[genes]` | `var` | `@meta.features` | `rowData` |
| measure `[<red>]` (per-dim stdev) | `uns["<red>"]["variance"]` | `DimReduc@stdev` | (reducedDim attr) |
| aux passthrough (`uns`/`@misc`) | `uns` | `@misc` | `metadata` |

`<red>` ∈ {pca, umap, tsne, harmony, …}; `<RED>` is its upper-cased Seurat key stem (pca→`PC`).

## Layer 2 — the convention table (what makes native tools choke)

**Seurat**
- a `DimReduc` `Key` **must end in `_`** (`PC_`, `UMAP_`) — Seurat errors otherwise.
- feature names: Seurat replaces `_` with `-`; the profile records the original mapping so it round-trips.
- `meta.data` rownames **must equal** `colnames(object)` (the cell names).
- counts are **genes × cells** (the transpose of AnnData), a `dgCMatrix`.
- the active identity must be a **factor**.

**AnnData / scanpy**
- `obs_names` (and `var_names`) must be **unique** and string — scanpy/anndata index on them.
- `rank_genes_groups`' `groupby` column must be a **categorical** dtype.
- `obsm` embeddings are conventionally `X_`-prefixed (`X_pca`) for scanpy plotting to find them.
- `.X` is float; integer counts live in `layers["counts"]`.

**SingleCellExperiment / scran-scater**
- the normalized assay is named **`logcounts`** (scran/scater look for it by name).
- `reducedDim` and `altExp` row/col-names must align with the parent's `colnames`.

## Layer 3 — recorded defaults for ambiguity

When the source under-determines the target, the choice is fixed and recorded, never silent:

- **Which measure becomes the "primary" (`X` / `@data`)?** `lognorm` if present, else `raw`. The chosen
  source location is recorded in the field's `provenance` so the inverse conversion restores it. (A raw
  measure with no lognorm sibling thus lands in `layers["counts"]`, leaving `X = None` — faithful, and the
  reason `lstar convert` reports the matrix in `counts`.)
- **A facet measured on a cell subset** → typed **partial coverage** (an `index` into the span axis), not
  a padded dense matrix and not a separate axis.
- **Anything with no target slot** → `ds.dropped` (visible in the conversion report), never discarded.

## Validation — the native-acceptance check

`lstar convert --check` (default on) closes the loop: it opens the produced object in its native library
and runs a canonical-ops smoke — scanpy `normalize_total`/`log1p`/`pca`/`rank_genes_groups`; Seurat
`NormalizeData`/`FindVariableFeatures`/`ScaleData`/`RunPCA`; scran/scater `logNormCounts`/`modelGeneVar`/
`runPCA` — plus the Layer-2 invariants above. The heavy analysis libraries are optional: absent → the
check degrades to *open + structural invariants* and reports `ops skipped`, never a hard failure.
`--strict` turns a check failure into a non-zero exit. This is the contract being machine-verified: not
"did the bytes round-trip" but "will the native toolchain accept what we produced".

> Scope: the smoke runs a *representative* set of canonical ops, not every downstream tool. The deeper
> guard against *silent* mis-placement (a tool that runs but on the wrong matrix) is Layer 1 + the `state`
> typing — which is the point of routing through L★ rather than converting format-to-format directly.

## Both backends honor this contract

The same role→slot mapping governs whether a conversion runs through the native package or lstar's
**package-free** codec (`--backend direct`): the h5py h5ad reader/writer, the base-R Seurat reader/writer,
and the base-R SCE *reader* produce byte-for-value the same L★ ↔ native mapping as
`anndata`/`SeuratObject`/`SingleCellExperiment` would. The Seurat package-free *writer* materializes the
object from a **pinned SeuratObject schema** (a recent `Assay5`/`DimReduc`/`Seurat` `setClass`, lifted
verbatim from source) with the S4 class identity forged to `SeuratObject`, so a real SeuratObject session
reconstructs and accepts it — verified by the native-acceptance check on every CI run, which is exactly
what keeps the pinned schema from drifting.

**Asymmetry, by design.** Reading a serialized object packagelessly only needs base R (`readRDS` +
slot-walk); *manufacturing* one needs the destination's class machinery. For Seurat that machinery is a
flat, forgeable `setClass` schema. For **SCE the write side stays native-only**: a valid
`SingleCellExperiment` is a deep `SummarizedExperiment` + `GRanges` (rowRanges) + internal-`DataFrame`
hierarchy whose nested validity invariants make a forged object impractical to keep correct — so
`--backend direct` to an SCE target reports that `SingleCellExperiment` is required.
