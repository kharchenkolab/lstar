# lstar ⇄ SingleCellExperiment

Everything about moving data between Bioconductor's **SingleCellExperiment** (SCE) — and its
`SummarizedExperiment` (SE) base — and lstar's L★ model: how to read, write, and convert it, the exact
representation mapping, the SCE/SE variants the reader handles, what survives a conversion, and the
asymmetry that makes SCE a **read-anywhere, write-native-only** target. SCE embodies the same
shared-object principle L★ generalizes (one object many packages read and write), so the mapping is
direct.

**On this page:** [Quick start](#quick-start) · [Install](#requirements--install) ·
[Reading](#reading-sce-into-l) · [Writing](#writing-sce-from-l) · [Converting](#converting) ·
[Representation](#the-l-representation-the-sce-profile) ·
[Preserved vs dropped](#what-is-preserved-vs-dropped) · [Versions](#versions--variants-recognized) ·
[Native validity](#native-validity--the---check-smoke) · [Collections](#collections) ·
[Examples](#worked-examples) · [Troubleshooting](#troubleshooting) · [API](#api--references)

Entry points: **R** `read_sce(sce)` / `write_sce(ds)` (in-memory), `read_sce_backed(h5ad)`
(disk-backed), and the `lstar convert` CLI. SCE lives in R; converting from AnnData (Python) routes
through a `.lstar.zarr` store.

## Quick start

```r
library(lstar)
ds  <- read_sce(sce)              # SingleCellExperiment -> L* dataset
sce2 <- write_sce(ds)             # L* dataset -> SingleCellExperiment
ds$dropped                        # what L* had no field for — recorded, not lost
```

```bash
lstar convert sample.h5ad sample.rds --to sce   # AnnData -> SingleCellExperiment (.rds)
lstar convert sample.rds sample.h5ad            # SCE/Seurat .rds -> AnnData (class sniffed)
lstar inspect sample.rds                        # print the L* structure of an .rds, no write
```

## Requirements & install

| route | with the native package | package-free (`--backend direct`) |
|---|---|---|
| SCE `.rds` → L★ store (**read**) | R + `SingleCellExperiment` + `S4Vectors` | **base R + `lstar`** (S4 slot-walk; no SCE) |
| L★ store → SCE `.rds` (**write**) | R + `SingleCellExperiment` + `S4Vectors` | **— native only** (see below) |

```r
install.packages("/path/to/lstar/R", repos = NULL, type = "source")
# Bioconductor (only for the formats you touch):
# BiocManager::install(c("SingleCellExperiment", "S4Vectors"))
```

**The read/write asymmetry — read it before relying on `--backend direct`.** Reading a serialized SCE
packagelessly needs only base R (`readRDS` + an S4 slot-walk). *Manufacturing* a valid SCE does not:
it is a deep `SummarizedExperiment` + `GRanges` (`rowRanges`) + internal-`DataFrame` hierarchy whose
nested validity invariants make a forged object impractical to keep correct. So the **SCE *writer* is
native-only** — `--backend direct` to an SCE target stops with a message that `SingleCellExperiment` is
required. (Seurat, whose schema is flat and forgeable, *can* be written package-free; SCE cannot.) The
analysis stack (`scran`/`scater`) is needed only for the optional `--check` smoke, not for conversion.

## Reading SCE into L★

`read_sce(sce)` turns a live `SingleCellExperiment` (or a plain `SummarizedExperiment`) into an
`lstar_dataset`. Assays become measures, `reducedDims` become embeddings (with their `rotation` as
loadings), `colData`/`rowData` columns become cell/gene fields, `altExps` become a second feature axis,
and `colPairs`/`rowPairs` become relations.

```r
ds <- read_sce(sce)
ds$profiles    # e.g. c("singlecellexperiment@0.1", "SingleCellExperiment@1.24.0")
ds$dropped     # e.g. "metadata/<name>", a nested colData column, an unmappable colPair
```

Cells with **`NULL` dimnames** (common in real SCEs keyed by a `Barcode` `colData` column rather than
`colnames`) get stable synthesized labels (the `Barcode` column if present, else `cell1…`), so the
cell axis is never empty.

`read_sce_backed(h5ad)` builds an SCE whose assay is an on-disk `HDF5Array` `DelayedMatrix` — the
matrix is never materialized (see [Converting](#converting)).

## Writing SCE from L★

`write_sce(ds)` materializes a `SingleCellExperiment` (requires `SingleCellExperiment` + `S4Vectors`):
measures become assays (transposed to **genes × cells**), embeddings become `reducedDims` (with a
matching `<coord>_loadings` field re-attached as the `rotation` attribute), cell/gene fields become
`colData`/`rowData`, a second feature axis becomes an `altExp`, and `relation` fields become
`colPairs`/`rowPairs`. The normalized assay is named **`logcounts`** so scran/scater find it by name.

## Converting

**The one command.** A `.rds` target defaults to Seurat; pass `--to sce` for SingleCellExperiment:

```bash
lstar convert x.h5ad y.rds --to sce     # AnnData -> SCE
lstar convert sce.rds out.h5ad          # SCE -> AnnData (the .rds class is sniffed: Seurat vs SCE)
lstar convert sce.rds out.lstar.zarr    # SCE -> L* store (keeps everything)
```
`--report` shows fields kept + `dropped`; `--check` (default on; `--strict`) opens the result in
scran/scater; `--backend` selects native vs direct (recall: SCE **write** is native-only).

**Same language (R, in memory).** SCE ↔ Seurat is two functions:

```r
sce <- write_sce(read_seurat(seurat_obj))   # Seurat -> SCE
so2 <- write_seurat(read_sce(sce))          # ... and back
```

**Across languages (`.h5ad` ↔ SCE).** The `.lstar.zarr` store is the bridge — Python writes it, R reads
it:

```r
sce <- write_sce(lstar_read("pbmc.lstar.zarr")); saveRDS(sce, "pbmc_sce.rds")
```

**Disk-backed (bounded memory).** A streamed `.h5ad` (from `lstar.convert_to_h5ad` in Python) opens as
an SCE whose assay is an on-disk `HDF5Array` `DelayedMatrix`:

```r
sce <- read_sce_backed("atlas.h5ad")   # the assay is a DelayedMatrix; the matrix is never read in
```
`HDF5Array` is an optional `Suggests`; absent → a clear error.

## The L★ representation (the SCE profile)

The mapping is bidirectional and field-preserving — the inverse view of the combined table in
[`../mapping.md`](../mapping.md).

**Identity → axes:**

| SCE | L★ axis | origin / role |
|---|---|---|
| `colnames(sce)` (or the `Barcode` colData column) | `cells` | observed / observation |
| `rownames(sce)` | `genes` | observed / feature |

**Field rules:**

| SCE slot | L★ field | role · span · state/subtype |
|---|---|---|
| assay `counts` | `counts` | measure · (cells, genes) · raw |
| assay `logcounts` | `X` | measure · (cells, genes) · lognorm |
| any other assay | its name | measure · (cells, genes) · (state unset) |
| `colData[k]` | `k` | numeric→`measure`, factor→`label` (+ induced factor axis) · (cells) |
| `rowData[k]` | `k` | numeric→`measure`, factor→`label` · (genes) |
| `reducedDim("<RED>")` | `<red>` (lower-cased, e.g. `PCA`→`pca`) | embedding · (cells, `<red>` coord axis) |
| `attr(reducedDim, "rotation")` | `<red>_loadings` | loading · (genes, `<red>`) |
| `altExp(name)` (e.g. ADT/spike-ins) | a feature axis `name`; assays → `name.<assay>` | measure · (cells, `name`) |
| `colPair(name)` | `colpair_name` | relation · (cells, cells) |
| `rowPair(name)` | `rowpair_name` | relation · (genes, genes) |
| `metadata(sce)` | — | recorded in `dropped` (free-form study-level list, not typed) |

`state` typing keeps `counts` (raw) and `logcounts` (lognorm) distinct, so a conversion places each in
the right slot of the target (Seurat `counts`/`data`, AnnData `layers['counts']`/`X`). PCA scores and
their loadings share one `pca` coordinate axis, so the `rotation` survives a round-trip. An `altExp`
is the SCE analogue of a Seurat second assay / an AnnData (MuData) modality — it becomes a **second
feature axis**, not a layer.

## What is preserved vs dropped

The shared-vocabulary core — `counts`/`logcounts`/scaled assays, `reducedDims` (scores + rotation),
`colData`/`rowData`, clusterings, `altExps` (a second modality), and `colPairs`/`rowPairs` graphs —
survives every conversion among SCE, Seurat, and AnnData. What has no L★ field is recorded in
`ds$dropped`, never silently lost:

- **`metadata(sce)`** — the free-form study-level list — is not typed; its names are recorded.
- a **nested `colData`/`rowData` column** (a `DataFrame`, `GRanges`, or list inside the table) can't be
  coerced to a vector field and is recorded (e.g. `colData/<col> (DataFrame)`). `Rle` (run-length)
  columns *are* unpacked and kept.
- an `altExp` whose name collides with an existing axis, or a `colPair`/`rowPair` whose dimensions
  don't match, is recorded rather than mis-placed.

**Keep the `.lstar.zarr` store and nothing is dropped** — the manifest describes only what a *native
target* couldn't carry.

## Versions & variants recognized

- **Full SCE vs. plain `SummarizedExperiment`.** `reducedDims`, `altExps`, and `colPairs`/`rowPairs`
  are SCE-only; the reader guards those accessors so a plain `SummarizedExperiment` (or an odd SCE
  subclass such as `ReprocessedFluidigmData`) reads without crashing — you get assays + col/rowData.
- **`Rle` columns** in `colData`/`rowData` are unpacked to plain vectors.
- **`NULL` dimnames** — cells keyed by a `Barcode` `colData` column get synthesized labels.
- **`SingleCellExperiment` version** is recorded as `SingleCellExperiment@<version>` in `ds$profiles`.
- The package-free reader additionally locates gene names from a `GRangesList` `rowRanges`
  (`rowRanges@partitioning@NAMES`) when the assay has no rownames.

## Native validity & the `--check` smoke

A round-trip preserves L★'s representation; it does not guarantee scran/scater accept the object.
`lstar convert --check` (default on) opens the produced SCE and runs a canonical-ops smoke —
`scater::logNormCounts` / `scran::modelGeneVar` / `scater::runPCA` — plus the SCE invariants:

- the normalized assay is named **`logcounts`** (scran/scater look for it by name).
- `reducedDim` and `altExp` row/col-names align with the parent's `colnames`.

`scran`/`scater` absent → the check degrades to *open + structural invariants*; `--strict` makes a
failure non-zero.

## Collections

SCE represents one experiment (one `cells × genes` matrix plus aligned slots), not a heterogeneous
multi-sample collection. lstar has **no SCE collection profile**: a multi-sample study is kept as an
L★ collection or carried in formats that represent one natively (Conos; a Seurat v5 split assay; a
flattened AnnData — see [`../conversions.md`](../conversions.md#collections-convert-too)).
`MultiAssayExperiment` (the Bioconductor multi-assay container) and `SpatialExperiment`
(`spatialCoords`) are **out of current scope** — `SpatialExperiment` support is planned; spatial
coordinates already round-trip through the AnnData/Seurat profiles as a `spatial` observed coordinate
axis.

## Worked examples

A worked SCE read/write is [§8 of `../examples.md`](../examples.md); the Seurat ↔ SCE in-memory
conversion is in [`../conversions.md`](../conversions.md) §"Same-language conversions"; and the
cross-language chain through SCE (AnnData → Seurat → SCE → … → AnnData, verifying the core survives) is
[`../../examples/roundtrip_xlang.sh`](../../examples/roundtrip_xlang.sh). SCE version handling and the
package-free read path are cross-validated in `conformance/sce_versions.sh`.

## Troubleshooting

- **`--backend direct` won't write SCE.** By design — a valid SCE can't be forged without
  `SingleCellExperiment`. Install it (and `S4Vectors`), or write to a `.lstar.zarr` store / Seurat
  instead.
- **scran/scater can't find normalized values.** The normalized assay must be named `logcounts`;
  `write_sce` does this, but if you renamed it, restore the name.
- **A `colData` column is gone after conversion.** It was a nested (`DataFrame`/`GRanges`/list) column —
  check `ds$dropped`. `Rle` columns are kept; nested ones can't be coerced to a vector field.
- **Cells came back as `cell1, cell2, …`.** The source SCE had `NULL` `colnames` and no `Barcode`
  column; provide cell names (or a `Barcode` colData column) for stable labels.

## API & references

- **R**: `read_sce(sce)`, `write_sce(ds)`, `read_sce_backed(h5ad, layer=, assay_name=)`;
  `lstar_read`/`lstar_write`.
- The deterministic role→slot contract: [`../mapping.md`](../mapping.md). The L★ model:
  [`../model.md`](../model.md). The conversion overview: [`../conversions.md`](../conversions.md).
  (SCE has no entry in the proposal's profile catalog yet — this page is the normative SCE profile
  reference, source-verified against `R/R/profile_sce.R`.)
