# lstar ⇄ Seurat

Everything about moving data between **Seurat** objects (the R/Seurat in-memory object, serialized as
`.rds`) and lstar's L★ model: how to read, write, and convert it, the exact representation mapping, the
many object versions the reader recognizes (legacy v2 → v5), what survives a conversion, and the
conventions Seurat's own tools require on write-back. Seurat maps onto L★ unusually closely — a
`DimReduc` natively bundles an embedding with its loadings (L★'s shared coordinate axis), and
`@commands` is a provenance log.

**On this page:** [Quick start](#quick-start) · [Install](#requirements--install) ·
[Reading](#reading-seurat-into-l) · [Writing](#writing-seurat-from-l) · [Converting](#converting) ·
[Representation](#the-l-representation-the-seurat-profile) ·
[Preserved vs dropped](#what-is-preserved-vs-dropped) · [Versions](#versions--variants-recognized) ·
[Native validity](#native-validity--the---check-smoke) · [Collections](#collections) ·
[Examples](#worked-examples) · [Troubleshooting](#troubleshooting) · [API](#api--references)

Entry points: **R** `read_seurat(so)` / `write_seurat(ds)` (in-memory), `read_seurat_backed(h5ad)`
(disk-backed), and the `lstar convert` CLI. Seurat lives in R; to convert from AnnData (Python) the
conversion routes through a `.lstar.zarr` store.

## Quick start

```r
library(lstar)
ds  <- read_seurat(seurat_obj)              # Seurat object -> L* dataset
so2 <- write_seurat(ds)                      # L* dataset -> Seurat object (v5)
ds$dropped                                   # what L* had no field for — recorded, not lost
```

```bash
lstar convert sample.rds sample.h5ad --report   # Seurat -> AnnData + a fidelity report
lstar convert sample.h5ad sample.rds            # AnnData -> Seurat (.rds), bridged through R
lstar inspect sample.rds                        # print the L* structure of an .rds, no write
```

## Requirements & install

| route | with the native package | package-free (`--backend direct`) |
|---|---|---|
| Seurat `.rds` → L★ store (**read**) | R + `SeuratObject` | **base R + `lstar`** (S4 slot-walk; no SeuratObject) |
| L★ store → Seurat `.rds` (**write**) | R + `SeuratObject` | **base R + `lstar`** (materializes from a pinned SeuratObject schema) |

```r
install.packages("/path/to/lstar/R", repos = NULL, type = "source")   # the R package
```

Seurat/SCE conversion legs of the CLI run a short `Rscript` bridge; set `LSTAR_RLIB` / `LSTAR_RSCRIPT`
if `lstar` isn't on the default R library path. Full Seurat is **not** needed for conversion — only the
optional `--check` smoke uses `Seurat`.

## Reading Seurat into L★

`read_seurat(so, assay = SeuratObject::DefaultAssay(so))` turns a live Seurat object into an
`lstar_dataset`. It first **detects the object's version** and dispatches accordingly (legacy v2
through v5 — see [Versions](#versions--variants-recognized)). Each assay's layers, reductions, graphs,
metadata, identities, and (for spatial objects) image coordinates are typed into L★ fields.

```r
ds <- read_seurat(seurat_obj)
ds$profiles      # e.g. c("seurat@0.1", "SeuratObject@5.0.2") — the detected variant
ds$dropped       # slots with no L* field (e.g. a ChromatinAssay fragment file)
```

`read_seurat_backed(h5ad)` builds a Seurat **v5** object whose counts live on disk as a BPCells
`IterableMatrix` — the matrix is never materialized (see [Converting](#converting)).

## Writing Seurat from L★

`write_seurat(ds)` materializes a Seurat **v5** object, returning each field to its native slot and
enforcing the conventions Seurat's own functions require (below). Counts are written **genes × cells**
(a `dgCMatrix`, the transpose of AnnData). Fields with no Seurat slot are written to `@misc` and listed
in `ds$dropped` (per-sample PCA from a collection, for instance, goes to `so@misc$lstar_dropped`). A
collection (per-sample axes) is materialized as a **split v5 assay** (see [Collections](#collections)).

## Converting

**The one command.** `lstar convert` detects both endpoints and routes through the L★ store:

```bash
lstar convert s.rds s.h5ad              # Seurat -> AnnData
lstar convert s.rds s.lstar.zarr        # Seurat -> L* store (keeps everything)
lstar convert x.h5ad y.rds --to sce     # land the .rds as SingleCellExperiment instead of Seurat
```
`--report` shows fields kept + `dropped`; `--check` (default on; `--strict` to gate) opens the result
in Seurat and runs a canonical-ops smoke; `--backend auto|native|direct` picks the codec.

**Same language (R, in memory).** Seurat ↔ SingleCellExperiment is two functions:

```r
sce <- write_sce(read_seurat(seurat_obj))   # Seurat -> SCE
so2 <- write_seurat(read_sce(sce))          # ... and back
```

**Across languages (`.h5ad` ↔ Seurat).** The `.lstar.zarr` store is the bridge — Python writes it, R
reads it (or vice versa), with no `reticulate`:

```r
so <- write_seurat(lstar_read("pbmc.lstar.zarr")); saveRDS(so, "pbmc.rds")
```

**Disk-backed (bounded memory).** A streamed `.h5ad` (from `lstar.convert_to_h5ad` in Python) opens as
a disk-backed Seurat v5 — counts stay on disk as BPCells, never materialized:

```r
so <- read_seurat_backed("atlas.h5ad")   # Seurat v5 ops stream off the on-disk IterableMatrix
```
On a 40,220 × 20,138 atlas this peaks at ~0.6 GB RSS (mostly loading R + packages) vs ~7 GB to read
the matrix in. `BPCells` is an optional `Suggests`; absent → a clear error.

## The L★ representation (the Seurat profile)

The mapping is bidirectional and field-preserving — the inverse view of the combined table in
[`../mapping.md`](../mapping.md).

**Identity → axes.** Cells are shared across assays; each assay contributes its own feature axis:

| Seurat | L★ axis | origin / role |
|---|---|---|
| `colnames(object)` | `cells` | observed / observation |
| `rownames(assay)` | `<assay>` (e.g. `RNA`, `ADT`) — one feature axis **per assay** | observed / feature |

**Field rules.** Per assay, the three layers map by `state`; reductions bundle three fields on a shared
coordinate axis:

| Seurat slot | L★ field | role · span · state/subtype |
|---|---|---|
| `assay$counts` | `counts` | measure · (cells, `<assay>`) · raw |
| `assay$data` | `X` | measure · (cells, `<assay>`) · lognorm |
| `assay$scale.data` (usually HVGs only) | `scale.data` | measure · (cells, `<assay>`) · scaled — **partial coverage** when over a feature subset |
| `obj@meta.data[k]` | `k` | by dtype: categorical→`label`, numeric→`measure` · (cells) |
| `Idents(obj)` | `ident` | label · (cells) · `active_ident` (active-vs-stored distinction preserved) |
| `assay@meta.features[k]` | `k` | by dtype · (`<assay>`) |
| `obj@reductions[k]` (`DimReduc`) | `k` + `k_loadings` + `k_stdev` | embedding (cells, `k`) · loading (`<assay>`, `k`) · measure (`k`) |
| `obj@graphs[k]` / `obj@neighbors[k]` | `k` / `nn_k` | relation · (cells, cells) · subtype guessed |
| a 2nd assay (ADT/ATAC) | a 2nd feature axis (multimodal) | measures over (cells, `<assay2>`) |
| `obj@images[fov]` (Visium/Slide-seq/FOV) | `spatial` | embedding · (cells, **observed** coord axis) · subtype=spatial (pixels → `dropped`) |
| `obj@commands` | — | recorded as `provenance` on the produced fields |
| `obj@misc[k]` | — | opaque catch-all → `dropped` |

The `DimReduc` row is the payoff: a direct `as.Seurat()`/`as.AnnData()` typically drops gene loadings,
but L★ keeps scores and loadings on one shared coordinate axis, so they ride through every conversion.
A **Signac `ChromatinAssay`** (scATAC) keeps its defining peak **genomic ranges** (as `seqnames` /
`start` / `end` feature fields, `subtype=genomic_pos`); its external fragment file has no cross-format
home and is recorded in `dropped`.

## What is preserved vs dropped

The shared-vocabulary core — raw/normalized/scaled expression, reductions (scores + loadings + stdev),
`meta.data`/`meta.features`, identities/clusterings, multimodal assays, graphs, and spatial centroids —
survives every conversion among Seurat, AnnData, and SCE. What has no slot — arity-3 tensors, trees,
fitted models, a `ChromatinAssay` fragment file, image pixels — is written to `@misc` and **listed in
`ds$dropped`**, never silently discarded. **Keep the L★ store and nothing drops** — `dropped` describes
only what a native target couldn't carry.

## Versions & variants recognized

The reader detects the object's vintage and adapts, so a collaborator on a different Seurat release
doesn't break the conversion:

- **Legacy v2** — the lowercase `seurat` S4 class that predates `Assay`/`SeuratObject` entirely — is
  detected and read through a dedicated path that pulls its fixed slots (`raw.data`/`data`/`scale.data`,
  the `dr` list of `dim.reduction`s, `meta.data`, `ident`, `snn`, the multimodal `assay` list) via
  `attr()`, so the ancient class need not be defined in the running R. Write-back emits a modern object
  (old→new is the point).
- **v3 / v4** — `Assay` with fixed `counts`/`data`/`scale.data` slots.
- **v5** — `Assay5` with named layers; a fallback covers `SeuratObject < 5`.
- **Per-assay subtype** is recorded: `SCTAssay` (SCTransform residuals), `ChromatinAssay` (Signac
  scATAC: peaks + genomic ranges).
- **A split v5 assay** (`split(assay, f = sample)`) is recognized as a **collection** (see below).
- **`scale.data` over the variable features only** is kept as typed **partial coverage**, not padded.
- **Spatial coordinates** live in `@images`, *not* `Reductions` — captured as a `spatial` observed
  coordinate axis (multi-section subsets use partial coverage); pixel images go to `dropped`.

The detected `<format>@<version>` is recorded in `ds$profiles`.

## Native validity & the `--check` smoke

A round-trip through L★ preserves L★'s representation; it does not guarantee Seurat's own tools accept
the object. `lstar convert --check` (default on) opens the produced Seurat object and runs a
canonical-ops smoke — `NormalizeData` / `FindVariableFeatures` / `ScaleData` / `RunPCA` — plus the
Seurat invariants the profile enforces on write-back:

- a `DimReduc` `Key` **must end in `_`** (`PC_`, `UMAP_`) — Seurat errors otherwise.
- feature names: Seurat replaces `_` with `-`; the profile records the original mapping so it
  round-trips.
- `meta.data` rownames **must equal** `colnames(object)` (the cell names).
- counts are **genes × cells**, a `dgCMatrix`.
- the active identity must be a **factor**.

`Seurat` absent → the check degrades to *open + structural invariants*; `--strict` makes a failure
non-zero.

## Collections

A multi-sample study is a **collection**, not one matrix. A **split v5 assay**
(`split(assay, f = sample)`) is exactly that — layers named `<root>.<sample>`, each covering only its
sample's cells — and `read_seurat` reconstructs it as an L★ collection (per-sample `cells.<s>` axes).
Conversely, **Seurat v5 preserves a collection** that other formats would flatten: a Conos integration
writes to a split v5 assay with no fabricated corrected matrix —

```r
so <- write_seurat(write_conos(con))   # -> Seurat v5 SPLIT assay: per-sample raw layers over the UNION
#   genes, the joint graph as Graphs(), the embedding as a DimReduc, clusters/sample in meta.data.
#   read_seurat(so) reads it straight back as an L* collection; per-sample PCA -> so@misc$lstar_dropped.
```
(AnnData, by contrast, flattens a collection to one union matrix — keep the L★ store for the full
per-sample structure.)

## Worked examples

Runnable scripts in [`../../examples/`](../../examples): `convert_h5ad_to_seurat.sh` (a complete
`.h5ad → Seurat` converter). The split-v5 collection round-trip is example
[7 in `../examples.md`](../examples.md). Conos↔Seurat is covered by `conformance/conos.sh`.

## Troubleshooting

- **`DimReduc` key error in Seurat after a conversion.** L★ enforces the `_`-terminated key on
  write-back; if you hand-build a `DimReduc`, end its `Key` in `_`.
- **Feature names changed `_` → `-`.** Seurat munges underscores in feature names; the profile records
  the original so a round-trip through L★ restores it.
- **`scale.data` covers fewer genes than `counts`.** Expected — Seurat keeps `scale.data` over the
  variable features only; L★ types it as partial coverage, not an error.
- **Something is missing.** Check `ds$dropped` / `so@misc$lstar_dropped`. Convert to a `.lstar.zarr`
  store to keep everything.
- **`--backend direct` write produced an object an old R rejects.** The package-free writer materializes
  from a pinned recent SeuratObject schema; read it in a current SeuratObject session (or use
  `--backend native`).

## API & references

- **R**: `read_seurat(so, assay=)`, `write_seurat(ds)`, `read_seurat_backed(h5ad, group=, assay=)`;
  `lstar_read`/`lstar_write`. **Collections**: `write_conos(con)` / `read_conos(ds)`.
- The deterministic role→slot contract: [`../mapping.md`](../mapping.md). The L★ model:
  [`../model.md`](../model.md). The conversion overview: [`../conversions.md`](../conversions.md).
  Normative profile rules: [`../../misc/Lstar_proposal.md`](../../misc/Lstar_proposal.md) Appendix B.3.
