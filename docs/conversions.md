# Converting between formats

The long-term aim of L★ is that tools read and write the `.lstar.zarr` format directly. In the short
term, the most immediately useful thing lstar does is act as **glue between the formats people already
use** — AnnData, Seurat, SingleCellExperiment, pagoda2, and Conos — across Python and R. This page is
the practical guide to that.

## The idea in one picture

Every supported format has a **profile**: a precise, bidirectional mapping between that format's native
object and the L★ model. A *reader* turns a native object into L★ fields; a *writer* turns L★ fields
back into a native object. Because L★ sits in the middle, you convert format **X** to format **Y** by
composing them:

```
   AnnData ──read_anndata──▶               ──write_seurat──▶ Seurat
   Seurat  ──read_seurat ──▶   L★ dataset  ──write_sce────▶ SingleCellExperiment
   SCE     ──read_sce    ──▶  (axes+fields) ──write_anndata▶ AnnData
   Conos   ──write_conos ──▶               ──lstar.write──▶ .lstar.zarr
```

`convert(X → Y) = write_Y(read_X(object))`. The L★ dataset in the middle is an ordinary value (or an
on-disk `.lstar.zarr` store), so the same conversion works in memory within one language, or across
languages by writing the store from one and reading it from the other.

**Why route through a model at all, instead of writing X→Y converters directly?** Two reasons, and
they are the selling point:

1. **A shared vocabulary makes the conversion *meaningful*, not just mechanical.** Every profile maps
   the same biological object to the same L★ signature — log-normalized expression becomes
   `data : measure, state=lognorm` whether it came from AnnData's `X`, a Seurat `data` layer, or
   pagoda2's `getExpressionBlock()`. Principal components become a `pca` embedding *together with* their
   gene `pca_loadings` sharing one coordinate axis. So a conversion preserves the *meaning* of each
   piece, and things a one-off converter typically drops — like PCA loadings on an AnnData→Seurat
   conversion — survive (see [What is preserved](#what-is-preserved-and-what-is-not)).
2. **Loss is reported, never silent.** Anything a target format has no place for is written to that
   format's sidecar (`uns` / `Misc` / `metadata`) and recorded in the dataset's `dropped` list — so you
   can see exactly what a conversion could not carry, instead of discovering it missing later.

## One command: the `lstar convert` CLI

For the common cases you don't need to write any glue — the CLI detects each format from its path,
routes through the L★ store (in-process for Python formats, via a short `Rscript` bridge for Seurat/SCE),
and reports what crossed:

```bash
lstar convert pbmc.h5ad pbmc.rds              # AnnData -> Seurat (.rds), bridged automatically
lstar convert atlas.h5ad atlas.lstar.zarr     # AnnData -> L★ store
lstar convert s.rds s.h5ad --report           # Seurat -> AnnData + the full fidelity report
lstar convert x.h5ad y.rds --to sce           # a .rds target as SingleCellExperiment (default is Seurat)
lstar inspect data.h5ad                        # read + report its L★ structure, no write
```

Two things make it more than a one-liner:

- **A fidelity report** (`--report`, `--report-json FILE`) lists every axis/field with its role, state,
  span, and `provenance`, and — crucially — **`dropped`**: what the target could not represent, made
  visible rather than silently lost.
- **A native-acceptance check** (`--check`, on by default; `--strict` to gate the exit code) opens the
  produced object in its *own* library and runs a canonical-ops smoke (scanpy / Seurat / scran-scater), so
  you know the native analysis tools will accept it — not just that the bytes round-tripped. The heavy
  analysis libraries are optional; absent → the check degrades to open + structural invariants. What lands
  where (and why it is deterministic) is the explicit contract in [`mapping.md`](mapping.md).

When the target is an L★ store, the write layout is tunable: `--zarr-format {v3,v2}` (v3 default),
`--compression {none,gzip,zstd}` (+ `--compression-level`), `--chunk-elems`, and `--shard-elems` (v3
sharding — many chunks into fewer, still byte-range-readable objects). `--viewer` additionally precomputes
the `viewer@0.1` navigators and writes them with the viewer's per-field compression layout (`--compress-primary`
also compresses the gene-major counts); see [`format.md`](format.md).

Format detection: `.h5ad`→AnnData, `.h5mu`→MuData, `.rds`→Seurat/SCE (sniffed), `.lstar.zarr`/`.zarr`→
store. Override with `--from`/`--to`. Seurat/SCE legs need R with the `lstar` package on the path (set
`LSTAR_RLIB`/`LSTAR_RSCRIPT` if it isn't on the default library path).

### What you need installed (and the package-free fallback)

A conversion needs the *format* library for each endpoint it touches — **not** the analysis packages
(scanpy / full Seurat / scran are only used by the optional `--check`). And where lstar has its own codec,
the format package itself is optional too: `--backend auto` (default) uses the native package when it's
importable and falls back to lstar's **package-free** codec otherwise.

| convert | with the native package | package-free fallback (`--backend direct`) |
|---|---|---|
| `.h5ad` ↔ store | `anndata` | **`h5py`** (lstar reads/writes the HDF5 encoding directly) |
| Seurat `.rds` ↔ store | R + `SeuratObject` | **base R + the `lstar` R package** (reads via S4 slot-walk; writes a pinned-schema object) |
| SCE `.rds` → store (read) | R + `SingleCellExperiment` | **base R + the `lstar` R package** (S4 slot-walk) |
| store → SCE `.rds` (write) | R + `SingleCellExperiment` | — (native only: a valid SCE needs the SummarizedExperiment / GRanges machinery) |
| `.h5mu` ↔ store | `mudata` | — (native only for now) |
| store ↔ store | — | — |

So `pip install lstar-sc h5py` converts `.h5ad` with no anndata, and a Seurat `.rds` reads/writes with just
base R + lstar (no SeuratObject). Force a path with `--backend native|direct`. When the package-free path
hits something only the native package can handle — an unrecognized on-disk version, or an external-pointer
/ `BPCells`-backed matrix — it stops with a clear message naming the package to install (lstar then uses it
automatically). The deterministic role→slot contract both paths honor is [`mapping.md`](mapping.md).

> The package-free codecs are kept honest by CI: every Tier-A path is cross-validated against the native
> package (forced `--backend direct`, then read back and value-compared) in `conformance/convert_cli.sh`.

## The functions

You only need a handful. The reader/writer for each format:

| Format | Read into L★ | Write from L★ | Language | Package |
|---|---|---|---|---|
| **AnnData** (`.h5ad`/`.zarr`) | `read_anndata(adata)` | `write_anndata(ds)` | Python | `lstar.profiles.anndata` |
| **Seurat** (legacy v2 → v5) | `read_seurat(so)` | `write_seurat(ds)` | R | `lstar` |
| **SingleCellExperiment** | `read_sce(sce)` | `write_sce(ds)` | R | `lstar` |
| **Conos** (a collection) | `write_conos(co)` → L★ | `read_conos(ds)` → Conos | R | `lstar` |
| **pagoda2** | *(via Conos members; standalone reader planned)* | | R | `lstar` |
| **L★ store** | `lstar.read(path)` / `lstar_read(path)` | `lstar.write(ds, path)` / `lstar_write(ds, path)` | Py / R | `lstar` |

The L★ dataset they hand back is the universal intermediate: in Python a `lstar.Dataset`, in R an
`lstar_dataset` (a list of `$axes` and `$fields`). Either can be written to a `.lstar.zarr` store and
read by any of the three languages.

> **Per-format guides.** For a complete, format-expert reference on each endpoint — every native slot
> and its L★ signature, version handling, the native conventions the format's own tools require, and
> troubleshooting — see [`formats/anndata.md`](formats/anndata.md),
> [`formats/seurat.md`](formats/seurat.md), and
> [`formats/singlecellexperiment.md`](formats/singlecellexperiment.md). This page is the cross-format
> overview; those are the deep dives.

## Same-language conversions (in memory)

When both formats live in the same language, conversion is a one-liner.

**R: Seurat → SingleCellExperiment** (and back), preserving expression, embeddings, and metadata:

```r
library(lstar)
ds  <- read_seurat(seurat_obj)   # Seurat object  -> L* (assays->measures, DimReducs->embeddings, meta.data->fields)
sce <- write_sce(ds)             # L*             -> SingleCellExperiment (measures->assays, embeddings->reducedDims)
# ...and the reverse direction is just the other two functions:
ds2 <- read_sce(sce)
so2 <- write_seurat(ds2)
```

**Python: AnnData ↔ AnnData through L★** (a faithful round-trip; useful to confirm nothing is lost, and
as the in-memory half of a cross-language conversion):

```python
from lstar.profiles.anndata import read_anndata, write_anndata
ds     = read_anndata(adata)     # AnnData -> L*: X/layers -> measures, obsm -> embeddings, obsp -> relations, .raw kept
adata2 = write_anndata(ds)       # L* -> AnnData: each field written back to the slot it came from
print(ds.dropped)                # e.g. ['uns/neighbors'] -- what L* had no field for (recorded, not lost)
```

## Cross-language conversions (the common case: h5ad ↔ Seurat)

AnnData lives in Python and Seurat in R, so a conversion crosses the language boundary. The
`.lstar.zarr` **store on disk is the bridge** — Python writes it, R reads it (or vice versa). No shared
memory, no `reticulate`, no re-implementation of either format's reader in the other language.

```bash
# (1) Python: AnnData -> L* store
python3 -c '
import anndata as ad, lstar
from lstar.profiles.anndata import read_anndata
lstar.write(read_anndata(ad.read_h5ad("pbmc.h5ad")), "pbmc.lstar.zarr")'

# (2) R: L* store -> Seurat object, saved as an .rds
Rscript -e '
library(lstar)
so <- write_seurat(lstar_read("pbmc.lstar.zarr"))
saveRDS(so, "pbmc.rds")'
```

That is a complete `h5ad → Seurat` converter. The runnable, commented version (which also reports what
survived) is [`examples/convert_h5ad_to_seurat.sh`](../examples/convert_h5ad_to_seurat.sh); a
longer chain that goes AnnData → Seurat → SCE → … → AnnData and verifies the core survives the whole
trip is [`examples/roundtrip_xlang.sh`](../examples/roundtrip_xlang.sh).

## Low-memory conversion (streaming)

The conversions above load the source matrix into memory. For a large atlas that can be many
gigabytes. lstar offers a **bounded-memory** path that, instead, reads the source matrices straight
from disk and streams them block-by-block — so peak memory is ~one block, not the whole matrix.
This works in **both directions** between an `.h5ad` file and an L\* store, because both are on-disk
formats.

**Into an L\* store** (`h5ad → .lstar.zarr`): read the `.h5ad` in backed mode (its matrices stay on
disk) and stream them out.

```python
import lstar
lstar.convert_anndata("atlas.h5ad", "atlas.lstar.zarr")   # peak memory ~ one block, not the matrix
#   equivalently, explicit form:
#   import anndata as ad
#   from lstar.profiles.anndata import read_anndata
#   lstar.write(read_anndata(ad.read_h5ad("atlas.h5ad", backed="r")), "atlas.lstar.zarr", stream=True)
```

**Out to an `.h5ad`** (`.lstar.zarr → h5ad`): read the store lazily (its measures stay on disk as
`LazyCSX` proxies) and stream each one into the h5ad's sparse groups. The small parts (`obs`/`var`/
`obsm`/graphs) go through anndata so their on-disk encoding is reused, not re-implemented.

```python
lstar.convert_to_h5ad("atlas.lstar.zarr", "atlas.h5ad")   # peak memory ~ one block
#   equivalently, on any L* dataset whose measures are lazy/backed sources:
#   from lstar.profiles.anndata import write_anndata_streamed
#   write_anndata_streamed(lstar.read("atlas.lstar.zarr", lazy=True), "atlas.h5ad")
```

On a real 40,220 × 20,138 dataset (77.6M nonzeros) the `h5ad → L\*` direction converts in **~110 MB**
of peak memory instead of **~1.3 GB** eager (≈12× less); a 40k × 5k store streams **out** to `.h5ad`
in **~32 MB** instead of ~120 MB (and faster, since nothing is fully materialized first). Each
measure keeps its native orientation (a CSR `X` stays CSR, a CSC layer stays CSC — no transpose),
and the result is byte-identical to the eager conversion. The same streaming applies to **L\* → L\***
recompression/re-chunking: `lstar.write(lstar.read(src, lazy=True), dst, stream=True)` rewrites a
store without ever holding it whole.

### Disk-backed native targets (bounded all the way to a usable object)

The streamed `.h5ad` above is itself a disk-backed representation, so the bounded path runs *all the
way into a native object* — you never hold the matrix, even at the destination:

| Target | Open the streamed `.h5ad` as | Backed by |
|---|---|---|
| **AnnData** (Python) | `anndata.read_h5ad(path, backed="r")` | h5ad on disk (`X` stays a `SparseDataset`) |
| **Seurat v5** (R) | `read_seurat_backed(path)` | **BPCells** on-disk `IterableMatrix` |
| **SingleCellExperiment** (R) | `read_sce_backed(path)` | **HDF5Array** `DelayedMatrix` |

```r
library(lstar)
# after: python -c 'import lstar; lstar.convert_to_h5ad("atlas.lstar.zarr", "atlas.h5ad")'
so  <- read_seurat_backed("atlas.h5ad")   # counts live on disk (BPCells); Seurat v5 ops stream off it
sce <- read_sce_backed("atlas.h5ad")      # assay is an on-disk HDF5Array DelayedMatrix
```

On the 40,220 × 20,138 Marrow atlas (1.3 GB `.h5ad`), building either disk-backed object peaks at
**~0.6 GB** of RSS — most of which is just loading R and the heavy packages — versus ~7 GB to read
the matrix in. The matrix is never materialized; the object is fully usable, the data byte-for-byte
identical. `BPCells` and `HDF5Array` are *optional* (R `Suggests`): each reader errors with a clear
message if its package is absent, so lstar's required footprint stays small.

Boundary: bounded-memory conversion works whenever **both ends are on-disk formats** — `.h5ad` ↔
`.lstar.zarr`, and the disk-backed native objects above. Converting into a *fully in-memory* object
(a Seurat `dgCMatrix` in an `.rds`, an in-memory AnnData `X`) still, by definition, needs the matrix
in RAM at the destination.

## What is preserved, and what is not

A conversion preserves a field when **both** formats have a place for the same object. The set that all
of AnnData, Seurat, and SCE share — the **shared-vocabulary core** — survives every conversion among
them:

| Object | L★ signature | AnnData | Seurat | SCE |
|---|---|---|---|---|
| raw counts | `counts` : measure, `state=raw` | `layers['counts']` | `counts` layer | `counts` assay |
| normalized expression | `X`/`data` : measure, `state=lognorm` | `X` | `data` layer | `logcounts` assay |
| scaled expression | `scaled` : measure, `state=scaled` | a layer | `scale.data` layer | an assay |
| PCA scores **+ loadings** | `pca` embedding + `pca_loadings` loading (shared axis) | `obsm['X_pca']` + `varm['PCs']` | a `DimReduc` (incl. loadings) | a `reducedDim` (+ rotation) |
| UMAP / t-SNE | `umap` / `tsne` : embedding | `obsm['X_umap']` | a `DimReduc` | a `reducedDim` |
| cell / gene metadata | columns over `(cells)` / `(genes)` | `obs` / `var` | `meta.data` / `meta.features` | `colData` / `rowData` |
| clustering / cell type | a `label` over `(cells)` | an `obs` column | `Idents` / a `meta.data` column | a `colData` column |
| a second modality (ADT/ATAC) | a measure over `(cells, proteins/peaks)` — a second **feature axis** | a MuData modality | a second assay | an `altExp` |
| cell–cell / gene–gene graph | a `relation` over `(cells,cells)` / `(genes,genes)` | `obsp` / `varp` | a `Graph` / `Neighbor` | `colPairs` / `rowPairs` |
| spatial coordinates | `spatial` : embedding over an **observed** coordinate axis (`subtype=spatial`) | `obsm['spatial']` | `so@images` (Visium/FOV/Slide-seq centroids) | — (SpatialExperiment `spatialCoords`: planned) |

The PCA-loadings row is the concrete payoff: a direct AnnData→Seurat conversion usually discards the
gene loadings, but because L★ keeps the scores and loadings on one shared `pca` axis, they ride through.
The same shared-axis trick carries **joint integration** outputs (WNN / MOFA+ / totalVI): the factor
scores (an embedding) and per-modality loadings share one factor axis, and the modality weights ride as
cell measures. A **modality measured on only some cells** rides over the shared `cells` axis as
**partial coverage** (an `index`), not a padded matrix or a separate axis.

**What a given target cannot hold is reported in `dropped`:**

- Multimodal, graphs, and partial coverage now have native slots in the relevant formats (rows above),
  so they survive — they are **no longer dropped**. What genuinely has no slot — **arity-3 tensors,
  trees, and fitted models** — is written to the format's sidecar (`uns` / `@misc` / `metadata`) and
  listed in `dropped`. Format-specific extras with no cross-format home (e.g. a Signac
  `ChromatinAssay`'s external **fragment file**) are recorded too; its **peak ranges** *are* kept (as
  `seqnames`/`start`/`end` feature fields).
- Through **AnnData**: anything in `uns` that isn't a recognized structure is recorded in `dropped`
  (re-surfaced under `adata.uns['lstar/dropped']` on write-back).
- **Keep the L★ store and nothing is dropped at all** — the `dropped` manifest only describes what a
  *native target* couldn't carry.

The rule of thumb: **convert through L★ and keep the L★ store** if you want nothing lost; convert all
the way to a native format and you keep its shared-vocabulary core plus a recorded manifest of what it
couldn't carry. If you keep the data in L★, even the off-vocabulary pieces are preserved.

## Versions are recognized, not assumed

The readers detect the version/variant of the object they're given and adapt, rather than assuming one
layout — so conversions don't break when a collaborator is on a different release:

- **Seurat**: a **very old v2** object — the lowercase `seurat` S4 class that predates `Assay`/
  `SeuratObject` entirely — is detected (class `seurat`, not `Seurat`) and read through a dedicated path
  that pulls its fixed slots (`raw.data`/`data`/`scale.data`, the `dr` list of `dim.reduction`s,
  `meta.data`, `ident`, `snn`, multimodal `assay` list) via `attr()`, so the ancient class need not be
  defined in the running R; write-back emits a modern object (old→new is the point). Then v3/v4
  (`Assay`, fixed `counts`/`data`/`scale.data` slots) vs. v5 (`Assay5`, layers); a
  fallback covers SeuratObject < 5. Per-assay class is recorded — `SCTAssay` (SCTransform residuals),
  `ChromatinAssay` (Signac scATAC: peaks + genomic ranges). A **split v5 assay**
  (`split(assay, f = sample)`) is recognized as a **collection**. A `scale.data` over the variable
  features only is kept as **partial coverage**. **Spatial coordinates** (Visium/FOV/Slide-seq) live in
  `so@images`, *not* in `Reductions` — they are captured as a `spatial` observed coordinate axis
  (mirroring the AnnData `obsm['spatial']` path; multi-section subsets use partial coverage), with the
  pixel images recorded in `dropped`.
- **SingleCellExperiment**: full SCE vs. a plain `SummarizedExperiment` (SCE-only accessors are
  guarded); S4 `Rle` columns are unpacked; cells keyed by a `Barcode` colData column (NULL dimnames)
  get synthesized labels.
- **pagoda2**: the modern `getRawCounts()` accessor vs. the legacy `$counts` slot.
- **AnnData**: the library version, the on-disk sparse layout (legacy `h5sparse_format` < 0.7 and modern
  `encoding-type`), and the `.raw` slot (kept on its own gene axis when its gene set differs).

The detected `<format>@<version>` is recorded in `ds.profiles`, so a downstream reader knows exactly
what produced the data.

## Collections convert too

A multi-sample study is a **collection**, not one matrix (see [the model](model.md)). The collection
profiles convert it as such:

```r
library(lstar)
ds <- write_conos(conos_object)   # a Conos object -> an L* collection:
#   a `samples` axis; per-sample cells.{s}/genes.{s} axes + counts.{s} measures;
#   a union `cells` axis with a `sample` label; the joint embedding, clusters, and integration graph.
lstar_write(ds, "study.lstar.zarr")   # reads back identically in Python / C++ / the browser
co  <- read_conos(ds)             # and back: rebuilds a live Conos (per-sample Pagoda2 + joint layer)
```

A split Seurat v5 object converts the same way through `read_seurat` (example
[7 in examples.md](examples.md#7-r-seurat-v5-split-collection)). The heterogeneity these preserve —
samples with different cells, even different gene sets — is exactly what concatenating into one matrix
would erase.

**A collection converts to the native collection formats, too** — and *without* a corrected expression
matrix, which a graph-based integration like Conos never computes (it integrates in graph space only):

```r
so <- write_seurat(ds)   # -> a Seurat v5 SPLIT assay: per-sample raw layers over the UNION genes,
#   the joint graph as Graphs(), the embedding as a DimReduc, clusters/sample in meta.data.
#   read_seurat(so) reads it straight back as an L* collection. Per-sample PCA -> so@misc$lstar_dropped.
```
```python
from lstar import read, write_anndata
a = write_anndata(read("study.lstar.zarr"))   # -> ONE AnnData (a single matrix, so this flattens):
#   X = raw joint counts; obs = sample + clusters; obsm["X_*"] = embedding; obsp = the graph (aliased to
#   `connectivities` + uns["neighbors"] so sc.tl.leiden / sc.tl.umap run with no extra prep).
```

Neither target needs a corrected matrix — Seurat v5 (un-integrated split) and scanpy (Harmony/BBKNN/scVI
leave `X` raw and put the integration in `obsm`/`obsp`) both store graph-based integration natively.
Fidelity asymmetry: **Seurat v5 preserves the collection** (split layers round-trip to a collection);
**AnnData flattens it** to one union matrix. Keep the L★ store to retain the full per-sample structure.
Covered by `conformance/conos.sh` (a real Conos object + a synthetic divergent-genes collection).

---

See the [worked examples](examples.md) for runnable code, and the normative profile rule catalog
(every native location and its L★ signature, for AnnData, Seurat, pagoda2, Conos, and cacoa) in
[`../misc/Lstar_proposal.md`](../misc/Lstar_proposal.md), Appendix B.
