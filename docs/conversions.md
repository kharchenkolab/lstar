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

## The functions

You only need a handful. The reader/writer for each format:

| Format | Read into L★ | Write from L★ | Language | Package |
|---|---|---|---|---|
| **AnnData** (`.h5ad`/`.zarr`) | `read_anndata(adata)` | `write_anndata(ds)` | Python | `lstar.profiles.anndata` |
| **Seurat** (v3/v4/v5) | `read_seurat(so)` | `write_seurat(ds)` | R | `lstar` |
| **SingleCellExperiment** | `read_sce(sce)` | `write_sce(ds)` | R | `lstar` |
| **Conos** (a collection) | `write_conos(co)` → L★ | *(read-back deferred)* | R | `lstar` |
| **pagoda2** | *(via Conos members; standalone reader planned)* | | R | `lstar` |
| **L★ store** | `lstar.read(path)` / `lstar_read(path)` | `lstar.write(ds, path)` / `lstar_write(ds, path)` | Py / R | `lstar` |

The L★ dataset they hand back is the universal intermediate: in Python a `lstar.Dataset`, in R an
`lstar_dataset` (a list of `$axes` and `$fields`). Either can be written to a `.lstar.zarr` store and
read by any of the three languages.

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

Boundary: bounded-memory conversion works whenever **both ends are on-disk formats** (`.h5ad` ↔
`.lstar.zarr` today). Converting all the way into an *in-memory* object — a Seurat `dgCMatrix` in an
`.rds`, an in-memory AnnData `X` — still needs the matrix in RAM at the destination, *unless* that
format is targeted as a *disk-backed* representation (AnnData `backed`, Seurat v5 **BPCells**, SCE
`HDF5Array`), which lets the streamed `.h5ad` above be opened without a full load. Those disk-backed
target adapters are the next step.

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

The PCA-loadings row is the concrete payoff: a direct AnnData→Seurat conversion usually discards the
gene loadings, but because L★ keeps the scores and loadings on one shared `pca` axis, they ride through.

**What a given target cannot hold is reported in `dropped`:**

- Through **Seurat or SCE**: neighbor/other **graphs** (`relation` fields), arity-3 tensors, trees, and
  fitted models have no native slot, so they are written to `@misc` / `metadata` and listed in
  `dropped`. (An AnnData ↔ L★ ↔ AnnData round-trip *does* keep `obsp` graphs, because AnnData has a
  place for them.)
- Through **AnnData**: anything in `uns` that isn't a recognized structure is recorded in `dropped`
  (re-surfaced under `adata.uns['lstar/dropped']` on write-back).

The rule of thumb: **convert through L★ and keep the L★ store** if you want nothing lost; convert all
the way to a native format and you keep its shared-vocabulary core plus a recorded manifest of what it
couldn't carry. If you keep the data in L★, even the off-vocabulary pieces are preserved.

## Versions are recognized, not assumed

The readers detect the version/variant of the object they're given and adapt, rather than assuming one
layout — so conversions don't break when a collaborator is on a different release:

- **Seurat**: v3/v4 (`Assay`, fixed `counts`/`data`/`scale.data` slots) vs. v5 (`Assay5`, layers); a
  fallback covers SeuratObject < 5. A **split v5 assay** (`split(assay, f = sample)`) is recognized as
  a **collection**.
- **pagoda2**: the modern `getRawCounts()` accessor vs. the legacy `$counts` slot.
- **AnnData**: the library version, and the `.raw` slot (kept on its own gene axis when its gene set
  differs from the main one).

The detected `<format>@<version>` is recorded in `ds.profiles`, so a downstream reader knows exactly
what produced the data.

## Collections convert too

A multi-sample study is a **collection**, not one matrix (see [the model](model.md)). The collection
profiles convert it as such:

```r
library(lstar)
ds <- write_conos(conos_object)   # a Conos object -> an L* collection:
#   a `samples` axis; per-sample cells.<s>/genes.<s> axes + counts.<s> measures;
#   a union `cells` axis with a `sample` label; the joint embedding, clusters, and integration graph.
lstar_write(ds, "study.lstar.zarr")   # reads back identically in Python / C++ / the browser
```

A split Seurat v5 object converts the same way through `read_seurat` (example
[7 in examples.md](examples.md#7-r-seurat-v5-split-collection)). The heterogeneity these preserve —
samples with different cells, even different gene sets — is exactly what concatenating into one matrix
would erase.

---

See the [worked examples](examples.md) for runnable code, and the normative profile rule catalog
(every native location and its L★ signature, for AnnData, Seurat, pagoda2, Conos, and cacoa) in
[`../misc/Lstar_proposal.md`](../misc/Lstar_proposal.md), Appendix B.
