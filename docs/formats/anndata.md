# lstar ⇄ AnnData

Everything about moving data between **AnnData** (the scverse/scanpy in-memory object and its
`.h5ad` / `.zarr` on-disk forms) and lstar's L★ model: how to read, write, convert, and stream it,
the exact representation mapping, what survives a conversion, and the version quirks the reader
handles. AnnData is lstar's most-developed profile and the only one with a full on-disk **streaming**
path, because `.h5ad` is itself a disk format.

**On this page:** [Quick start](#quick-start) · [Install](#requirements--install) ·
[Reading](#reading-anndata-into-l) · [Writing](#writing-anndata-from-l) ·
[Converting](#converting) · [Representation](#the-l-representation-the-anndata-profile) ·
[Preserved vs dropped](#what-is-preserved-vs-dropped) · [Versions](#versions--variants-recognized) ·
[Native validity](#native-validity--the---check-smoke) · [Collections](#collections) ·
[Examples](#worked-examples) · [Troubleshooting](#troubleshooting) · [API](#api--references)

Entry points: **Python** `read_anndata` / `write_anndata` (in-memory), `convert_anndata` /
`convert_to_h5ad` (streamed, on-disk), and the `lstar convert` CLI. AnnData lives in Python; to land
the result in an R object (Seurat/SCE) the conversion routes through a `.lstar.zarr` store.

## Quick start

```bash
lstar convert sample.h5ad sample.lstar.zarr --report   # AnnData -> L* store + a fidelity report
lstar convert sample.h5ad sample.rds                   # AnnData -> Seurat (.rds), bridged through R
lstar inspect sample.h5ad                              # print the L* structure of an .h5ad, no write
```

```python
from lstar.profiles.anndata import read_anndata, write_anndata
ds     = read_anndata(adata)     # AnnData object -> L* dataset
adata2 = write_anndata(ds)       # L* dataset -> AnnData object (each field returned to its slot)
print(ds.dropped)                # what L* had no field for (e.g. unrecognized uns entries) — recorded, not lost
```

## Requirements & install

| route | with the native package | package-free (`--backend direct`) |
|---|---|---|
| `.h5ad` ↔ L★ store | `anndata` | **`h5py`** — lstar reads/writes the HDF5 encoding directly |
| in-memory `read_anndata`/`write_anndata` | `anndata` | — (operates on a live `AnnData`) |

```bash
pip install lstar-sc                 # the Python package; add h5py for the package-free .h5ad codec
```

The analysis stack (`scanpy`) is **not** needed for conversion — only for the optional `--check`
native-acceptance smoke (below). So `pip install lstar-sc h5py` converts `.h5ad` with neither anndata
nor scanpy installed.

## Reading AnnData into L★

`read_anndata(adata, kind="sample")` turns a live `AnnData` into an `lstar.Dataset`. It does **not**
copy the matrices eagerly when the object is `backed`: `X`, each `layers[...]`, and `.raw.X` are held
as on-disk streaming sources (`LazyCSX`) so a multi-gigabyte atlas opens in megabytes.

```python
import anndata as ad
from lstar.profiles.anndata import read_anndata

ds = read_anndata(ad.read_h5ad("atlas.h5ad", backed="r"))   # matrices stay on disk (lazy)
ds.field("X").values        # a LazyCSX proxy until you touch it
ds.profiles                 # e.g. ['anndata@0.10.8'] — the detected source version
```

The reader **detects the object's version/layout and adapts** (see [Versions](#versions--variants-recognized)):
the anndata version, the on-disk sparse encoding, and whether a `.raw` slot is present go into
`ds.profiles`/`ds.provenance`, so a downstream reader knows exactly what produced the data.

## Writing AnnData from L★

`write_anndata(ds)` materializes an `AnnData`, returning each field to the slot it came from (recorded
in `provenance`), and honoring the conventions scanpy/anndata expect (unique string `obs_names`,
`X_`-prefixed `obsm` keys, categorical `obs` columns for `groupby`, float `X` with integer counts in
`layers["counts"]` — see [Native validity](#native-validity--the---check-smoke)). Fields L★ carries
that AnnData has no slot for are written under `uns["lstar/dropped"]` and listed in `ds.dropped`.

For large datasets, `write_anndata_streamed(ds, "out.h5ad")` (and the convenience
`convert_to_h5ad(store, "out.h5ad")`) stream each measure block-by-block into the `.h5ad`'s sparse
groups rather than building the matrix in memory.

## Converting

**The one command.** `lstar convert` detects both endpoints from their paths and routes through the
L★ store:

```bash
lstar convert pbmc.h5ad pbmc.rds              # -> Seurat (.rds); --to sce for SingleCellExperiment
lstar convert pbmc.h5ad pbmc.lstar.zarr       # -> L* store (keeps everything; nothing dropped)
lstar convert s.rds s.h5ad --report           # Seurat/SCE -> AnnData + the full fidelity report
```
`--report`/`--report-json` lists every field with role/state/span/provenance and the `dropped` set;
`--check` (default on; `--strict` to gate the exit code) opens the result in scanpy and runs a
canonical-ops smoke; `--backend auto|native|direct` picks the codec.

**Same language (Python, in memory).** `write_anndata(read_anndata(adata))` is a faithful round-trip —
useful to confirm nothing is lost.

**Across languages (the common case: `.h5ad` ↔ Seurat/SCE).** AnnData is Python, Seurat/SCE are R, so
the `.lstar.zarr` store on disk is the bridge — no `reticulate`, no shared memory:

```python
import anndata as ad, lstar
from lstar.profiles.anndata import read_anndata
lstar.write(read_anndata(ad.read_h5ad("pbmc.h5ad")), "pbmc.lstar.zarr")   # Python writes the store
```
```r
so <- lstar::write_seurat(lstar::lstar_read("pbmc.lstar.zarr")); saveRDS(so, "pbmc.rds")  # R reads it
```

**Streaming (bounded memory).** Because both `.h5ad` and `.lstar.zarr` are on-disk, conversion between
them runs in ~one block of memory, in **both** directions; each matrix keeps its native orientation
(a CSR `X` stays CSR, a CSC layer stays CSC — no transpose), and the result is byte-identical to the
eager path:

```python
import lstar
lstar.convert_anndata("atlas.h5ad", "atlas.lstar.zarr")   # h5ad -> store, peak ~ one block
lstar.convert_to_h5ad("atlas.lstar.zarr", "atlas.h5ad")   # store -> h5ad, peak ~ one block
```
On a 40,220 × 20,138 atlas (77.6M nonzeros) the `h5ad → L★` direction peaks at ~110 MB vs ~1.3 GB
eager (≈12×). The streamed `.h5ad` then opens disk-backed all the way into a native object:
`anndata.read_h5ad(path, backed="r")` (Python), or in R `read_seurat_backed(path)` (BPCells) /
`read_sce_backed(path)` (HDF5Array).

## The L★ representation (the AnnData profile)

The mapping is **bidirectional and field-preserving**: each slot has a fixed L★ signature, so meaning
survives even where a one-off `as.Seurat()` would drop it (PCA gene loadings are the canonical
example). This is the inverse view of the combined table in [`../mapping.md`](../mapping.md).

**Identity → axes.** The two indices become the dataset's primary axes:

| AnnData index | L★ axis | origin / role |
|---|---|---|
| `.obs` index (`obs_names`) | `cells` | observed / observation |
| `.var` index (`var_names`) | `genes` | observed / feature |

**Field rules.** Every slot maps to a typed field over those axes:

| AnnData slot | L★ field | role · span · state/subtype |
|---|---|---|
| `X` | `X` (the "primary" measure) | measure · (cells, genes) · `state` from `uns['lstar/state']` else inferred |
| `layers[k]` | `k` | measure · (cells, genes) · state guessed from the name (`counts`→raw, `logcounts`/`data`→lognorm, `scaled`→scaled) |
| `.raw.X` | `raw` (on its **own** gene axis if the gene set differs) | measure · (cells, raw-genes) · raw |
| `obs[k]` | `k` | by dtype: categorical→`label`, numeric→`measure` · (cells) |
| `var[k]` | `k` | by dtype · (genes) |
| `obsm["X_<red>"]` | `<red>` (e.g. `pca`, `umap`) | embedding · (cells, `<red>` coord axis) |
| `obsm["spatial"]` | `spatial` | embedding · (cells, **observed** coord axis) · subtype=spatial |
| `varm["<RED>s"]` | `<red>_loadings` (paired to the matching `obsm`) | loading · (genes, `<red>`) |
| `obsp[k]` | `k` | relation · (cells, cells) · subtype guessed (distances→similarity) |
| `varp[k]` | `k` | relation · (genes, genes) |
| `uns[k]` (unrecognized) | — | recorded in `dropped` (read-only catch-all; see below) |

The **state** typing is what makes routing safe across formats: a measure tagged `state=lognorm`
lands in Seurat's `data` slot / SCE's `logcounts`, never confused with raw counts. PCA scores
(`obsm['X_pca']`) and their gene loadings (`varm['PCs']`) are tied onto **one shared `pca`
coordinate axis**, so they convert together. A modality measured on a subset of cells rides over the
shared `cells` axis as **partial coverage** (an `index`), not a padded matrix.

## What is preserved vs dropped

The **shared-vocabulary core** — raw/normalized/scaled expression, PCA (scores + loadings) / UMAP /
t-SNE, `obs`/`var` metadata, clusterings, cell–cell/gene–gene graphs (`obsp`/`varp`), a second
modality, and spatial coordinates — survives every conversion among AnnData, Seurat, and SCE.

What has no cross-format slot is written to the target's sidecar and **listed in `ds.dropped`**, never
silently lost. For AnnData specifically, **`uns`** is the main source of drops: any `uns` entry that
isn't a recognized structure (a neighbors graph, an `X_*` embedding companion) is recorded and
re-surfaced under `adata.uns["lstar/dropped"]` on write-back. Arity-3 tensors, trees, and fitted
models have no L★ field and go the same way. **Keep the `.lstar.zarr` store and nothing is dropped at
all** — the `dropped` manifest only describes what a *native target* couldn't carry.

## Versions & variants recognized

The reader does not assume one layout — it detects and adapts, so conversions don't break across
collaborator releases:

- **anndata library version** — recorded as `anndata@<version>` in `ds.profiles`.
- **on-disk sparse encoding** — legacy `h5sparse_format` (anndata < 0.7) and the modern
  `encoding-type` attribute are both read.
- **`.raw`** — kept on its own gene axis when its gene set differs from `.var` (older pipelines stash
  pre-HVG counts there). The modern (`raw/X`) vs legacy (`raw.X`) on-disk location is handled.
- **`uns['neighbors']`** — recent anndata migrates `distances`/`connectivities` into `.obsp`; the
  reader follows either layout.

## Native validity & the `--check` smoke

A `X → store → X` round-trip proves L★ kept *its own* representation; it does not prove scanpy will
accept the object. `lstar convert --check` (default on) closes that gap: it opens the produced
AnnData and runs a canonical-ops smoke — `sc.pp.normalize_total` / `log1p` / `pca` /
`rank_genes_groups` — plus the AnnData invariants the profile enforces:

- `obs_names` / `var_names` are **unique** strings (scanpy indexes on them).
- a `rank_genes_groups` `groupby` column is a **categorical** dtype.
- `obsm` embeddings are `X_`-prefixed (`X_pca`) so scanpy's plotting finds them.
- `.X` is float; integer counts live in `layers["counts"]`.

`scanpy` absent → the check degrades to *open + structural invariants* and reports `ops skipped`;
`--strict` turns a failure into a non-zero exit.

## Collections

A multi-sample study is a **collection**, not one matrix. AnnData represents one matrix, so a
collection (e.g. a Conos integration) **flattens** to a single AnnData over the union of genes — a
documented, lossy projection:

```python
from lstar import read, write_anndata
a = write_anndata(read("study.lstar.zarr"))   # X = raw joint counts; obs = sample + clusters;
#   obsm["X_*"] = joint embedding; obsp = the integration graph (aliased to `connectivities` +
#   uns["neighbors"] so sc.tl.leiden / sc.tl.umap run with no extra prep).
```
This is faithful to how scanpy stores graph-based integration (Harmony/BBKNN/scVI leave `X` raw and
put the result in `obsm`/`obsp`) — but the per-sample structure (divergent gene sets) is unioned away.
**Seurat v5 preserves a collection** (split layers); **AnnData flattens it.** Keep the L★ store to
retain the full per-sample structure.

## Worked examples

Runnable scripts in [`../../examples/`](../../examples): `convert_h5ad_to_seurat.sh` (a complete
`.h5ad → Seurat` converter that reports what survived) and `roundtrip_xlang.sh` (AnnData → Seurat →
SCE → … → AnnData, verifying the core survives the whole chain). The Python round-trip and lazy
streaming are walked in [`../examples.md`](../examples.md) §2–3.

## Troubleshooting

- **`X` is `None` after a conversion.** The source had only raw counts and no log-normalized layer;
  the raw measure correctly lands in `layers["counts"]`, leaving `X` empty (faithful — normalize to
  populate `X`).
- **A `.var`/`.obs` column came back as object/string, not categorical.** `read_anndata` types numeric
  columns as `measure` and only factor/categorical ones as `label`; recast in scanpy if you need a
  categorical `groupby`.
- **Something is missing.** Check `ds.dropped` / `adata.uns["lstar/dropped"]` — it enumerates exactly
  what the target couldn't hold. Convert to a `.lstar.zarr` store instead to keep everything.
- **`--backend direct` failed on a feature.** The package-free h5py codec stops with a message naming
  what only `anndata` can handle; rerun with `--backend native` (or `auto`).

## API & references

- **Python**: `read_anndata(adata)`, `write_anndata(ds)`, `convert_anndata(h5ad, store)`,
  `convert_to_h5ad(store, h5ad)`, `write_anndata_streamed(ds, h5ad)`; `lstar.read`/`lstar.write`,
  `lstar.read(path, lazy=True)`, `stream_col_stats`.
- **R (disk-backed AnnData targets)**: `read_seurat_backed(h5ad)`, `read_sce_backed(h5ad)`.
- The deterministic role→slot contract: [`../mapping.md`](../mapping.md). The L★ model:
  [`../model.md`](../model.md). The conversion overview: [`../conversions.md`](../conversions.md).
  Normative profile rules: [`../../misc/Lstar_proposal.md`](../../misc/Lstar_proposal.md) Appendix B.2.
