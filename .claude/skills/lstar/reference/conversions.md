# Format conversion (reference)

lstar's near-term value: glue between single-cell formats. `convert(X → Y) = write_Y(read_X(obj))`,
with an L★ dataset (in memory) or a `.lstar.zarr` store (on disk, bridges languages) as the universal
intermediate.

## Readers / writers

| Format | read → L★ | L★ → write | Lang | Import |
|---|---|---|---|---|
| AnnData (`.h5ad`/`.zarr`) | `read_anndata(adata)` | `write_anndata(ds)` | Py | `from lstar.profiles.anndata import read_anndata, write_anndata` |
| Seurat (legacy v2 → v5) | `read_seurat(so[, assay])` | `write_seurat(ds)` | R | `library(lstar)` |
| SingleCellExperiment | `read_sce(sce)` | `write_sce(ds)` | R | `library(lstar)` |
| Conos (collection) | `write_conos(co)` → L★ | *(read-back deferred)* | R | `library(lstar)` |
| L★ store | `lstar.read(p)` / `lstar_read(p)` | `lstar.write(ds,p)` / `lstar_write(ds,p)` | Py/R | — |

Intermediate value: Python `lstar.Dataset`; R `lstar_dataset` (`$axes`, `$fields`). Either → a
`.lstar.zarr` store readable by Py/R/C++/JS.

## Conversion matrix

- **Same language (in memory):** R Seurat↔SCE (`write_sce(read_seurat(so))`, `write_seurat(read_sce(sce))`);
  Py AnnData↔AnnData (round-trip / fixed point).
- **Cross language (bridge via store):** AnnData(Py) ↔ Seurat/SCE(R) — write the store on one side,
  read it on the other. No reticulate; the store is the contract.
- **Chains:** AnnData → Seurat → SCE → … → AnnData; the shared-vocabulary core survives any length
  (see `examples/roundtrip_xlang.sh`).

## What survives (shared-vocabulary core) vs. what drops

Survives every AnnData↔Seurat↔SCE conversion: `counts`(raw), `X`/`data`(lognorm), `scaled`,
`pca` **+ `pca_loadings`** (one shared axis — loadings survive, which direct converters drop),
`umap`/`tsne`, cell/gene metadata, `label`s (clusters/celltype).

Recorded in `ds.dropped` (written to the target's sidecar, never silent):
- through **Seurat/SCE**: `relation` fields (neighbor graphs), arity-3 tensors, trees, models — no slot.
- through **AnnData**: unrecognized `uns` entries (re-surfaced under `adata.uns['lstar/dropped']`).
- AnnData↔L★↔AnnData keeps `obsp` graphs + `.raw` (AnnData has slots for them).

Rule: keep the L★ store to lose nothing; convert to a native format to keep its core + a `dropped` manifest.

## Mapping specifics (for accurate docs/answers)

- **anndata** (`python/src/lstar/profiles/anndata.py`): `X`→`X` measure; `.raw.X`→`raw` (own
  `genes_raw` axis if divergent); `layers[k]`→measures; `obs/var[k]`→label|measure; `obsm[k]`→embedding
  (`X_pca`→`pca`); `varm[k]`→`<coord>_loadings`; `obsp/varp`→relations; `uns`→dropped. Provenance records
  the exact native slot for exact write-back.
- **seurat** (`R/R/profile_seurat.R`): assay layers `counts`/`data`/`scale.data` ↔ measures
  `counts`/`X`/`scale.data` (states raw/lognorm/scaled), **transposed** (Seurat is genes×cells);
  `meta.data`→cell fields; `DimReduc`→embedding (+ `<rn>_loadings`); split v5 assay→collection.
  No graph slot (relations drop).
- **sce** (`R/R/profile_sce.R`): assays↔measures (`counts`→raw, `logcounts`→`X`/lognorm), transposed;
  `colData`/`rowData`↔cell/gene fields; `reducedDims`↔embeddings (+ rotation loadings). No graph slot.
- **conos** (`R/R/profile_conos.R`): collection — `samples` axis, per-sample `cells.<s>`/`genes.<s>` +
  `counts.<s>` + `pca.<s>`, union `cells`, `sample` design label, joint embedding/clusters/graph.

## Version recognition (don't assume one layout)
Seurat legacy **v2** (pre-`Assay` lowercase `seurat` class, slots via `attr()`) / v3/v4 `Assay` vs v5
`Assay5` (+ `GetAssayData` fallback for SeuratObject<5; split v5 →
collection); pagoda2 `getRawCounts()` vs `$counts`; AnnData lib version + `.raw`. Detected
`<format>@<version>` recorded in `ds.profiles`; unrepresentable → `dropped`.

## Examples
`examples/convert_h5ad_to_seurat.sh` (commented h5ad→Seurat .rds, reports preserved/dropped),
`examples/cross_language_demo.sh`, `examples/roundtrip_xlang.sh` (multi-format chain), `examples/conos_collection_demo.R`.
Conformance: `conformance/cross_format.sh` (AnnData↔Seurat↔SCE core preserved). Full guide:
`docs/conversions.md`; normative profile catalog: `misc/Lstar_proposal.md` Appendix B.
