# Format conversion (reference)

lstar's near-term value: glue between single-cell formats. `convert(X → Y) = write_Y(read_X(obj))`,
with an L★ dataset (in memory) or a `.lstar.zarr` store (on disk, bridges languages) as the universal
intermediate.

## The CLI (`lstar convert`) — the quickest path

```bash
lstar convert a.h5ad a.rds                 # AnnData -> Seurat; detect by extension, bridge Py↔R via store
lstar convert a.h5ad a.lstar.zarr          # -> store  (--to sce for SCE; --from/--to to override detection)
lstar convert a.rds  a.h5ad --report       # + fidelity report: every field + provenance + `dropped`
lstar inspect a.h5ad --report-json r.json  # read + structured report, no write
```
- Detection: `.h5ad`→anndata, `.h5mu`→mudata, `.rds`→Seurat/SCE (sniffed), `.lstar.zarr`/`.zarr`→store.
- `--check` (default on; `--strict` gates exit code): opens the result in its native library + runs a
  canonical-ops smoke (scanpy `pca`/`rank_genes_groups`; Seurat `RunPCA`; scran `modelGeneVar`) — proves
  native tools accept it, not just round-trip. Heavy libs optional → degrades to open + structural.
- Seurat/SCE legs need R + the `lstar` package (`LSTAR_RLIB`/`LSTAR_RSCRIPT` if not on the default path).
- **`--backend auto|native|direct`** — the package-free fallback. `auto` (default) uses the format's native
  package when present, else lstar's own codec, so the domain packages aren't *required*. Without them:
  `.h5ad`↔store (read+write) needs only `h5py`; **Seurat `.rds`**↔store (read+write) and **SCE `.rds`**→store
  (**read**) need only base R + the `lstar` R package (no SeuratObject/SingleCellExperiment — readers walk
  S4 slots via `attr()`, the Seurat writer builds a pinned-schema object). Native-only: **SCE write** (a
  valid SCE needs the SummarizedExperiment/GRanges machinery) and `.h5mu`. `--backend direct` forces the
  codec; at a wall (unknown version, `BPCells`-backed matrix) it raises a clear "install X" error. Analysis
  packages (scanpy/Seurat/scran) are only for `--check`, never for converting.
- Entry point: `python -m lstar convert …` or the `lstar` console script. Code: `python/src/lstar/cli.py`
  (+ `_native_check.py`, `profiles/anndata_direct.py`); the deterministic role→slot contract is `docs/mapping.md`.

## Readers / writers

| Format | read → L★ | L★ → write | Lang | Import |
|---|---|---|---|---|
| AnnData (`.h5ad`/`.zarr`) | `read_anndata(adata)` | `write_anndata(ds)` | Py | `from lstar.profiles.anndata import read_anndata, write_anndata` |
| Seurat (legacy v2 → v5) | `read_seurat(so[, assay])` | `write_seurat(ds)` | R | `library(lstar)` |
| SingleCellExperiment | `read_sce(sce)` | `write_sce(ds)` | R | `library(lstar)` |
| Conos (collection) | `write_conos(co)` → L★ | `read_conos(ds)` → Conos | R | `library(lstar)` |
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
- through **Seurat**: arity-3 tensors, trees, models — no slot. (Cell-cell `relation` graphs **are** kept,
  as Seurat `Graphs()`.)
- through **SCE**: `relation` fields (neighbor graphs), arity-3 tensors, trees, models — no slot.
- through **AnnData**: unrecognized `uns` entries (re-surfaced under `adata.uns['lstar/dropped']`).
- AnnData↔L★↔AnnData keeps `obsp` graphs + `.raw` (AnnData has slots for them).

Rule: keep the L★ store to lose nothing; convert to a native format to keep its core + a `dropped` manifest.

## Collections (Conos and other multi-sample sets)

A **collection** (`kind="collection"`, e.g. from `write_conos`) is per-sample raw counts over divergent
`genes.<s>` axes + a **joint layer over the union cells** — a graph, an embedding, clustering(s). A Conos
object integrates in **graph space only**: there is *no corrected/batch-aligned expression matrix* (it's
expensive to compute). That's fine — neither native collection format needs one:

- **`write_seurat(collection)` → Seurat v5 split assay** (R). Per-sample counts become split layers
  (`counts.<s>`) over the **union** gene set (a gene absent in a sample is 0 in its layer); the joint graph
  → `Graphs()`, the joint embedding → a `DimReduc`, clustering(s)/`sample` → `meta.data`. No integrated
  assay is fabricated — an un-integrated-but-jointly-analyzed v5 object is exactly this. **`read_seurat`
  reads it back as a `collection`** (per-sample structure preserved). Per-sample `pca.<s>` (sample-local
  rotations, not a joint reduction) → `dropped` (recorded in `so@misc$lstar_dropped`).
- **`write_anndata(collection)` → one AnnData** (Py). A single AnnData *is* one matrix, so this **flattens**:
  `X` = raw joint counts (== conos `getJointCountMatrix()`), `obs` = `sample` + clustering, `obsm["X_*"]` =
  embedding, `obsp` = the graph (also aliased to `connectivities` + `uns['neighbors']` so `sc.tl.leiden` /
  `sc.tl.umap` run with no extra prep — exactly how scanpy stores a graph-integrated dataset). Reading it
  back gives a `sample`, not a collection. Disjoint (cross-species) gene sets union into one sparse matrix —
  Seurat v5's split keeps them honestly separate, so prefer Seurat there.

Fidelity asymmetry: **Seurat v5 preserves the collection** (split layers round-trip); **AnnData flattens it**
(one matrix). Keep the L★ store to retain the full per-sample structure. Tested in `conformance/conos.sh`
(real Conos `small_panel` + a synthetic divergent-genes collection).

## Mapping specifics (for accurate docs/answers)

Full deterministic role→slot contract (and the per-format conventions native tools require — Seurat
`_`-terminated `DimReduc` keys, scanpy categorical `groupby`, SCE `logcounts` name): **`docs/mapping.md`**.
Verify with `lstar convert --check` (native-acceptance: open + canonical-ops smoke). Per-profile summary:

- **anndata** (`python/src/lstar/profiles/anndata.py`): `X`→`X` measure; `.raw.X`→`raw` (own
  `genes_raw` axis if divergent); `layers[k]`→measures; `obs/var[k]`→label|measure; `obsm[k]`→embedding
  (`X_pca`→`pca`); `varm[k]`→`<coord>_loadings`; `obsp/varp`→relations; `uns`→dropped. Provenance records
  the exact native slot for exact write-back.
- **seurat** (`R/R/profile_seurat.R`): assay layers `counts`/`data`/`scale.data` ↔ measures
  `counts`/`X`/`scale.data` (states raw/lognorm/scaled), **transposed** (Seurat is genes×cells);
  `meta.data`→cell fields; `DimReduc`→embedding (+ `<rn>_loadings`); cell-cell `relation`↔`Graphs()`;
  split v5 assay→collection, **and** a `collection`→v5 split assay over the union genes (`write_seurat`).
- **sce** (`R/R/profile_sce.R`): assays↔measures (`counts`→raw, `logcounts`→`X`/lognorm), transposed;
  `colData`/`rowData`↔cell/gene fields; `reducedDims`↔embeddings (+ rotation loadings). No graph slot.
- **conos** (`R/R/profile_conos.R`): collection — `samples` axis, per-sample `cells.<s>`/`genes.<s>` +
  `counts.<s>` + `pca.<s>`, union `cells`, `sample` design label, joint embedding/clusters/graph.
  `write_conos(co)` imports; `read_conos(ds)` rebuilds a live Conos (per-sample Pagoda2 + joint layer).

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
