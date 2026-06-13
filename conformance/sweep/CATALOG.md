# Corpus catalog — what each real dataset exercises, and its CI synthetic

A working map of the **local real corpus** (curated loaders + the breadth sweeps) and how each real
dataset maps to a **synthetic CI** stand-in. The rule: real data is local-only; every structural feature
a real dataset carries must be represented by a synthetic CI fixture (distilled, not the real bytes).
This file is the checklist for "is this real structure represented synthetically?" — the
`test_synth_faithful` guard enforces it for the curated loaders.

Legend for **CI synthetic**: the generator/fixture that stands in for the real dataset in CI; `—` = not
represented (gap) · `(structure)` = the structure is synthesized though no single 1:1 generator exists.

## Curated loaders (`python/tests/corpus.py` ↔ `python/tests/synth.py`)

| real dataset | source/format | exercises | CI synthetic |
|---|---|---|---|
| `pbmc68k_reduced` | scanpy `.h5ad` | obs categoricals (bulk_labels/phase/louvain), `*_colors`, `uns['pca']`, neighbors `OverloadedDict`, `.raw` (same gene set), **logreg** one-vs-rest DE | `synth.pbmc68k_like()` (real scanpy pipeline on synthetic counts) |
| `pbmc3k_processed` | scanpy `.h5ad` | louvain/colors/pca, `.raw` over a **divergent** gene set (HVG subset), DE substrate | `synth.pbmc3k_like()` (subset_hvg=True) |
| `pbmc3k_with_de` | scanpy `rank_genes_groups` | **t-test/wilcoxon** full DE bundle; **pairwise** (reference-group) kept verbatim | DE run on the synthetic pbmc3k |
| `pancreas_velocity` | scVelo `.h5ad` | `spliced`/`unspliced` layers, `clusters`+colors, `uns['velocity_graph']` (+`_neg`) cell×cell | `synth.velocity()` |
| `citeseq_matrices` / `citeseq_mudata` | real minipbcite subsample (`.mtx`) | RNA+ADT → genes/proteins feature axes (shared by Py MuData + R Seurat/SCE) | `synth.citeseq_*` + `synth.write_citeseq_mtx` (one generator, both languages) |
| `minipbcite` | downloaded `.h5mu` (411 cells) | global `celltype` factor, per-mod PCA/uns, global joint embedding | `synth.citeseq_mudata_annotated()` |
| `marrow_backed` | TMS Senis Marrow `.h5ad` (1.2 GB, local) | realistic-size (40 k×20 k) **backed** read; many real factor axes; color/pca promotion | — (size; CI exercises the backed proxy on small synthetic) |

## SCE breadth sweep — Bioconductor `scRNAseq` (`sweep_scrnaseq.R`)

**56/61 PASS, 0 profile bugs** (after 3 sweep-caught fixes). Representative structures and their CI stand-ins:

| real structure (example datasets) | exercises | CI synthetic |
|---|---|---|
| plain assays + colData/rowData (most) | counts/logcounts, factor col/rowData | `sce_versions.sh` cases |
| NULL dimnames (`BachMammaryData`, `ErnstSpermatogenesisData`) | cells keyed by a `Barcode` colData column, no colnames | (label-synthesis path; sce_versions could add a NULL-dimnames case — **gap**) |
| S4 `Rle` / nested `DataFrame`/`GRanges` cols (`BuettnerESCData`, `BunisHSPCData`, `DarmanisBrainData`) | Rle unpacking; uncoercible nested cols recorded | (**gap** — not in a synthetic fixture; sweep-only) |
| plain `SummarizedExperiment` (`ReprocessedFluidigmData`) | SCE-only accessors guarded | (**gap** — sweep-only) |
| `reducedDims` + `altExps` (various) | embeddings (+rotation→loadings); ADT/ERCC altExps → feature axes | `sce_versions.sh` `+reducedDims`, `+altExps` |
| `colPairs`/`rowPairs` | cell-cell / gene-gene graphs → relations | `sce_versions.sh` `+colPairs/rowPairs` |

> Note the three sweep-caught SCE structures (NULL dimnames, S4 Rle/nested cols, SE-not-SCE) are
> **handled** by the profile but **not yet represented by a synthetic CI fixture** — candidates to add so
> CI guards against regressions, not just the local sweep.

## Seurat breadth sweep — SeuratData (`sweep_seurat.R`) — **10/10 PASS**

| real object | classes | exercises | CI synthetic |
|---|---|---|---|
| `pbmc3k` / `pbmc3k.final` | Assay (v4) | RNA; **HVG-subset PCA loadings** (2000/13714) | `seurat_versions.sh` v3 Assay + `HVG-subset loadings` |
| `cbmc`, `bmcite` | Assay+Assay | RNA+ADT CITE-seq (bmcite: spca loadings 2000/17009) | `seurat_versions.sh` `multimodal RNA+ADT (synth)` |
| `thp1.eccite` | Assay ×4 | **4-modality** ECCITE-seq (RNA+ADT+HTO+GDO) | (multimodal path covers N assays; a 4-assay synthetic case is a **gap**) |
| `ifnb`, `panc8`, `pbmcsca`, `hcabm40k` | Assay | integration / multi-tech / multi-method | (collection path; `seurat_versions.sh` v5-split covers collection) |
| `celegans.embryo` | Assay | non-human RNA | (encoding-only; covered) |
| `pbmcMultiome` (`pbmc.atac`/`pbmc.rna`) | **ChromatinAssay** + Assay | RNA+ATAC; **108 k peak ranges** (seqnames/start/end) + fragments | (ChromatinAssay has **no synthetic fixture** — needs Signac; sweep-only — **gap**) |
| Azimuth `*ref` (11) | SCTAssay refs | large reference atlases | load-skip (need Azimuth loader); **not represented** |

## scVelo velocity / trajectory (`sweep_velocity.py`) — **4/4 PASS**

Local-only `.h5ad` cached under `testdata/velocity/` (download via an isolated scvelo venv).

| real dataset | source/format | exercises | CI synthetic |
|---|---|---|---|
| `dentategyrus` | scVelo `.h5ad` (mouse hippocampus) | `ambiguous`/`spliced`/`unspliced` layers, `clusters`+colors | `synth.velocity()` / `pancreas_velocity_small.h5ad` |
| `gastrulation_erythroid` | scVelo `.h5ad` (mouse gastrulation) | spliced/unspliced; many factor axes (stage/celltype/sequencing.batch) | `synth.velocity()` (layers) — **extra factor breadth is sweep-only** |
| `bonemarrow` | scVelo `.h5ad` (**human**) | spliced/unspliced; non-mouse organism | `synth.velocity()` (encoding-only; covered) |
| `pancreas` (full) | scVelo `.h5ad` (mouse pancreas) | spliced/unspliced; the full-size source of the CI subsample | `synth.velocity()` + the committed subsample fixture |

> The committed CI fixture `pancreas_velocity_small.h5ad` (with `uns['velocity_graph']`) stands in for the
> velocity structure; the sweep additionally confirms the layers survive on 4 real datasets across 2 organisms.

## CITE-seq RNA+ADT — 10x public `.h5` via Seurat (`sweep_citeseq_10x.R`) — **3/3 PASS**

| real dataset | source/format | exercises | CI synthetic |
|---|---|---|---|
| `5k_pbmc_protein_v3` | 10x `.h5` → RNA `Assay5` + ADT `Assay` | RNA 33538 + ADT **32** antibodies → genes/proteins feature axes | `seurat_versions.sh` `multimodal RNA+ADT (synth)` + `citeseq` mtx fixture |
| `pbmc_1k_protein_v3` | 10x `.h5` (smaller) | RNA + ADT 17; small-cell-count CITE-seq | same |
| `malt_10k_protein_v3` | 10x `.h5` (**MALT lymphoma**) | non-PBMC tissue CITE-seq | same (tissue is encoding-irrelevant) |

## Multiome RNA+ATAC — 10x public `.h5` (`sweep_multiome_10x.R` Seurat + `sweep_mudata.py` MuData) — **PASS**

| real dataset | source/format | exercises | CI synthetic |
|---|---|---|---|
| `pbmc_granulocyte_sorted_3k` | 10x ARC `.h5` → RNA `Assay5` + Signac `ChromatinAssay` | 36601 genes + **98319** peak ranges (chr:start-end → seqnames/start/end) | (ChromatinAssay — **no synthetic fixture**, needs Signac; sweep-only — **gap**) |
| `human_brain_3k` | 10x ARC `.h5` (**brain**) → ChromatinAssay | 134030 peaks; different peak set + tissue | same **gap** |
| `pbmc_multiome_3k.h5mu` | built `.h5mu` → MuData `peaks` axis | ATAC modality → a **`peaks` feature axis from MuData** (not Signac); `(cells, peaks)` measure | (MuData `peaks` axis — **no synthetic fixture**; `test_mudata.py` covers genes/proteins only — **gap**) |

## MuData `.h5mu` breadth (`sweep_mudata.py`) — **3/3 PASS (after the bug-5 fix)**

| real dataset | source/format | exercises | CI synthetic |
|---|---|---|---|
| `minipbcite` | downloaded `.h5mu` (411 cells) | global `celltype` factor, **per-modality PCA loadings of different dim** (rna 50 / prot 31) — *the bug-5 collision*, global joint embedding | `synth.citeseq_mudata_annotated()` + `test_mudata.py` MOFA case — **but the different-dimension per-modality PCA collision has no synthetic fixture (gap)** |
| `5k_pbmc_protein.h5mu` | built from 10x CITE-seq `.h5` | RNA+ADT modalities → genes/proteins; large real CITE-seq | `synth.citeseq_mudata()` |
| `pbmc_multiome_3k.h5mu` | built from 10x multiome `.h5` | RNA+ATAC → genes/**peaks**; peak coords as feature labels | (MuData `peaks` axis — **gap**, see above) |

## Integration read AS collections (`sweep_integration.R`) — **3/3 PASS**

Real SeuratData integration objects split by their sample column into a heterogeneous L* collection.

| real dataset | split by | exercises | CI synthetic |
|---|---|---|---|
| `ifnb` | `stim` (2) | per-sample axes + counts; union cells + `sample` label; round-trips | `collection.sh` synthetic 2-sample collection |
| `panc8` | `dataset` (8, 5 techs) | 8 heterogeneous per-sample gene sets | `collection.sh` (covers the shape; 8-sample real scale is sweep-only) |
| `pbmcsca` | `Method` (9) | 9 samples; largest (31021 cells) | `collection.sh` (shape) |

> These exercise the **same collection model the Conos profile builds** (per-sample `cells.<s>`/`genes.<s>`
> + `samples` axis + union `cells`), reached from Seurat integration objects rather than a Conos `.rds` —
> i.e. a second, format-independent producer of the collection shape, all round-tripping.

## Conos / pagoda2 (`sweep_conos.R`)

| real object | exercises | CI synthetic |
|---|---|---|
| `conI.rds` (4 samples), `con.rds` (8 samples) | collection (multi-sample, heterogeneous), joint graph, per-sample embeddings | `collection.sh` synthetic collection |
| real pagoda2 object | — | **gap** (no real object yet; mock pagoda2 covers the viewer schema) |

## Version variety (graceful recognition)

| real form | exercises | CI synthetic |
|---|---|---|
| `dentategyrus_anndata011.h5ad` (anndata 0.11; encoding-version 0.1.0 / dataframe 0.2.0) vs the legacy pre-0.7 `.h5ad` | both read by anndata 0.8 + `read_anndata` AND by streaming `convert_anndata` → identical L* | `test_versions.py` (synthetic old/new); the real cross-version read is sweep-only |

## Gaps this catalog surfaces (synthetic representation missing for a handled real structure)

These real structures **work** in the profile (sweep-verified) but have **no synthetic CI fixture**, so CI
can't guard them against regression — prioritized to add as small synthetic cases:

1. SCE **NULL dimnames** · **S4 Rle / nested cols** · **SE-not-SCE** — add to `sce_versions.sh`.
2. Seurat **ChromatinAssay** (synthetic peaks + ranges, gated on Signac in the local tier) and a
   **4-assay** (ECCITE-style) case.
3. **Backed/realistic-size** promotion on a small synthetic backed store (the Marrow analogue).
4. **MuData `peaks` feature axis** (ATAC modality from a `.h5mu`, not Signac) — `test_mudata.py` covers
   genes/proteins only; add a synthetic RNA+ATAC MuData case (peak-coord feature labels + `(cells, peaks)`).
5. **MuData per-modality PCA loadings of *different* dimensionality** (the bug-5 collision: rna 50 PCs /
   prot 31 PCs both keyed `varm['PCs']`) — `test_mudata.py`'s MOFA case shares one *equal-length* factor
   axis, so it does NOT guard the unequal-length namespacing fix. Add a synthetic two-modality case whose
   `varm['PCs']` differ in width to lock in the fix against regression.
