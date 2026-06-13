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

## Spatial — AnnData (`sweep_spatial.py`) — **13/13 PASS**

Conceptual spatial support (Tier 3): `obsm['spatial']` → observed coordinate axis (subtype `spatial`),
round-trips to `obsm['spatial']`; `uns['spatial']` (images/scalefactors) stays in the lossless passthrough
(deferred from typing, but byte-for-byte preserved). Local-only `.h5ad` under `testdata/spatial/`.

| real dataset | source/format | exercises | CI synthetic |
|---|---|---|---|
| `V1_*` Visium (9) | scanpy `visium_sge` `.h5ad` | `obsm['spatial']` → observed coord axis (subtype `spatial`); `uns['spatial']` scalefactors+images survive in passthrough; human+mouse; breast/heart/lymph/kidney/brain | (**gap** — no synthetic spatial fixture: a synth AnnData with `obsm['spatial']` + a `uns['spatial']` scalefactors/image stub would lock in the observed-coord-axis + passthrough behavior) |
| `Targeted_/Parent_…Cerebellum` | Visium targeted-vs-parent pair | panel restriction (1186 vs 36601 genes) over the same spots/coords | (**gap** — same spatial fixture; the panel-restriction is just a gene-axis subset) |
| `V1_Mouse_Brain_Sagittal_Posterior` (+ `_Section_2`) | Visium multi-section pair | two slices of one experiment = a spatial **collection** (each its own object/coords) | `collection.sh` (collection shape) + the spatial **gap** for the per-slice coords |
| `sq_merfish` | squidpy MERFISH `.h5ad` | imaging-based; **`spatial` + `spatial3d`** (3-D coords; 3d kept but untyped-as-spatial); 150-gene panel | (**gap** — spatial fixture; a *second/3-D* coordinate set untyped-as-spatial is a recorded behavior, not yet synth) |
| `sq_seqfish` | squidpy seqFISH `.h5ad` | imaging-based mouse embryo; `spatial` + `X_umap` | (spatial **gap**) |
| `sq_slideseqv2` | squidpy Slide-seqV2 `.h5ad` | bead array (41k); **spatial neighbor graph** (`spatial_connectivities`/`_distances` obsp → relation, round-trips); `deconvolution_results` (composition kept as embedding) | spatial-graph → `obsp`-relation path is covered by graph fixtures; spatial coords are the **gap** |
| `sq_imc` | squidpy IMC `.h5ad` | **protein** spatial (34-marker panel), no `uns['spatial']` | (spatial **gap**; protein-as-genes-axis is encoding-irrelevant here) |

> **Deferred-and-recorded (correct), not lost:** Visium tissue images + scalefactors survive verbatim in
> the lossless passthrough (verified identical after `write_anndata`); the spatial neighbor graph survives
> as a cells×cells relation. **Not silently lost** — `ds.dropped` is empty for all 13.

## Spatial — Seurat (`sweep_spatial.R`) — 5/5 load, **coords captured (fix #6)**

| real dataset | classes | exercises | CI synthetic |
|---|---|---|---|
| `stxBrain` anterior1/2 + posterior1/2 | `Assay5` + `VisiumV2` image | 4 Visium sections = a multi-section collection; coords in `so@images` (`GetTissueCoordinates` + `ScaleFactors`) | `seurat_versions.sh` spatial-coords (FOV) case — coords → `spatial` axis, pixels → `dropped` (fix #6) |
| `ssHippo` | `Assay` + `SlideSeq` image | high-res Slide-seqV2 (53,173 beads); coords in `so@images` | covered by the same FOV synthetic case |

> **Fixed** (REPORT bug #6, `14b0225`): `read_seurat` now mirrors the AnnData path — `so@images`
> coordinates (`GetTissueCoordinates()`) become a `spatial` **observed coordinate axis** (subtype
> `spatial`; multi-section subsets use partial coverage), and the pixel images are recorded in
> `ds$dropped` (no longer silent). The CI synthetic is the `CreateFOV` spatial-coords case in
> `seurat_versions.sh`; images themselves stay deferred, as on the AnnData side.

## Perturbation (`sweep_perturbation.py`) — **3/3 PASS**

scPerturb harmonized `.h5ad` (Zenodo 7041848). Exercises FACTOR-AXIS induction at scale (categorical
perturbation/guide → derived factor axis) + bounded-memory backed-read round-trip. Local-only under
`testdata/perturbation/`.

| real dataset | source/format | exercises | CI synthetic |
|---|---|---|---|
| `NormanWeissman2019_filtered` | scPerturb `.h5ad` (111k cells) | **high-cardinality factor induction**: `perturbation`=**237** (single+combo CRISPRa) + `guide_id`=290; `nperts` combo arity; backed read → `write(stream=True)` | `synth.py` factor-axis path covers induction; a **237-level / combinatorial** factor axis at this scale is **sweep-only** (gap — a synth obs categorical with hundreds of levels would guard induction at cardinality) |
| `SrivatsanTrapnell2020_sciplex2` | scPerturb `.h5ad` (24k cells) | **dose as ordered factor** (`dose_value`=8) × drug (`perturbation`=5); also surfaced finding #7 (categorical `var['ensembl_id']` → 58,302-level degenerate factor axis) | the identifier-induction case is now guarded + covered by `test_induce.py`'s near-unique case (fix #7); ordered-dose factor synth fixture still a **gap** |
| `DatlingerBock2021` | scPerturb `.h5ad` (39k cells) | combo 2nd-guide column (`perturbation_2`); 384-level `sample`; scifi-RNA-seq guide layout | (induction covered generically; multi-perturbation-column layout is sweep-only) |

> Induction is dtype-driven (a categorical column auto-induces); these objects confirm it holds for
> hundreds–thousands of levels and that the heavy count matrix round-trips with bounded memory from a
> backed read. The **degenerate-factor-from-identifier** behavior (finding #7) is now **guarded** in
> `model.py` (`14b0225`): auto-induction skips a categorical whose levels are near-unique relative to its
> span, so an identifier column no longer mints a giant factor axis (an explicit `induce()` still works).

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
6. **Spatial coordinate axis + `uns['spatial']` passthrough** — no synthetic AnnData fixture carries
   `obsm['spatial']` + a `uns['spatial']` (scalefactors/image stub). Add a small synth case to guard the
   observed-coord-axis (subtype `spatial`) + lossless-passthrough behavior the 13 real spatial objects pass.
7. **Seurat `so@images` (Visium/Slide-seq coords)** — currently a profile *gap* (REPORT bug #6: coords
   silently dropped), not just a fixture gap. Once the Seurat profile captures `GetTissueCoordinates()` into
   a `spatial` coord axis, add a synthetic Seurat spatial object to guard it.
8. **High-cardinality + combinatorial factor induction** — Norman2019's `perturbation` (237 single+combo
   levels) / `guide_id` (290) are sweep-only; a synth obs categorical with hundreds of levels (incl. combo
   "A+B" labels) would guard induction at cardinality. Pairs with finding #7 (degenerate factor from a
   categorical *identifier* column — decide a cardinality heuristic in model.py before adding a fixture).
