# Corpus acquisition catalog — real-data sources to grow the local test corpus

*Discovery + cataloging pass (2026-06-13). Prioritized sources for the **LOCAL** test corpus and the
breadth sweeps (`conformance/sweep/`). **None of this is committed to or downloaded by github CI** — CI
runs the synthetic-but-faithful fixtures (`python/tests/synth.py`) + small committed fixtures. This file
catalogs *where real data lives and exactly how to fetch it* so a follow-up agent or the user can
populate the gitignored `testdata/` dir. See `python/tests/CORPUS.md` (curated loaders),
`conformance/sweep/CATALOG.md` (breadth sweeps), and `misc/format_coverage.md` (gap analysis).*

## Why expand (the complaint)

The corpus today is thin and lopsided toward one shape (one `cells × genes` matrix + obs/var + a couple
of embeddings + DE):

| format | what we have now |
|---|---|
| AnnData | pbmc68k_reduced, pbmc3k(_processed), TMS-Marrow (local 1.2 GB), pancreas-velocity, paul15, burczynski06 |
| MuData | minipbcite (411 cells, RNA+ADT) + a tiny committed CITE-seq `.mtx` fixture |
| Seurat | ~10 SeuratData datasets installed (lazy) + pbmcMultiome ATAC; Azimuth refs unloaded |
| SCE | the 61 Bioconductor `scRNAseq` datasets (56 pass) |
| Conos | 2 local `.rds` (conI 4-sample, con 8-sample) |

**Systematically under-represented scenarios:** multimodal beyond one CITE-seq object, *all* of spatial,
trajectory/velocity beyond one pancreas, integration/collection benchmarks (the Conos model's home turf),
perturbation (zero), cell-cell communication & eQTL (arity-3 tensors, zero), and the long tail of
**format/version variety** (old `.h5ad`, Seurat v2/v3/v5 serialized, `.loom`, `.h5Seurat`, TileDB-SOMA /
cellxgene census). The catalog below targets each, with **several candidates per scenario**.

## What is already installed here (acquisition is cheap for these)

- **Python** (system 3.8): `scanpy 1.9.1`, `anndata 0.8.0`, `mudata 0.2.3`, `loompy 3.0.8`. **Missing**
  (need a venv from `/tmp/python`, the standalone py≥3.10): `scvelo`, `squidpy`, `spatialdata`/
  `spatialdata-io`, `cellxgene-census` / `tiledbsoma`, `pertpy`, `muon`.
- **R 4.4.1** in `.Rlib`: `SeuratData` (21-dataset manifest, several installed lazy), `scRNAseq`,
  `Signac`, `pagoda2`, `alabaster.*`, `zellkonverter`-adjacent. **Missing** (need install): `Seurat`,
  `SingleCellExperiment` runtime, `conos`, `SpatialExperiment`, `spatialLIBD`, `STexampleData`,
  `TENxVisiumData`, `zellkonverter`, `HDF5Array`, `BPCells`.
- **Local on-disk real data:** `~/p21/pagoda2/misc/conos_two_sample_integration/` (conos 2-sample +
  two pagoda2 `.rds`), `~/cacoa/age/tab.muris/` (TMS Marrow as `.h5ad`, **`.h5seurat`**, and
  `.rds` — three serializations of one dataset = a free format-variety case), `~/p2/examples/` (pagoda2
  objects incl. an `atac` dir), `~/p21/lstar/testdata/` (minipbcite, pancreas velocity).

---

## Tier A — Multimodal (one shared `cells` axis + multiple feature axes)

The L★ multimodal shape is **one shared cells axis + N feature axes** (`genes`/`proteins`/`peaks`/…).
We have exactly one real CITE-seq object end-to-end; we need ADT-at-scale, RNA+ATAC (multiome), and the
"more than two modalities" cases (TEA/ASAP/ECCITE).

### CITE-seq (RNA + ADT)
| dataset | source / method | size | organism/tissue | exercises / gap filled |
|---|---|---|---|---|
| **cbmc** (Stoeckius CITE-seq) | `SeuratData::InstallData("cbmc"); data("cbmc")` | 8,617 cells, RNA+ADT (13 ADT) | human cord blood | **already installed lazy**; real 2-Assay CITE-seq, the canonical ADT case. Promote from sweep-only to a curated R loader |
| **bmcite** | `SeuratData::InstallData("bmcite")` | 30,672 cells, RNA+25 ADT, wnnUMAP | human bone marrow | larger ADT panel + **WNN joint embedding** (a cross-modal embedding over the shared cells axis) + spca loadings on 2000/17009 (subset coverage) |
| **pbmc5k CITE (10x)** | 10x public: `5k_pbmc_protein_v3_filtered_feature_bc_matrix.h5` from 10xgenomics.com/datasets | ~5k cells, RNA + ~32 ADT | human PBMC | a **raw 10x feature-barcode `.h5`** with a `Gene Expression` + `Antibody Capture` feature-type split — exercises the importer's feature-type partitioning, not a pre-split object |
| **minipbcite** (have it) | already in `testdata/` | 411 cells | human PBMC | keep as the small `.h5mu` anchor |

### 10x Multiome (RNA + ATAC, shared cells, peaks axis)
| dataset | source / method | size | organism/tissue | exercises / gap filled |
|---|---|---|---|---|
| **pbmcMultiome** | `SeuratData::InstallData("pbmcMultiome")` → `pbmc.rna` (Assay) + `pbmc.atac` (**ChromatinAssay**) | 11,909 cells | human PBMC | **already installed lazy**; 108k peak ranges (seqnames/start/end) + fragment file ref — the only ChromatinAssay we have. Promote to a curated loader |
| **10k PBMC multiome (10x)** | 10x public `pbmc_granulocyte_sorted_10k` filtered `.h5` + ATAC fragments | ~10k cells, RNA + ~100k peaks | human PBMC | a **raw multiome `.h5`** (feature-type split RNA/peaks) + an external fragments TSV — the pre-Seurat importer path |
| **Chen 2019 SNARE-seq** / **SHARE-seq** | scglue / GEO `.h5ad` (paired RNA+ATAC) | varies | mouse | a non-10x paired RNA+ATAC encoding (different peak-axis provenance) |

### >2 modalities: TEA-seq / ASAP-seq / ECCITE-seq (3–4 feature axes)
| dataset | source / method | size | organism/tissue | exercises / gap filled |
|---|---|---|---|---|
| **thp1.eccite** | `SeuratData::InstallData("thp1.eccite")` | RNA + ADT + **HTO** + **GDO** (4 assays) | human THP-1 (ECCITE-seq CRISPR) | **already installed lazy**; the **4-modality** case (also a perturbation screen — see Tier E). No synthetic 4-assay fixture exists → top gap |
| **TEA-seq** (Swanson 2021) | GEO GSE158013 / scPerturb-adjacent `.h5` | ~10k cells, RNA+ADT+ATAC (tri-modal) | human PBMC | the only **RNA+protein+ATAC** trimodal shape (three distinct feature axes at once) |
| **ASAP-seq** (Mimitou 2021) | GEO GSE156478 | ADT + ATAC (protein+chromatin, no RNA) | human | a multimodal object with **no RNA modality** — tests that `genes` isn't assumed present |

---

## Tier B — Spatial (its own tier; coords now, images deferred)

Spatial is entirely absent from the corpus. The near-term L★ target is **coordinates as a 2-D embedding
+ spatial metadata + (for imaging-based platforms) per-molecule / per-cell feature axes**; images,
polygons, and multiscale rasters are deferred (`misc/format_coverage.md` Tier 3). Cover one example per
platform family so the coord/units/scalefactor handling is exercised broadly.

### Sequencing-based (spot/bead arrays)
| dataset | source / method | size | platform | exercises / gap filled |
|---|---|---|---|---|
| **stxBrain** | `SeuratData::InstallData("stxBrain")` — 4 sections (anterior/posterior ×2) | 12,167 spots | **10x Visium** (mouse brain) | **in the SeuratData manifest already**; a Visium SpatialObject with `obsm['spatial']` + scalefactors + a **multi-section collection** (4 paired slices = a spatial *collection*) |
| **stxKidney** | `SeuratData::InstallData("stxKidney")` | 1,438 spots | Visium (mouse kidney) | small Visium case |
| **ssHippo** | `SeuratData::InstallData("ssHippo")` | bead-level | **Slide-seqV2** (mouse hippocampus) | high-resolution bead array (no spot grid) — different coord density than Visium |
| **scanpy `visium_sge`** | `sc.datasets.visium_sge(sample_id=...)` — **27 real 10x samples** (breast/heart/lymph node/brain/kidney, human+mouse, incl. targeted-vs-parent pairs) | ~1–4k spots each | 10x Visium | **already in installed scanpy** — a whole Visium sub-sweep for free; targeted-vs-parent pairs add a panel-restriction case |
| **squidpy `slideseqv2()` / `visium_hne_adata()` / `visium_fluo_adata()`** | `sq.datasets.*()` (needs squidpy) | preprocessed | Slide-seqV2, Visium | pre-annotated AnnData with `obsm['spatial']` + image keys in `uns['spatial']` |
| **spatialLIBD** (Bioc) | `spatialLIBD::fetch_data("spe")` | 12 samples, ~48k spots | Visium (human DLPFC) | a **SpatialExperiment** (R side) — 12-sample spatial collection with manual layer annotations; exercises `read_sce`/SpatialExperiment coord handling |
| **STexampleData / TENxVisiumData** (Bioc) | `STexampleData::Visium_humanDLPFC()`, `TENxVisiumData::HumanBreastCancer()` etc. | varies | Visium (multi-tissue) | the canonical SpatialExperiment fixtures; TENxVisiumData has single- and **multi-section** experiments across organisms |
| **Visium-HD** | 10x public (`Visium_HD_Human_Colon_Cancer`) via `spatialdata-io.visium_hd()` or Space Ranger `.h5` | millions of 2µm bins (binned to 8/16µm) | Visium HD | **bin-resolution** (huge spot count) — a size + multi-resolution-binning stressor |
| **Slide-seqV2 (squidpy/Stickels)** | as above | ~40k beads | mouse | bead array scale |
| **Stereo-seq** | STOmics `.gef`/`.h5ad` via `spatialdata-io.stereoseq()` | very large, nanoscale | Stereo-seq (mouse/axolotl) | nanoball field-of-view coords; another distinct coord provenance |

### Imaging-based (single-molecule / single-cell; subcellular)
| dataset | source / method | size | platform | exercises / gap filled |
|---|---|---|---|---|
| **squidpy `merfish()`** | `sq.datasets.merfish()` (Moffitt hypothalamus) | ~73k cells | **MERFISH** | imaging-based AnnData: targeted gene panel (~150 genes), cell-level coords + a **per-molecule** table — the small-feature-axis spatial case |
| **squidpy `seqfish()`** | `sq.datasets.seqfish()` (Lohoff embryo) | ~57k cells | **seqFISH** | mouse embryo, 3D-ish coords |
| **Xenium** | 10x public (`Xenium_V1_human_breast`) via `spatialdata-io.xenium()` or `Seurat::LoadXenium()` | ~160k cells, 280-gene panel | **10x Xenium** | subcellular transcript table + cell/nucleus boundaries; the modern in-situ standard |
| **CosMx / SMI** | NanoString public via `spatialdata-io.cosmx()` or `Seurat::LoadNanostring()` | ~100k cells, 1000-plex | **NanoString CosMx** | a different imaging vendor's coord/QC layout |
| **MERFISH brain (Allen/Vizgen)** | Vizgen MERSCOPE public via `spatialdata-io.merscope()` | ~500k–1M cells | MERSCOPE | atlas-scale imaging spatial |
| **squidpy `imc()` / `mibitof()` / `four_i()`** | `sq.datasets.imc()` etc. | small | IMC / MIBI-TOF / 4i | **protein** (not RNA) spatial — a `proteins` feature axis over spatial coords; antibody panels |
| **STexampleData seqFISH/MERFISH/SlideSeq** | `STexampleData` loaders | small | multi-platform | SpatialExperiment-side imaging fixtures |

> **SpatialData (`.zarr`) route:** `spatialdata-io` readers (`xenium`, `visium_hd`, `merscope`,
> `cosmx`, `stereoseq`, `curio`, …) → write to SpatialData `.zarr`; published example `.zarr` stores are
> downloadable from the [SpatialData datasets page](https://spatialdata.scverse.org/en/latest/tutorials/notebooks/datasets/README.html)
> (e.g. `xenium_rep1_io`, `visium_hd_3.0.0_io`, `visium_associated_xenium_io`, `merfish`, `mibitof`,
> `mouse_liver`). These are the **native-Zarr spatial** comparison point for the L★ store.

---

## Tier C — Trajectory / RNA velocity

We have one velocity object (pancreas). scVelo ships several developmental datasets with
spliced/unspliced layers — each a different lineage topology (linear, branching, cyclic).

| dataset | source / method | size | organism/tissue | exercises / gap filled |
|---|---|---|---|---|
| **scv pancreas** (have small) | `scv.datasets.pancreas()` | 3,696 cells | mouse pancreas endocrine | branching endocrinogenesis; already the velocity anchor |
| **scv dentategyrus** | `scv.datasets.dentategyrus()` | ~3k cells | mouse dentate gyrus | neurogenesis branching; a second velocity topology |
| **scv dentategyrus_lamanno** | `scv.datasets.dentategyrus_lamanno()` | ~18k cells | mouse | the original La Manno velocyto loom-derived object → **`.loom` provenance** |
| **scv gastrulation / gastrulation_erythroid / gastrulation_e75** | `scv.datasets.gastrulation*()` | ~10–140k cells | mouse embryo | atlas-scale velocity; erythroid is a clean linear trajectory |
| **scv forebrain** | `scv.datasets.forebrain()` | ~1,720 cells | human wk10 embryo | human glutamatergic neuronal lineage; small |
| **scv bonemarrow** | `scv.datasets.bonemarrow()` | ~5,780 cells | human bone marrow | hematopoiesis |
| **CellRank pancreas/lung** | `cellrank` example data (builds on scVelo) | varies | mouse | adds a **transition matrix** (cell×cell `relation`) + macrostates + fate probabilities (group axis) |
| **Monocle3 cell_dataset** | `monocle3` example (`cds` / `.rds`) or its `.h5ad` exports | varies | — | a `cell_data_set` (SCE-subclass) carrying a **principal graph** / pseudotime — a different trajectory serialization |
| **dyngen/dyntoy synthetic** | `dyntoy` (R) | tiny | — | known-ground-truth branching for a topology-labeled `factor` axis |

> velocity exercises: `spliced`/`unspliced`/`Ms`/`Mu`/`velocity` **layers** (measures over cells×genes),
> `var['fit_*']` (per-gene measures), `uns['velocity_graph']`/`_neg` (cell×cell **relations** — note:
> in `uns`, not `obsp`), latent-time (cell measure). Multiple topologies stress the group/lineage axis.

---

## Tier D — Integration / collections / cross-species (the Conos / collection model)

This is the design's home turf ([[feedback-collection-not-tensor]]) and currently thinnest on *real,
benchmark-grade* examples. These exercise: a `samples` axis + per-sample `cells.<s>`/`genes.<s>` +
union `cells` axis + joint embedding + joint graph (`relation`), and **batch as a first-class factor**.

| dataset | source / method | size | organism/tissue | exercises / gap filled |
|---|---|---|---|---|
| **scIB pancreas** | figshare `https://figshare.com/ndownloader/files/24539828` (`.h5ad`) | ~16k cells, **8 technologies** | human pancreas | the canonical **batch-integration benchmark**: a known batch (`tech`) factor + harmonized cell-type labels; the cross-technology collection |
| **scIB immune (human)** | figshare 12420968 (`.h5ad`) | 33,506 cells, 10 donors, 5 studies | human bone marrow + blood | multi-study immune integration; 16 harmonized cell types |
| **scIB immune (human+mouse)** | figshare 12420968 | — | human + **mouse** | **cross-species** integration (orthology-mapped genes axis) |
| **scIB lung atlas** | figshare `https://figshare.com/ndownloader/files/24539942` (`.h5ad`) | 32,426 cells, 16 batches, 2 techs | human lung | batch + technology factors |
| **panc8** | `SeuratData::InstallData("panc8")` | 14,892 cells, 5 techs (SMARTSeq2/Fluidigm/CelSeq/CelSeq2/inDrops) | human pancreatic islets | **already in manifest**; the Seurat-side 8-tech-style integration collection (a split object) |
| **ifnb** | `SeuratData::InstallData("ifnb")` | 13,999 cells, CTRL vs STIM | human PBMC | **already installed lazy**; a 2-condition integration benchmark (stimulation as a factor) |
| **pbmcsca** | `SeuratData::InstallData("pbmcsca")` | 31,021 cells, **7 technologies** | human PBMC | **already in manifest**; HCA cross-technology benchmark — 7-way method comparison |
| **hcabm40k** | `SeuratData::InstallData("hcabm40k")` | 40,000 cells, 8 donors | human bone marrow | a clean donor-collection at 40k |
| **Tabula Muris Senis** (have Marrow local) | local `~/cacoa/age/tab.muris/`; full atlas from figshare/cellxgene | Marrow 40k local; full ~350k | mouse multi-tissue/age | **a collection flattened to one matrix** (use AS a collection: tissue+age+mouse-id factors); already the realistic-size anchor |
| **conos conI / con / acon** | local `~/p21/pagoda2/misc/.../*.rds` (conI 4-sample, con 8-sample; acon 8.7 GB) | 4–8 samples | — | **the native collection object** — keep as the Conos round-trip anchor; acon is the large stressor |
| **Tabula Sapiens** | figshare / cellxgene `.h5ad` | ~500k cells, ~24 organs | human multi-organ | a huge human atlas collection (organ/donor/method factors) — backed-read + collection at scale |

---

## Tier E — Perturbation (Perturb-seq / CRISPR / drug response)

**Zero coverage today.** These add a **perturbation/guide factor axis** (often a categorical with a
control level), guide-assignment metadata, and (for Perturb-seq) a guide×cell relation. The **scPerturb**
resource is the harmonized one-stop shop: 44 datasets, uniform `.h5ad`, harmonized `obs` keys
(`perturbation`, `perturbation_type`, `ncounts`, …).

| dataset | source / method | size | type | exercises / gap filled |
|---|---|---|---|---|
| **scPerturb collection** | scperturb.org / Zenodo RNA `10.5281/zenodo.7041848`, ATAC `10.5281/zenodo.7058381` (`.h5ad`) | 44 datasets | CRISPR + drug + epigenome | the **bulk perturbation sweep** — harmonized `obs` so one loader covers many; pick a small + a large |
| **Norman 2019** | scPerturb `NormanWeissman2019_filtered.h5ad` | ~70k cells, 105 single + combo CRISPRi | Perturb-seq (K562) | **genetic-interaction** combos → a perturbation factor with multi-guide levels |
| **Replogle 2022** | scPerturb `ReplogleWeissman2022_*` | >1.1M cells (genome-scale) | Perturb-seq (K562/RPE1) | **atlas-scale** perturbation; backed-read stressor + thousands of perturbation levels |
| **Papalexi 2021** (ECCITE) | scPerturb / `pbmc.eccite` | ~20k cells | ECCITE-seq (CRISPR + ADT + HTO) | perturbation **× multimodal** (also Tier A) |
| **thp1.eccite** | `SeuratData` (installed) | — | ECCITE-seq | already-local perturbation+4-modality (see Tier A) |
| **sci-Plex 3** | figshare `22122701` / scPerturb (`.h5ad`) | ~650k cells, 188 drugs × doses | drug-response | **dose** as an ordered factor + drug factor; high-cardinality treatment axis |
| **Frangieh 2021** | scPerturb `FrangiehIzar2021_*` | ~218k cells | Perturb-CITE-seq (melanoma) | CRISPR × CITE-seq (RNA+ADT+guide) |
| **Datlinger 2021** | scPerturb `DatlingerBock2021` | — | scifi-RNA-seq CRISPR | another guide-assignment layout |
| **pertpy built-ins** | `pertpy.data.*` (e.g. `norman_2019`, `papalexi_2021`, `sciplex3_*`, `dixit_2016`) | varies | mixed | programmatic loaders for many of the above (needs pertpy) |

---

## Tier F — Cell–cell communication & eQTL (arity-3 tensor cases)

These exercise relations/fields with **arity > 2** — the (source, target, ligand-receptor) and
(gene, SNP, cell-type) shapes that don't fit cells×genes. Mostly tabular results, not count matrices;
good stressors for the field/relation model beyond the 2-axis norm.

### Cell–cell communication (CCC)
| dataset | source / method | size | exercises / gap filled |
|---|---|---|---|
| **LIANA output** | `liana` (R/Py) on any annotated object → a tidy result table (source, target, ligand_complex, receptor_complex, scores) | small | a **(source_celltype, target_celltype, interaction)** arity-3 result — model as a relation over two `factor` axes + an interaction axis with score measures |
| **CellPhoneDB output** | run CellPhoneDB → `means.txt` / `pvalues.txt` (interaction × cell-type-pair) | small | the classic interaction×pair matrix; a 2-D table indexed by a composite (A|B) pair axis |
| **CellChat object** | `CellChat` `.rds` (`net`/`netP` arrays) | small | a **3-D array** (sender × receiver × pathway) — a true arity-3 tensor `field` |
| **Tensor-cell2cell** | `cell2cell` tensor output | small | an explicit **(sample × LR × sender × receiver)** 4-D tensor — the highest-arity case |
| **liana test data** | `liana::liana_test()` / `liana.testing` sample object | tiny | a ready-made small annotated object to generate CCC output on |

### sc-eQTL
| dataset | source / method | size | exercises / gap filled |
|---|---|---|---|
| **OneK1K** | onek1k.org / GEO GSE196830 (eQTL summary stats + the 1.27M-cell `.h5ad`) | 1.27M cells, 982 donors, 14 cell types | a **(gene × SNP × cell-type)** eQTL table — arity-3 field; plus a donor-collection atlas |
| **sc-eQTLGen** | eqtlgen.org (consortium summary stats) | large | cell-type-specific eQTL summary statistics tables |
| **eQTL Catalogue slice** | ebi.ac.uk/eqtl (per-study summary stats, `.tsv`) | per-gene | gene–variant association tables as a `relation` over a `genes` axis + a `variants` axis |

---

## Tier G — Format / version variety (the silent-corruption surface)

Different serializations of *the same logical data* are where round-trip bugs hide
([[feedback-graceful-version-recognition]]). Acquire one small example of each on-disk format/version.

### AnnData `.h5ad` schema versions
| version | source / method | exercises / gap filled |
|---|---|---|
| **anndata <0.7 (legacy)** | old GEO supplementary `.h5ad` (e.g. pre-2020 deposits) or write with `anndata==0.6` | pre-`encoding-type` layout: no `categories` group, raw groups, old sparse layout — importer must detect & degrade |
| **anndata 0.7** | `anndata==0.7` `write` of pbmc3k | the `encoding-type`/`encoding-version` transition |
| **anndata 0.8** | current installed (have several) | baseline |
| **anndata ≥0.10** | write with anndata 0.10+ (`CSRDataset` rename, new sparse) | the 0.10 sparse-class rename we already hit in CI; nullable extension dtypes land here |
| **`.h5ad` from scanpy `read_visium`** | `sc.datasets.visium_sge()` | spatial-flavored `.h5ad` (`uns['spatial']` scalefactors) |

### Seurat / SeuratObject serialized `.rds` versions
| version | source / method | exercises / gap filled |
|---|---|---|
| **v2 Seurat (`seurat` S4)** | old published `.rds` (pre-2018, e.g. archived pbmc3k v2) | the oldest class layout (slots differ entirely) — graceful-recognition stressor |
| **v3/v4 `Assay`** | `pbmc3k.final` (installed), `cbmc`, most SeuratData | the dominant published class |
| **v5 `Assay5`** | `UpdateSeuratObject` or a v5-authored object | layered assay (`Assay5 ⊄ Assay`) |
| **v5 split** | a v5 object with split layers (per-sample) | the collection-as-layers case |
| **`SCTAssay`** | any SeuratData run through SCTransform | residuals/corrected-counts slots |
| **`ChromatinAssay`** | `pbmcMultiome` `pbmc.atac` (installed) | Signac peak ranges + fragments |
| **`.h5Seurat`** | local `~/cacoa/age/tab.muris/Marrow.h5seurat` | the **SeuratDisk HDF5** serialization (vs `.rds`) — a different on-disk Seurat format we already have locally |

### SingleCellExperiment / SummarizedExperiment variants
| variant | source / method | exercises / gap filled |
|---|---|---|
| **plain SCE** | the 61 `scRNAseq` datasets (swept) | baseline — already covered |
| **SE-not-SCE** | `scRNAseq::ReprocessedFluidigmData()` | a `SummarizedExperiment` lacking SCE accessors (already handled; needs a CI fixture) |
| **SCE with altExps** | `scRNAseq` CITE datasets | ADT/ERCC as altExps → feature axes |
| **SpatialExperiment** | spatialLIBD / STexampleData | the spatial SCE subclass (Tier B) |
| **`cell_data_set` (Monocle3)** | monocle3 `.rds` | an SCE subclass with a principal graph (Tier C) |

### Other interchange formats
| format | source / method | exercises / gap filled |
|---|---|---|
| **`.loom`** | `loompy` (installed 3.0.8); scVelo `dentategyrus_lamanno`, velocyto outputs, or `sc.read_loom` | the **loom** layout (row=genes, col=cells, attrs) — a transposed-orientation importer path; we have loompy locally |
| **`.h5mu`** | minipbcite (have); any MuData write | multimodal HDF5 (covered, keep) |
| **Conos / pagoda2 `.rds`** | local objects (have) | the native collection + the p2 viewer object |
| **TileDB-SOMA / cellxgene census** | see Tier H | the columnar/Arrow-backed atlas format |

---

## Tier H — cellxgene census (TileDB-SOMA): a programmatic firehose for size + diversity

The single biggest **programmatic** source — tens of millions of cells, hundreds of datasets, fully
queryable. Pull **diverse, controlled slices** rather than whole datasets; each slice is real, current,
and arrives as an `AnnData`. Needs `cellxgene-census` (+ `tiledbsoma`) in a py≥3.10 venv.

```python
import cellxgene_census
with cellxgene_census.open_soma(census_version="stable") as census:   # pin a dated version for reproducibility
    # 1) tiny tissue/cell-type slice (a few hundred cells)
    adata = cellxgene_census.get_anndata(
        census, organism="Homo sapiens",
        obs_value_filter="tissue_general == 'lung' and cell_type == 'B cell' and disease == 'COVID-19'",
        column_names={"obs": ["assay","cell_type","tissue","disease","sex","dataset_id"]})
    # 2) row-range slice for a deterministic small pull:  obs_coords=slice(0, 1000)
    # 3) a whole curated dataset by dataset_id (see census_datasets table) → atlas-scale, read backed
```

| slice recipe | size | exercises / gap filled |
|---|---|---|
| tissue × cell_type × disease filter | hundreds–thousands | tiny real slices on demand; many organisms/tissues/assays → factor-axis diversity |
| `obs_coords=slice(0, N)` | exactly N | deterministic small pulls for fixtures |
| a single `dataset_id` | 10k–1M+ | a real curated dataset, **backed** read (size stressor) |
| `census["census_data"]["homo_sapiens"]` X stream | tens of millions | the out-of-core / `stream_col_stats` ceiling |
| **CELLxGENE-hosted embeddings** | per-dataset | precomputed scVI/Geneformer embeddings → embedding fields from a non-PCA source |

> Census also exposes **mouse** (cross-species) and a `census_datasets` table to discover by
> collection/tissue/assay. It is the recommended way to hit "millions of cells" and "many tissues"
> without curating dozens of downloads.

---

## Tier I — Size variety (tiny → atlas)

| tier | example | source | use |
|---|---|---|---|
| **tiny** (10²–10³) | scanpy `krumsiek11`/`toggleswitch`/`blobs`; STexampleData small; LIANA test obj; census `slice(0,500)` | installed / cheap | fast unit fixtures, edge cases |
| **small** (10³–10⁴) | pbmc3k, cbmc, ifnb, panc8, scv pancreas | installed | the working corpus |
| **medium** (10⁴–10⁵) | pbmcsca (31k), bmcite (31k), TMS Marrow (40k, local), Norman (70k), MERFISH (73k) | installed/download | realistic-size streaming/threading |
| **large** (10⁵–10⁶) | sci-Plex (650k), Tabula Sapiens (500k), OneK1K (1.27M) | download/census | backed read; collection gather |
| **atlas** (10⁶–10⁷) | Replogle (1.1M), cellxgene census streams (tens of M) | census | the out-of-core ceiling; bounded-memory reductions |

---

## Prioritized top-20 "acquire next"

Ranked by (gap severity × acquisition ease). Everything marked **installed** is a near-free win
(already in `.Rlib`/scanpy — just promote from sweep-only to a curated, version-logging loader).

| # | scenario | dataset | method | ~size | gap filled |
|---|---|---|---|---|---|
| 1 | Multimodal RNA+ATAC | **pbmcMultiome** | `SeuratData` (installed) | 11.9k | only real ChromatinAssay / multiome; no synthetic fixture exists |
| 2 | Multimodal CITE-seq | **cbmc** + **bmcite** | `SeuratData` (installed) | 8.6k / 30.7k | real 2-Assay ADT + WNN joint embedding; ADT at scale |
| 3 | Multimodal (4-way) | **thp1.eccite** | `SeuratData` (installed) | — | 4-modality (RNA+ADT+HTO+GDO) + perturbation; biggest multimodal gap |
| 4 | Spatial Visium | **stxBrain** (4 sections) | `SeuratData` | 12.2k | first spatial object; multi-section spatial collection |
| 5 | Spatial (sweep) | **scanpy `visium_sge`** (27 samples) | `sc.datasets.visium_sge()` (installed) | 1–4k ea | a whole Visium sub-sweep for free |
| 6 | Spatial imaging | **squidpy merfish + seqfish** | `sq.datasets.*()` (needs squidpy) | 57–73k | imaging-based, targeted panel, molecule table |
| 7 | Integration benchmark | **scIB pancreas (8-tech)** | figshare ndownloader 24539828 | 16k | canonical batch-integration; known `tech` factor |
| 8 | Integration / collection | **panc8** + **pbmcsca** | `SeuratData` | 14.9k / 31k | multi-technology collections (Seurat side) |
| 9 | Perturbation | **scPerturb Norman2019** | Zenodo 7041848 (`.h5ad`) | 70k | first Perturb-seq; combo-guide factor axis |
| 10 | Trajectory/velocity | **scv dentategyrus + gastrulation_erythroid** | `scv.datasets.*()` | 3k / 10k | second + linear velocity topology |
| 11 | Census firehose | **cellxgene-census slices** | `cellxgene_census.get_anndata(...)` | any | programmatic diversity + atlas scale; many tissues/organisms |
| 12 | Format variety | **TMS Marrow `.h5seurat`** | local `~/cacoa/age/tab.muris/Marrow.h5seurat` | 40k | SeuratDisk HDF5 format (have it locally, unused) |
| 13 | Format variety | **`.loom`** | `loompy` (installed) / `scv.dentategyrus_lamanno` | ~18k | transposed loom orientation |
| 14 | Spatial Slide-seq | **ssHippo** | `SeuratData` | bead-level | high-res bead array (vs spot grid) |
| 15 | Spatial SCE | **spatialLIBD (12 samples)** | `spatialLIBD::fetch_data("spe")` | 48k | SpatialExperiment 12-sample collection (R coord path) |
| 16 | CCC (arity-3) | **LIANA / CellChat output** | run `liana`/`CellChat` on cbmc | small | first arity-3 (source×target×interaction) tensor |
| 17 | Format version | **old `.h5ad` (<0.7) + ≥0.10** | GEO legacy + anndata 0.10 write | small | schema-version recognition span |
| 18 | Spatial Xenium | **10x Xenium breast** | `spatialdata-io.xenium()` / `Seurat::LoadXenium()` | 160k | modern subcellular in-situ; transcript table |
| 19 | Perturbation (atlas) | **scPerturb Replogle2022** | Zenodo 7041848 | 1.1M | atlas-scale perturbation; backed-read stressor |
| 20 | Cross-species | **scIB immune human+mouse** | figshare 12420968 | — | cross-species integration (orthology gene axis) |

### Fastest wins (already installed — promote to curated loaders, ~0 download)
pbmcMultiome, cbmc, bmcite, thp1.eccite, panc8, pbmcsca, ifnb, hcabm40k, celegans.embryo, stxBrain,
stxKidney, ssHippo (all SeuratData); `sc.datasets.visium_sge` (27 Visium); the local
`Marrow.h5seurat` / TMS `.rds`; `loompy`. These alone fill multimodal, spatial-Visium, integration, and
two format-variety gaps **without a single new download**, just by lifting them from sweep-only into
curated, version-logging corpus loaders.

---

## Notes / constraints honored

- **Local-only.** Nothing here is added to CI; CI keeps running synthetic-faithful fixtures. This
  catalogs *sources*, it does not download large data (per the task).
- **Collections stay collections** ([[feedback-collection-not-tensor]]): scIB/panc8/pbmcsca/conos/TMS
  are to be loaded **as** multi-sample collections (batch/tech/donor factors + per-sample axes), never
  flattened to one aligned matrix (TMS-Marrow's single-matrix form is the explicit exception, labeled).
- **Graceful version recognition** ([[feedback-graceful-version-recognition]]): Tier G exists precisely
  to feed the importers structurally different serializations of the same logical object.
- **Light deps** ([[feedback-lstar-light-deps]]): squidpy/spatialdata/scvelo/pertpy/cellxgene-census
  are acquisition-side tools for *generating* local fixtures, not lstar runtime deps.
</content>
</invoke>
