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

## Conos / pagoda2 (`sweep_conos.R`)

| real object | exercises | CI synthetic |
|---|---|---|
| `conI.rds` (4 samples), `con.rds` (8 samples) | collection (multi-sample, heterogeneous), joint graph, per-sample embeddings | `collection.sh` synthetic collection |
| real pagoda2 object | — | **gap** (no real object yet; mock pagoda2 covers the viewer schema) |

## Gaps this catalog surfaces (synthetic representation missing for a handled real structure)

These real structures **work** in the profile (sweep-verified) but have **no synthetic CI fixture**, so CI
can't guard them against regression — prioritized to add as small synthetic cases:

1. SCE **NULL dimnames** · **S4 Rle / nested cols** · **SE-not-SCE** — add to `sce_versions.sh`.
2. Seurat **ChromatinAssay** (synthetic peaks + ranges, gated on Signac in the local tier) and a
   **4-assay** (ECCITE-style) case.
3. **Backed/realistic-size** promotion on a small synthetic backed store (the Marrow analogue).
