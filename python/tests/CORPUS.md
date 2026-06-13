# Test corpus — two tiers: real (local) + synthetic-faithful (CI)

The profile tests are grounded in **real** datasets (real tool output, real structures), not
hand-fabricated simulations. But the corpus runs in **two tiers**:

- **Locally (default):** the loaders in [`corpus.py`](corpus.py) fetch/derive **real** data, cached in
  the gitignored `testdata/` dir; the large atlas is local-only. This is the honest grounding, and the
  broad real breadth lives in [`../../conformance/sweep/`](../../conformance/sweep). Real data has
  repeatedly caught bugs fabricated tests could not (see "Bugs caught" below).
- **In CI (`LSTAR_SYNTHETIC_CORPUS=1`):** each loader returns a **synthetic-but-faithful** stand-in from
  [`synth.py`](synth.py) — synthetic counts pushed through the *real* scanpy/mudata pipeline, so the same
  library code produces the same structures (categoricals, `*_colors`, `uns['pca']`, the `neighbors`
  OverloadedDict, `rank_genes_groups` dtypes, RNA+ADT modalities, `velocity_graph`). **No real datasets
  are committed to or downloaded by github** (they're large/slow). Keeping the synthetic fixtures
  structurally representative of the real corpus is the explicit contract; the local real runs verify it.

The table below describes the **real** (local) sources; the CI column names the synthetic stand-in.

## Datasets

| loader | real source (local) | CI synthetic stand-in | what it grounds |
|---|---|---|---|
| `pbmc68k_reduced()` | `sc.datasets.pbmc68k_reduced()` (~few MB, cached) | `synth.pbmc68k_like()` | obs categoricals (bulk_labels/phase/louvain) + `.raw` (same gene set as X); `*_colors`; `uns['pca']`; a **one-vs-rest** `rank_genes_groups` (method=**logreg** → names+scores only); the `neighbors` `OverloadedDict` |
| `pbmc3k_processed()` | `sc.datasets.pbmc3k_processed()` (~23 MB, cached) | `synth.pbmc3k_like()` | louvain/colors/pca + `.raw` over a **divergent** gene set (X=HVG subset); the substrate for DE variants |
| `pbmc3k_with_de(method, reference)` | scanpy `rank_genes_groups` on pbmc3k | (runs on the synthetic pbmc3k) | **t-test/wilcoxon** DE (full names/scores/logfoldchanges/pvals/pvals_adj) and **pairwise** (`reference=<group>`) |
| `pancreas_velocity()` | scVelo pipeline on `scv.datasets.pancreas()` (subsampled, cached; needs scvelo) | `synth.velocity()` | `spliced`/`unspliced` layers, `clusters` categoricals + colors, and `uns['velocity_graph']` (+ `_neg`) cell×cell |
| `citeseq_matrices()` / `citeseq_mudata()` | subsampled from real `minipbcite` (80 cells × 27 genes + 29 proteins) | `synth.citeseq_*` | RNA+ADT multimodal → genes/proteins feature axes; shared by the Python MuData test and the R Seurat/SCE tests (one generator, both languages) |
| `minipbcite()` | downloaded real CITE-seq MuData (411 cells, RNA+ADT) | `synth.citeseq_mudata_annotated()` | global `celltype` factor axis, per-mod PCA/uns, global joint embedding |
| `marrow_backed()` | local TMS Senis droplet Marrow h5ad (1.2 GB, **local-only**) | — (skips in CI) | realistic-size (40220 × 20138) read-backed; many real annotations → factor axes; promotions on a real atlas. *Skips with a note where absent.* |

## Key sourcing facts (verified)

- **scVelo's velocity graph is in `uns['velocity_graph']` / `velocity_graph_neg`** (cell × cell sparse),
  **not `obsp`** — the kNN graph (`distances`/`connectivities`) is in `obsp`, but the *velocity* graph
  is in `uns`. The profile types these uns square matrices as `relation`s (else they'd be dropped).
- **The raw scVelo pancreas h5ad is pre-velocity** (spliced/unspliced only); running the pipeline
  (`filter_and_normalize → moments → velocity → velocity_graph`) produces the graph. Locally,
  `pancreas_velocity()` regenerates it this way from `scv.datasets.pancreas()` and caches under
  `testdata/`; in CI, `synth.velocity()` constructs the same structure (spliced/unspliced layers + a
  cell×cell `uns['velocity_graph']`) without scvelo. Nothing is committed.
- **No standard public dataset reliably ships pandas nullable extension dtypes** (Int64/boolean/string);
  the documented path (used in `test_nullable`) is to derive nullable columns from a real dataset's real
  values + a real missingness rationale.
- Real `varm['PCs']` carries **NaN** (undefined loadings); NaN is a valid float that round-trips through
  the dense encoding (comparisons use `equal_nan=True`).
- AnnData allows the **same column name in `obs` and `var`** (e.g. `n_counts`); L* disambiguates the
  second (`n_counts.genes`) so neither is lost, with provenance preserving the native location.

## Version tracking (R objects)

"A Seurat object" is several structurally-different classes; the profile records the version of **every
object it reads** in `ds$profiles` and the corpus tests log it as they go:
- `object@<v>` — the object's serialized SeuratObject version (a v3/v4 object loaded under Seurat 5 still
  reports its own; `UpdateSeuratObject` does **not** promote `Assay`→`Assay5`).
- `assay@<name>:<class>` — **per assay** (multimodal mixes classes), e.g. `assay@RNA:Assay5`,
  `assay@ADT:Assay`, `assay@SCT:SCTAssay`. The assay *class* (not the object version) is the real
  version indicator (`Assay5` ⊄ `Assay`; `SCTAssay`/`ChromatinAssay` ⊂ `Assay` — branch SCT/Chromatin
  *before* the generic `Assay`).

| object | version recorded | covers |
|---|---|---|
| constructed v3 | `assay@RNA:Assay` | v3/v4 Assay path |
| constructed v5 | `assay@RNA:Assay5` | v5 layered path |
| v5 split | `assay@RNA:Assay5` (split layers) | integration / collection |
| SCT | `assay@SCT:SCTAssay` | SCTransform residuals |
| multimodal | `assay@RNA:Assay5 assay@ADT:Assay5` | per-modality feature spaces |
| real pbmc3k.final | `object@5.4.0 assay@RNA:Assay` | published v4 (updated) |
| real cbmc | `assay@RNA:Assay assay@ADT:Assay` | published CITE-seq |

## Bugs caught by real data (that fabricated tests missed)

1. `copy.deepcopy` on `uns` exploded on the `neighbors` `OverloadedDict` (cyclic back-reference to the
   AnnData/obsp) → spine-copy + leaf-reference instead.
2. Restoring `uns['neighbors']` crashed anndata's overloaded setter (obsp graphs injected by the getter
   became `None` through the passthrough) → capture from the raw backing dict.
3. `obs`/`var` same-name columns collided in the field namespace → silent data loss → disambiguate.
4. scVelo's `velocity_graph` (in `uns`) was being dropped → typed as a cell-cell relation.
