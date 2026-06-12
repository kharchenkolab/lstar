# Test corpus — real datasets the suite is grounded in

The profile tests are grounded in **real** datasets (real tool output, real structures), not
hand-fabricated simulations. Loaders live in [`corpus.py`](corpus.py); small datasets are fetched on
demand and cached in the gitignored `testdata/` dir, the velocity fixture is committed under
`fixtures/`, and the large atlas is local-only. Real data has repeatedly caught bugs fabricated tests
could not (see "Bugs caught" below).

## Datasets

| loader | source | size | what it grounds |
|---|---|---|---|
| `pbmc68k_reduced()` | `sc.datasets.pbmc68k_reduced()` | ~few MB, cached | real obs categoricals (bulk_labels/phase/louvain) + `.raw` counts over a divergent gene set; real `*_colors`; real `uns['pca']`; a real **one-vs-rest** `rank_genes_groups` (method=**logreg** → names+scores only); the `neighbors` `OverloadedDict` |
| `pbmc3k_processed()` | `sc.datasets.pbmc3k_processed()` | ~23 MB, cached | real louvain/colors/pca; the substrate for real DE variants |
| `pbmc3k_with_de(method, reference)` | scanpy `rank_genes_groups` on pbmc3k | computed | real **t-test/wilcoxon** DE (full names/scores/logfoldchanges/pvals/pvals_adj) and real **pairwise** (`reference=<group>`) |
| `pancreas_velocity()` | committed `fixtures/pancreas_velocity_small.h5ad` (~0.8 MB) | committed | **real scVelo output** (150 cells × 250 genes, subsampled): real `spliced`/`unspliced` layers, `clusters` categoricals + colors, and a real `uns['velocity_graph']` (+ `_neg`) |
| `marrow_backed()` | local TMS Senis droplet Marrow h5ad | 1.2 GB, **local-only** | realistic-size (40220 × 20138) read-backed; many real annotations → factor axes; promotions on a real atlas. *Skips with a note where absent (CI grounds on the downloadable corpus).* |

## Key sourcing facts (verified)

- **scVelo's velocity graph is in `uns['velocity_graph']` / `velocity_graph_neg`** (cell × cell sparse),
  **not `obsp`** — the kNN graph (`distances`/`connectivities`) is in `obsp`, but the *velocity* graph
  is in `uns`. The profile types these uns square matrices as `relation`s (else they'd be dropped).
- **The raw scVelo pancreas h5ad is pre-velocity** (spliced/unspliced only) and in **old (<0.7) on-disk
  format**; running the pipeline (`filter_and_normalize → moments → velocity → velocity_graph`) produces
  the graph. The committed fixture was generated this way (isolated venv; see commit history) and
  re-assembled in the anndata-0.8 main env for compatibility.
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
