# Coverage tracker — formats, versions, cases, examples

A living matrix of what the conversion profiles must cover, with the **number of examples** backing each
row. Legend: ✓ covered · ◐ partial (works but a sub-case is recorded-not-typed) · ✗ gap (planned).
**Real** = a real published/atlas object (local or downloaded); **synth** = a constructed real-CLASS
fixture (the library's own constructors, deterministic, CI-able). See `python/tests/CORPUS.md` (Python
corpus) and `conformance/{seurat,sce}_versions.sh` + `real_corpus_r.sh` (R corpus).

## Example counts (current)

| corpus | real examples | synthetic fixtures |
|---|---|---|
| AnnData (.h5ad) | **4** (pbmc68k_reduced, pbmc3k_processed, TMS-Marrow[local], pancreas-velocity) + micropatterns[local] | DE/nullable derived from real values |
| MuData (.h5mu) | **2** (minipbcite CITE-seq[download]; citeseq fixture[committed]) | 0 |
| Seurat | **3** (pbmc3k.final v4[local], cbmc CITE-seq[local], citeseq fixture[committed]) | **4** (v3, v5, v5-split, SCT) |

> **Shared real CITE-seq fixture** (`python/tests/fixtures/citeseq/*.mtx`, 80 cells × 27 genes + 29
> proteins, subsampled from real minipbcite, ~27 KB) — drives **both** the MuData and Seurat RNA+ADT
> multimodal tests on the *same real* data. RNA+ADT is abundant; nothing here is simulated.
| SingleCellExperiment | **2** (ZeiselBrain[local], citeseq fixture[committed]) | structures (reducedDims/factors/metadata) on the real base |
| Conos / pagoda2 | 0 | 2 (mock collection, mock pagoda2) |

## AnnData (.h5ad)

| case | status | examples |
|---|---|---|
| encodings: dense / csr / csc / utf8 | ✓ | all real + synth |
| categorical (ordered + `-1` missing) | ✓ | real pbmc68k/pbmc3k/Marrow categoricals |
| nullable Int64 / boolean / string (values+mask) | ✓ | derived from real pbmc3k values (no public h5ad ships them) |
| numpy structured arrays (rank_genes_groups) | ✓ | real scanpy DE |
| versions: 0.7 legacy h5sparse / 0.8 / ≥0.10 | ✓ | test_versions + CI(latest) + local(0.8); backed-class rename handled |
| X / raw(divergent genes) / layers | ✓ | pbmc68k .raw; pancreas spliced/unspliced layers |
| obs/var (categorical/numeric/nullable/bool) | ✓ | real; **obs/var same-name collision** handled (`n_counts`/`n_counts.genes`) |
| obsm/varm (pca/umap/loadings) | ✓ | real; NaN in loadings handled |
| obsp/varp graphs (relations) | ✓ | real distances/connectivities |
| uns: params / `*_colors` / pca-variance | ✓ promoted | real pbmc68k/Marrow |
| uns: rank_genes_groups (t-test/wilcoxon=full; logreg=score-only; pairwise=passthrough) | ✓ | real pbmc3k/pbmc68k DE variants |
| uns: velocity_graph (cell×cell in uns) | ✓ | real scVelo pancreas fixture |
| uns: neighbors OverloadedDict / generic tail | ✓ passthrough | real (caught deepcopy + restore bugs) |
| backed mode (bounded-memory) | ✓ | Marrow backed + backed_targets |
| spatial obsm['spatial'] + uns['spatial'] images | ✗ separate tier | (Visium/Xenium — deferred) |

## MuData (.h5mu) — multimodal

| case | status | examples |
|---|---|---|
| modalities → canonical feature axes (RNA→genes, ADT→proteins, ATAC→peaks) | ✓ | minipbcite (RNA+ADT), synth |
| per-modality X / layers / var | ✓ | real + synth |
| per-modality obsm (own PCA/UMAP) | ✓ | minipbcite |
| global obs (categoricals → factor axes, `<mod>:` prefix) | ✓ | minipbcite celltype/leiden/louvain/leiden_wnn |
| global obsm (WNN / MOFA) | ✓ | minipbcite X_wnn_umap / X_mofa |
| per-modality + global uns | ✓ passthrough | minipbcite |
| obsmap/varmap: aligned cells | ✓ | minipbcite (all cells in both) |
| obsmap/varmap: **partial overlap** (0=absent) | ◐ per-mod `cells.<mod>` axis | no real partial example yet |
| RNA+ATAC multiome | ◐ | no real `.h5mu` multiome example yet (path = same as RNA+ADT) |

## Seurat

| case | status | examples |
|---|---|---|
| `Assay` (v3/v4) | ✓ | synth v3 + real pbmc3k.final |
| `Assay5` (v5) | ✓ | synth v5 |
| `SCTAssay` (residuals + SCTModel.list) | ◐ (data typed; SCTModel recorded) | synth SCT |
| `ChromatinAssay` (Signac scATAC) | ✗ (Signac not installed) | — (ranges/fragments → recorded) |
| v5 **split**/integration (per-sample layers → collection) | ✓ | synth split; write-back re-splits |
| multimodal RNA+ADT / RNA+ATAC | ✓ | synth RNA+ADT + real cbmc |
| reductions: embeddings / loadings / **stdev** | ✓ | synth + real |
| feature.loadings / scale.data over **HVG subset** | ◐ recorded (partial coverage) | real pbmc3k.final (2000/13714) |
| graphs (dgCMatrix) / `Neighbor` (nn.idx/dist) | ◐ graphs ✓; Neighbor ✗ | real graphs |
| meta.data factors / `Idents` (active) | ✓ | synth + real |
| images (VisiumV1/V2/FOV) / `@commands` | ✗ separate tier / not typed | — |
| version tracking (per-assay class + object version) | ✓ | recorded in `profiles`, logged by corpus |

## SingleCellExperiment

| case | status | examples |
|---|---|---|
| assays (counts/logcounts/multiple) | ✓ | synth + real |
| reducedDims (+ `rotation` attr → loadings) | ✓ | synth +reducedDims |
| **altExps** (ADT / spike-ins → feature axes) | ✓ | real citeseq ADT + real ZeiselBrain ERCC/repeat |
| colData / rowData factors | ✓ | synth + real |
| metadata (free-form) / colPairs | ◐ metadata recorded; colPairs ✗ | synth |

## Conos / pagoda2 / cross-cutting

| case | status | examples |
|---|---|---|
| collection (multi-sample, heterogeneous) | ✓ | synth collection + Seurat-v5-split |
| pagoda2 viewer schema + DE/pseudobulk over factor axis | ✓ | mock pagoda2 |
| real conos / pagoda2 objects | ✗ | (no real example yet) |
| factor-axis induction · nullable masks · lossless passthrough · DE/pseudobulk bundles | ✓ | across the above |
| cross-language Py ↔ C++ ↔ R ↔ JS/WASM (encodings + induced_by) | ✓ | conformance/*.sh + js.sh |

## Top gaps (prioritized)

1. **Spatial** (Visium/Xenium/CosMx) — its own tier (images, coordinate frames, molecules). Deferred per decision.
2. **Signac ChromatinAssay** — genomic ranges/fragments (needs Signac); peaks-as-assay already works.
3. **MuData partial-overlap** + **RNA+ATAC multiome** real example — source a `.h5mu` multiome.
4. **Seurat `Neighbor`** (nn.idx/dist) → relations; **`@commands`** provenance.
5. **Faithful partial coverage** (subset `scale.data`/loadings, SCT residuals) instead of recorded-as-dropped.
6. **Real conos/pagoda2** objects in the local corpus.
