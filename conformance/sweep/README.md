# Corpus sweep — run the profiles against *many* real datasets

A handful of hand-picked examples is a demo, not a corpus. These harnesses sweep the profiles across
**tens to hundreds** of real datasets from the major repositories + local objects, recording pass/fail +
quirks, so we find the long tail of real-world structure (not just what we thought to construct).

The sweep corpus is **local-only** (datasets are large/many; cached in `testdata/` and Bioconductor/
SeuratData caches) — these are NOT in CI. CI runs the small committed fixtures + the curated downloads
(see `python/tests/CORPUS.md`, `SUPPORT.md`). The sweep is what you run locally to stay honest.

| harness | source | scope |
|---|---|---|
| `sweep_scrnaseq.R` / `sweep_scrnaseq_driver.sh` | Bioconductor `scRNAseq` (61 real SCEs) | `read_sce`→`write_sce`→`lstar_write`; subprocess-isolated per dataset |
| `sweep_anndata.py` | scanpy.datasets + local atlases + **scVelo velocity** + **10x CITE-seq** | `read_anndata`→`validate`→`write` |
| `sweep_velocity.py` | scVelo real velocity/trajectory `.h5ad` (`testdata/velocity/`) | `read_anndata`→`validate`→`write` (spliced/unspliced layers) |
| `sweep_mudata.py` | real `.h5mu` (minipbcite + 10x CITE-seq + 10x multiome, `testdata/mudata_examples/`) | `read_mudata`→`write`→`write_mudata`; subprocess-isolated per dataset |
| `sweep_citeseq_10x.R` | 10x public CITE-seq `.h5` (5k/1k PBMC, 10k MALT) → Seurat RNA+ADT | `read_seurat`→`write_seurat`→`lstar_write` |
| `sweep_multiome_10x.R` | 10x public multiome `.h5` (PBMC 3k, brain 3k) → Seurat RNA + Signac ChromatinAssay | `read_seurat`→`write_seurat`→`lstar_write` |
| `sweep_integration.R` | SeuratData integration (ifnb / panc8 / pbmcsca) **read AS collections** | split-by-sample → per-sample axes → `lstar_write`→`lstar_read` |
| `sweep_conos.R` | local real `Conos` objects (`*.rds`) | `write_conos`→`lstar_write` |
| `sweep_seurat.R` / `install_and_sweep_seurat.R` | SeuratData (per-package enumeration) | `read_seurat`→`write_seurat` (installs are heavy) |
| `sweep_seurat_refs.R` | SeuratData Azimuth `*ref` atlases | records load-skip (need Azimuth loader) |

How to acquire the local datasets (all cached under the **gitignored** `testdata/`, never committed):
- **scVelo velocity** (`testdata/velocity/*.h5ad`): an isolated venv with scvelo, e.g.
  `python3 -m venv /tmp/scvelo_venv && /tmp/scvelo_venv/bin/pip install scvelo` then
  `scv.datasets.dentategyrus(file_path=...)` etc. (dentategyrus, gastrulation_erythroid, bonemarrow, pancreas).
- **10x CITE-seq / multiome** (`testdata/citeseq_10x/*.h5`, `testdata/multiome_10x/*.h5`): `curl` the public
  `*_filtered_feature_bc_matrix.h5` from `cf.10xgenomics.com` (URLs in this dir's git log / the catalog).
- **`.h5mu`** (`testdata/mudata_examples/*.h5mu`): built from the 10x `.h5` by splitting `feature_types`
  into modalities (`mudata.MuData({"rna":..., "prot"/"atac":...}).write(...)`); minipbcite auto-downloads.

Run e.g. `Rscript conformance/sweep/sweep_scrnaseq.R` or `python3 conformance/sweep/sweep_mudata.py`
(writes `/tmp/sweep_*.tsv`); aggregate into `REPORT.md`. A **FAIL/VALIDATE-ERR** = a profile bug to fix;
a **LOADERR/SKIP** = the dataset itself needs an extra package or is too big (not our bug — recorded
separately).
