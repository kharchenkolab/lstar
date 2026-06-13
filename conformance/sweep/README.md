# Corpus sweep — run the profiles against *many* real datasets

A handful of hand-picked examples is a demo, not a corpus. These harnesses sweep the profiles across
**tens to hundreds** of real datasets from the major repositories + local objects, recording pass/fail +
quirks, so we find the long tail of real-world structure (not just what we thought to construct).

The sweep corpus is **local-only** (datasets are large/many; cached in `testdata/` and Bioconductor/
SeuratData caches) — these are NOT in CI. CI runs the small committed fixtures + the curated downloads
(see `python/tests/CORPUS.md`, `SUPPORT.md`). The sweep is what you run locally to stay honest.

> **One entry point:** `bash conformance/sweep/retest_local.sh` runs the faithfulness guard + every
> cached-data sweep and ends with a TRIAGE of what to update. The recurring directive — when to run it,
> how to read the outcomes, how to update the synthetics and the docs — is **[`RETEST.md`](RETEST.md)**.
> The table + acquisition notes below are the reference the orchestrator and directive draw on.

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
| `sweep_spatial.py` (+ `sweep_spatial_fetch*.py`) | 10x Visium (`visium_sge`, 9 samples) + squidpy imaging (merfish/seqfish/slideseqv2/imc) | `read_anndata`→`validate`→`write`→`write_anndata`; asserts the `spatial` observed coord axis + `uns['spatial']` passthrough round-trip; subprocess-isolated per dataset |
| `sweep_spatial.R` | SeuratData stxBrain (4 Visium sections) + ssHippo (Slide-seqV2) | `read_seurat`→`lstar_write`→`write_seurat`; checks `so@images` coords survive as a `spatial` observed coord axis (now captured — REPORT bug #6 fixed in `14b0225`) |
| `sweep_perturbation.py` | scPerturb `.h5ad` (Datlinger / sciPlex2 / Norman2019) | backed `read_anndata`→`validate`→`write(stream=True)`; checks perturbation/guide categoricals induce factor axes (Norman: 237/290 levels); subprocess-isolated per dataset |

How to acquire the local datasets (all cached under the **gitignored** `testdata/`, never committed):
- **scVelo velocity** (`testdata/velocity/*.h5ad`): an isolated venv with scvelo, e.g.
  `python3 -m venv /tmp/scvelo_venv && /tmp/scvelo_venv/bin/pip install scvelo` then
  `scv.datasets.dentategyrus(file_path=...)` etc. (dentategyrus, gastrulation_erythroid, bonemarrow, pancreas).
- **10x CITE-seq / multiome** (`testdata/citeseq_10x/*.h5`, `testdata/multiome_10x/*.h5`): `curl` the public
  `*_filtered_feature_bc_matrix.h5` from `cf.10xgenomics.com` (URLs in this dir's git log / the catalog).
- **`.h5mu`** (`testdata/mudata_examples/*.h5mu`): built from the 10x `.h5` by splitting `feature_types`
  into modalities (`mudata.MuData({"rna":..., "prot"/"atac":...}).write(...)`); minipbcite auto-downloads.
- **Spatial Visium** (`testdata/spatial/V1_*.h5ad`, `Targeted_/Parent_*.h5ad`): `python3
  conformance/sweep/sweep_spatial_fetch.py` — `sc.datasets.visium_sge(sample_id=...)` (installed scanpy)
  downloads each Space Ranger output (filtered `.h5` + `spatial/` dir) and caches the assembled AnnData
  (`obsm['spatial']` + `uns['spatial']`) as `.h5ad`. A breadth spread of 9 of the 27 samples (tissues ×
  organisms + a multi-section pair + a targeted-vs-parent pair).
- **Spatial imaging** (`testdata/spatial/sq_*.h5ad`): needs a squidpy side-venv (squidpy ≠ an lstar dep) —
  `python3 -m venv /tmp/sq && /tmp/sq/bin/pip install squidpy`, then `/tmp/sq/bin/python3
  conformance/sweep/sweep_spatial_fetch_sq.py` caches `sq.datasets.{merfish,seqfish,slideseqv2,imc}()` as
  `.h5ad` (written by anndata 0.11 → a newer schema, a free format-variety bonus). The sweep itself reads
  them with lstar's system Python.
- **Spatial Seurat** (`stxBrain` + `ssHippo`): `R_LIBS=.Rlib Rscript conformance/sweep/sweep_spatial.R`
  auto-installs the SeuratData packages into `.Rlib` (local-only) on first run.
- **Perturbation** (`testdata/perturbation/*.h5ad`): `curl` the scPerturb harmonized `.h5ad` from Zenodo
  10.5281/zenodo.7041848 (record `7041849`), e.g.
  `curl -L -o testdata/perturbation/NormanWeissman2019_filtered.h5ad
  https://zenodo.org/api/records/7041849/files/NormanWeissman2019_filtered.h5ad/content`. The sweep uses
  Datlinger (34 MB) / sciPlex2 (139 MB) / Norman2019 (667 MB) — a varied-size spread; the genome-scale
  Replogle gwps (8.8 GB) is intentionally skipped.

Run them all at once with `bash conformance/sweep/retest_local.sh` (tolerant; covers whatever is cached;
see [`RETEST.md`](RETEST.md)), or one at a time, e.g. `Rscript conformance/sweep/sweep_scrnaseq.R` or
`python3 conformance/sweep/sweep_mudata.py` (each writes `/tmp/sweep_*.tsv`); aggregate into `REPORT.md`.
A **FAIL/VALIDATE-ERR** = a profile bug to fix; a **LOADERR/SKIP** = the dataset itself needs an extra
package or is too big (not our bug — recorded separately).
