# Corpus sweep — run the profiles against *many* real datasets

A handful of hand-picked examples is a demo, not a corpus. These harnesses sweep the profiles across
**tens to hundreds** of real datasets from the major repositories + local objects, recording pass/fail +
quirks, so we find the long tail of real-world structure (not just what we thought to construct).

The sweep corpus is **local-only** (datasets are large/many; cached in `testdata/` and Bioconductor/
SeuratData caches) — these are NOT in CI. CI runs the small committed fixtures + the curated downloads
(see `python/tests/CORPUS.md`, `COVERAGE.md`). The sweep is what you run locally to stay honest.

| harness | source | scope |
|---|---|---|
| `sweep_scrnaseq.R` | Bioconductor `scRNAseq` (61 real SCEs) | `read_sce`→`write_sce`→`lstar_write`; the first ~45, skipping >80k cells / missing-dep loaders |
| `sweep_anndata.py` | scanpy.datasets + local atlases | `read_anndata`→`validate`→`write` |
| `sweep_conos.R` | local real `Conos` objects (`*.rds`) | `write_conos`→`lstar_write` |
| `sweep_seurat.R` | SeuratData (24 real) | `read_seurat`→`write_seurat` (installs are heavy) |

Run e.g. `Rscript conformance/sweep/sweep_scrnaseq.R` (writes `/tmp/sweep_*.tsv`); aggregate into
`REPORT.md`. A **FAIL** = a profile bug to fix; a **LOADERR/SKIP** = the dataset itself needs an extra
package or is too big (not our bug — recorded separately).
