# Local full-corpus retest — the recurring directive

This is the standing recipe for the periodic "prove the profiles against the **large local corpus**,
keep the **synthetic CI fixtures faithful**, and **sync the docs**" loop. It is run **rarely** (before a
release, after a profile change, when an upstream format library bumps, or when we add a corpus) — not in
CI. CI runs only the small committed fixtures + curated downloads (synthetic-only,
`LSTAR_SYNTHETIC_CORPUS=1`); this is the local tier that keeps that synthetic tier honest.

## TL;DR

```bash
bash conformance/sweep/retest_local.sh            # guard synthetics + sweep all cached corpora + triage
RETEST_HEAVY=1 bash conformance/sweep/retest_local.sh   # + the GB-scale SeuratData sweeps
RETEST_GATE=1  bash conformance/sweep/retest_local.sh   # + the full run.sh conformance gate (real loaders) first
```

The orchestrator is **tolerant**: a missing dataset / absent package / single crashing dataset degrades
to `SKIP`/`LOADERR`/`CRASH` for that item only, never aborts. It edits nothing — it ends with a **TRIAGE**
naming exactly what to update. It exits non-zero only on a real regression (a profile `FAIL`/`VALIDATE-ERR`).

## The two tiers (why this exists)

- **Synthetic tier (CI):** `python/tests/synth.py` generates synthetic counts and runs the *real*
  scanpy/mudata/Seurat pipelines over them, so the fixtures have real *structure* on synthetic *values* —
  offline, fast, committable. Gated by `conformance/run.sh` under `LSTAR_SYNTHETIC_CORPUS=1`.
- **Real tier (local, this directive):** the actual published objects (Bioconductor `scRNAseq`, SeuratData,
  scanpy/squidpy datasets, scPerturb, 10x `.h5`, real Conos `.rds`), cached under the gitignored
  `testdata/` + the `.Rlib`/Bioconductor/SeuratData caches. This is what catches the long tail.

The contract that ties them: **every structure the real corpus carries must be represented in a synthetic
fixture.** `python/tests/test_synth_faithful.py` enforces it automatically (real signature ⊆ synthetic
signature, over field/axis roles, subtypes, feature axes). When it fails, a synthetic has drifted — that's
the signal to update `synth.py`, not a profile bug.

## The loop (what the orchestrator runs, and what you do with each outcome)

1. **Acquire / refresh** the local datasets you want covered — see [`README.md`](README.md) ("How to
   acquire the local datasets"). Everything lands under `testdata/` (never committed). Skip what you don't
   have; the orchestrator covers whatever is cached.

2. **Faithfulness guard** (`[1]` in the orchestrator). 
   - **GAP** (`synthetic missing <role/axis/subtype>`) → a real structure left the synthetic fixtures.
     **Update `python/tests/synth.py`** to reproduce that structure (add the layer/obsm/uns/categorical
     that induces it), re-run the guard until clean. This is the **"update synthetics"** step.

3. **Corpus sweeps** (`[2]`) → one `/tmp/sweep_<modality>.tsv` per modality, each row a dataset with a
   `status`:
   - **FAIL / VALIDATE-ERR** → a **profile bug on real data**. Fix the profile
     (`python/src/lstar/profiles/*.py` or `R/R/profile_*.R`), **add a fixture that reproduces it** (a CI
     case in `conformance/*.sh` and/or a `synth.py` structure), then re-sweep. This is how every sweep-found
     bug (#1–#7 in REPORT.md) was closed.
   - **LOADERR / SKIP / CRASH / TIMEOUT** → the *dataset* needs a package, is absent, or is too big — **not
     our bug**. Record it in the REPORT tally; don't chase it.
   - **PASS on a structure no synthetic covers yet** → add that structure to `synth.py` (and a CI fixture)
     so it's gated going forward — then it's no longer sweep-only.

4. **Docs sync** (`[3]` triage). After any change, refresh the three coverage docs so they never lie:
   - `conformance/sweep/REPORT.md` — the **tally** (per-modality PASS/FAIL counts; each finding, flipped
     from "gap" to "**Fixed** (`<commit>`)" once closed; relabel any pre-fix output snapshot honestly).
   - `conformance/sweep/CATALOG.md` — the **corpus map** (dataset → what it exercises → which CI synthetic
     stands in for it; add a row for a newly-covered case).
   - `SUPPORT.md` — the **user-facing matrix** (`status` / `real` / `CI (synth)` columns; flip a `◐`/`✗`/gap
     to `✓` when closed, and leave `real` empty honestly if a genuine object is unobtainable).
   Also re-grep for stale enumerations (`grep -rn "v3/v4/v5" docs/ SUPPORT.md README.md .claude/skills`),
   the way the v2 pass did — a newly-supported version/format must be added everywhere coverage is listed.

5. **Re-run CI green** after the changes (the synthetic tier must still pass): push and confirm the
   `conformance` workflow is green — the synthetic fixtures you updated are what CI exercises.

## Conventions worth keeping (so this stays reliable)

- **Subprocess isolation** for big/fragile sweeps: one dataset per timeout-guarded child process, so a
  segfault/hang/OOM fails only that dataset (`sweep_scrnaseq_driver.sh` is the template; the Python sweeps
  isolate per-file internally). Never let one bad object kill the sweep.
- **Status vocabulary** (the `status` column every sweep emits): `PASS` · `FAIL` (round-trip/profile error)
  · `VALIDATE-ERR` (store fails `validate`) · `LOADERR` (dataset needs a missing package) · `SKIP` (absent)
  · `CRASH`/`TIMEOUT` (isolated). Triage keys on these.
- **Large R programs run from a temp file, never `Rscript -e '...'`** — a big `-e` string overflows R's
  ~8 KB command-line buffer and is silently ignored (R prints `WARNING: '-e ...'` and runs nothing),
  which then trips `set -o pipefail` at a trailing `grep` for a baffling failure. Quoted heredoc → temp
  file, `$RLIB` via `commandArgs`, `</dev/null` (see `seurat_versions.sh` / `seurat_v2.sh`).
- **Nothing real is committed.** `testdata/` is gitignored; fixtures in the repo are synthetic or tiny
  curated subsets. A genuine object that can't be obtained (e.g. an ancient-Seurat build that won't
  compile) leaves the `real` cell empty — say so, don't fabricate.
