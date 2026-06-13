# Testing & CI coverage — improvement plan

## Decisions (review) & execution order

Per review, the scope is locked and ordered. Principle reaffirmed: **all real data is local; nothing real
on github** — distill the format's lstar requirements into a *small synthetic* CI test instead.

1. **Structural hole: JS/WASM in CI** — build it.
2. **Structural hole: synthetic-faithfulness guard** — build it.
3. **Corpus catalog** — keep a systematic working file (`conformance/sweep/CATALOG.md`) mapping each real
   dataset → what it exercises → its CI synthetic stand-in.
4. **Partial-coverage `index` arrays** — implement typed partial-coverage in lstar core + tests
   (unblocks faithful partial-overlap multiome; was the cross-repo blocker).
5. **Joint-method storage shapes** (WNN / MOFA+ / totalVI) — not our job to *run* them, but it **is** our
   job to cover how their results are stored/shaped. Add examples + tests.
6. **Spatial — conceptual only** — named axes (a `spatial` coordinate axis) as recommended; **no** serious
   image/vendor support. Small synthetic test; note in SUPPORT.md.
7. **Old serialized-object versions** — install old packages on the side (local), implement support **iff
   the version is cleanly recognized**; note version support in SUPPORT.md. No real data on github.
8. **Version-matrix CI job** + coverage/native-R-validate/fuzz (the cheap durable items).

Everything real stays local (U1–U8 disposition): spatial real = deferred; Azimuth/old-objects/backed/big
atlases = local installs on the side; joint methods = we test storage shape only.

---


Analysis of the current test surface (the two-tier model, `SUPPORT.md`, the sweep `REPORT.md`, the
conformance suite, and the multimodal proposals) and a prioritized plan to harden it. Each item is
tagged **value** / **effort** and whether it belongs in **CI** (synthetic, fast, hermetic) or the
**local** real-corpus tier. The last section flags what is *unreasonable* to test in CI or even locally.

## Where we are

- **CI (synthetic-only, hermetic):** two jobs — `python-cpp` (Python + C++ + Py↔C++) and
  `r-cross-format` (R profiles + Py↔R↔C++). Fixtures come from `synth.py` (synthetic counts through the
  *real* scanpy/mudata/Seurat pipelines). No real data committed or downloaded.
- **Local real corpus:** `corpus.py` loaders + `conformance/sweep/` — real breadth (scRNAseq **56/61**,
  SeuratData lazy **10/10**, multiome, real Conos). This is what catches the long tail; it has caught
  **5 real profile bugs** no synthetic fixture had.
- **Cross-language conformance:** `conformance/*.sh` covers Py/C++/R for encodings, induction, nullable,
  aux, DE, chunked, kernels, block reader; `js.sh` covers JS/WASM.

## The two structural holes (highest leverage)

### H1 — JS/WASM is not in CI at all
`js.sh` (the browser/WASM reader/writer/view + WASM kernels) runs **only locally** (`run.sh`); neither CI
job invokes it. A JS/WASM regression ships unnoticed. **Fix:** add a `js-wasm` CI job — set up & **cache**
emsdk, build the WASM kernels, generate a store with the Python writer, run `js.sh`'s kernel/reader/view
tests. *Value: high · Effort: low-moderate (emsdk install is the only real cost, and it caches).* **CI.**

### H2 — no automated "synthetic faithfully represents real" guard
The core contract of the two-tier model — synthetic CI fixtures must structurally mirror the real corpus
— is enforced only by human judgment. Nothing fails when they drift (e.g. a new scanpy version changes
`rank_genes_groups` dtypes). **Fix:** a `test_synth_faithful.py` that, *when real data is available*
(skips in CI), reads both the synthetic stand-in and its real counterpart through the profile and asserts
the **structural signatures match** — same set of {field roles, axis roles, encodings, promoted-uns keys,
DE-bundle shape, feature-axis names} — comparing *structure, not values*. Run it in the local tier and as
a pre-release gate. *Value: high · Effort: moderate.* **Local (gates the contract).**

## CI hardening (synthetic, reasonable)

| # | item | value | effort | notes |
|---|---|---|---|---|
| C1 | **JS/WASM CI job** (H1) | high | low-mod | the single biggest hole |
| C2 | **Dependency version matrix** | med | low | CI pins one stack; add a job pinning `anndata==0.8`+`zarr<3` beside latest so version-recognition regressions surface (real old *serialized* objects stay local — see U6) |
| C3 | **Coverage measurement** | med | low | `pytest-cov` (Py) + `covr` (R) → artifact/summary; makes untested branches visible |
| C4 | **Native R `validate()`** | med | mod | R conformance currently round-trips to Python to validate; a native R validator self-checks the R writer earlier and lets `real_corpus_r.sh` assert locally |
| C5 | **Property/fuzz tests for encodings** | med | mod | randomized dense/CSR/CSC/categorical/nullable round-trips (hypothesis in Py) catch edge cases curated fixtures miss (empty axes, all-missing mask, 0-nnz, int64→BigInt boundaries) |
| C6 | **Fail-loud on unexpected SKIPs** | low | low | CI steps that SKIP optional backends print SKIP but pass; add a CI-only assertion that the *expected* synthetic cases ran (so a silent SKIP of a case we *can* test in CI is caught) |

## Coverage expansion (close typeable gaps + new tiers)

| # | item | tier | notes |
|---|---|---|---|
| E1 | **`scale.data` over HVG subset** (Seurat) | profile | currently `◐` recorded. Re-route like loadings — a measure over a `<assay>_scaled_features` subset axis — but mark the axis so `write_seurat` rebuilds it as the assay's `scale.data` layer, **not** a spurious assay. Needs the "subset-axis measure ≠ new modality" guard. |
| E2 | **Seurat `@commands` provenance** | profile | analysis history → capture in the aux passthrough (like `uns`), not dropped. Low value but completes "never lose". |
| E3 | **SCT residuals / `SCTModel`** | profile | `◐` recorded; faithfully typing the residual scale-factors is the remaining SCT work. |
| E4 | **Spatial tier** (own milestone) | both | Visium/Xenium/CosMx/Slide-seq. Design the profile: a `spatial` coordinate axis, image arrays, molecule tables. CI gets a **small synthetic** Visium-like fixture (coords + tiny image + a few molecules); real vendor objects are local-only. SeuratData `stxBrain`/`stxKidney`/`ssHippo` are already installed for the local side. |
| E5 | **MuData partial-overlap** | profile/CI | the `cells.<mod>` path exists but is untested; add a **constructed** partial-overlap MuData fixture now (don't wait for a real one). Full fidelity is gated on lstar typed `index` arrays (cross-repo — see U7). |
| E6 | **Real `.h5mu` multiome, real pagoda2 object** | local corpus | data gaps, not code; source and add to the sweep. |
| E7 | **Azimuth reference atlases** | local corpus | wire the Azimuth loader for SCTAssay/large-reference breadth (local only — see U3). |
| E8 | **Periodic sweep re-run** | local | repos drift; document a cadence (and ideally a local cron) to re-run `conformance/sweep/` and refresh `REPORT.md`, so new upstream structures are caught. |

## Prioritized order

1. **C1 (JS in CI)** and **H2 (faithfulness guard)** — close the two structural holes.
2. **E1 + E2** — finish the typeable Seurat gaps (consistency with the four just closed).
3. **C2 + C3 + C4** — version matrix, coverage, native R validate (cheap, durable).
4. **C5** — encoding fuzzing.
5. **E5** — constructed MuData partial-overlap.
6. **E4 (spatial tier)** — a milestone of its own; do after the above.
7. **E6/E7/E8** — corpus expansion + cadence.

## Unreasonable to test in CI (or even locally) — flagging per your ask

- **U1 — Real spatial imaging at vendor scale** (Xenium/CosMx/MERSCOPE/Visium-HD): multi-GB image
  tensors and vendor-specific on-disk formats. *Plan:* mock a **small** synthetic Visium-like fixture in
  CI; the real vendor-format matrix and real image sizes stay local/manual. Hosting them in CI is
  unreasonable.
- **U2 — Full `Seurat` umbrella on the GitHub runner** (for `SCTransform`): it installs as a binary but
  **fails to load** on the minimal runner (missing system libs for the heavy deps). Debugging the runner
  image is a time-sink with low payoff. *Recommendation:* keep SCTransform coverage **local-only**
  (it runs wherever full Seurat loads), and don't over-invest in making the runner load it.
- **U3 — Azimuth reference atlases in CI**: need the `Azimuth` package + GB-scale reference downloads +
  a finicky disk-dataset loader. Local-only.
- **U4 — Disk-backed real backends at scale** (BPCells / HDF5Array / backed AnnData on a real atlas): CI
  does a small **smoke**; genuine bounded-memory behavior on multi-GB data needs a real machine →
  local-only by design.
- **U5 — Multi-GB real atlases / large Conos** (TMS Marrow 1.2 GB, `acon.rds` 8.7 GB): size precludes
  CI; local-only (already the design).
- **U6 — Genuinely old *serialized* objects** (Seurat v2 `.rds`, anndata <0.7 h5ad): reading real
  ancient files needs ancient package versions that don't co-install cleanly with current ones. We test
  version **recognition** via the libraries' own constructors; reading real legacy files is a
  manual/local check, not CI.
- **U7 — lstar typed partial-coverage `index` arrays**: spec'd but unimplemented in lstar core, so a
  *fully faithful* partial-overlap multiome round-trip can't be CI-gated yet — it currently falls back to
  the per-modality `cells.<mod>` axis. This is a cross-repo dependency, not a test-harness problem.
- **U8 — Joint-method algorithms themselves** (totalVI / MOFA+ / WNN): not lstar's job. We test the
  **storage shape** of their outputs (embedding + per-facet loadings + modality weights), not the
  algorithms; reproducing the methods in CI is out of scope and unreasonable.
