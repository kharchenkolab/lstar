# Cross-language parity contract

lstar ships the same capabilities on four surfaces — **C++** (`core/`, the header-only `libstar`),
**Python** (`python/`, pybind11), **R** (`R/`, cpp11), and **JS/WASM** (`js/`, embind). The project's
core aim (see `.claude/CLAUDE.md`) is that these behave *identically*, ideally by sharing one
implementation rather than maintaining four that drift. This document is the contract that keeps them
from drifting, and the checklist for anyone adding or changing a shared feature.

## The rule

**One recipe, N thin bindings.** A shared computation lives in exactly one place:

- **Numeric / algorithmic logic → the C++ core.** Kernels (`col_sum_by_group`, `markers_one_vs_rest`,
  `overdispersion`), orientation flips (`csc_to_csr`, `csr_to_csc`), and the viewer cell reorder
  (`hilbert_index`, `cell_order_pos`, `viewer_cell_order`) are written once in `core/include/lstar/lstar.hpp`
  and bound into Python (`_accel.cpp`), R (`lstar_cpp.cpp`), and WASM (`lstar_wasm.cpp`). A binding may
  keep a numpy/Matrix fallback for the no-accelerator install, but it MUST match the core (guarded by a
  cross-impl conformance test).
- **Policy / heuristics that can't live in a kernel → a single declarative source + a linter.** The
  viewer's grouping-detection policy (the preferred-name list, min/max group counts) is canonical in
  `conformance/viewer_policy.json`; each surface carries a copy (`_PREFERRED_GROUPINGS`,
  `.VIEWER_PREFERRED_GROUPINGS`, `js/core/policy.ts`) that `conformance/policy_linter.py` asserts equal.
- **Contract shapes → the format spec.** Field names, roles, spans/orientation, encodings, and states
  are fixed by `docs/format.md` (e.g. the `viewer@0.1` table) and enforced by `validate()`.

Do **not** re-implement a shared computation per language "because it's small". The bugs this contract
exists to prevent were all small: a CSR check that threw instead of normalizing, a reorder stubbed to
identity on one surface, a category list sorted on two surfaces and first-seen on the third, an int cast
that truncated a float measure on one surface only.

## Checklist for a shared-feature change

1. **Land the logic once** — in the core (bound everywhere) or as a single-sourced policy + linter entry.
2. **Normalize inputs at the boundary**, don't reject them. If a surface can't consume an input shape a
   sibling accepts (e.g. a CSR matrix), convert it via the shared kernel — never throw where another
   surface would succeed.
3. **Add/extend a cross-surface conformance leg** over the feature's full input space. The viewer legs
   live in `conformance/viewer*.sh`:
   - `viewer.sh` — Python vs native-R prep on a synthetic base (all fields, incl. `counts_cellmajor` + `_order`).
   - `viewer_js.sh` — lstar's own `extend-viewer.ts` vs Python on **CSC + CSR + competing-groupings** inputs.
   - `viewer_corpus.sh` — convert curated corpus datasets (real local / synthetic-faithful in CI) and
     cross-check every surface. This is where realistic structure (non-alphabetical categoricals,
     several competing groupings, real embeddings, native CSR) exercises what synthetic fixtures miss.
     Includes a synthetic `dense_primary` dataset so the dense measure-read path is covered on every
     surface (the corpus is otherwise all-sparse).
   - `encoding_invariance.test.ts` — the live-viewer compute (`LstarView.colStats`/`subsampleDE`) gives
     identical results for a measure written `{dense, csc, csr}`. Metamorphic: encoding must not change
     the answer, so no reference is needed.
   - `policy_linter.py` — the single-sourced policy constants match across surfaces.

   **Assert the artifacts, not just "no error".** Where a failure degrades silently — a caught exception
   that falls back rather than raising — assert the outputs *exist* (`od_score`/`stats_*`/`markers_*`/
   `counts_cellmajor*`), or a skipped step reads as a pass.
4. **Cover the input-space axes**, not just the happy path. **Encoding is a first-class axis:** any path
   that reads a field must be exercised across every encoding it can receive, on every surface that
   implements it — never only the encoding the corpus happens to carry. The measure encoding is
   `{csc, csr, dense}`: a dense measure (e.g. an SCE `logcounts` assay or a scaled/dense AnnData `X`) is
   stored at `/values`, a sparse one at `/data`+`/indices`+`/indptr`, so a measure-read path must accept
   all three on every surface. Also: label encoding `{categorical, utf8}` incl. non-alphabetical order;
   groupings `{single, several competing}`; embedding `{present, absent}`; basis `{raw, lognorm}`;
   coverage `{full, partial}`; nullable `{mask, none}`.
5. **Wire the leg into CI** — `.github/workflows/ci.yml` (r-cross-format for Py+R, js-wasm for Py+JS) and
   `conformance/run.sh` (the full local run, incl. the real corpus). A cross-surface check that isn't in
   CI, or that skips silently, does not count as coverage.
6. **Update the docs** when the change affects compatibility or the public contract (`docs/format.md`,
   the per-format guides), per the doc discipline in `.claude/CLAUDE.md`.

## Why the harness is shaped this way

Divergences hide for two structural reasons: the fixtures don't cover the input space (a CSC-only
synthetic never exercises CSR; an all-sparse corpus never exercises a dense measure; single-grouping never
exercises detection), and the cross-surface check ignores what actually diverges (the fields
`counts_cellmajor` / `_order`, or — for a path that fails silently — whether the output exists at all). So
the contract is: **the fixture must span the input space (field encoding included), and the check must
compare every field and assert it exists.** When a new divergence class appears, the fix is not just the
code — it's the fixture axis + the assertion that would have caught it.

## Cross-surface scope & known asymmetries

Not everything is uniform across surfaces, and some asymmetries are deliberate. This is the authoritative
list (from the parity audit) so a real gap is never confused with an intentional one.

**By design (surface-native capabilities):**
- **Format IO profiles.** AnnData / MuData live in Python; Seurat / SingleCellExperiment / Conos / pagoda2
  live in R. Each ecosystem's object model is native to one language; all funnel into the *same* L* core
  representation (`cross_format.sh` proves the chain), so this is not a core-representation divergence.
- **Interactive query API.** `LstarView` / Crossfilter / live `colStats` / `scalarToRGBA` are JS-only; the
  Python/R `view()` delegate to the pagoda3 viewer. The precompute-once half (`extend_for_viewer`) is
  uniform across Python/R/JS and fully conformance-covered.
- **`validate()` is the canonical structural validator (Python).** R/JS/C++ have no separate validator;
  instead every surface's *output* store is Python-validated in CI (the viewer / cross-format legs do this).
- **Content-based `state` inference is AnnData-scoped.** The AnnData reader (native *and* direct, now in
  agreement) infers raw/lognorm/scaled from content; the R Seurat/SCE profiles infer `state` from the slot
  name (reliable for canonical `counts`/`data`/`logcounts`/`scale.data` slots).
- **Single-file `.lstar.zarr.zip` is uniform (read + write) on all four surfaces** — not a tier. The
  STORED constraint (see docs/format.md §Packaging) makes every reader a seek + copy with no zip-layer
  codec, so the *implementations* differ by each surface's **role, not capability**: JS reads a *hosted*
  zip by issuing an HTTP `Range` into the archive (`ZipStore` / `httpZipSource` — the reason STORED
  matters); C++/R read a *local* zip by extracting its STORED entries to a temp dir (a copy, no
  decompression) and running the normal reader; Python uses zarr's `ZipStore`. Writers all pack STORED
  (ZIP64-aware) and reject a DEFLATE `.zip` on read. Enforced by `zip.sh` / `zip_r.sh` / `zip_js.sh`.

**Scoped / follow-up (a real gap, intentionally deferred — do not silently widen):**
- **Streaming reductions directly from a `.zip`** are not wired: the bounded-memory `stream_col_stats` /
  block readers run over a *directory* store, and C++/R read a `.zip` by extracting it first (so a huge
  single-file store is materialized to a temp dir before streaming). JS's `ZipStore` *is* range-capable
  but isn't threaded into the streaming reducers yet. Fine for the current use (stream over a directory;
  load/host a single file); revisit if bounded-memory streaming over a huge single-file store is needed.
- **DE analysis API (`pseudobulk`, `collection_pseudobulk`, `de_bundle`, `de_factors`) is Python-only**, but
  its per-(group,gene) reduction now routes through the shared `col_sum_by_group` core kernel (was a numpy
  per-group loop that duplicated it). The interactive selection-DE in the pagoda3 viewer is a *separate* JS
  implementation that should likewise call the WASM `subsampleDeRank`/`colSumByGroup` kernels rather than
  reimplement the reduction. Port the DE *API* to R if it must be cross-surface; the *kernel* is now shared.
- **Depth-normalized streaming reducers (`stream_col_stats` depth/population args, streamed pseudobulk)**
  are R/C++-only (R is the pagoda2 host); Python's `stream_col_stats` is lognorm-only.
- **`subsample_de_rank`** kernel is bound on Python/R/WASM but the live selection-DE is the JS viewer's own
  (`LstarView.subsampleDE`); the kernel is available for callers who want cross-surface-identical ranking.
- **`uncertainty`** round-trips through Python and the C++ core, but is not threaded through the R bridge
  (rare, AnnData-specific).
- **Compression codec.** gzip is the portable codec (all surfaces read it); C++/R also read zlib; JS reads
  gzip only. Python can write any numcodecs codec, but a non-gzip store is not portable — prefer gzip.

**Metadata that is non-normative** (excluded from the "byte-identical store" contract): a viewer field's
`provenance` stamp carries surface-specific detail (`method`/`curve`/`grid`), so `cmd_equiv` compares the
*data* fields (stats/markers/od/`counts_cellmajor[_order]`) and not provenance. The `od_score` tolerance is
looser than the rest because Python derives its per-gene variance naively while the core uses a stable
centered form (they agree through the F-test; unify the od variance to tighten it).
