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
   - `policy_linter.py` — the single-sourced policy constants match across surfaces.
4. **Cover the input-space axes**, not just the happy path: counts encoding `{csc, csr}`; label encoding
   `{categorical, utf8}` incl. non-alphabetical order; groupings `{single, several competing}`; embedding
   `{present, absent}`; basis `{raw, lognorm}`.
5. **Wire the leg into CI** — `.github/workflows/ci.yml` (r-cross-format for Py+R, js-wasm for Py+JS) and
   `conformance/run.sh` (the full local run, incl. the real corpus). A cross-surface check that isn't in
   CI, or that skips silently, does not count as coverage.
6. **Update the docs** when the change affects compatibility or the public contract (`docs/format.md`,
   the per-format guides), per the doc discipline in `.claude/CLAUDE.md`.

## Why the harness is shaped this way

Divergences hid for two structural reasons, both now closed: the fixtures didn't cover the input space
(a CSC-only synthetic never exercised CSR; single-grouping never exercised detection), and the
cross-surface equality check ignored the fields that diverged (`counts_cellmajor` / `_order`). So the
contract is: **the fixture must span the input space, and the equality check must compare every field.**
When a new divergence class appears, the fix is not just the code — it's the fixture axis + the field
comparison that would have caught it.
