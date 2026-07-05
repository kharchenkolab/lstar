# lstar (L★)

## What lstar is

L★ is a general model for single-cell/spatial omics data, built from **axes** (the entities you
index by — cells, genes, samples, clusters) and **fields** (typed data over them — counts,
embeddings, loadings, graphs, labels, designs). Because everything is just axes and fields, one
small model spans cases a fixed `cells × genes` container strains on: multi-sample (even
cross-species) integrations kept as a *collection* of heterogeneous samples, CITE-seq with a second
feature axis, case-control cohorts carrying a statistical *design*.

Its most immediate use is **lossless, explicit conversion** between the formats people already use —
AnnData, Seurat, SingleCellExperiment, pagoda/conos — routing each through one shared-vocabulary L★
store, preserving the meaning of every piece and *reporting* (never silently dropping) whatever a
target can't hold. It reads/writes a portable Zarr-based format and streams heavy operations in
bounded memory, so multi-gigabyte datasets convert on a laptop.

## Components (one model, four surfaces)

| Surface | Path | Role |
|---|---|---|
| C++ core | `core/` | `libstar` — header-only: the model, chunked+gzip Zarr IO, the fast kernels |
| Python | `python/` | the `lstar` package on zarr-python, binding the optional compiled C++ accelerator |
| R | `R/` | the `lstar` package; the format profiles (Seurat, SCE, Conos) live here |
| Browser/Node | `js/` | a TypeScript reader (zarrita) + the kernels compiled to WebAssembly, for viewers |

`conformance/` holds the shared round-trip / cross-format / cross-language test suite.

## Scope vs pagoda3

lstar owns the shared substrate — the model, the on-disk format, all four surfaces, and every reusable
building block: Zarr IO, the fast kernels, and the JS/WASM **reader, store backends, and
viewer-extension prep**. pagoda3 is a *downstream app/viewer* that consumes these; anything reusable by
another consumer belongs in lstar (on whatever surfaces need it), not pagoda3.

## Core aim: cross-language consistency

**We aim to implement consistent logic and feature sets across C++, Python, R, and JS/WASM — ideally
by minimizing code duplication** (share the C++ core; bind it from Python and R; compile it to WASM
for the browser) rather than reimplementing the same logic four times and letting the versions drift.
When adding or changing a feature, treat all four surfaces as one system: keep the behavior, the API
shape, and the results identical across languages, and extend the `conformance/` suite so any drift
is caught. **See `docs/parity.md`** for the concrete contract and the checklist to follow (one recipe /
N thin bindings; single-sourced policy + linter; a cross-surface conformance leg spanning the input
space, wired into CI).

## Testing
Roundtrip and other conversions with a corpus of tests datasets provides a powerful test. See python/tests/corpus.py and corpus data on mendel.
CI in the repo is kept light, relying on synthetic data (don't check in real datasets into the repo - if you need to add a tests case, figure out how to generate appropriate synthetic)
Periodically run a cross-surface **parity/duplication audit** (method + prior findings in `docs/parity.md` and `misc/parity_audit.md`): fan out read-only agents over the C++/Python/R/JS surfaces to hunt feature-parity gaps, per-language reimplementations of shared logic, and behavioral/policy divergences. The conformance legs + `conformance/policy_linter.py` catch known drift; new *classes* of divergence need a fresh agent sweep (not currently automatable).
Test a reader/store backend at the layer a real consumer uses — **decoded field *values* read back through the reader** (`openLstar` / `read`), compared across *every* backend — not merely that it opens or that raw `get()` bytes match. "Opens" ≠ "reads data correctly": a backend that returns nothing for every chunk still opens (metadata is separate) and yields non-empty axes/fields while silently zeroing all data. See `js/test/store_backends.test.ts` (and `misc/zipstore_bug_postmortem.md`).

## Docs
`docs/` holds **coherent, published documents** that describe the software *as it currently is* — its functionality and current state, not how it got there. Each has a specific audience and goal — most are **end-user-facing** (bioinformaticians who *use* the package: the format spec, the model, recipes); a few (e.g. `docs/parity.md`) are a **contributor** contract. Before touching one: identify its audience and goal, read enough of it to know its structure and narrative, then add only what serves *that* goal for *that* audience — an edit must fit the document, not sit beside it as a patch. Write for the reader (a working bioinformatician for user docs), not for yourself.

**These are documents, not journals or changelogs.** Bug post-mortems, testing-methodology lessons, audit narratives, progress reports, "how we found/fixed X", and other debug/historical notes do **not** go in `docs/` — put a terse, durable principle in this file (CLAUDE.md) and any detailed write-up in `misc/` (gitignored local dev notes) instead. On a real change, first decide whether it's even in-scope for a given document (minor bugfixes usually aren't); a change to a supported format/version or a CLI/API change **is** in-scope and must be reflected — and must also update the R `NEWS.md` (regenerate `man/*.Rd` if roxygen changed) and any downstream **skills/recipes or handoff notes that document lstar's usage** (ABA recipes, the pagoda3 handoff), which drift silently otherwise.
