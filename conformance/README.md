# lstar conformance suite

Each `*.sh` here is a self-contained cross-language gate: it builds a store in one binding and proves
another binding reads it back with values intact. `run.sh` is the master runner; `sweep/` is the local
real-corpus loop (`sweep/RETEST.md`). Most scripts are synthetic and run in CI; the real-data ones skip
cleanly when their datasets/packages are absent.

## Origin coverage — the standing rule

**A symmetric interchange format must be tested from every *producer*, not just every *consumer*.** L★
lets Python, R, C++, and JS all *write* a store and all *read* one. A test of the shape

```
Py → store → R → store → Py        # "round-trips through R"
```

*looks* symmetric but the payload **always originates in Python**. It exercises the Python writer, the R
reader, the Python reader — and the R writer **only on data R just read back from a Python store** (R as
*rewriter*, never as *author*). Any writer code path reachable only from data a language *constructs
natively* is invisible to it.

That blind spot shipped a real bug. Field `provenance` round-tripped `Py → R → Py` green for months
because R read it as a JSON *string* and wrote it back as a string — internally consistent. The latent
`is.character()`-only writer branch was never hit, because R-read provenance was always a string. It broke
only when provenance **originated in R as a named list** (pagoda2's case), which that branch silently
dropped. The old assertion was `is.character(p)` — it asserted the *wrong* native type and passed anyway.

Two consequences, both now rules:

1. **Producer-origin coverage.** For every rich feature (categorical, mask, partial coverage, aux,
   arity-n, provenance, induced axes, collections), there must be a leg where **each binding authors the
   feature in its native form from scratch** — not just rewrites a store another language wrote. The
   R-authored legs use `structure(list(...), class="lstar_dataset")` (see `categorical.sh` for the
   template); the JS-authored leg is `js/test/writer_make.ts` (emits every encoding).

2. **Cross-reader agreement over same-language round-trip.** When language X writes a store, assert that a
   **different** language Y reads *equal* values — never that X reads its own store back. A
   same-language round-trip can be self-consistently wrong (read + write bugs in one language cancel). A
   second reader breaks that symmetry. Assert the **contracted native type + exact nested content**
   (named list in R, dict in Python, the actual values), not mere survival/presence.

The C++ core is the engine under the R binding, so **every R-authored leg also exercises the C++ writer**,
and every "C++/R reads the Python store" leg exercises the C++ reader — C++ as producer/consumer is
covered transitively (no standalone C++-authoring harness needed; the binding is a thin shell over the
same `lstar::write`).

### Anti-patterns (smell tests)

- An assertion that only checks a value *exists* / has the right *length* / is non-empty — tighten it to
  the native type + exact content. `is.character(p)` passing is the canonical smell.
- A feature exercised only from the language that first implemented it. Ask: *who authored this payload?*
  If the answer is always "Python", the other writers are untested for it.
- Trusting `X → store → X`. Prefer `X → store → {Y, Z}` (fan-out to other readers) **and** `{Y,Z} → X`.

## Coverage matrix (rich-feature cross-language tests)

| test | feature | Py-authored | R-authored | JS-authored |
|---|---|:---:|:---:|:---:|
| `categorical.sh` | factor codes / categories / ordered / −1 | ✓ | ✓ | ✓ (writer) |
| `nullable.sh` | uint8 validity mask | ✓ | ✓ | ✓ (writer) |
| `partial.sh` | partial coverage over a derived axis | ✓ | ✓ | ✓ (writer) |
| `arity3.sh` | n-D dense tensor | ✓ | ✓ | — |
| `aux.sh` | passthrough subtree (uns / @misc) | ✓ | ✓ (mutates a leaf) | ✓ (writer) |
| `provenance.sh` | field provenance (recipe / facet) | ✓ | ✓ (case 6 originates in R) | — |
| `induce.sh` | induced factor-axis link | ✓ | ✓ | ✓ (writer) |
| `collection.sh` | collection-of-samples (not a tensor) | reads | ✓ (R origin) | — |
| `collection_true.sh` | `collection_from` heterogeneity (divergent + disjoint/cross-species genes) | ✓ | ✓ | — |
| `conos.sh` | Conos graph-only collection -> Seurat v5 split + AnnData (no corrected matrix) | reads | ✓ (real Conos + R synth) | — |
| `de.sh` | DE bundle (`rank_genes_groups`) | ✓ | reads | — | one-directional by design |
| `js.sh` (writer) | every encoding, chunked + gzip | reads | — | ✓ |

JS authors via `writer_make.ts` → `writer_crossread.py` (Python cross-reads, value-equal). The `de.sh`
bundle is intentionally Python-origin (it *is* a scanpy result); the R side only consumes it.

## Adding a new rich-feature test — checklist

1. Python authors the feature from scratch → a **different** reader (R via the C++ core, and/or Python)
   asserts native type + exact content.
2. **R authors the same feature from scratch** (`structure(list(...), class="lstar_dataset")`) → Python
   cross-reads it value-equal. (This is the leg that catches writer-drop bugs.)
3. If the JS writer can emit it, add it to `js/test/writer_make.ts` + assertions to `writer_crossread.py`.
4. Tag the test header `# Origin coverage: Py-authored ✓ | R-authored ✓ …` and add a row above.

## Running locally

```bash
bash conformance/run.sh                 # build core + install R pkg + Python/cross-impl/cross-format
bash conformance/<name>.sh              # a single gate (needs .Rlib + python/src on PYTHONPATH)
LSTAR_EMCC_PYTHON=/path/to/py3.10 bash conformance/js.sh   # JS/WASM (needs emsdk + js/node_modules)
```
