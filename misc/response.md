# Response to the pagoda2.1 suggestions (lstar side)

Re: `misc/suggestions.md` (S1–S5). Written 2026-06-13 after implementing. Bottom line: **S1, S3, S4a,
S5, and S2 are done and tested**; **S4b is a deliberate "don't merge yet — generalize the signature
instead."** Two of your "confirm whether it already works" items flipped on inspection, noted below.

Everything below round-trips across Python, C++, and R (the path that matters for pagoda2), and is
guarded by CI conformance unless stated.

---

## S1 — recipe params as provenance — **DONE, but it was broken** (not "already works")

You asked us to confirm provenance survives a round-trip. It did **not**: the R cpp11 binding silently
dropped field `provenance` — Python preserved it, but `pagoda2 → lstar(R) → pagoda2` returned `NULL`. So
the recipe would have been lost on exactly your path. **Fixed**: provenance is now carried across the R
boundary as an **opaque JSON string** (the same robust trick the `aux` passthrough uses), preserving
arbitrary nesting. Guarded by `conformance/provenance.sh` (Py→store→R→store→Py, recipe intact).

**Use it like this** (R): `ds$fields[["counts"]]$provenance` is a JSON string;
`jsonlite::fromJSON(...)` to read, `jsonlite::toJSON(recipe, auto_unbox=TRUE)` to set. Python sees a
normal dict.

**One correction to your ask — params vs. vectors.** Put the *small* recipe scalars in provenance
(`model`, `depthScale`, `log_base`, `winsor_caps`). Do **not** put the *large precomputed per-axis
vectors* there (the CLR per-row divisor, length n_cells; the per-column IDF, length n_genes): provenance
is a JSON metadata dict read on every store open — serializing a 40k-element array into it is a
bloat-and-correctness mistake. Store those as **arity-1 fields** — a `measure` over `cells` / `genes`
with `subtype="recipe_scalar"`. They then stream, chunk, and compress like any other field, and the
kernel reads them as the `r[]` / `c[]` vectors (see S4b). Net: recipe = a few scalars in provenance + (if
you want byte-exact reproduction) two arity-1 fields.

**S1′ (typed recipe vocabulary).** Not done — opaque provenance round-trips faithfully today, which is
enough to be lossless. Typing `plain`/`clr`/`tfidf` only buys cross-kernel reuse (the WASM viewer
applying CLR without pagoda2 in the loop); revisit alongside S4b if/when that's wanted.

## S2 — collection-level streaming grouped reducer — **DONE (Python), this is the form**

Implemented `lstar.collection_pseudobulk(ds, factor, field="counts", lognorm=False)` — the collection
generalization of `pseudobulk`. **Form (build your Conos path against this):**

- `factor` — a categorical `label` over the **union** `cells` axis (your joint clustering). Its
  categories are the output groups (the induced rule-2 factor axis).
- `field` — the per-sample measure root; it finds every `<field>.<s>` measure over `(cells.<s>, genes…)`.
- It **walks the per-sample measures one at a time** (the joint matrix is never materialized — bounded
  by the largest sample, not the union), maps each sample's cells to their union positions to get each
  row's group, and **accumulates in float64**. Genes are **label-aligned** across samples (a gene absent
  from a sample contributes 0).
- Output: `pb.<factor>.{mean,frac}` measures over `(factor, genes)`, with `provenance.collection=True`.

Verified streamed == a materialized reference for both shared and *differing* per-sample gene sets
(`python/tests/test_collection_reduce.py`). This is the **primitive extension**, not a "markers bundle"
— exactly as you framed it.

*Scope note:* v1 materializes each sample block (bounded memory relative to the union, the win you
wanted). A truly *block-streaming* version (each sample itself read in column blocks from a backed
source) is the natural next step when you point this at on-disk Conos stores — it slots onto the same
signature; say the word and we'll add it on the C++ side behind `stream_col_stats`/grouped-sum.

## S3 — partial `index` over a derived union axis — **DONE (already worked; now guarded)**

Confirmed: a partial-coverage measure whose `index_axis` is a **derived** `kind="collection"` union axis
validates and round-trips — `index` is keyed by axis *name*, with no observed-vs-derived assumption. Your
facet membership mask maps 1:1 onto it. Locked in cross-language by `conformance/partial.sh`
(Py→R→Py over a derived union axis). Nothing for you to work around.

## S4a — determinism contract — **DONE (taken, as you expected)**

Stated in `docs/principles.md §5` and **tested exactly** (`== 0`, not within tolerance): the C++ kernel's
`mean`/`var`/`nnz` are **bit-identical across 1/2/4/8 threads** (float64 accumulation, column-parallel,
no cross-thread reduction) — `conformance/stream_reduce.sh` + `python/tests/test_determinism.py`. You can
rely on a facet's summaries being identical whichever library computes them — the precondition for S4b.

## S4b — one shared core — **decision: don't merge; generalize the signature**

Your scalars+enum boundary is right (and it's the same insight from the multimodal review — CLR's divisor
behaves like `depth[row]`, IDF like a per-column cap; pure per-entry functions of precomputed scalars,
which is simultaneously what makes the kernel thread-invariant *and* keeps "CLR" out of the core). But a
hard codebase **merge** couples pagoda2's hot path to libstar's release cycle for low immediate gain
(`misc2.cpp` works today). 

**What we'll do instead:** generalize `libstar`'s reducer toward the `(raw block, r[] per-row scalar,
c[] per-col scalar, transform enum t)` signature — which the fused `depth+log1p` reducer already half-is
(it takes optional `depth` + `depthScale`). Exposed as a bounded core primitive, pagoda2 can call it
*optionally* (header-only) and feed it the `r[]`/`c[]` it computed (your CLR divisor / IDF, stored as the
S1 `recipe_scalar` fields), with the core never learning "CLR." **Keep the repos separate.** Revisit a
true merge only on a concrete benchmark, and per your decision rule, the first kernel that can't fit
scalars+enum is the signal to keep it in `misc2.cpp`. We haven't built the generalized signature yet —
it's the obvious next lstar kernel task and is gated on a real consumer (you), so let's sequence it with
your Conos work.

**On the "DE/markers crossed the line" caution — we half-agree.** lstar's DE is mostly *interchange
typing* (it represents scanpy's already-computed `rank_genes_groups` as a `(factor,genes)` bundle — it
does not run a test); `pseudobulk`/grouped-sum is a *general primitive*. So the core hasn't really
crossed into method-specific *compute*. The line we'll hold: **never implement specific statistical tests
in the core** — keep the bundle typing (interchange) and the general reductions; anything test-specific
stays in the user package.

## S5 — facet-set provenance on joint products — **DONE (gated on S1, now unblocked)**

Convention: a joint product's embedding records the contributing feature axes in
`provenance["input_axes"]` (a list of feature-axis names). The MuData profile **auto-populates** it by
inference — the feature axes carrying a `loading` over the embedding's factor axis (MOFA: `["genes",
"proteins"]`) — and it round-trips (now that S1 makes provenance survive R). **Your side:** when you
write `reductions[["WNN"]]` with `facets=c("RNA","ADT")` to lstar, set the same key
(`provenance$input_axes <- c("genes","proteins")`); it will round-trip losslessly.

---

## Sequencing for the Conos work you're about to start

1. **Use today:** `collection_pseudobulk` (S2 form above), partial `index` over the union axis (S3),
   provenance for the recipe (S1 — small params in provenance, large vectors as `recipe_scalar` fields),
   `input_axes` on joint reductions (S5), and the determinism contract (S4a) for reproducible summaries.
2. **Coordinate before building:** the **block-streaming** collection reducer (S2 scope note) and the
   **generalized `(r[], c[], enum)` kernel** (S4b) — both are lstar-side work we'll do *with* you, gated
   on your concrete Conos store + a benchmark, so we don't build them speculatively.
3. **Not doing:** merging the codebases (S4b), computing normalization in lstar (the retracted virtual
   `recipe` field), or adding statistical tests to the core.
