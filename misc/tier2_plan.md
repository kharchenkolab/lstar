# Tier 2 (per-group results) — implementation plan + risk projection

*Draft for review, companion to `format_coverage.md` and `induction_design.md`. Tier 2 = analysis
results computed **per group** (differential expression, pseudobulk, PAGA…) represented as ordinary
fields over an induced **factor axis**, using the `categorical` (codes-into-axis) encoding. The point of
this doc is the **risk projection** the model asked for — does the treatment cause descriptor
explosions, costly remappings, or translation ambiguities?*

## 0. What it concretely builds

1. **`categorical` encoding** (foundation, shared with Tier-1 dtype fidelity and the spatial molecule
   dictionary): a label stored as integer codes + a categories axis + `ordered` + `-1` missing.
2. **`induce()` + `factor` role + canonical identity** (per `induction_design.md`): a categorical field
   induces its factor axis at read time, data-driven.
3. **Per-group result fields** over `(factor, genes)` — DE/markers + pseudobulk sufficient statistics.
4. *(secondary)* PAGA (relations over the factor axis), dendrograms, `@commands`/params provenance.

### The precedent already in the tree (and what it reveals)

`R/R/profile_pagoda2.R:67-86` already emits, **per clustering**: one `groups_<g>` axis + **5 fields**
(`stats_*_{sum,sumsq,nexpr}`, `markers_*_{lfc,padj}`). Two things it reveals up front:
- It's **dense** (K×G matrices), not sparse — computed markers cover all genes.
- It is **internally inconsistent in orientation**: stats are `(groups, genes)` but markers are
  `(genes, groups)` (line 81 vs 84). That inconsistency, *inside one profile*, is concrete proof that
  "(group, gene) vs (gene, group)" must be **canonicalized** or it will leak ambiguity everywhere.
  It also tags the axis `role="feature"` — the mislabel `induction_design.md` fixes with `factor`.

## 1. Implementation phases

- **P1 — `categorical` encoding** in `model.py`/`zarr_io.py` + the C++/R/JS readers (codes, categories,
  `ordered`, `-1`). This is the foundation; it also fixes the Tier-1 categorical round-trip bug.
- **P2 — `induce()` + `factor` role + canonical identity + `validate()` check**; retrofit the profiles
  to induce categorical columns and set `induced_by` (data-driven, per §3 of the induction note).
- **P3 — DE + pseudobulk read/write** for AnnData (`uns['rank_genes_groups']`), generalize the pagoda2
  pseudobulk path; capture Seurat marker tables where present.
- **P4 — secondary** (PAGA/dendrogram/provenance), only if wanted.

## 2. Risk projection (the actual question)

### R1 — Descriptor explosion: **real but bounded; mitigable.**
Sources, by severity:
- **Pairwise / per-condition DE is O(K²)** (or O(K·conditions)). A 20-cluster all-pairs DE is 190
  comparisons. *But the explosion already exists in the source* — the user put 190 results in `uns`; we
  only reflect them. The risk is that *our* mapping makes it **worse** (5 loose top-level fields per
  comparison → ~950 fields). **Mitigation:** (a) model one DE result as a **bundle** — companion arrays
  (lfc + padj + score) sharing one `(factor, genes)` field group, not 4-5 independent top-level fields;
  (b) **don't auto-type pairwise/conditional DE** — keep it lossless *passthrough* (the long tail), type
  only the common **one-vs-rest** case. That caps the typed surface at ~1 bundle per clustering.
- **Multiple clusterings/resolutions** (leiden@0.5/1.0/2.0, cell_type, …) — N clusterings × 1 DE bundle +
  pseudobulk. Linear, ~5-8 fields each. Manageable; each is a real result the user computed.
- **Derivable `names`.** scanpy's `names` (per-group ranked gene list) is *recoverable by argsort* of the
  scores. **Don't store it** — reconstruct the ranking on write. Saves a field per group and a duplicate.
Verdict: real risk, fully contained by **bundle + don't-type-pairwise + drop names**. Typed Tier-2
surface ≈ (1 DE bundle + 1 pseudobulk bundle + 1 factor axis) per clustering.

### R2 — Factor-axis proliferation: **non-issue (a manifest entry, not storage).**
Every categorical column (sample, condition, batch, donor, sex, phase, 10 annotations…) induces a factor
axis — 10-20 of them in a real study, most spanned only by their own inducing label. *But the factor
axis **is** the `categories` array of the categorical encoding, which we store regardless.* Naming it an
axis adds a manifest entry + a role, **zero extra bytes**. So this is clutter, not cost. (Optional: only
*name* a factor axis when a non-self field spans it — but since it's free, eager is fine.) Verdict: not a
real problem; worth stating so it isn't mistaken for one.

### R3 — Time-consuming remappings: **not a performance risk.**
On **read** from AnnData, scores/lfc/pvals are ordered by each group's *ranking* (aligned to `names`),
not by gene index — so building a gene-keyed `(factor, genes)` field requires a **scatter** keyed by
`names[g][i] → gene_index`. That is O(nnz) one-time: ~K·G ≈ 20·20k = 4·10⁵ hashmap lookups — milliseconds.
On **write**, regenerate `names` by argsort per group — same cost. Verdict: cheap; the concern here is
*correctness* (the scatter must use the right gene axis, R4b), not speed. The only wrinkle is **tie
order** in argsort not exactly reproducing scanpy's `names` for equal scores — cosmetic.

### R4 — Translation ambiguities: **the real risk cluster.** Five distinct ones:
- **(a) DE comparison semantics.** `uns['rank_genes_groups'].params['reference']` is `'rest'`
  (one-vs-rest) **or** a specific group (one-vs-one). A `(factor, genes)` field implicitly means
  one-vs-rest; a reference-group DE is a *different* quantity. **Must record `reference` on the DE bundle**
  (and `method`, `groupby`). Pairwise has no `(factor, genes)` home at all → passthrough (R1).
- **(b) Which gene axis?** DE may be computed on `use_raw=True` (the `.raw` gene set) or a HVG subset, not
  the main `genes` axis. The `(factor, genes)` field must span **the gene axis the DE was actually tested
  on** — which can be `genes_raw` or an HVG axis. Get this wrong and the scatter (R3) misaligns silently.
  **Must resolve the DE's gene axis from `params['use_raw']`/the tested gene set**, not assume `genes`.
- **(c) Storage asymmetry across formats.** AnnData stores DE in `uns`; **Seurat and SCE do not** (markers
  are a detached data.frame / a metadata list). So DE **round-trips only through AnnData**. AnnData→L*→
  Seurat puts DE in `@misc` (opaque) or the loss manifest. This is **inherent, not fixable** — "support
  DE" is an *asymmetric* capability, and must be reported via `dropped`, not silently. Also, a Seurat
  marker df's `cluster` column links to *which* identity class is **ambiguous** when detached from the
  object — so Seurat→L* DE capture is best-effort.
- **(d) Method-dependent stat presence.** wilcoxon → scores+pvals; logreg → scores only; t-test →
  different. **Not every companion (lfc/score/pval/padj) exists.** The bundle's fields must be
  **optional**, never assumed-present. Record `method`.
- **(e) Density mismatch.** scanpy DE is dense (all genes ranked); Seurat markers are *filtered*
  (significant only) — sparse, with a filter threshold that is itself a result descriptor absent from
  AnnData. The sparse-or-dense `(factor, genes)` measure holds both, but an AnnData→L*→Seurat→L* trip can
  **change density** (the Seurat filter drops genes). Flag: "markers" is not one well-defined object
  across ecosystems.
Verdict: this is where Tier 2 actually bites. None are show-stoppers, but they force the DE representation
to carry **descriptors (`reference`, `method`, `groupby`, gene-axis) and be optional/sparse-tolerant**,
and they make Seurat/SCE DE a capture-not-bridge.

### R5 — Internal-shape canonicalization: **must decide, and the precedent is already wrong.**
Two candidate canonical shapes for DE/markers: **(A)** sparse-or-dense measures over `(factor, genes)`
(matrix-aligned, what the viewer reduces), vs **(B)** a long-form *tidy* table over a derived `de_rows`
axis with `group`/`gene`/`lfc`/`padj` columns (how DE is consumed/sorted). The pagoda2 precedent uses (A)
but **inconsistently orients it** (groups×genes for stats, genes×groups for markers). **Decision needed:**
pick (A) as canonical with a fixed orientation (recommend `(factor, genes)` to match the measure
convention) + carry the descriptors from R4, and offer (B) as a derived *view* (a `get_markers()` that
flattens to a tidy frame). Settle the orientation **once** or it leaks everywhere (the precedent proves it).

### R6 — New encoding across four languages: **real implementation cost, not a design risk.**
`categorical` (codes + categories + ordered + `-1`) must read identically in Python/R/C++/JS, like
`utf8` did. The C++/JS readers are the long pole. Bounded, known work; flagged so it's budgeted.

### R7 — Axis/field name sharing: **confirm, minor.**
Field `leiden` (label over cells) and axis `leiden` (its categories) share a name (separate namespaces:
`axes/` vs `fields/` zarr groups, and separate `ds.axes`/`ds.fields` dicts — so legal). Caution: any
tooling assuming a globally-unique name must be checked; `validate()` and cross-refs must treat a span
entry as an *axis* reference unambiguously.

## 2b. Alignment — how values map to cells/genes (vs pagoda2's named vectors)

lstar keeps labels **once, on the axis**; a field over `(cells)` is *positional* against the `cells`
axis labels (value at position *i* ↔ `cells.labels[i]`). pagoda2/conos instead carry `names` on every
vector and realize the cell-mapping per-vector (which is what lets a factor be mismatched/reordered/
partial). lstar **resolves that name-alignment once at ingest** (align the factor's names → `cells` axis
positions) and stores a positional field — same name-correct mapping, but one authoritative labelling,
so fields cannot disagree on order. A **mismatched/partial factor** (covers a cell subset) becomes a
full-length `categorical` with `-1` (missing) for uncovered cells — the *same* sentinel as NA-categorical
fidelity, so "doesn't cover all cells" and "has NAs" are one mechanism (exactly an AnnData `obs` column
with NaN gaps). The only case `-1`-padding handles poorly — a factor over a *tiny* fraction of a huge
atlas — is the spec'd **partial `coverage`** (a field over a subset axis + a `membership` relation back
to `cells`), which is **spec'd-but-not-implemented** and an add-on *if* faithful sparse-subset factors
are needed rather than padded ones.

## 3. Decisions to settle before coding

1. **Canonical DE shape + orientation:** `(factor, genes)` sparse-or-dense measures (A) as canonical,
   fixed orientation, + a tidy derived view (B). Fix the orientation the pagoda2 precedent got
   inconsistent.
2. **DE bundle carries descriptors** (`reference`, `method`, `groupby`) and is **optional per stat**.
3. **DE spans the gene axis it was tested on** (`genes`/`genes_raw`/HVG), resolved from `params`.
4. **Type only one-vs-rest; passthrough pairwise/conditional** (avoid the O(K²) explosion).
5. **Don't store the DE *ranking* array** (`rank_genes_groups['names']` — the per-group gene order;
   reconstruct by argsort of the scores). *(This is about the DE result's ordering only — not axis labels
   or the names-based value→cell/gene alignment, which lstar keeps; see "Alignment" below.)*
6. **Seurat/SCE DE = capture-to-`@misc`/loss-manifest**, not a round-trip — documented asymmetry.

## 4. Recommended scope (the problem-free core)

**Build:** `categorical` encoding → `induce()`/`factor` → **one-vs-rest DE as a `(factor, genes)` bundle
with descriptors** + **pseudobulk** (generalize the existing libstar `col_sum_by_group` path) → tidy
derived view. AnnData read/write; Seurat/SCE capture + loss-record.
**Defer / passthrough:** pairwise & conditional DE, PAGA/dendrogram, provenance.
This core hits the highest-value gap (AnnData DE, currently fully dropped), reuses the kernel we just
optimized, feeds the pagoda2.1 viewer directly, and **structurally avoids every explosion** (R1/R2) and
**confines the ambiguities** (R4) to recorded descriptors + an honest capture-not-bridge story for
Seurat/SCE.

**Bottom line:** the approach is sound and won't blow up *if* we (i) bundle rather than emit loose
fields, (ii) type only one-vs-rest and passthrough the O(K²) tail, (iii) drop the derivable `names`, and
(iv) make the DE bundle carry `reference`/`method`/gene-axis descriptors. The genuine, unavoidable cost
is the cross-format **asymmetry** (DE round-trips only via AnnData) — which is a fact about the formats,
not our model, and is handled by recording it rather than pretending symmetry.
