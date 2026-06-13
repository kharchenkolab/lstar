# Tier 2 — implementation plan

*Build plan for `format_coverage.md` Tier 2, per the design in `induction_design.md` and the risk
projection in `tier2_plan.md`. Sequenced P1→P3 with a correctness gate (conformance + tests) between
phases. P1 is the heavy, cross-language piece; P2/P3 are mostly Python + R profile work.*

## On-disk layout (the contract all four languages implement)

A categorical label `leiden` over `cells`, sharing categories with its induced factor axis `leiden`:

```
fields/leiden/.zattrs lstar = { role:"label", span:["cells"], encoding:"categorical",
                                categories:"leiden", ordered:false }   # categories = factor-axis name
fields/leiden/codes           = int32 zarr array over cells; value k -> categories[k]; -1 = missing
axes/leiden/.zattrs   lstar = { kind:"axis", role:"factor", origin:"derived", induced_by:"leiden" }
axes/leiden/labels            = utf8 (the categories), stored once, as any axis
```

The categories live **once**, as the factor axis's labels; the field stores only `codes` + the axis
name + `ordered`. This is the single mechanism behind categorical fidelity, factor axes, and (with
`categories` pointing at an *existing* axis like `genes`) the spatial molecule dictionary.

**The one in-memory API decision:** a lightweight `Categorical` value wrapper (`codes:int[]`,
`categories` (= the factor axis labels, by reference), `ordered:bool`; `np.asarray()` → decoded labels),
pandas-free. Profiles convert pandas `Categorical` / R `factor` ↔ this. (Alternative: add
`categories`/`ordered` fields to `Field` and keep `values`=codes — rejected: makes `values` ambiguous.)

---

## P1 — the `categorical` encoding (foundation; all four languages + conformance)

Delivers Tier-1 dtype fidelity (ordered + `-1` missing) immediately, and the factor-axis substrate.

- **Python model** (`model.py`): a `Categorical` value class; `add_field` recognizes a pandas
  `Categorical` (or a `(codes, categories, ordered)`) → encoding `categorical`; `_infer_role`→label.
  Plain string arrays stay `utf8` until `induce()` (P2).
- **Python IO** (`zarr_io.py`): `_write_values` branch for `categorical` → write `codes` (int32, `-1`
  fill) + meta `categories`/`ordered`; `_read_values`/`_lazy_values` → build `Categorical` (categories
  from the referenced axis labels). Mirrors the existing `utf8` path.
- **C++ core** (`lstar.hpp`): `Field` gains a categorical case (`codes` NdArray + `categories` axis
  name + `ordered`); `read()`/`write()` handle it (codes array + the categories are the referenced
  axis's labels, already read). Re-vendor to `R/inst/include`.
- **R** (`lstar_cpp.cpp`, `lstar.R`): R `factor`/`ordered` ↔ `categorical` natively (factor = 1-based
  codes + `levels` + `ordered`); map to 0-based codes + categories axis on write, rebuild `factor` on
  read.
- **JS** (`js/core/reader.ts`): read `categorical` → decode codes via the categories axis (for the
  viewer); `view.ts` metadata path uses it for crossfilter.
- **Validate** (`validate.py`): a `categorical` field's `categories` axis must exist and its length must
  bound the codes; `-1` only where allowed.
- **Conformance + tests:** new `conformance/categorical.sh` (a categorical with `ordered=TRUE` + a
  missing `-1` round-trips Python↔C++↔R, labels/order/missing preserved) + `python/tests/test_categorical.py`;
  add to `run.sh`/CI. **Gate:** byte-identical round-trip across languages before P2.
- **Docs:** `format.md` (the encoding) + `model.md` (label = codes-into-axis).

## P2 — `induce()` + `factor` role + canonical identity

- **Model** (`model.py`): `Dataset.induce(field) -> axis_name` implementing the three rules
  (categorical label → `factor` axis with categories=field categories, name=field name, `induced_by`
  set, collision check per `induction_design.md §4`; embedding → coordinate axis; union → union axis).
  Add `"factor"` to the role vocabulary. **Eager auto-induce**: `add_field` of a `Categorical` label
  induces its factor axis automatically (data-driven).
- **Validate**: induced axis labels == labels derived from `induced_by`; `factor` role on label-induced
  axes; name-collision rule (reuse on identical labels, error on different).
- **Retrofit profiles** to induce + set `induced_by`, replacing hand-rolled `add_axis`:
  - `profiles/anndata.py`: `obs`/`var` categoricals → `categorical` fields + induce (replaces the
    string-coercion in `_by_dtype_series`); `_coord_axis`/`_pair_coord_axis` set `induced_by`.
  - `R/R/profile_pagoda2.R`: `groups_<g>` → `induce()` of the clustering label; **drop `role="feature"`
    → `factor`**, and **fix the orientation** (canonical `(factor, genes)` — currently stats are
    `(groups,genes)` but markers `(genes,groups)`).
  - `R/R/profile_seurat.R` (Idents + factor meta.data), `R/R/profile_sce.R`, `R/R/profile_conos.R`:
    factors → categorical + induce; record `induced_by`.
- **Tests:** induction round-trips; two inducers with identical labels → one shared axis (alignment);
  `validate` catches a drifted induced axis. **Gate** before P3.

## P3 — DE + pseudobulk (the per-group results)

**Representation (settled):** a DE result is a **bundle** of optional measures over **`(factor, genes)`**
(fixed orientation), sparse-or-dense, sharing one sparsity pattern, with field meta carrying
`reference`/`method`/`groupby`:

```
fields/de.<factor>.lfc   role=measure span=[<factor>, <gene-axis>]  subtype=de  meta={reference, method, groupby}
fields/de.<factor>.padj  role=measure span=[<factor>, <gene-axis>]  subtype=de   (optional)
fields/de.<factor>.score role=measure ...                                        (optional)
fields/pb.<factor>.mean  role=measure span=[<factor>, <gene-axis>]  subtype=pseudobulk
fields/pb.<factor>.frac  role=measure ...   (+ sufficient stats sum/sumsq/nexpr as today)
```

- **AnnData read** (`anndata.py`): parse `uns['rank_genes_groups']`; **only when
  `params['reference']=='rest'`** (one-vs-rest) → build the `(factor, genes)` bundle: resolve the
  **gene axis from `params['use_raw']`** (`genes` vs `genes_raw`/HVG), scatter each group's
  score/lfc/pval/padj from ranked order to gene-keyed positions via `names`. **Pairwise/reference-group
  → passthrough** (verbatim in the lossless `uns` tail, not typed). Don't store `names`.
- **AnnData write**: gather the bundle → per-group ranked arrays (argsort scores → regenerate `names`),
  emit `uns['rank_genes_groups']` structured arrays + `params`.
- **Pseudobulk**: lift the `lstar_cpp_col_sum_by_group` path out of `profile_pagoda2.R` into a shared
  helper emitting canonical-orientation `(factor, genes)` stats; profiles call it.
- **Seurat/SCE** (`profile_seurat.R`/`profile_sce.R`): **capture** a marker data.frame if present →
  bundle; **export** → `@misc`/`metadata` + record in `dropped` (no native slot — the documented
  asymmetry).
- **Tidy view**: `markers(ds, factor)` (Py) / `lstar_markers()` (R) → long-form `(group, gene, lfc,
  padj)` frame for ergonomics.
- **Tests/conformance:** AnnData DE read→L*→write reproduces `rank_genes_groups` on the representable
  (one-vs-rest) content; pseudobulk matches the in-memory reduction; orientation is canonical;
  pairwise DE survives as passthrough.

## P4 — deferred (list, don't build now)

PAGA (`uns['paga']` → relations over the factor axis + diffmap/pseudotime), dendrograms, `@commands`/
param provenance. Each is "ordinary fields over the existing factor axis" — no new machinery, just
profile plumbing — so they slot in later without rework.

---

## Sequencing, effort, and the gates

```
P1 categorical encoding ──gate(cross-lang round-trip)──▶ P2 induce()/factor ──gate(induction+validate)──▶ P3 DE/pseudobulk ──▶ (P4)
   heaviest: Py+C+++R+JS+conformance        Python model + retrofit profiles      Python+R profiles + tidy view
```

- **P1 is the long pole** — a new encoding across four readers + conformance, comparable to adding
  `utf8`. P2/P3 are mostly Python `model.py` + the R/Python profiles.
- **Risk containment is structural** (per `tier2_plan.md`): bundle (not loose fields), type only
  one-vs-rest, drop the ranking `names`, descriptors on the DE bundle, `-1`-padding for subset factors.
- **Decisions still to confirm before P3** (from `tier2_plan.md §3` + the alignment note): canonical DE
  orientation `(factor, genes)` ✔ proposed; one-vs-rest-typed / pairwise-passthrough ✔; gene-axis from
  `params` ✔; Seurat/SCE = capture-not-bridge ✔; `-1`-padding default vs implement partial `coverage`
  (defaulting to `-1`-padding; `coverage` only if a sparse-subset use case appears).
- **What stays out of scope:** partial `coverage` (unless needed), pairwise/conditional DE typing, P4.
