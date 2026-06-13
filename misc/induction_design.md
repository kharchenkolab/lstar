# Induction in L★ — current state and an `induce()` design

*Draft for review, companion to `format_coverage.md`. Concerns the model primitive that lets per-group
results (DE, PAGA), coordinate spaces (PCA/UMAP/physical x,y), and collection union axes all be
*ordinary fields* rather than special cases.*

## 1. Current state — induction is specified, not implemented

`docs/model.md` ("Three induction rules") makes induction the model's central economy:
1. **Coordinate axis** — an `embedding`'s columns *are* an axis (the 50 PCA columns are the `pca` axis,
   shared by cell scores and gene loadings). Physical x/y/z is the same with an *observed* coord axis.
2. **Group axis** — the categories of a `label` *are* a derived axis (the 9 `leiden` clusters are the
   `leiden` axis), so per-group results are fields over it.
3. **Union axis** — the union of several axes with `membership` back-relations (collection `cells` over
   samples; metacells over cells; cells over transcript points).

The **code does none of this automatically:**
- `Dataset.add_field` (`model.py:72`) infers a field's span by length-matching to **existing** axes; an
  unmatched dimension raises ("pass span=[...]"). It never creates an axis. A `label` is recognized
  (`_infer_role`) but its category axis is not made.
- `Axis.induced_by` (`model.py:24`) exists and round-trips (`zarr_io.py:50,98`) but is **set by nothing**
  — always `None`.
- Profiles hand-roll axis creation, **inconsistently**: anndata `_coord_axis`/`_pair_coord_axis` invent
  synthetic labels (`pca0…`, `role="coordinate"`); the pagoda2 viewer makes `groups_<grouping>` with the
  real category names but `role="feature"`; conos makes per-sample axes by hand. None set `induced_by`.
- Identity is **name-only**; induced axes are **always materialized** (labels stored); there is **no
  `induce()`**.

**Consequence:** the model's economy is real on paper but reinvented per profile, with no guarantee that
independent results over the same clustering land on the *same* axis (they align only if two profiles
happen to choose the same name *and* labels). The hook (`induced_by`) is there; the mechanism is not.

## 2. The unifying claim

> **A `label` is "integer codes into an induced axis."** Make that one thing first-class and three
> currently-separate problems collapse into it:
> - **categorical fidelity** (Tier 1) — `ordered` and the `-1` missing sentinel are just properties of
>   the codes + the categories axis;
> - **group axes** (Tier 2) — the categories axis a label codes into *is* the group axis DE/PAGA span;
> - **molecule-scale dictionaries** (spatial) — a `gene` label over 10⁸ molecules stored as codes into
>   the existing `genes` axis (not 10⁸ UTF-8 strings).

Same primitive (`label = codes → axis`), three payoffs. Induction is therefore not a Tier-2 feature; it
is the load-bearing model primitive, and getting it right once pays off across categoricals, group
results, designs (a `design` induces a contrast axis), embeddings, collections, and spatial frames.

## 3. `induce()` semantics

A single operation that turns an inducing field into its derived axis, consistently, recording the link.

```
induce(field) -> axis_name                               # Dataset.induce (Py) / lstar_induce (R)
```

| Inducing field (role) | Induced axis | Labels | Axis role | Shared by |
|---|---|---|---|---|
| `label` over X (categorical) | factor axis | the field's **categories**, in category order | **`factor`** (new) | per-group measures/relations |
| `embedding` over (X, C) | coordinate axis `C` | column keys (`pca0…` or given) | `coordinate` | scores + `loading` (genes×C) |
| `design` over samples | contrast axis | coefficient/contrast names | `coordinate`/`design` | DE-over-contrasts |
| union of axes | union axis | union of member labels | `observation` (usually) | `membership` relations to members |

`induce()` (a) computes the canonical labels, (b) creates the axis **iff** an axis with the same
*canonical identity* (§4) doesn't already exist — else returns the existing one, (c) sets
`axis.induced_by = field.name` (and a `rule` tag), (d) assigns a consistent role. The model gains the
dead `induced_by` hook as live, checkable metadata.

**Induction is data-driven, not declared.** A profile is a *static* spec: it encodes the format's
structure and the *rules* ("a `meta.data`/`obs` column → a label/measure over cells; an embedding → a
coordinate axis"), but it cannot know a study's actual columns — those are arbitrary. So induction fires
**at read time, keyed on the data's type**: a profile that ingests a column which *is* categorical/
ordinal induces its factor axis and stores it as codes (§6) *because of the dtype*; a real-valued column
is a `measure` with no axis. A study's bespoke `treatment_arm`/`severity` columns induce their axes with
no profile change. Direct `add_field` follows the same rule: the value's type/role drives induction.
`induce()` is the shared implementation the profiles (and `add_field`) call — not a step a user lists.

The induced categorical axis gets **`role="factor"`** (not the overloaded `feature`, and not a name
prefix) so a reader can tell "these 9 are clusters" from "these 20k are genes." Factor axes are exactly
the ones whose per-group result fields are typically **sparse** (filtered DE/markers over factor×gene).

## 4. Canonical identity (so independent results align)

The point of induction is that *two separate computations over the same clustering produce fields over
the same axis*. Define an induced axis's **identity = (rule, inducing-field name, ordered label set)**:

- **Canonical name = the inducing field's name, bare** (the axis induced by the `leiden` label is
  `leiden`; the field and axis share a name because they are the same concept seen two ways — a
  cells→category map and the category set. Axes and fields are separate namespaces, so this is legal and
  intended).
- **Reuse, don't duplicate:** `induce()` on a label whose canonical axis already exists *with the same
  ordered labels* returns it. Same labels ⇒ same axis ⇒ DE from scanpy and markers from pagoda2 land on
  one axis and align.
- **Collisions are explicit, never silent.** If an axis of that name already exists with *different*
  labels (a real clash — e.g. a label literally named `genes`, or two unrelated clusterings a profile
  named the same), it is an error the writer resolves (rename the field, or suffix the axis) — never a
  silent merge of two different label sets onto one axis.
- `validate()` enforces the invariant: an induced axis's labels **equal** the labels derived from its
  `induced_by` field. This is what makes induction *checkable* rather than conventional.

## 5. Materialize vs derive — resolved by *sharing the categories*

The question: does an induced axis **store** its labels (materialized) or get **computed on read** from
the inducing field (derived)?
- *Materialize* — self-describing, trivial cross-language reads, addressable without the inducing field;
  but the labels duplicate the field's categories and can **drift**.
- *Derive* — single source of truth, no drift, minimal store; but the reader must recompute (and for a
  label, scan it for unique categories).

**The tension dissolves once a label is stored as codes + categories (the categorical encoding, §6):**
the categories are stored **once** (as part of the label field), and the **group axis simply *is* that
categories array** — not a copy (so no duplication, no drift, no scan) and not recomputed (so no extra
reader logic). It is *shared*, not materialized-twice or derived-by-scan.

Recommendation:
- **Group/contrast axes:** the inducing categorical label stores its categories once; the axis points at
  them. (Storage: the categories live in `axes/<name>/labels`; the label field stores only codes +
  `induced_by` to that axis. One array, shared.)
- **Coordinate axes** (synthetic `pca0…`): materialize — they're tiny and the labels are conventional.
- **Union axes:** materialize — computing the union needs all members anyway.
- In all cases the axis carries `induced_by`, so `validate()` can confirm consistency and a reader can
  regenerate if it ever needs to.

Net: **store the labels once (with the inducing field), reference them as an axis** — derive's no-drift
guarantee with materialize's self-description.

## 6. The `label = codes-into-axis` encoding (the concrete addition)

Today labels serialize as `utf8` (one string per element). Add a label encoding that stores **integer
codes into a referenced axis**:

```
field 'leiden' : role=label, span=(cells), encoding=categorical
  codes        : int over cells          # value k = k-th label of the categories axis; -1 = missing
  categories   -> axis 'leiden'          # the induced group axis (ordered)
  ordered      : bool                     # category order is meaningful
```

- When `categories` is **induced from this field's own values**, you get the **group axis** (Tier 2) and
  exact **categorical fidelity** (ordered + `-1`), with the categories stored once.
- When `categories` references an **existing axis** (e.g. `genes`), you get the **dictionary encoding**
  for the molecule `gene` label — codes into `genes`, no per-row strings (spatial, 10⁸ rows).

So `categorical` is the single encoding behind Tier-1 categorical round-trip, Tier-2 group axes, and the
spatial molecule dictionary. `utf8` remains for free-text/identifier labels with no meaningful category
set.

## 7. What to implement (scoped)

1. `Dataset.induce()` / `lstar_induce()` + a `group` axis role; populate `induced_by`.
2. The `categorical` label encoding (codes + categories-axis ref + `ordered` + `-1` missing), in
   `model.py`/`zarr_io.py` and the C++/R readers.
3. `validate()` check: induced axis labels == labels derived from `induced_by`.
4. Retrofit profiles to call `induce()` and set `induced_by` (anndata `_coord_axis`/`_pair_coord_axis`,
   anndata `obs`/`var` categoricals, pagoda2 `groups_*`, conos per-sample/union axes).
5. A canonical naming convention for induced axes.

## 8. Resolved decisions (review pass 1)

- **Induction is data-driven, not declared.** A profile is a static format+rules spec and cannot
  enumerate a study's columns; induction fires at read time keyed on the data's *type* (categorical/
  ordinal → factor axis + codes; embedding → coordinate axis; real → measure). `induce()` is the shared
  implementation profiles and `add_field` call — not a step a user lists. *(was: auto vs explicit)*
- **Bare axis names** (`leiden`, not `groups_leiden`), with **explicit collision resolution** (§4):
  reuse on identical labels, error on a name clash with different labels — never a silent merge.
- **A `factor` axis role** for label-induced (categorical) axes — distinct from `feature`/`coordinate`,
  self-describing; the factor axes are the sparse-result-bearing ones.
- **Ordering:** respect an ordered factor/categorical's category order verbatim; a plain string vector
  with no inherent order → **first-appearance** order.
- **Materialize coordinate & union axes** (tiny, computed-from-members), and **share** factor-axis labels
  with the inducing categorical label (§5) — confirmed.

## 9. Still open

- Whether `add_field` should induce **eagerly** (on add) or **lazily** (on first reference / on write) —
  both are "data-driven"; the question is just *when* within a build session. Lean: eager, since the
  axis is a deterministic function of the field and `induced_by` makes it traceable.
- The exact **`factor` role name** (`factor` vs `category`) and whether the contrast axis from a `design`
  is its own role or also `factor`.
