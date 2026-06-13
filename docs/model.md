# The L★ model

This is a precise, worked description of the L★ data model. It builds on the normative proposal in
[`../misc/Lstar_proposal.md`](../misc/Lstar_proposal.md) (Part 2 — the core model, and Appendix A —
the Zarr schema); read that for the full specification and the formal profile catalog. Here the goal
is to explain *what each construct is and why*, with enough examples that you can read and write an
L★ dataset by hand. Where lstar's current implementation covers only part of the spec, a
**Implemented today** note says so.

## The whole model in one sentence

A dataset is a set of **axes** — the entities you index by — and **fields** — typed data defined over
one or more of those axes. Everything else (embeddings, graphs, trees, factor models, spatial
geometry, fitted models, multi-sample collections, experimental designs) is a *convention* expressed
with axes and fields, not a separate kind of object.

That economy is the whole point. A fixed schema gives you a few named slots (`obs`, `obsm`, `obsp`,
`layers`, …); when a method produces something that doesn't fit a slot, it lands in `uns`/`misc` as an
opaque blob. In L★ a new method adds a *role* and maybe an *axis* — never a schema slot — so the
long tail of analyses (RNA-velocity graphs, gene-regulatory networks, cell–cell communication
tensors, fate probabilities, fitted models) all have a first-class, queryable home.

## Axes — the entities you index by

An **axis** is a named, ordered set of element *labels*, with provenance. The two ordinary axes are
`cells` (barcodes) and `genes` (symbols); a study also has `samples`; analysis introduces more
(`pca`, `umap`, a `leiden` cluster axis, …).

Two properties carry weight later:

- **Observed vs. derived.** An axis is *observed* if the experiment gives it (cells, genes, peaks,
  proteins, transcripts, spots, donors) or *derived* if analysis produces it (clusters, metacells,
  neighborhoods, terminal states, factors, tree nodes). A **coordinate axis** is a special derived
  case whose elements are the *dimensions of a space* — the 50 columns of a PCA, the 2 axes of a UMAP.
- **Label-keyed.** An axis is identified by its labels, not by integer position. A field refers to
  elements by label — which is what lets a field cover only *some* of an axis (partial coverage), or
  carry labels not currently present (a superset). Access is a join on labels; an uncovered element
  is *absent*, not silently filled with `NA`.
- **Namespaced.** Each axis name is local to the dataset or sample that defines it. Sample A's `cells`
  and sample B's `cells` are *distinct* axes. This is exactly what lets a multi-sample study stay a
  collection of heterogeneous parts instead of being concatenated onto one shared axis (see
  [Collections](#collections-a-collection-is-not-a-tensor)).

"Observation" and "feature" are not kinds of axis — they're a *labeling a profile applies* to the two
axes of a measure: in `cells × genes`, cells are the observations and genes the features. Which axis
plays which role is a profile choice.

> **Implemented today.** lstar models `origin` (`observed`/`derived`) and `role`
> (`observation`/`feature`/`coordinate`/…). Per-sample namespaced axes are used by the collection
> profiles (e.g. `cells.GSM5746259`). **Partial coverage is implemented**: a field may cover a subset of
> a span axis via an `index` of positions into it (see Fields, below) — a modality measured on only some
> cells (10x multiome barcode whitelists, CITE-seq dropout) lands over the shared `cells` axis with an
> index, not a separate axis or an NA-padded matrix.

## Fields — typed data over axes

A **field** is typed data over an ordered tuple of axes — its **span**. The arity of the span is the
shape of the thing:

- **arity 1** — a vector: a metadata column, pseudotime, a cluster label.
- **arity 2** — a matrix or relation: counts, PCA loadings, a kNN graph, proportions.
- **arity 3** — a tensor: a `group × group × lr_pair` communication tensor; an `celltype × gene ×
  variant` eQTL table.

A field carries these attributes — but **only `values` is required**; the rest is the *resolved view*,
inferred from the data and overridden only to unlock a role-specific behavior:

| attribute | what it is |
|---|---|
| **span** | the axes it ranges over, e.g. `["cells","genes"]`, `["cells","cells"]` (a relation) |
| **role** | a semantic type tag — what kind of object this is and what invariants apply (below) |
| **coverage / index** | `full` (default), or `partial` — the field covers only a subset of one span axis (`index_axis`), keyed by an `index` of integer positions into it (one per value row along that axis) |
| **values + encoding** | the data, *stored* (`dense`/`csr`/`csc`/`coo`/`ragged`/`raster`) or *virtual* (`recipe`) |
| **state** | for measures: `raw`/`lognorm`/`scaled`/… so a method can refuse an invalid operation |
| **uncertainty** | an optional companion field of the same span (a p-value, a variance) |
| **provenance** | method, parameters, seed, version, and the input fields/axes it was built from |

Two consequences of "only `values` is required" are worth internalizing:

- **Inference does the work.** Hand the binding a character vector and it becomes a `label`; a numeric
  vector a `measure`; a sparse matrix the size of `cells × genes` a `measure` over those axes. You
  specify `state=lognorm` only so that differential-expression guards apply, or force an integer
  vector to be a `label` rather than a `measure` — not to make storage work. A field with a missing or
  wrong role still stores and reads back faithfully; **the role enables behavior, it does not gate
  storage.**
- **Partial factors are normal.** A `celltype` defined on 8,000 of 12,438 cells is *partial* — the
  rest are unannotated, not `NA`. A factor referencing 2,000 barcodes absent from the object keeps
  them as out-of-axis labels that resolve on a later join. Nothing is forced into a dense, NA-padded
  table.

> **Implemented today.** Encodings `dense`, `csc`, `csr`, `coo`, and `utf8` (strings, lstar's encoding
> for `label`/string fields), each optionally with a nullable validity `mask`. `state`, `subtype`,
> relation flags, and `provenance` are stored. **Partial coverage** (`coverage="partial"` + an `index`
> into `index_axis`) is implemented and round-trips across Python, C++, and R. `recipe` (virtual
> fields), `ragged`, and `raster` are specified but not yet implemented.

## Roles — what a field *means*

A **role** declares what kind of object a field is and what invariants apply. A reader that doesn't
recognize a role still sees *a typed field over known axes* — so unknown roles degrade gracefully
rather than breaking. The core vocabulary:

| role | span shape | examples |
|---|---|---|
| **measure** | `(obs, feature)` or arity-1 | counts, normalized/scaled expression, QC metrics, marker stats |
| **embedding** | `(axis, coordinate-axis)` | PCA/UMAP/tSNE scores; spatial x/y/z; CellRank fate probabilities |
| **loading** | `(feature, coordinate-axis)` | PCA loadings; MOFA / cNMF / SCENIC weights (shares the embedding's axis) |
| **relation** | `(axisA, axisB)` | kNN/SNN graph, scVelo transitions, a SCENIC gene network, an orthology map, a deconvolution. Flags `directed`/`weighted`, `subtype` ∈ {similarity, transition, spatial, interaction, regulatory, correspondence, membership} |
| **label** | `(axis) → categories` | a `leiden` clustering, a cell-type annotation. *Induces a group axis* equal to its categories |
| **sequence** | `(axis, position)` ragged | TCR/BCR CDR3 strings, CRISPR-edited barcodes |
| **design** | `(axis)` | a cacoa linear-model design over `samples` (formula + contrast). *Induces a coefficient axis* |
| **transform** | a fitted model | a trained scVI/CellTypist/UMAP projector — *applied* to data, re-appliable to new cells |

The vocabulary is **open**: a new method introduces a role (or a namespaced `x-myorg:foo` term), and
the model doesn't change. The long-tail example from the proposal — RNA velocity — is just four
ordinary fields, no schema slot:

```text
velocity        role=measure    over (cells, genes)   a displacement on the gene basis
velocity_graph  role=relation   over (cells, cells)   directed, subtype=transition
fate            role=embedding  over (cells, terminal_state)   a derived axis of end states
grn             role=relation   over (genes, genes)   directed, subtype=regulatory
```

> **Implemented today.** Profiles emit `measure`, `embedding`, `loading`, `relation`, `label`, and
> `design` (the per-cell `sample` partition). `sequence` and `transform` are specified; `models/`
> (fitted transforms) are not yet written.

## Three induction rules (why the model stays small)

Most of L★'s economy comes from three rules that *induce* a derived axis from a field, so that
downstream results are ordinary fields rather than special cases:

1. **Coordinate axis.** An `embedding` places an axis into a coordinate space; the columns it produces
   *are* an axis. The 50 columns of a PCA **are** the `pca` axis — shared by the cell scores
   (`embedding` over `cells × pca`) and the gene loadings (`loading` over `genes × pca`). One rule
   unifies PCA, UMAP, MOFA, cNMF, SCENIC, CellRank. Physical x/y/z is the same construction with an
   *observed* coordinate axis and no loading.
2. **Group axis.** The categories of a `label` *are* a derived axis — the 9 categories of `leiden`
   become the `leiden` axis — so per-group results (markers, PAGA connectivity, composition,
   communication) are ordinary fields over it. This is why "results that aren't per-cell" need no
   special slot.
3. **Union axis.** The union of several axes, with `relation`/`membership` fields back to each: the
   joint `cells` of a collection over its samples; metacells over cells; cells over transcript points.

> **Implemented today (group axis).** Rule 2 is live: a `label` stored in the `categorical` encoding
> (integer codes + an ordered category set) **induces** a derived **`factor` axis** whose labels *are*
> its categories. `Dataset.induce(field)` (Python) does this — and `add_field` of a categorical fires it
> **automatically** (induction is data-driven: it keys on the value's type, not on a profile listing
> it). Identity is canonical — a label's bare name + its ordered label set — so independent results over
> the same clustering land on **one** axis and align; a name clash with *different* labels is an error,
> never a silent merge. The axis carries `induced_by` back to its field, and `validate()` checks the
> axis labels still equal the field's categories, so induction is *checkable*, not merely conventional.
> Coordinate (rule 1) and union (rule 3) axes are still created by the profiles directly. See
> [`format.md`](format.md) for the on-disk `categorical` encoding.

A few **named patterns** are bundles of the above, not new constructs:

- A **table** (`obs`, `var`, `cellMeta`, `sampleMeta`) is the set of arity-1 fields over one axis,
  presented as a data frame. It's a *view*: each column remains an individual field with its own role,
  state, coverage, and provenance.
- A **tree/hierarchy** (a cluster dendrogram, a lineage) is a derived node axis + a `relation` parent
  field + a `membership` leaf-map.
- **Geometry**: cell/spot **positions** are an `embedding` over an *observed* coordinate axis (position
  is an attribute of the entities), whereas an **image** is a `raster` field over its own `(y, x,
  channel)` base axes (position is the index). `transform` fields relate coordinate frames.

## Provenance — a dataset is a provenance graph

Every field and axis records the inputs it was built from, so a dataset *is* a provenance graph:
`umap` ← `knn` ← `pca` ← `data` ← `counts`. That makes results reproducible, lets a guard refuse an
operation on a measure of the wrong `state` (clustering raw counts), and gives an automated agent an
auditable lineage.

## A whole dataset, listed

A complete single-cell object is a short, language-independent listing of its axes and fields — this
is the canonical example from the proposal (§2.5):

```text
Sample "pbmc_donorA"
  axes
    cells    12,438 barcodes            [observed]
    genes    18,200 symbols             [observed]
    pca      50 dims                    [derived: induced by the 'pca' embedding]
    umap     2 dims                     [derived: induced by the 'umap' embedding]
    leiden   9 clusters                 [derived: induced by the 'leiden' label]
  fields
    counts         measure    cells × genes    state=raw      sparse
    data           measure    cells × genes    state=lognorm  sparse
    n_umi          measure    cells            (QC)
    pct_mt         measure    cells            (QC)
    pca            embedding  cells × pca      induces the 'pca' latent axis
    pca_loadings   loading    genes × pca      shares the 'pca' axis
    umap           embedding  cells × umap     induces the 'umap' axis
    knn            relation   cells × cells    similarity, weighted
    leiden         label      cells → leiden   induces the 'leiden' factor axis (categorical)
    markers.lfc    measure    leiden × genes   on the 'leiden' axis (canonical factor-first)
    markers.padj   measure    leiden × genes   + uncertainty
    dendrogram     hierarchy  over leiden
```

Notice there are no special slots — markers (a per-cluster result) and the dendrogram (a tree) are
just fields over the induced `leiden` axis. See
[`examples.md` §1](examples.md#1-build-write-read-validate-python) for this built in code.

## Collections — a collection is not a tensor

A `kind="collection"` dataset represents *heterogeneous samples* — possibly different cells **and**
different gene sets, even different species — **without** concatenating them into one matrix. Because
axes are namespaced, each member keeps its own `cells`/`genes`; the joint analysis lives over a
*derived union* `cells` axis. The canonical shape (the conos profile, proposal Appendix B.5):

```text
Collection "BM_integration"   (2 human, 2 mouse)
  members (each a sample with its OWN axes)
    human_1   cells 9,102 × genes 19,400
    mouse_1   cells 7,330 × genes 21,800     (a DIFFERENT gene axis)
  axes
    samples   4
    cells     union = 31,934   ids = (sample, cell)   [derived]
    joint     12 clusters                              [derived]
  fields
    sampleMeta   table     over samples     { species, donor, tissue }
    inclusion    relation  cells → sample.cells          which union cell came from which sample
    orthology    relation  genes(human) × genes(mouse)   cross-species feature map
    knn          relation  cells × cells                 the joint integration graph
    umap         embedding cells × umap                  joint layout
    joint        label     cells → joint                 clusters spanning all samples
```

No concatenated matrix exists: `featureVector("CD8A")` *gathers* per sample — human cells return CD8A,
mouse cells the orthologous Cd8a, absent elsewhere. Alignment within a sample is legitimate (its facets
share cells); across samples you keep a collection joined by a graph. A single Pagoda2 object is the
`sample` unit; a Conos object is the `collection` unit.

> **Implemented today.** The R profiles ingest a **Conos** object and a split **Seurat v5** assay as
> collections: a `samples` axis, per-sample `cells.{s}`/`genes.{s}` axes and `counts.{s}` measures, a
> union `cells` axis, a `sample` design label, and the joint embedding / clusters / graph. The
> orthology and inclusion relations and the four packaging modes (3.2) are specified but not yet
> emitted. See [`examples.md` §9](examples.md#9-r-a-conos-collection).

## Validation invariants

A conforming store satisfies (what `lstar.validate` checks today): every field `span` references
existing axes; a field's shape matches its axes' lengths — or, for **partial coverage**, `len(index)`
along `index_axis` (with `index` in range, and `index_axis` one of the span axes); a `relation` (csc/csr
2-D) spans exactly two axes; a 1-arity `label` is `utf8` and as long as its axis; an induced factor
axis's labels match its inducing field's categories. Roles/states outside the core registry are
*warnings*, not errors — the vocabulary is open.

---

Next: the on-disk [format spec](format.md) (Appendix A realized), and the worked
[examples](examples.md). The full normative text — including the profile rule catalog for AnnData,
Seurat, pagoda2, Conos, and cacoa — is in [`../misc/Lstar_proposal.md`](../misc/Lstar_proposal.md).
