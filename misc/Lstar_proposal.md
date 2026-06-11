# L★ — a data model and interchange format for single-cell and spatial omics

*Proposal, working draft. (The model is named **L★**; where the star is inconvenient to type, "L*" is the same name. Filename `L0star_contract.md` is legacy; a rename to `Lstar_proposal.md` is reasonable.)*

## Abstract

L★ is a proposed common model and on-disk interchange format for single-cell and spatial omics datasets. It rests on two constructs — **axes** (the entities one indexes by) and **fields** (typed data defined over axes) — and represents embeddings, graphs, trees, factor models, spatial geometry, fitted models, multi-sample collections, and experimental designs as conventions over those two constructs rather than as fixed slots. Existing ecosystem formats (AnnData/MuData, Seurat, the pagoda-verse, SOMA, SingleCellExperiment) are expressed as **profiles**: precise, versioned, bidirectional mappings between a native object and L★. The format is realized on Zarr. The central premise of this proposal is that a small, general core, together with well-governed profiles and a shared vocabulary, can provide lossless interchange and a natural home for the long tail of analyses that current fixed schemas relegate to `uns`/`misc`. That premise holds only if the profile and vocabulary standardization (Part 5) is treated as seriously as the model itself.

The document is normative in Parts 2–4 (the model, the serialization, the profile mechanism) and discursive in Parts 1 and 5; the hardened Zarr schema and the profile catalog are Appendices A and B. Names and the exact vocabulary are draft.

---

# Part 1 — Motivation

## 1.1 The problem

Every mature single-cell container fixes a *schema*: one `cells × genes` matrix plus a fixed set of named slots (AnnData's `obs`, `var`, `obsm`, `varm`, `obsp`, `layers`, `uns`; Seurat's assays / reductions / graphs; SingleCellExperiment's assays / `reducedDims` / `colData`). For the routine workflow this works well — the slot names are, in effect, the vocabulary. Two situations strain it.

**The long tail of analyses has no natural home.** Off the main path, structures are stored as opaque entries in `uns`/`misc` or require bespoke objects:

- **RNA velocity** — the velocity vectors fit in `layers`, but the velocity *graph* becomes an opaque `uns['velocity_graph']`, and CellRank's fate probabilities are placed in `obsm` as though they were an embedding.
- **Gene regulatory networks** — a TF→target graph defined on the *gene* axis has no modeled location, so SCENIC distributes a separate `.loom`.
- **Cell–cell communication** — a `cluster × cluster × ligand-receptor` tensor is defined on the *cell-type* axis, which the schema does not name.
- **Fitted models** — a trained scVI, scArches, or foundation model is not data; it is a function applied to data.

**Heterogeneous multi-sample data cannot be a single tensor.** Four samples — two human, two mouse — cannot be concatenated into one matrix, because they do not share a gene set; yet a single-tensor schema requires it. Integration methods (conos, and multi-sample analysis generally) keep samples as a *collection* linked by a graph, which the fixed schemas cannot represent natively.

## 1.2 The approach

The approach is not to abandon schemas but to factor them differently. Each structure above is, at base, the same two things: a set of labeled **axes** (cells, genes, clusters, samples, factors, transcripts, terminal states, ligand-receptor pairs, …) and typed **fields** defined over some of those axes. L★ fixes a small, general schema of exactly those two constructs, and expresses the ecosystem-specific structure (which matrix is primary, what the slots are named) as **profiles** layered on top. The familiar formats become profiles; the long tail and heterogeneous multi-sample data become ordinary fields and axes; all of it is serialized on Zarr. The velocity example, in L★, is a set of first-class, queryable, provenance-bearing fields, with no schema change:

```
velocity        field  role=vector       over (cells, genes)            a displacement on the gene basis
velocity_graph  field  role=relation     over (cells, cells)            directed, transition
fate            field  role=probability  over (cells, terminal_state)   a derived axis of end states
grn             field  role=relation     over (genes, genes)            directed, regulatory
```

The objective is an ecosystem-neutral interchange standard: a conforming dataset can be read in any language, converted losslessly between ecosystem formats where they overlap, and converted with a recorded, explicit boundary where they do not.

## 1.3 Scope and non-goals

**In scope:** the core model (Part 2); a normative Zarr serialization (Part 3 and Appendix A); a profile mechanism with precise, bidirectional mappings for the major ecosystems (Part 4 and Appendix B); and the standardization and governance required for profiles to interoperate (Part 5).

**Out of scope:** a new storage engine (L★ uses Zarr/OME-NGFF, as AnnData-zarr and SpatialData do); a particular in-memory API (each language binds idiomatically — conformance is defined by the on-disk schema and the profiles, not by an API); a compute or analysis framework; and replacing AnnData, Seurat, or SOMA, which instead become profiles. L★ is a model, a format, and a profile system, not a library.

## 1.4 Relation to existing standards

- **AnnData / MuData (scverse).** The de facto single-cell container and a Zarr format, well suited to the single-tensor case and (MuData) to shared-cell modalities. L★ generalizes to heterogeneous collections and to arbitrary axes and structures, and treats AnnData/MuData as profiles. *Why not simply add more AnnData slots?* Because the additions do not converge: a slot for velocity graphs, then one for gene-regulatory networks, then one for cell-type-pair communication tensors, then one for cross-species cell-type mappings — each is a structure over a different pair of axes (`cells × cells`, `genes × genes`, `celltype × celltype`, `celltypeA × celltypeB`). The axis-and-field construction expresses all of them uniformly, and represents AnnData itself as one profile over it.
- **SpatialData (scverse).** SpatialData defines an elements model for spatial data — images, label masks, shapes, and points, related by coordinate transforms. L★ can express these as fields and axes (2.3), which allows a SpatialData object and other spatial formats (Visium, Xenium, Giotto) to be represented in, and bridged through, a single model rather than converted pairwise.
- **TileDB-SOMA / the SOMA specification.** TileDB-SOMA implements SOMA: a way of storing single-cell data as on-disk arrays that can be queried directly from cloud storage without loading the whole dataset, together with a programming interface for that access. Its emphasis is storage and scalable access rather than a conceptual interchange model. L★ is complementary — it defines the model and a Zarr realization, and could in principle be stored on a SOMA backend.
- **Bioconductor SingleCellExperiment / MultiAssayExperiment.** These embody the same principle L★ generalizes: a shared object that many independent packages read and write, so that tools interoperate without pairwise converters. The Bioconductor community has the deepest experience making this work at scale — hundreds of interoperating packages built on a common class — and therefore both relevant expertise and a stake in a cross-language version of the idea. L★ extends the shared-object approach beyond a single language and to the structures (multi-sample collections, networks, designs, fitted models) that SCE/MAE do not natively cover; a Bioconductor profile (Appendix B) is a natural point of collaboration.
- **OME-NGFF.** The imaging Zarr conventions that L★ reuses for raster encodings.

L★ does not compete on storage technology. It proposes a common conceptual model and a profile system so that these ecosystems interoperate, and gives the off-main-path structures a first-class representation.

## 1.5 Benefits

Each benefit below is realized only if the standardization of Part 5 holds; that condition is stated where relevant.

- **Lossless, mechanical interchange.** Because a field names its axes, export to a native format is a deterministic mapping and import is its inverse; a sequence such as "read an h5ad, normalize in pagoda, integrate in conos, test in cacoa, write back an h5ad" is a chain of profile mappings over one logical object, and the parts that a given target cannot hold are *reported* rather than lost (4.3). This requires the shared core vocabulary of 5.1.
- **A representation for the long tail.** Velocity graphs, regulatory networks, communication tensors, fate probabilities, fitted models, and lineage trees are first-class typed fields with provenance, rather than opaque `uns` entries. A new method adds a role and an axis, not a schema slot, so the format does not change as methods do.
- **Heterogeneous collections without a special case.** Because axes are namespaced, a cross-species or mixed-modality collection is an ordinary object — the conos example in Appendix B is a concrete case that concatenation cannot represent at all.
- **Provenance and validity by construction.** The field-to-input provenance graph (2.4) makes any result reproducible and allows a function to refuse an invalid operation (for example, differential expression on a measure whose `state` is `raw`).
- **Packaging that scales from a workstation to an atlas.** One model serves a standalone sample, a lazily fetched remote collection, and a portable single-file bundle, with lossless conversion among them and with samples shared across collections without duplication (3.2).
- **A substrate suited to automated analysis.** A small, stable, self-describing model with typed fields and explicit provenance is well suited to inspection and safe modification by analysis agents ([`design_for_ai.md`](./design_for_ai.md)).

In summary, the routine case remains simple (two axes and a dozen fields), while the long-tail, multi-sample, spatial, and cohort cases are expressed in the same vocabulary — a property that a fixed schema cannot provide and that only standardized profiles can deliver.

---

# Part 2 — The core model (normative)

An **axis** is a set of entities one indexes by — the cells in a sample, the genes, a set of clusters, the samples in a study. A **field** is data defined over one or more axes. Axes are the entities; fields are everything measured or computed about them. The most ordinary object — a count matrix with barcodes and symbols, plus a cluster label — is exactly two axes (`cells`, `genes`) and three fields (`counts` over `(cells, genes)`, `cluster` over `(cells)`); a UMAP adds a field over `(cells, 2 layout dimensions)`; markers add a field over `(genes, clusters)`.

## 2.1 Axes

An axis is a named, ordered set of element labels, with provenance. By origin it is **observed** (given by the experiment: cells, genes, peaks, proteins, transcripts, spots, bins, donors, variants, image pixel/position axes) or **derived** (produced by analysis: clusters, metacells, neighborhoods, regions, clones, terminal states, ligand-receptor pairs, coefficients, factors/programs/regulons, tree nodes).

A **coordinate (dimension) axis** is a special case whose elements are the components of a space into which other axes are placed by an `embedding` (2.3). It is **observed** for physical x/y/z (a metric space) or **derived/latent** for reductions (principal components, UMAP dimensions, factors); latent coordinate axes differ only in that they typically carry a feature-side `loading`.

An axis has two properties that the rest of the model relies on.

- It is **label-keyed**: an axis is identified by its element labels (barcodes, gene symbols, cluster names), not by integer position. A field therefore refers to elements by label, which is what allows a field to cover only some of an axis's elements, or to refer to labels that are not currently present (2.2).
- It is **namespaced**: each axis name is local to the dataset or sample that defines it — sample A's `cells` axis and sample B's `cells` axis are distinct axes. Combining samples therefore does not force them onto a shared axis, which is what allows a multi-sample dataset to remain a collection of heterogeneous parts rather than a single concatenated matrix (3.2, Appendix B).

("Observation" and "feature" are not kinds of axis but a labeling that a profile applies to the two axes of a measure: in a `cells × genes` matrix, cells are the observations — the independent units of which there are many — and genes are the features measured on each. Which axis plays which role is a profile choice.)

## 2.2 Fields

A field is typed data over an ordered tuple of axes. It has the following attributes:

- **span** — the axes it ranges over. arity 1 is a vector (a metadata column, pseudotime, a label); arity 2 is a matrix or relation (counts, loadings, a graph, genotype, proportions); arity 3 is a tensor (a communication `group × group × lr_pair`; an eQTL `celltype × gene × variant`).
- **role** — a semantic type tag declaring what kind of object the field is, and what invariants apply (the vocabulary is in 2.3).
- **index / coverage** — the field is keyed by labels and need not cover its whole axis: it may be defined on a subset of an axis (partial coverage) or carry labels not currently in the axis (a superset). Access is a join on labels; uncovered elements are *absent*, not silently filled with NA.
- **values** — the data, in some **encoding** (dense, sparse, edge-list, ragged, raster), and either *stored* or *virtual* (`recipe`): a virtual field is computed on demand from a recorded recipe over its provenance inputs — for example, normalized `data` kept as a recipe over `counts` and library sizes rather than as a materialized matrix (this is precisely pagoda2.1's matrix-views).
- **state** (for measures) — raw, lognorm, scaled, clr, …, so that methods can refuse invalid operations.
- **uncertainty** (optional) — a companion field of the same span; for example, scVI provides a variance alongside each expression value, and a cacoa shift carries a p-value and z-score.
- **provenance** — method, parameters, seed, version, time, and references to the input fields and axes.

**Only `values` is required.** The full list of attributes is the *resolved* view, not a form one fills in: a binding given a vector infers the rest — a factor or character vector becomes a `label`, an ordered factor an ordinal label, a numeric vector a `measure`, a named vector shorter than the axis a partial-coverage field. One specifies more only to enable a role-specific behavior (marking a numeric `state=lognorm` so that differential-expression guards apply; forcing an integer vector to be a `label` rather than a measure). A field with an unspecified or mis-inferred role still stores and reads back faithfully; the role enables behavior, it does not gate storage.

**Partial and incompletely specified factors are normal.** A `celltype` defined on 8,000 of 12,438 cells is partial — the remaining cells are unannotated, not NA; a factor referencing 2,000 barcodes absent from the object retains them as out-of-axis labels that resolve on a later join and are otherwise ignored. Nothing is forced into a dense, NA-padded table.

## 2.3 Roles, induced axes, and named patterns

A **role** is a semantic type tag on a field declaring what kind of object it is and what invariants apply. A consumer that does not recognize a role still sees a typed field over known axes. The recommended core roles (the vocabulary is versioned and extensible — Part 5):

- **`measure`** — a quantity over two axes (e.g. `cells × genes`): counts, normalized expression, antibody intensity, pathway activity, spot-by-celltype proportions. Carries a `state`.
- **`embedding`** — coordinates placing one axis into a coordinate space (`axis × coordinate-axis`): PCA and UMAP scores (latent dimensions); spatial x/y/z (observed); CellRank fate probabilities (`cells × terminal-state`).
- **`loading`** — feature contributions to a coordinate axis (`feature × coordinate-axis`): PCA loadings; MOFA, cNMF, and SCENIC weights. A loading shares its coordinate axis with the matching embedding.
- **`relation`** — edges between two axes (`axisA × axisB`): a kNN graph, conos alignments, an scVelo transition matrix, a SCENIC gene network, an orthology map, a deconvolution. Carries flags `directed` and `weighted`, a `subtype` (similarity, transition, spatial, interaction, regulatory, correspondence, membership), and an optional constraint (simplex, partition).
- **`label`** — a categorical assignment of an axis (`axis → categories`): a `leiden` clustering, a cell-type annotation. A label induces a group axis equal to its set of categories.
- **`sequence`** — variable-length per-element data (ragged): TCR/BCR CDR3 strings, CRISPR-edited barcodes.
- **`design`** — an experimental-design specification over an axis (a formula and contrast): a cacoa linear-model design over `samples`. A design induces a coefficient axis.
- **`transform`** — a fitted, applicable model, described next.

Some artifacts produced during analysis are not data but **fitted models** — a trained scVI network, a CellTypist classifier, a fitted UMAP projector. L★ represents these with the `transform` role. Unlike every other field, which is read, a `transform` is applied to data to produce new fields, and (the property that distinguishes it) can be re-applied to *new* cells. Its stored values are the learned parameters (neural-network weights, a PC rotation, per-gene kinetic rates — themselves fields, often large and held out-of-core) together with an apply contract naming what it consumes and produces (`scVI: cells × genes → cells × latent`; `CellTypist: cells × genes → a label`). Storing the model, rather than only its output, is what allows a new sample to be placed into an existing reference space (the reference-to-query step in scArches and Azimuth) and what records which parameters produced a given embedding.

**Three induction rules** account for most of the model's economy:

1. **Coordinate axis.** An `embedding` places an axis into a coordinate space; when the dimensions are derived and carry a feature-side `loading`, the shared axis is latent. Concretely, the 50 columns produced by PCA *are* the `pca` axis, shared by the cell scores (`embedding`) and the gene loadings (`loading`). The single rule unifies PCA, UMAP, MOFA, cNMF, SCENIC, and CellRank macrostates. Physical x/y/z is the same construction with an observed coordinate axis and no loading.
2. **Group axis.** The categories of a `label` *are* a derived axis — the 9 categories of a `leiden` label become the `leiden` axis — and per-group results (markers, PAGA connectivity, communication, composition) are ordinary fields over it. This is why results that are not per-cell are not a special case.
3. **Union axis.** The union of several axes, with `relation`/`membership` inclusion fields back to each: the joint `cells` of a collection over its samples; metacells over cells; cells over transcript points; donors over cells.

A few **named patterns** are bundles of the above, not new constructs:

- A **table** is a named bundle of the arity-1 fields over one axis, presented as a data frame — `obs`, `var`, `cellMeta`, and `sampleMeta` are tables over the cells, genes, and samples axes. It is a view: the columns remain individual fields with their own role, state, coverage, and provenance (uniform in their axis, not in coverage).
- A **tree/hierarchy** is a derived node axis, a `relation` parent field, and a `membership` leaf-map (a cluster dendrogram, a lineage phylogeny, a cell-type tree).
- **Geometry** distinguishes two cases: cell or spot **positions** are an `embedding` field over an observed coordinate axis (position is an attribute of the entities one indexes by), whereas an **image** is a raster field over its own `(y, x, channel)` base axes (here position is the index). `transform` fields relate coordinate frames.

## 2.4 Provenance

Because every field and axis records the inputs from which it was built, a dataset is a **provenance graph** — `umap` records that it was built from `knn`, which was built from `pca`, which was built from the normalized `data`, which was built from `counts`. This makes results reproducible, allows a guard to refuse an operation on a measure of the wrong `state` (clustering raw counts, for example), and provides an auditable lineage for automated analysis (1.5, [`design_for_ai.md`](./design_for_ai.md)).

## 2.5 An example dataset

A whole single-cell object is a short, language-independent listing of its axes and fields:

```
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
    leiden         label      cells → leiden   induces the 'leiden' group axis
    markers.lfc    measure    genes × leiden   on the 'leiden' axis
    markers.padj   measure    genes × leiden   + uncertainty
    dendrogram     hierarchy  over leiden
```

---

# Part 3 — Serialization and packaging (normative)

## 3.1 The store

On disk a dataset is the two registries serialized: a Zarr group with child groups `axes/`, `fields/`, and `models/`. All L★ metadata is placed under an `"lstar"` key in each group's `.zattrs`, so that it does not collide with native metadata, and consolidated metadata (`.zmetadata`) is required for access over HTTP. A field's `span` references axes by name, so placement and reading are deterministic lookups rather than inference from array shape — a `cells × cells` relation and a `cells × genes` measure are unambiguous even when the two counts are equal. The `state`, `coverage`, relation flags, and `uncertainty` companion are the only semantics a generic reader must honor; an unrecognized `role` or `subtype` degrades to "a typed field over known axes." The exact group, array, and attribute schema is given in Appendix A.

## 3.2 Packaging

Because a collection references its samples by identity (axes are namespaced), packaging is decoupled from logical structure:

```
1) Standalone sample — one self-contained store, a valid object on its own:
   human_1.lstar.zarr/

2) Collection, bundled — one store physically containing its members (portable; one .zip):
   BM_integration.lstar.zarr/
     samples/human_1/ …          members embedded
     axes/ fields/ models/        collection-level (the union axis, joint fields, links, results)

3) Collection, by reference — a thin overlay; members live wherever they are:
   BM_integration.lstar.zarr/
     .zattrs lstar.members: [ {id: human_1, uri: ./human_1.lstar.zarr,        hash: …},
                              {id: mouse_1, uri: s3://atlas/mouse_1.lstar.zarr, hash: …} ]
     axes/ fields/ models/        collection-level only (no member bytes)

4) Mixed — some members embedded, some referenced.
```

A sample store is always openable on its own; the same sample may be referenced by several collections without duplication; a referenced collection loads lazily over HTTP. A referenced collection can be bundled (members inlined) or a bundled one split (members extracted) without loss. The single obligation: a reference pins an `id` and a content `hash`, and the provenance graph records the inputs of each joint field, so that a changed referent renders the dependent fields detectably stale.

Language bindings are non-normative: each language exposes L★ idiomatically (Python over the Zarr store, R, JavaScript in a viewer). Conformance is defined by the schema and the profiles, not by any API.

---

# Part 4 — Profiles (the load-bearing specification)

This is where L★ connects to existing software. A **profile** is a precise, versioned, bidirectional mapping between a native format and L★. Because L★ is the more general model, a profile is most naturally read as a mapping *from* the native structure *into* L★ (a reader); the writer is the same mapping applied in reverse, lossy wherever a rule is read-only.

## 4.1 What a profile is

A profile is a versioned document with two parts: an *identity* (which native axes become the L★ observation and feature axes) and an ordered list of *rules*, each mapping a native location to an L★ field signature:

```
NATIVE_LOCATION   <dir>   field "NAME" : ROLE over (AXES) [attrs]      # notes
```

where `<dir>` is `<->` (bidirectional), `->` (read-only; lossy on write), or `<-` (write-only). For example, two rules from the `anndata` profile:

```
.X            <->  field "X"   : measure   over (cells, genes)  state=…   # the primary matrix
.obsm[{k}]    <->  field "{k}" : embedding over (cells, dim_axis(k))      # X_pca→pca, X_umap→umap, spatial
```

The full rule syntax, the helper functions, and the five reference profiles — `anndata`, `seurat`, and `pagoda2` (sample-level) and `conos` and `cacoa` (collection-level), each with a worked instance — constitute the normative bridge specification in Appendix B. A profile is what makes a given AnnData, Seurat, or conos object become L★ axes and fields, and back.

## 4.2 Readers, writers, converters, and conformance

For a given profile `P` — for example `P = anndata` — a conforming implementation provides three operations:

- a **reader** `read_P`, which applies the profile's rules forward to turn a native object into an L★ dataset;
- a **writer** `write_P`, which applies the bidirectional rules in reverse; any field matched only by a read-only rule, or by no rule, is written to the profile's sidecar (`uns`/`misc`) and recorded in an `lstar/dropped` manifest, so that loss is never silent;
- a **converter** `convert(P → Q) = write_Q(read_P(·))`, which moves data between two formats and reports what survives.

Conformance is established by two round-trip tests that every profile must publish: (a) `write_P(read_P(native))` equals `native` on the profile's declared lossless subset; and (b) `read_P(write_P(dataset))` preserves every L★ field the target can hold. A conformance level (core, extended, spatial) states which roles and capabilities a reader and writer implement.

## 4.3 A worked conversion

Consider converting an AnnData object to a Seurat object through L★. The input `.h5ad` contains: a log-normalized `X`; raw `layers['counts']`; `obsm['X_pca']` with `varm['PCs']`; `obsm['X_umap']`; an `obsp['connectivities']` neighbor graph; an `obs['leiden']` clustering; and an `uns['velocity_graph']` from scVelo.

Reading it with the `anndata` profile produces these L★ fields:

- `data` and `counts` — measures over `(cells, genes)`, states lognorm and raw;
- `pca` (embedding) and `pca_loadings` (loading), sharing the `pca` coordinate axis;
- `umap` — an embedding;
- `knn` — a relation, subtype similarity;
- `leiden` — a label over cells;
- `velocity_graph` — a relation over `(cells, cells)`, directed, subtype transition.

Writing the result with the `seurat` profile places the shared-vocabulary fields faithfully:

- `counts` and `data` → assay layers;
- `pca` → a `DimReduc`, including its loadings (which a direct h5ad-to-Seurat conversion typically discards);
- `umap` → a `DimReduc`;
- `knn` → a graph;
- `leiden` → `Idents`.

The one field outside the shared vocabulary — `velocity_graph`, a directed `cells × cells` relation for which Seurat has no slot — is written to `@misc` and recorded in the `lstar/dropped` manifest. Two conclusions follow: every field in the common vocabulary survives the round trip, and the single field that cannot is reported rather than silently lost. Had the data remained in L★, the velocity graph would have been preserved as well.

---

# Part 5 — Standardization, versioning, and governance

This section concerns the conditions under which L★ delivers interoperability rather than merely relocating incompatibility. L★ deliberately shifts the standardization weight from the schema onto the profiles and the vocabulary. A thin, general model helps only if the conventions are shared: if two writers both emit valid L★ but one names normalized expression `data` with `state=lognorm` while the other names it `normalized` with `state=log1p`, or one stores clusters as a `label` named `leiden` while another stores them as a `measure` named `cluster_id`, conversion fails even though both files parse. The incompatibility is not removed but relocated from the schema to the vocabulary. For this reason the standardization and governance described below are central to the proposal, not ancillary.

## 5.1 What must be standardized

The guiding principle is to maximize what is *derivable*, so that less must be agreed. A field's role and named axes already carry most of its meaning (a `relation over (genes, genes)` is a gene network without further declaration). What still requires agreement:

1. **The core specification** — the model (Part 2) and the Zarr schema (Appendix A). Small, versioned, rarely changed.
2. **The role registry** — the canonical roles, their invariants, and controlled vocabularies for `state`, relation `subtype`, and `encoding`. Extensible through namespaced terms (`x-myorg:foo`), so that a new method does not require a specification change; a documented path promotes a widely used extension into the core.
3. **Per-profile naming** — within the `anndata` profile, that `X_umap` is an embedding named `umap`, `obsp['distances']` is a `similarity` relation, and so on. Most concrete disagreement lives here.
4. **Cross-profile vocabulary alignment** — the central requirement. For a conversion between two formats to be meaningful, the same concept must map to the same L★ signature from both sides: "normalized expression" must become `measure, state=lognorm` whether it originates in AnnData, Seurat, or pagoda2. This is achieved by a small **shared core vocabulary** into which all profiles map the common objects.

Items 1 and 2 should be kept small; the effort belongs in 3 and 4, since interoperability depends on them. We have compiled draft, best-guess mappings for the major formats — AnnData, Seurat, the pagoda-verse, and SingleCellExperiment/MultiAssayExperiment — covering most of the structures in common single-cell use; these define a starting shared vocabulary (illustrated below). The subsequent work is to converge these mappings over time, through the conformance suite (5.3) and an open process (5.4), rather than to mandate them at the outset.

A starter shared core vocabulary (v0, illustrative). Each common object maps to one canonical L★ signature that every profile targets:

- raw counts → `counts` : measure, `state=raw`
- normalized expression → `data` : measure, `state=lognorm`
- scaled expression → `scaled` : measure, `state=scaled`
- principal components → `pca` (embedding) with `pca_loadings` (loading), sharing the `pca` axis
- UMAP / t-SNE → `umap` / `tsne` : embedding
- neighbor graph → `knn` : relation, `subtype=similarity`
- clustering → a label named for the method (e.g. `leiden`) over `cells`
- cell-type annotation → `celltype` : label over `cells`
- markers → `markers.<grouping>` (`.lfc`, `.padj`, …) : measure over `(genes, group)`, with uncertainty
- cell / gene / sample metadata → columns of `cellMeta` / `geneMeta` / `sampleMeta`
- sample partition → `sample` : label over the union `cells`
- spatial position → `spatial` : embedding over an observed coordinate axis

An h5ad `X` (log-normalized), a Seurat `data` layer, and pagoda2's `getExpressionBlock()` all map to `data : measure, state=lognorm`; convergence means every profile agreeing on signatures of this kind.

## 5.2 Versioning

Three independently versioned components, each declared in every dataset and profile:

- **`spec_version`** (the core model and Zarr schema), semantic versioning; a dataset records the version under which it was written.
- the **role/vocabulary registry version**, semantic versioning; a new core role or controlled-vocabulary term increments it, while a namespaced extension does not.
- a **`profile_version`** per profile, declaring the `spec_version` and registry version it targets. A dataset's root records which profile wrote it (`profiles: ["anndata@0.3", …]`), so that a reader can interpret idiosyncrasies and know which conformance suite applied.

Compatibility rule: a reader supporting `spec_version` *M.x* must read any *M.y* dataset, treating unrecognized roles as opaque typed fields; a profile change that alters a mapping is breaking and increments the profile major version.

## 5.3 Conformance testing

Vocabulary agreements that are not tested erode. The enforcement mechanism is a shared **conformance suite**: a corpus of reference datasets, the round-trip tests of 4.2, and cross-profile coverage reports (for instance, `convert(anndata → seurat)` must preserve the shared-vocabulary core). A profile conforms at level L if and only if it passes the suite at level L. The suite, rather than prose, is the authoritative definition of the standard, as reference implementations were for AnnData and MuData.

## 5.4 Governance

The most difficult aspect, and the most likely point of failure, is social rather than technical. The arrangements that have succeeded for comparable standards are a small specification, reference implementations in the major languages, a conformance suite, and a neutral steward operating a lightweight proposal process with an explicit quality bar (the scverse model for AnnData, MuData, and SpatialData; the SOMA consortium; OBO and OME for ontologies and imaging). For L★ specifically:

- Provide the core specification, reference reader/writer implementations in Python (over Zarr) and R, and the conformance suite before requesting adoption; the reference implementations, not the text, constitute the standard.
- Implement the first three profiles (`anndata`, `seurat`, `pagoda2`) and demonstrate a lossless `anndata`↔`seurat` round trip on the shared-vocabulary core. This demonstration is the principal evidence of value; without it, adoption is unlikely.
- Place the role registry and the shared core vocabulary under a public, versioned proposal process with a small steering group, hosted by an existing neutral body where possible (scverse is a natural candidate, given that AnnData, MuData, and SpatialData are already maintained there) rather than by a single laboratory.
- Accept namespaced extensions freely; promote an extension into the core only when two independent implementations and suite coverage exist, which keeps the core small.

The honest assessment: if the community does not converge on items 3 and 4 of 5.1, L★ produces several dialects and offers no advantage over the present situation. The mitigations are those that have worked elsewhere — a small core, reference implementations, a conformance suite, and neutral governance — but they are necessary rather than sufficient, and they require adoption that L★ does not yet have.

## 5.5 Costs and limits

Stated plainly, so that the proposal can be evaluated rather than oversold:

- **Indirection.** Profile mapping occurs at input/output boundaries, not during computation; reading L★ on Zarr is comparable to reading AnnData-zarr, and a binding can expose native-speed array access. The abstraction is nonetheless additional software to build, maintain, and learn.
- **Dependence on the vocabulary.** Without convergence on 5.1(3)–(4), the lossless-interchange benefit does not materialize; the vocabulary above is a draft, not a ratified standard.
- **Not everything is naturally axis-and-field-shaped.** Continuous or unbounded index spaces (a parameter sweep, a continuous trajectory coordinate) fit only awkwardly as derived axes, and deeply nested structures rely on the collection-of-collections recursion rather than a flat representation. None of this is disqualifying, but a proposal that omitted it would be incomplete.

---

# Appendix A — Zarr serialization schema (normative)

A dataset is a Zarr group with three child groups — `axes/`, `fields/`, `models/` — and a root attribute block. All L★ metadata is placed under an `"lstar"` key in each group's `.zattrs`. Consolidated metadata (a single `.zmetadata`) is required for access over HTTP.

```
<root>/                      .zattrs: { "lstar": { "spec_version": "0.1",
                                                   "kind": "sample" | "collection",
                                                   "profiles": ["anndata@0.1", ...],   # written-by
                                                   "axes": [...], "fields": [...] } }   # index (also in .zmetadata)

  axes/<axis>/               .zattrs: { "lstar": { "kind": "axis",
                                                   "origin": "observed" | "derived",
                                                   "role": "observation" | "feature" | "coordinate" | null,  # profile hint
                                                   "induced_by": "<field>" | null,
                                                   "provenance": {...} } }
              labels         1-D array of element labels (str or int).  For a coordinate axis, the dimension names.

  fields/<field>/            .zattrs: { "lstar": { "kind": "field",
                                                   "role": "<role>",
                                                   "span": ["<axis>", ...],          # ordered
                                                   "state": "<state>" | null,        # measures
                                                   "encoding": "dense"|"csr"|"csc"|"coo"|"edge_list"|"ragged"|"raster"|"recipe",
                                                   "coverage": "full" | "partial",
                                                   "directed": bool, "weighted": bool, "subtype": "<...>",   # relations
                                                   "uncertainty": "<field>" | null,
                                                   "provenance": {...} } }
              # value arrays, by encoding:
              dense:      values                        # N-D array matching span
              csr/csc:    data, indices, indptr         # + lstar.shape
              coo/edge:   row|source, col|target, weight
              ragged:     values, offsets               # sequence
              raster:     a multiscale group (OME-NGFF) over the image base axes
              recipe:     no value arrays — lstar.recipe = { op, inputs:[<field>...], params }   # virtual; computed on read
              # partial coverage adds, per non-full axis:
              index/<axis>                              # the covered labels (a subset of axes/<axis>/labels)

  models/<model>/            .zattrs: { "lstar": { "kind": "model", "framework": "...",
                                                   "apply": { "in": ["<axis>",...], "out": ["<axis>",...] },
                                                   "weights": "<uri>" | inline, "provenance": {...} } }
```

A field's `span` references axes by name, so placement and reading are deterministic. The `state`, `coverage`, `uncertainty`, and relation flags are the only semantics a generic reader must honor; an unrecognized `role` or `subtype` degrades to a typed field over known axes. Collections add `samples/<id>/` sub-stores or `lstar.members` references (3.2).

---

# Appendix B — Profile mechanism and catalog (normative)

## B.1 The profile rule syntax

A profile is a versioned document declaring an identity (which native axes become the L★ observation and feature axes) and an ordered list of rules. Each rule has the form:

```
NATIVE_LOCATION   <dir>   field "NAME" : ROLE over (AXES) [attrs]      # notes
```

- `<dir>` is `<->` (bidirectional), `->` (read-only: native→L★, lossy on write), or `<-` (write-only).
- `NATIVE_LOCATION` is an exact path or slot pattern; `{k}` captures a key reused on the right.
- The right-hand side is the L★ field signature (2.2); helper functions are defined per profile: `byDtype(k)` (categorical→`label`, numeric→`measure`), `guess(k)` (state or subtype from a key, e.g. `counts`→raw, `distances`→similarity), `pair(k)` (match `varm[k]` to its `obsm`), and `dim_axis(k)` (the coordinate axis a key induces).

A profile must publish: the rule table; the helper definitions; a round-trip conformance result (4.2); and a `spec_version` and `profile_version` (Part 5). Rules are evaluated top to bottom and the first match applies; an unmatched native location falls to the profile's catch-all (typically `-> uns/misc`, recorded). This syntax is the contract a reader and writer implement; the mappings below are the specification, not paraphrase.

## B.2 profile `anndata` (scverse) — source: AnnData (.h5ad / .zarr)

```
identity:
  .obs index  := axis cells   origin=observed  role=observation
  .var index  := axis genes   origin=observed  role=feature
rules:
  .X                <->  field "X"        : measure  over (cells, genes)  state=(uns['lstar/state'] or "unknown")   # primary
  .layers[{k}]      <->  field "{k}"      : measure  over (cells, genes)  state=guess(k)
  .obs[{k}]         <->  field "{k}"      : byDtype(k) over (cells)        # one column of the cellMeta table
  .var[{k}]         <->  field "{k}"      : byDtype(k) over (genes)
  .obsm[{k}]        <->  field "{k}"      : embedding over (cells, dim_axis(k))     # X_pca→pca, X_umap→umap, spatial→observed
  .varm[{k}]        <->  field "{k}"      : loading   over (genes, dim_axis(pair(k)))
  .obsp[{k}]        <->  field "{k}"      : relation  over (cells, cells)  subtype=guess(k)
  .varp[{k}]        <->  field "{k}"      : relation  over (genes, genes)
  .uns[{k}]         ->   field "{k}"      : opaque                          # read-only catch-all (lossy on write)
write-back:
  any L★ field with no matching rule  <-  .uns["lstar/<name>"]   as a structured blob (recorded)
extensions:
  a 2nd feature axis                  <->  profile mudata  (mod[name])
  spatial points/shapes/images        <->  profile spatialdata
```

Fit: the single-sample core round-trips losslessly. Loss: group-axis results, arity-3 tensors, trees, and models survive only as `uns` entries.

## B.3 profile `seurat` — source: Seurat v5 object

Seurat maps onto L★ more closely than AnnData in two respects: a `DimReduc` natively bundles an embedding with its loadings (L★'s shared coordinate axis), and `@commands` is a provenance log.

```
identity:
  colnames(obj)            := axis cells   origin=observed  role=observation
  rownames(assay)          := axis <assay> origin=observed  role=feature        # one feature axis per assay
rules:
  assay$layers[counts|data|scale.data]  <->  field : measure over (cells,<assay>)  state=(raw|lognorm|scaled)
  obj@meta.data[{k}]                    <->  field "{k}" : byDtype(k) over (cells)
  Idents(obj)                           <->  field "ident" : label over (cells)
  assay@meta.features[{k}]              <->  field "{k}" : byDtype(k) over (<assay>)
  obj@reductions[{k}]  (DimReduc)       <->  { field "{k}"          : embedding over (cells, dim_axis(k)),
                                               field "{k}_loadings"  : loading   over (<assay>, dim_axis(k)),
                                               field "{k}_stdev"     : measure   over (dim_axis(k)) }
  obj@graphs[{k}] | obj@neighbors[{k}]  <->  field "{k}" : relation over (cells, cells)  subtype=guess(k)
  a 2nd assay (e.g. ADT)                <->  a 2nd feature axis (multimodal)
  obj@images[{fov}]  (Centroids/Segmentation/Molecules)  <->  geometry over the fov (spatialdata-style)
  obj@commands                          ->   provenance on the produced fields
  obj@misc[{k}]                         ->   field "{k}" : opaque    # catch-all (lossy)
```

Fit: reductions and provenance map exactly. Loss: group-axis results and trees fall to `Misc`.

## B.4 profile `pagoda2` — source: Pagoda2 R6 object

```
identity:
  rownames           := axis cells   origin=observed  role=observation
  colnames           := axis genes   origin=observed  role=feature
rules:
  $getRawCounts()                <->  field "counts" : measure over (cells,genes) state=raw   # canonical store
  $getExpressionBlock()          <->  field "data"   : measure over (cells,genes) state=lognorm  encoding=recipe  # virtual
  $cellMeta[{k}]                 <->  field "{k}" : byDtype(k) over (cells)
  $geneMeta[{k}]                 <->  field "{k}" : byDtype(k) over (genes)        # e.g. dispersion m, qv
  $reductions$PCA                <->  field "pca" : embedding over (cells, dim_axis(pca))
  $reductions$PCA (rotation)     <->  field "pca_loadings" : loading over (genes, dim_axis(pca))
  $embeddings[{k}]               <->  field "{k}" : embedding over (cells, dim_axis(k))
  $graphs[{k}]                   <->  field "{k}" : relation over (cells,cells) subtype=guess(k)
  $clusters/$clusterings[{k}]    <->  field "{k}" : label over (cells)            # + defaultGrouping
  $markerResults[{g}]            <->  field "markers.{g}" : measure over (genes, <g>)   # on the induced group axis
  $history                       ->   provenance
```

Fit: pagoda2.1 already exposes most of this; adopting the profile is largely renaming.

## B.5 profile `conos` (collection) — source: Conos R6 object

The reference collection: a set of heterogeneous samples linked by a joint graph, not a concatenated matrix.

```
identity:
  names(con$samples)             := axis samples   origin=observed
  union of member cells          := axis cells      origin=derived (union; ids = (sample, cell))
rules:
  con$samples[{id}]              <->  member sample  (read with its own sample profile: pagoda2 | seurat)
  con$getDatasetPerCell()        <->  field "sample" : label over (cells)        # the partition
  con$graph                      <->  field "knn"    : relation over (cells, cells)  joint
  con$embedding                  <->  field "umap"   : embedding over (cells, umap)  joint
  con$clusters[{t}]$groups       <->  field "{t}"    : label over (cells)            joint
  con$pairs / alignment results  <->  field "alignments" : relation over (cellsA, cellsB)  correspondence
  getJointCountMatrix()/getGeneExpression()  <->  gather (no concatenated matrix; partial across samples)
  (sample-level covariates)      <-   field over (samples)   # not native to conos; supplied as a `samples` table
```

A worked instance — four samples, two species:

```
Collection "BM_integration"   (2 human, 2 mouse)
  members (each a sample, namespaced, with its OWN axes)
    human_1   cells 9,102 × genes 19,400
    human_2   cells 8,517 × genes 19,400
    mouse_1   cells 7,330 × genes 21,800     (a DIFFERENT gene axis)
    mouse_2   cells 6,985 × genes 21,800
  axes
    samples   4
    cells     union = 31,934   ids = (sample, cell)   [derived]
    joint     12 clusters                             [derived]
  fields
    sampleMeta   table     over samples     { species, donor, tissue }
    inclusion    relation  cells → sample.cells
    alignments   relation  cellsA × cellsB              conos pairwise anchors
    orthology    relation  genes(human) × genes(mouse)  cross-species feature map
    knn          relation  cells × cells                joint graph
    umap         embedding cells × umap                 joint layout
    joint        label     cells → joint                clusters spanning all samples
```

Fit: this is conos restated as axes and fields. No concatenated matrix exists; `featureVector("CD8A")` gathers per sample (human cells return CD8A, mouse cells the orthologous Cd8a, or NA). The only additions made first-class are the `samples` table and the `orthology` map.

## B.6 profile `cacoa` (collection with a design) — source: Cacoa R6 object (`dev_lm`, v0.5.0)

Case-control analysis is a collection with a design. cacoa `dev_lm` generalized the binary case/control comparison to a full linear model.

```
identity:
  cao$data.object                := an underlying collection (a conos / sample list; read with profile conos)
  cao$sample.meta (rows)         := axis samples
rules:
  cao$sample.meta[{k}]           <->  field "{k}" : byDtype(k) over (samples)
  { cao$formula, cao$contrast,
    cao$model.matrices, cao$block.vars }  <->  field "design" : design over (samples)
                                               # induces axis coefficient;
                                               #   model_matrix : measure over (samples, coefficient)
                                               #   core/nuisance: label over (coefficient)
                                               #   blocks       : label over (samples)
  cao$cell.groups                <->  field "celltype" : label over (cells)
  cao$sample.per.cell            <->  field "sample" : label over (cells)
  cao$test.results[composition]  <->  field "composition" : measure over (celltype, samples)  + uncertainty (pval, z)
  cao$test.results[expr.shift]   <->  field "expr_shift"  : measure over (genes, celltype)     + uncertainty
  cluster-free per-cell shifts   <->  fields over (cells)
  cao$embedding                  <->  field "umap" : embedding over (cells, umap)
```

A worked instance — twelve samples, two conditions, donor-paired, adjusting for sex and age:

```
Collection "MS_cohort"
  axes
    samples      12
    cells        union of all sample cells     [derived]
    celltype     14                            [derived: joint annotation]
    coefficient  { conditionMS, sexM, age }    [derived: from the design]
  fields
    sampleMeta     table      over samples     { condition, sex, age, donor, batch }
    design         design     over samples     formula ~ condition + sex + age ;  contrast (condition: MS vs ctrl)
                     model_matrix   measure  samples × coefficient
                     core/nuisance  label    over coefficient
                     blocks         label    over samples            (permute within donor)
    annotation     label      cells → celltype
    composition    measure    celltype × samples    + uncertainty (stat / pval / z)
    expr_shift     measure    genes × celltype       + uncertainty
```

Fit: the entire linear-model design and the permutation results map onto a `design` (inducing a coefficient axis) and result fields, demonstrating that the `design` role is necessary.

## B.7 Readers, writers, converters, conformance

See Part 4.2 for the `read_P` / `write_P` / `convert(P → Q)` contracts and the two round-trip conformance tests every profile must publish.

---

# Appendix C — Design rationale

The model was reached in three passes. **Collapse:** beginning from approximately ten structure types (matrix, graph, tree, embedding, model, geometry, design, vector, sequence, correspondence), one observes that they overlap — a graph is a sparse matrix over (axis, axis); an embedding is a matrix over (axis, a coordinate axis); loadings are the same axis from the feature side; a link is a matrix over two cross-sample axes; a tree is a node axis with a parent relation and a leaf map — and so they reduce to {axis, field} with a role tag. **Test and extend:** applied to the single-cell and spatial scenario survey, derived/latent axes dissolve loadings, factors, and programs; group-axis induction dissolves markers and communication; union axes preserve collection heterogeneity; arity-3 fields handle communication and eQTL. The additions required are the roles `vector`, `sequence`, and `transform`, an optional `uncertainty` companion, label-keyed coverage, and provenance references. **Trim:** recognizing that sample, collection, and AnnData are profiles rather than primitives; keeping the role vocabulary open; observing that the provenance graph follows without additional construction. Rejected alternatives: a fixed enumerated taxonomy of structures (too rigid for the long tail); a single globally aligned tensor (the over-normalization the proposal rejects); and modeling fitted models or designs out of band (they fit as `transform` and `design` fields).

---

# Appendix D — Open questions

- **The core/extension boundary** — which roles and capabilities are core versus extension tiers (spatial geometry, sequences, genetics), advertised through a capability list.
- **Arity-n encoding** — sparse storage and access for arity-3 fields (communication, eQTL).
- **`embedding` versus a dedicated `coordinate` role** — physical coordinates are currently an `embedding` over an observed coordinate axis; whether to introduce a separate `coordinate` role (with metric operations and frame transforms) is open.
- **The shared core vocabulary (5.1)** — a v0 starter is drafted; ratifying and testing its membership, through the conformance suite, is the single most consequential remaining task before wide adoption.
- **Naming** — `L★`, *axis*/*field*, the role names, and the file extension (`*.lstar.zarr`) are provisional.
