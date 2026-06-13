# Format coverage — gaps & recommendations (Seurat / AnnData / spatial)

*Draft for review. Synthesizes a research pass over the Seurat (satijalab.org/seurat, SeuratObject
reference, Signac), scverse (anndata IO spec, sc-best-practices, MuData), and spatial (SpatialData,
squidpy, OME-NGFF, 10x Visium/Xenium, Vizgen, CosMx) docs, grounded against lstar's current profiles.*

## Where we are

The AnnData and Seurat profiles cover the **shared core** well — expression matrices/layers → measures,
obs/var → label/measure fields, obsm → embeddings, varm → loadings, obsp/varp → relations, DimReducs →
embedding+loading, graphs → relations, v3/v4-vs-v5 assay recognition, split-v5 → collection. Conos/SCE/
pagoda2 round-trip their cores; Conos correctly keeps per-sample heterogeneity.

Three things are **systematically dropped or absent** today:
1. **Everything in `uns` / `@misc` / `@commands`** — recorded as a name in `dropped`, *not even losslessly
   preserved*. This is where most analysis results live (DE, neighbor params, colors, PCA variance, …).
2. **Multimodal** (MuData, Seurat multi-assay) — only one feature (`genes`) axis is modeled.
3. **Spatial** beyond `obsm['spatial']` read as a plain 2D embedding — no images, polygons, molecules,
   or coordinate frames.

And several **encodings/roles are spec'd but unimplemented**: `ragged`, `raster`, `recipe` encodings;
`sequence`, `transform` roles; partial-axis coverage; uncertainty companion fields. Spatial is the
forcing function for most of these.

A baseline correctness issue cuts across everything: **dtype fidelity** — categoricals are coerced to
utf8 (losing `ordered` and the `-1` missing sentinel), and pandas nullable Int/boolean/string extension
dtypes (values+mask) aren't handled. These are silent round-trip corruptions, independent of any new
feature.

---

## Tier 1 — baseline fidelity + cheap structured results (do first; low risk, high value)

- **Categorical fidelity.** Preserve `ordered` and the `-1`/code-missing sentinel exactly; model a label
  as codes-into-a-categories-axis + an `ordered` flag. *The single most common silent corruption.*
- **Nullable / extension dtypes** (`nullable-integer/boolean/string` = values+mask). Carry an explicit
  null-mask rather than coercing to float-NaN (also keeps integer-ness — aligns with memory-lean dtype).
- **Color palettes** (`uns['<key>_colors']`, Seurat per-ident colors). Capture as a small field bound to
  the *category order* of the matching categorical; store the binding so reordering re-permutes colors.
- **Small structured `uns` param dicts** → typed records instead of dropped: `log1p.base` (gates
  re-normalization — semantically load-bearing), `pca.variance`/`variance_ratio` (a measure over the
  PCA-component axis — same axis as `X_pca`/`PCs`), `neighbors.params`/`umap`/`leiden` params.
- **Seurat: multiple stored idents.** Map *all* factor meta.data columns to labels, flag which is active
  (`active.ident`) — the active-vs-stored distinction is currently lost.
- **Seurat: `Neighbor` objects** (nn.idx/nn.dist KNN edge lists) → relations (index+distance), distinct
  from the symmetric `Graph`; plus DimReduc `stdev` (measure over dim axis) and `feature.loadings.projected`.
- **RNA-velocity layers** (`spliced`/`unspliced`/`Ms`/`Mu`/`velocity`) → measures (free given layer
  handling); `velocity_graph`/`_neg` → relations; `var['fit_*']` → measures over genes.
- **Make the `uns`/`@misc` long-tail passthrough _lossless_.** Upgrade `dropped` from "record the name"
  to "preserve the nested dict-of-elements verbatim and reproduce it on write." This is foundational:
  it makes round-trips safe *before* we type anything, and lets us promote recognized structures
  (Tiers 1–2) out of the tail incrementally.

## Tier 2 — analysis results need a derived "group/cluster" axis

The flagship gap. Build **one derived group axis** (clusters/cell-types) and reuse it:
- **Differential expression** — AnnData `uns['rank_genes_groups']` (per-group structured arrays of
  names/scores/lfc/pvals/pvals_adj) and Seurat marker tables. Honest shape is **ragged**: each group
  ranks a possibly-different subset/order of genes (`names` is a per-group permutation). Model as a
  group axis + measures over (group, gene) where dense, or ragged per-group gene-rankings where not.
  **Recommend exposing both a queryable typed view *and* lossless passthrough.** (Matches the
  collection-not-tensor instinct — don't force a dense group×gene tensor.)
- **PAGA / trajectory** — `uns['paga'].connectivities`(+tree) = relations over that same group axis;
  `obsm['X_diffmap']` embedding + `uns['diffmap_evals']` + `obs['dpt_pseudotime']` + `uns['iroot']`.
- **Dendrograms** (`uns['dendrogram_*']`) — linkage + ordering over the group axis.
- **Provenance** — Seurat `@commands` (name, timestamp, assay.used, params per run) and AnnData param
  dicts → `transform`/provenance records, so analysis history survives a round-trip.

## Tier 3 — multimodal (MuData / Seurat multi-assay)

This *is* lstar's "collection, not tensor" worldview applied to modalities:
- Each modality → its own feature axis (genes / proteins(ADT) / peaks(ATAC)), mapped by the existing
  per-format profile; a **shared cell axis** the modalities map into.
- MuData `obsmap`/`varmap` (1-based, 0=absent) → **alignment relations** (modality-local ↔ global index),
  preserving partial overlap rather than intersecting. Treat per-modality + maps as ground truth, the
  global `obs`/`var` as a derived cache (muon `.update()` semantics); preserve the `modality:` prefix.
- MuData-level `obsm` (WNN/totalVI/MOFA) → embeddings over the shared cell axis.
- Seurat CITE-seq (RNA+ADT assays) and 10x multiome (RNA + ChromatinAssay) are the same pattern.
- *Larger lift, real demand.* The `samples`/collection machinery already models most of it.

## Tier 4 — spatial (needs new encodings + a coordinate-frame concept)

**Feasible now (no new machinery):**
- Positions (spots/beads/centroids, `obsm['spatial']`) → `embedding` over an **observed coordinate axis**
  (`physical` = x,y[,z]); the only real change is *recording the unit (micron/pixel) and frame name*
  instead of an anonymous 2D embedding. Visium spots / Xenium centroids / Slide-seq beads all fit.
- Spatial neighborhood graphs → `relation(cells,cells)`, `subtype=spatial` (already in the vocabulary).
- Multi-FOV / multi-sample → **collections** (per-FOV namespaced axes + transforms to a global frame).

**Needs a spec'd-but-unimplemented encoding built:**
- **Molecules/transcripts** → a new observed **`molecules` axis** + ordinary fields: `position`
  (embedding over molecules×physical), `gene` (label, **dictionary-encoded** — millions of rows),
  `assignment` (relation molecules→cells; `overlaps_nucleus` flag), `qv` (measure). Points are
  *rectangular* (not ragged) but need **chunked/columnar (Parquet-like)** storage at 10⁷–10⁸ rows.
- **Cell/nucleus boundary polygons** → **`ragged`** field + a dedicated **`geometry` role** (offsets
  segment a flat vertex array by cell; closed rings). *Design fork:* pure-ragged stays light but can't
  express holes/multipolygons; embedding GeoParquet is faithful but a heavy geo dependency (light-deps
  tension). Visium spots / approximated cells are just centroid + `radius` (no geometry encoding needed).
- **Images** (H&E, DAPI, IF) → the spec'd **`raster`** encoding built on **OME-NGFF** (own y/x/channel
  base axes; pass the `multiscales` block through verbatim — adopt OME-NGFF metadata, don't reinvent).
  Support **referencing** an external OME-NGFF/OME-TIFF by URI as well as embedding (Xenium/Vizgen images
  are huge). Make raster an **optional/Suggests-tier** capability with graceful degradation.
- **Segmentation label masks** → a `raster` whose integer values cross-reference the `cells` axis.

**The one genuinely missing *concept*: named coordinate systems + frame-to-frame transforms.**
SpatialData/OME-NGFF make every element map to ≥1 named coordinate system (default `"global"`) with
composable transforms (Identity/Translation/Scale/Affine/Rotation/Sequence). lstar has no way to say
"this embedding and this image are registered to the same physical space." Today's `transform` role
means a *fitted model* (scVI/UMAP projector), not an *affine frame map* — different in kind. Without
this, you can store positions and images but not that they're co-registered — which is what makes
spatial "spatial." **Recommendation: adopt the SpatialData/OME-NGFF coordinate-system + transform model
rather than reinvent it** (a `frame` attribute on coordinate axes/embeddings + `transform` fields with
`subtype=coordinate` carrying the affine). Visium `scalefactors` → a Scale; CosMx `fov_positions` →
per-FOV Translations.

---

## Cross-cutting design decisions (the review agenda)

1. **Lossless passthrough first.** Upgrade `uns`/`@misc` handling from name-only `dropped` to verbatim,
   reversible preservation — independent of typing anything. Biggest safety win per unit effort.
2. **One derived group/cluster axis**, reused by DE, PAGA, markers, dendrograms. Build it once.
3. **DE results: typed/queryable view _and_ lossless passthrough** — and accept the ragged shape rather
   than forcing group×gene dense.
4. **dtype fidelity (categorical `ordered`/missing, nullable) is a baseline fix**, not a feature — do it
   regardless of the rest.
5. **Spatial coordinate frames: adopt SpatialData/OME-NGFF** (named coordinate systems + composable
   transforms) — the one new *concept*. Decide scope: full adoption vs a lighter `frame` attribute.
6. **Geometry encoding: ragged-only (light) vs GeoParquet (faithful/heavy)** — the central light-deps
   tension for polygons. Likely: ragged rings core + opaque GeoParquet passthrough for full fidelity.
7. **Raster = optional capability tier** (Suggests + graceful absence), reference-by-URI by default for
   large images; don't transcode OME-TIFF↔NGFF.
8. **Implement `ragged` + `raster` encodings and the `sequence`/`transform` roles** — spatial (and TCR/
   BCR sequences) is the forcing function; they're already in the spec.
9. **Profile granularity:** one `spatialdata` profile vs per-vendor importers. Graceful version detection
   argues for *some* per-vendor (e.g. Vizgen boundaries went HDF5→Parquet at v232; VisiumV1→V2 moved
   coordinates into FOV) — same principle as the Assay/Assay5 and modern/legacy-sparse handling.
10. **Multimodal scope:** is MuData/multi-assay in-scope now, or after the `uns`/group-axis/dtype work?
    It's a clean fit to the collection model but a larger lift.

---

## Condensed "where each artifact lives" catalog

**AnnData / scverse.** On-disk: every element carries `encoding-type`/`encoding-version` (anndata,
array, csr/csc_matrix, dataframe, categorical, string-array, nullable-{integer,boolean,string},
numeric-scalar, dict, awkward-array). `uns` is a recursive `dict`: neighbors(params + obsp pointers),
pca(variance/variance_ratio), rank_genes_groups(per-group structured arrays + params), `*_colors`,
dendrogram, paga, umap/leiden/log1p/hvg params. Velocity → layers + obsp graphs + var fit params.
MuData(.h5mu): `/mod/<name>` = full nested anndata + obsmap/varmap integer alignment + `mod-order`.

**Seurat.** `assays` (Assay v3/4 counts/data/scale.data vs Assay5 named layers; SCTAssay +
SCTModel.list per-gene NB params; ChromatinAssay peaks + GRanges/fragments/motifs/links/annotation);
`reductions` (DimReduc: cell.embeddings/feature.loadings/projected/stdev/jackstraw); `graphs`(dgCMatrix)
+ `neighbors`(nn.idx/dist); `meta.data`(cell) + assay `meta.features`(gene); `active.ident` + stored
idents; `images` (VisiumV1/V2; FOV → Centroids/Segmentation polygons/Molecules; per-FOV frames);
`commands`(provenance), `misc`/`tools`; v5 BPCells on-disk layers.

**Spatial primitives → frameworks.** points(SpatialData Points/Parquet · obsm['spatial'] · Centroids ·
Xenium cells.parquet); rasters(SpatialData Images xarray/OME-NGFF multiscale · uns['spatial'] hires/
lowres · VisiumV1 image · Xenium morphology.ome.tif); labels/masks(SpatialData Labels/OME-NGFF);
shapes/polygons(SpatialData Shapes/GeoParquet · Seurat Segmentation · Xenium *_boundaries.parquet);
molecules(SpatialData Points · Seurat Molecules · Xenium transcripts.parquet); transforms(SpatialData/
OME-NGFF coordinateTransformations · Visium scalefactors · CosMx fov_positions); coordinate
systems(SpatialData named extrinsic frames — *no lstar analogue today*).
