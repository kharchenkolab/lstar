# L★ — format & language support

What lstar can convert, read, and write today — by **format** (AnnData, MuData, Seurat,
SingleCellExperiment, Conos/pagoda2) and by **language** (Python, C++, R, JS/WASM) — with, for each
case, whether it is exercised by a **real** dataset and by the **synthetic CI** fixture that stands in
for it. This is the place to check "does lstar handle *my* object?" and the place we track where the
coverage still has holes.

> Status: early development. The matrices below reflect what round-trips and validates in the test
> suite, not a frozen API.

**Legend** — `✓` supported · `◐` partial (works, but a sub-case is *recorded as a known loss* rather
than fully typed — never silently dropped) · `✗` not yet · `—` not applicable.

## How testing works (so the columns mean something)

lstar is tested in **two tiers**, and a gap can live in either one — that distinction is tracked here on
purpose:

- **Real corpus (local).** The profiles run against *real* published objects — downloaded small
  datasets, large local atlases, and a breadth **sweep** over whole repositories (Bioconductor
  `scRNAseq` ≈ 61 SCEs, SeuratData, scanpy.datasets, real Conos objects). Real data is what catches the
  long tail; the sweep has caught **5 real profile bugs no fabricated fixture had** (SCE NULL dimnames,
  S4 `Rle` columns, SE-not-SCE accessors, and — most recently — silent loss of Signac `ChromatinAssay`
  peak ranges). See [`conformance/sweep/REPORT.md`](conformance/sweep/REPORT.md) and
  [`python/tests/CORPUS.md`](python/tests/CORPUS.md).
- **Synthetic CI (github).** CI runs **synthetic-only** — no real datasets are committed or downloaded
  (they're large/slow). The synthetic fixtures aren't hand-fabricated: they push **synthetic counts
  through the real scanpy/mudata/Seurat pipelines**, so the *same library code* produces the real
  structures (categoricals, `*_colors`, `uns['pca']`, the neighbors `OverloadedDict`,
  `rank_genes_groups` dtypes, RNA+ADT modalities, `velocity_graph`). Keeping these faithful to the real
  corpus is an explicit contract; the local real runs verify it.

So in the per-format tables, **Real** = a real object exercises the case; **CI** = the synthetic fixture
does. A case can therefore have three kinds of gap: the **profile** doesn't handle it (`✗`/`◐`), there's
**no real example** in the corpus, or the **synthetic fixture doesn't represent it** yet.

---

## Language / interface support

One C++ core underlies Python, R, and (via WebAssembly) the browser; R adds the format profiles. The
portable Zarr store is the interchange — every language reads and writes the *same* store, verified by
the cross-language conformance suite ([`conformance/`](conformance/)). Stores are **Zarr v3 by default**
(a per-node `zarr.json` + inline consolidated metadata), with legacy **v2** (`.zarray`/`.zgroup` +
a consolidated `.zmetadata`) available on request; all four surfaces read and write both formats.

| capability | Python | C++ (`libstar`) | R | JS / WASM |
|---|:---:|:---:|:---:|:---:|
| Read a store | ✓ | ✓ | ✓ | ✓ |
| Write a store | ✓ | ✓ | ✓ | ✓ |
| Zarr format **v2 + v3** (v3 default) | ✓ | ✓ | ✓ | ✓ |
| **v3 sharding** (`sharding_indexed`; byte-range-readable) | ✓ | ✓ | ✓ | ✓ |
| `validate()` (model/consistency checks) | ✓ | — | via Py/C++ | — |
| Encodings: dense / CSR / CSC / COO | ✓ | ✓ | ✓ | ✓ (read; write dense/CSC/CSR) |
| UTF-8 + categorical (codes/levels/ordered/`-1` missing) | ✓ | ✓ | ✓ | ✓ |
| Nullable validity mask (Int/bool/string) | ✓ | ✓ | ✓ | ✓ |
| Factor-axis **induction** (`induced_by` round-trip) | ✓ | ✓ | ✓ | ✓ |
| **Partial coverage** (a field on a subset of a span axis, via an `index`) | ✓ | ✓ | ✓ | ✓ |
| **Arity-3 fields** (CCC `sender×receiver×lr_pair`, eQTL `celltype×gene×variant`) | ✓ | ✓ | ✓ | — |
| Chunking + compression (**gzip / zstd**) | ✓ | ✓ (zstd needs libzstd) | ✓ (zstd needs libzstd) | ✓ |
| Per-field write layout (per-array codec / chunk / shard) | ✓ | ✓ | ✓ | ✓ |
| Byte-range read of a **compressed** array (chunk-granular; no whole-array decode) | ✓ | ✓ | ✓ | ✓ |
| Lazy / partial / over-network reads | ✓ | ✓ | ✓ | ✓ |
| Bounded-memory **streaming write** | ✓ | ✓ | ✓ | — |
| Disk-backed targets (don't materialize) | backed AnnData | — | BPCells / HDF5Array | — |
| Kernels: col mean/var · fused depth+log1p · grouped sum | ✓ | ✓ | ✓ | ✓ (colMeanVar, log1p) |
| **Deterministic reductions** (bit-identical across thread counts) | ✓ | ✓ | ✓ | — |
| Block reader (read a gene/feature subset) | ✓ | ✓ | ✓ | — |
| `provenance` round-trip (method params / normalization recipe) | ✓ | ✓ | ✓ | — |
| DE / markers / pseudobulk bundles · **collection** pseudobulk | ✓ | — | ✓ | — |
| **Format profiles** | AnnData, MuData | — | Seurat, SCE, Conos, pagoda2 | — |

Notes: the C++ core is the model + Zarr IO + kernels (no `validate`/profiles — those live in the
language packages). The JS/WASM surface **is the libstar C++ core compiled to WebAssembly** (the same
reader R/Python/C++ use — no separate JS Zarr reimplementation), so it reads v2 **and** v3 through one
recipe. It's the **viewer** data layer: read **and write** the store — every encoding both directions
(CSC/dense/categorical+factor/mask/partial/aux), **v2 or v3, chunked, sharded, gzip- or zstd-compressed**
(encode via a WASM writer module built from the same core) — and run the view kernels in-browser; it has
no format profiles or DE bundles. The writer also `addToStore`s derived fields onto an existing (e.g.
Python-written) store. A **compressed** array stays byte-range-readable: the reader decodes only the
covering chunk(s), so compression and sub-chunk reads coexist (the viewer's `viewer@0.1` prep uses this —
its stores are compressed per field by default: gene-major counts raw for instant gene-coloring,
`counts_cellmajor` zstd chunked+sharded, other navigators zstd single-chunk). Cross-language round-trips
covered by `conformance/{categorical,induce,nullable,partial,arity3,aux,provenance,collection,chunked,
read_block,stream_reduce,fused_view,js}.sh`. These test each binding as a **producer**, not just a
consumer: every rich encoding is *authored from scratch* in Python, R, and (for the core encodings) JS,
then **cross-read by a different language** asserting native type + exact values — see
[`conformance/README.md`](conformance/README.md) ("Origin coverage"). (The JS-write → Python/C++-read leg
is in `js.sh`.)

---

## AnnData (`.h5ad`) — Python

**Versions:** legacy < 0.7 (`h5sparse_format` on-disk) · 0.7+ (`encoding-type`) · 0.8 · ≥0.10 (incl. the
`CSRDataset`/`CSCDataset` backed-class rename) — all ✓. The legacy vs modern on-disk sparse layout is
distilled into a synthetic test (`test_legacy_format`, written with h5py — no old anndata needed), and CI
runs the core suite on **old anndata 0.8** (the `pinned-old` matrix leg) as well as latest. Real
grounding: `pbmc68k_reduced`, `pbmc3k_processed`, a local Tabula-Muris-Senis Marrow atlas (40 k × 20 k,
read backed), real scVelo output.

| case | status | real | CI (synth) | notes |
|---|:---:|:---:|:---:|---|
| `X` / `.raw` (divergent gene set) / `layers` | ✓ | ✓ | ✓ | pbmc68k `.raw`; pancreas spliced/unspliced |
| encodings dense / CSR / CSC / UTF-8 | ✓ | ✓ | ✓ | |
| categorical obs/var (ordered, `-1` missing) | ✓ | ✓ | ✓ | real louvain/phase/bulk_labels |
| nullable `Int64` / `boolean` / `string` (+mask) | ✓ | ✓ | ✓ | derived from real values (no public h5ad ships them) |
| obs/var **same-name** collision (`n_counts`) | ✓ | ✓ | — | disambiguated, neither lost |
| `obsm`/`varm` (pca/umap/loadings, incl. NaN) | ✓ | ✓ | ✓ | |
| `obsp`/`varp` graphs → relations | ✓ | ✓ | ✓ | distances/connectivities |
| `uns` params / `*_colors` / pca variance → typed | ✓ | ✓ | ✓ | promoted out of the passthrough, bound to axes |
| `uns['rank_genes_groups']` (t-test/wilcoxon/logreg/pairwise) | ✓ | ✓ | ✓ | one-vs-rest typed; pairwise kept verbatim |
| `uns['velocity_graph']` (scVelo, cell×cell in `uns`) | ✓ | ✓ | ✓ | typed as a relation (not dropped) |
| `uns` neighbors `OverloadedDict` / generic tail | ✓ | ✓ | ✓ | lossless passthrough |
| backed / bounded-memory convert | ✓ | ✓ | ◐ | real Marrow atlas (local); CI exercises the backed proxy |
| spatial coords `obsm['spatial']` → named `spatial` coordinate axis (conceptual) | ✓ | — | ✓ | observed coordinate axis; round-trips to `obsm['spatial']` |
| spatial images `uns['spatial']` / vendor frames / molecules | ✗ | — | — | deferred to a spatial tier; kept in the passthrough (not lost) |

## MuData (`.h5mu`) — Python (multimodal)

Modalities map to **canonical feature axes over one shared `cells` axis** (rna→`genes`, prot/adt→
`proteins`, atac→`peaks`) — the same shape as a Seurat multi-assay object. Real grounding: `minipbcite`
CITE-seq (downloaded); CI: synthetic RNA+ADT through the real mudata/scanpy pipeline.

| case | status | real | CI (synth) | notes |
|---|:---:|:---:|:---:|---|
| modalities → feature axes (genes/proteins/peaks) | ✓ | ✓ | ✓ | RNA+ADT |
| per-modality `X` / `layers` / `var` / `obsm` / `uns` | ✓ | ✓ | ✓ | own PCA/UMAP per modality |
| global `obs` categoricals → factor axes (`<mod>:` prefix) | ✓ | ✓ | ✓ | celltype/leiden |
| global `obsm` (WNN / MOFA joint embedding) | ✓ | ✓ | ✓ | |
| **joint-method shapes** — WNN modality weights (cell measures) + joint graph (`obsp`→relation) | ✓ | ✓ | ✓ | we type the *output shape*, not the algorithm |
| **MOFA/totalVI** — shared factor axis: scores (`obsm`) + per-modality loadings (`mod.varm`) | ✓ | ✓ | ✓ | one factor axis carries scores embedding + per-mod loadings (lstar induction) |
| `obsmap`/`varmap` aligned cells | ✓ | ✓ | ✓ | |
| **partial-overlap** modalities (a modality on a cell subset) | ✓ | — | ✓ | **typed partial coverage**: an `index` into the shared `cells` axis (not a `cells.<mod>` axis, not padded); round-trips Py↔C++↔R + back to MuData on the subset |
| RNA+ATAC multiome `.h5mu` | ◐ | — | — | path = RNA+ADT + partial coverage; no real `.h5mu` multiome sourced (real ATAC is covered via Seurat pbmcMultiome) |

## Seurat — R

**Classes/versions** (a "Seurat object" is several structurally-different classes; each is recorded in
`ds$profiles`): `Assay` (v3/v4) ✓ · `Assay5` (v5) ✓ · v5 **split**/integration → read as a *collection*,
re-split on write ✓ · `SCTAssay` ◐ · `ChromatinAssay` (Signac scATAC) ✓. Real grounding: a **10/10
PASS** SeuratData sweep (RNA, RNA+ADT, 4-modality ECCITE-seq, integration, HVG-subset loadings) + real
`pbmc3k.final`, `cbmc`, `pbmcMultiome`.

| case | status | real | CI (synth) | notes |
|---|:---:|:---:|:---:|---|
| `Assay` (v3/v4) | ✓ | ✓ | ✓ | real pbmc3k.final |
| `Assay5` (v5) + layers | ✓ | ✓ | ✓ | |
| v5 **split** (per-sample layers → collection) | ✓ | ✓ | ✓ | write-back re-splits |
| multimodal RNA+ADT / RNA+ATAC (other assays → feature axes) | ✓ | ✓ | ✓ | real cbmc, bmcite, thp1.eccite (RNA+ADT+HTO+GDO) |
| `ChromatinAssay` peaks + **genomic ranges** (seqnames/start/end) | ✓ | ✓ | — | real pbmcMultiome (108 k peaks); fragments **recorded** |
| `SCTAssay` (residuals + `SCTModel`) | ◐ | ✓ | ◐ | data typed; SCTModel recorded. CI runs SCT only when full Seurat loads |
| reductions: embeddings / loadings / **stdev** | ✓ | ✓ | ✓ | |
| loadings over **HVG subset** → subset feature axis | ✓ | ✓ | — | real pbmc3k.final (2000/13714) round-trips exactly; was dropped |
| `scale.data` over **HVG subset** → partial-coverage measure | ✓ | ✓ | — | typed as a partial measure over (cells, genes) keyed by a gene index; round-trips exactly (uses the new `index` machinery) |
| graphs (dgCMatrix) → relations | ✓ | ✓ | — | |
| `Neighbor` (nn.idx/dist) → weighted cell-cell relation | ✓ | — | ✓ | distance-weighted; was dropped |
| `meta.data` factors / active `Idents` | ✓ | ✓ | ✓ | active identity captured + restored |
| version tracking (per-assay class + object version) | ✓ | ✓ | ✓ | `assay@RNA:Assay5`, `object@5.4.0`, … |
| very old serialized objects (Seurat v2, pre-`Assay`) | ✓ | — | ✓ | dedicated **read path** for the legacy lowercase `seurat` S4 class (slots read via `attr()`, so the ancient class needn't be defined); old→new on write-back. CI fixture is built from Seurat 2.3.4's **authoritative** class defs (slot-exact vs the source) and read back with the class *undefined* — the real scenario. No `real` cell: a genuine full ancient-Seurat install is infeasible on the available toolchain (its shiny/plotly/hdf5r stack won't compile on R 3.6.2) |
| spatial coords (`so@images` Visium/FOV/Slide-seq → `spatial` axis) | ✓ | — | ✓ | mirrors the AnnData path; multi-section uses partial coverage; was a silent loss |
| image pixels / `@commands` | ✗ | — | — | deferred spatial tier / provenance not typed (pixels recorded in `dropped`) |

> Heavy backends (full `Seurat` umbrella for SCTransform, `Signac` for ChromatinAssay, `BPCells`) are
> optional: absent → the case degrades to a recorded SKIP, never a hard failure. CI loads the lightweight
> `SeuratObject` (the full `Seurat` package can install-but-fail-to-load on a minimal runner); the SCT
> path runs automatically wherever full Seurat loads (locally, and in CI once that's resolved).

## SingleCellExperiment — R

Also handles a plain **`SummarizedExperiment`** (no SCE-only accessors). Real grounding: a sweep over
Bioconductor `scRNAseq` — **56/61 PASS** (the 5 non-passes are missing loader packages / a dataset
segfault, not profile bugs) after 3 sweep-caught fixes; plus real `ZeiselBrain`.

| case | status | real | CI (synth) | notes |
|---|:---:|:---:|:---:|---|
| assays (counts / logcounts / multiple) | ✓ | ✓ | ✓ | |
| `reducedDims` (+ `rotation` → loadings) | ✓ | ✓ | ✓ | |
| **altExps** (ADT / spike-ins → feature axes) | ✓ | ✓ | ✓ | real ZeiselBrain ERCC; synthetic ADT |
| `colData` / `rowData` factors → factor axes | ✓ | ✓ | ✓ | |
| S4 `Rle` columns / nested `DataFrame`/`GRanges` | ✓ | ✓ | — | Rle unpacked; uncoercible nested cols **recorded** (sweep-caught) |
| **NULL dimnames** (cells keyed by `Barcode` colData) | ✓ | ✓ | — | labels synthesized (sweep-caught) |
| plain `SummarizedExperiment` (no reducedDims/altExps) | ✓ | ✓ | — | guarded accessors (sweep-caught) |
| `metadata` (free-form) | ◐ | ✓ | ✓ | recorded |
| `colPairs` / `rowPairs` (cell-cell / gene-gene graphs) | ✓ | — | ✓ | → relations over (cells,cells)/(genes,genes); round-trip nnz+values |

## Conos / pagoda2 — R

| case | status | real | CI (synth) | notes |
|---|:---:|:---:|:---:|---|
| collection (multi-sample, heterogeneous, even cross-species) | ✓ | ✓ | ✓ | real `conI.rds`/`con.rds` (4 & 8 samples); synthetic collection |
| joint graph / per-sample embeddings | ✓ | ✓ | ✓ | auto-upgrades an old igraph graph |
| `read_conos(ds)` → live Conos (round-trip) | ✓ | ✓ | — | rebuilds per-sample Pagoda2 + joint graph/embedding/clusters; `conformance/conos.sh` |
| collection → **Seurat v5** split assay | ✓ | ✓ | ✓ | per-sample layers over union genes + `Graphs()` + `DimReduc` + meta; reads back as a collection; no corrected matrix |
| collection → **AnnData** (one matrix; flattens) | ✓ | ✓ | ✓ | X = raw joint counts, graph → `obsp`/`connectivities` (+ `uns['neighbors']`), embedding → `obsm`; per-sample PCA dropped |
| pagoda2 viewer schema + DE/pseudobulk over a factor axis | ✓ | — | ✓ | mock pagoda2 (no real pagoda2 object in the corpus yet) |
| real pagoda2 object | ✗ | — | — | none sourced yet |

---

## Known gaps & roadmap

Prioritized, with **where** the gap is (profile / real corpus / synthetic fixture):

1. **Spatial tier** (Visium / Xenium / CosMx / Slide-seq) — *deliberately deferred beyond the concept*.
   Spatial **coordinates** are supported now (a named observed `spatial` coordinate axis); the deferred
   part is **images, vendor coordinate frames, and molecule tables** (kept in the passthrough, not lost).
   Affects Seurat images and the SeuratData `stx*`/`ssHippo` datasets.
2. **A real `.h5mu` multiome** — *corpus gap*: partial-overlap is now implemented (typed partial coverage
   via `index`, tested on a constructed fixture); what's missing is a *real* partial-overlap / multiome
   `.h5mu` in the local corpus (real ATAC is already covered via the Seurat `pbmcMultiome` sweep).
3. **Seurat `@commands`** provenance (analysis history) — *profile gap*: not typed (the `Neighbor` nn
   graph is now typed as a relation).
4. **Faithful partial coverage** — *profile refinement*: subset PCA **loadings are now typed** (over a
   `<reduction>_features` subset axis); subset `scale.data` and SCT residuals are still *recorded as
   losses* (`◐`) rather than typed.
5. **Real pagoda2 object** in the corpus, and a **real ATAC `.h5mu`** — *corpus gaps*.
6. **Azimuth reference atlases** (SeuratData `*ref`) — *corpus/tooling gap*: need the Azimuth loader to
   read them (they're disk datasets); would add SCTAssay breadth.

## Where the evidence lives

- [`python/tests/CORPUS.md`](python/tests/CORPUS.md) — the real datasets each test is grounded in (and
  their synthetic CI stand-ins).
- [`conformance/sweep/REPORT.md`](conformance/sweep/REPORT.md) — the breadth sweep over whole
  repositories, and every bug it caught.
- [`conformance/`](conformance/) — the cross-language / cross-format round-trip suite that backs the
  language matrix; [`conformance/README.md`](conformance/README.md) documents the **origin-coverage** rule
  (test every producer, not just consumers) and the per-test coverage matrix.
- [`docs/mapping.md`](docs/mapping.md) — the deterministic role→slot conversion contract (what lands
  where in each target, and why), and the **native-acceptance** check (`lstar convert --check`) that opens
  the produced object in its own library and runs a canonical-ops smoke — verifying native tools accept
  it, not just that bytes round-tripped (`conformance/convert_cli.sh`). The same contract governs the
  **package-free** backends (`--backend direct`): `.h5ad` ↔ store via `h5py` (no anndata) and Seurat
  `.rds` ↔ store via base R (no SeuratObject) — see the dependency table in `docs/conversions.md`.
