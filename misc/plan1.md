# lstar — implementation plan (plan1)

Working plan for the `lstar` package: a lightweight, fast "glue" library implementing the L★ model
(see [`Lstar_proposal.md`](./Lstar_proposal.md)) — a uniform in-memory representation, a Zarr
on-disk format, and bidirectional conversion to/from the common single-cell formats (AnnData,
Seurat, the pagoda-verse). Bindings in **R**, **Python**, and **C++**. Intended to be the shared
interchange/representation layer used by pagoda2.1 and our other packages, and usable by anyone else.

Items marked **[DECISION]** are the major forks to settle before/while building; the recommendation
is given but should be confirmed.

---

## 1. Goals and principles

- **Lightweight.** Small dependency footprint; native-format support is *optional* (loaded only when
  the host package is present), so `lstar` never drags AnnData/Seurat/scanpy as hard dependencies.
- **Fast.** The hot paths — sparse orientation conversion (CSR↔CSC↔COO), heterogeneous "gather"
  across collection samples, recipe (virtual-field) evaluation, chunk codec — should be C/C++ speed.
- **One canonical format.** The L★ Zarr store is the single on-disk representation; everything is a
  conversion to/from it.
- **Spec-driven, conformance-tested.** The spec and a shared conformance suite (round-trip tests) are
  the source of truth; multiple language implementations are kept honest by the suite, not by sharing
  every line of code.
- **Reusable.** pagoda2.1, conos, cacoa, and the viewer all read/write the same store.

---

## 2. Architecture **[DECISION — the central one]**

Two facts constrain the design: (a) reading a *live* native object is inherently language-bound (a
Seurat S4 object only exists in R; an AnnData object only in Python), so the **profile/ingest layer
is necessarily per-language**; (b) the *core* — the in-memory axis/field model, the Zarr store IO,
and the translation primitives — can be shared or reimplemented.

**Recommended: a shared C++ core for R & C++, with Python native on `zarr-python`.**

```
        ┌───────────────────────── spec + conformance suite (source of truth) ─────────────────────────┐
        │                                                                                               │
  Python package  (native, on zarr-python)            R package        C++ programs (pagoda2 src, …)
   - Dataset/Axis/Field model in Python                │  cpp11 bindings │   link libstar directly
   - L★ zarr IO via zarr-python (+ fsspec/dask)        └──────┬──────────┘
   - anndata profile (read/write live AnnData)                │
                                                       ┌───────┴────────────────────────────────┐
                                                       │  libstar  (C++ core)                     │
                                                       │   - axis/field/model registry            │
                                                       │   - L★ zarr IO (minimal, dep-light)      │
                                                       │   - translation primitives (CSR/CSC/COO, │
                                                       │     gather, recipe eval, codecs)         │
                                                       └──────────────────────────────────────────┘
   R-side profiles (seurat, pagoda2) read live R objects → hand matrices/metadata to libstar.
```

Rationale: Python gets the mature Zarr ecosystem (sharding, codecs, fsspec remote, dask laziness) for
free and native AnnData interop; R and C++ get one fast, self-contained engine without depending on a
weak R-Zarr library or on Python; the two Zarr IO implementations (zarr-python and libstar) are kept
consistent by the conformance suite. C++ is a *first-class binding target* (pagoda2's `src/` can use
libstar directly), which the user asked for.

Alternatives:

- **(B) Single C++ core for all three, Python bound to libstar.** Maximal consistency (one IO impl),
  but Python loses zarr-python/dask/fsspec and gains build weight (a C++ extension where a pure-Python
  package would do). Not recommended for a library meant to feel native to Python users.
- **(C) All-native now, C++ core later.** Python on zarr-python, R on `Rarr`/`pizzarr`, share only the
  spec+suite; add libstar once profiling justifies it. Lightest to start and fastest to a working
  demo, but R-Zarr write/v3/sharding support is uneven, and "C++ bindings" is deferred.

**Sub-decision — C++ Zarr IO** (only if a C++ core is built): avoid TensorStore (excellent but
Bazel-built → painful to vendor into a CRAN R package). Options: roll a minimal Zarr-v3 reader/writer
(JSON metadata + chunk blobs + a small codec set: raw/gzip/blosc; sharding is a thin add) against
`blosc2` + `zlib` + a header-only JSON lib (`nlohmann/json`); or vendor `z5`. **Lean: roll-minimal**,
because our store layout is narrow and well-defined, and it keeps packaging light.

**Sub-decision — core language.** The user's stack is C++ (Rcpp across pagoda2/sccore/conos) and asked
for C++ bindings, so **C++**. *Noted alternative:* Rust (`zarrs` crate is the best Zarr lib today,
binds cleanly to Python via PyO3 and R via extendr) would give better Zarr support and easier
packaging, at the cost of a toolchain our ecosystem doesn't use and slightly more friction calling it
from C++. Recommend C++ to match the ecosystem unless we want to invest in Rust deliberately.

---

## 3. Repository layout (monorepo)

```
lstar/
  README.md  LICENSE                         # MIT or BSD-3-Clause (open standard → permissive)  [DECISION: license]
  spec/                                       # the normative spec, curated from the proposal
    lstar-spec.md                             # core model + Zarr schema
    roles.md                                  # the role registry + controlled vocabularies
    vocabulary.md                             # the shared core vocabulary (v0)
    profiles/{anndata,seurat,pagoda2,conos,cacoa}.md
  core/                                       # libstar (C++)            [present iff arch A/B]
    include/lstar/*.hpp   src/*.cpp   CMakeLists.txt
  python/                                     # `lstar` Python package (native on zarr-python)
    pyproject.toml   src/lstar/{model,zarr_io,profiles/,vocab}.py
  R/                                          # `lstar` R package (cpp11 → libstar)
    DESCRIPTION  NAMESPACE  R/*.R  src/ (vendored libstar + cpp11 glue)
  conformance/                                # shared, language-agnostic test assets
    datasets/                                 # reference L★ stores + reference native files
    suite/                                    # round-trip + cross-profile coverage tests + runner
  misc/  Lstar_proposal.md  plan1.md
  .github/workflows/   ci.yml                 # build+test python & R; run conformance suite
```

---

## 4. The core model & API (each binding mirrors this)

In-memory objects: `Dataset` (holds `axes` and `fields` registries, optional `models`), `Axis`
(labels + kind + provenance), `Field` (span, role, encoding, state, coverage, uncertainty,
provenance, lazy values). Fields are **lazy by default** — values are zarr-backed and materialized on
access (Python: dask/zarr; C++: chunk-on-demand) so collections and >1M-cell data stay out-of-core.

Accessor surface (names final per language idiom; this is the contract):

```
open(path) / write(ds, path)                  # L★ zarr store IO (local or HTTP)
ds.axes() / ds.fields(axis=…) / ds.supports(cap)
ds.cellMeta() / featureMeta() / sampleMeta()  # table views (bundles of arity-1 fields)
ds.featureVector(name)                         # gather per-sample in a collection (partial-aware)
ds.reduction(name) / embedding(name) / grouping(name) / relation(name) / result(name) / model(name)
add_field(name, values, role=…, span=…, …)     # only `values` required; role/span/coverage inferred
```

Translation primitives in libstar (the speed-critical bits, shared): sparse `csr↔csc↔coo`, dense↔sparse,
`gather`/`scatter` across a union axis (heterogeneous collections), `recipe` evaluation (e.g. apply a
normalization recipe to raw counts on read), and the Zarr chunk codec.

---

## 5. Profiles & format support

Profiles are the per-language readers/writers (Appendix B of the proposal is the spec). Two ingest
modes: **(i) live object** (language-bound: a Seurat S4 in R, an AnnData in Python) and **(ii) native
file** (h5ad, h5Seurat, rds — cross-language in principle via HDF5).

- **anndata** (Python, M2): read/write live `AnnData`; also h5ad/zarr files. Optional dep on `anndata`.
- **seurat** (R, M3): read/write live `Seurat` (v5); `DimReduc` ⇄ shared latent axis, `@commands` ⇄
  provenance. Optional dep on `SeuratObject`.
- **pagoda2** (R, M3): read/write `Pagoda2`; mostly renaming given pagoda2.1's accessors.
- **singlecellexperiment / multiassayexperiment** (R/Bioconductor, M3–M5): the shared-class profile;
  `assays`→measures, `reducedDims`→embeddings, `colData`/`rowData`→arity-1 fields, `altExp`→a second
  feature axis, MAE→a collection. Key for Bioconductor-community adoption.
- **conos / cacoa** (R, M5): collection-level profiles.
- **mudata / spatialdata** (later): extensions.

Each profile ships its rule table + the two round-trip conformance results (proposal 4.2).

---

## 6. Zarr realization

Per Appendix A: `axes/`, `fields/`, `models/` groups; L★ metadata under an `"lstar"` `.zattrs` key;
consolidated `.zmetadata`; Zarr **v3** with **sharding** (few objects for HTTP) and the encoding set
`dense | csr | csc | coo | edge_list | ragged | raster | recipe`. Packaging: standalone sample,
bundled collection, or by-reference collection (members pinned by id+hash). The store **is** the
viewer format (the app reads a projection of it).

---

## 7. Conformance suite (built early, M0–M4)

The enforcement mechanism from proposal Part 5.3. A corpus of reference datasets + native files, and
two tests per profile: `write_P(read_P(native)) == native` on the lossless subset, and
`read_P(write_P(ds))` preserves what the target can hold; plus cross-profile coverage
(`convert(anndata→seurat)` preserves the shared-vocabulary core). Run in CI against both the Python
and R/C++ implementations — this is what keeps the two IO impls from drifting.

---

## 8. Dependencies & packaging (the "lightweight" budget)

- **Python**: `zarr>=3`, `numpy`, `scipy` (sparse). Optional extras: `[anndata]`, `[h5ad]`. Pure-Python
  package (no compiled extension) under arch A. PyPI.
- **R**: `cpp11`, `Matrix`. `SeuratObject`, `pagoda2`, `conos` in **Suggests** (profiles load only if
  present). Vendored libstar in `src/` with `blosc2`+`zlib`+`nlohmann_json`. Aim CRAN-installable.
- **C++/libstar**: `blosc2`, `zlib`, `nlohmann/json` (header). CMake. Installable standalone.

---

## 9. Milestones

- **M0 — scaffold.** Repo, license, CI skeleton, spec extracted into `spec/`, conformance harness stub.
- **M1 — L★ zarr core (Python).** `Dataset`/`Axis`/`Field`; read/write a store round-trip; lazy fields.
- **M2 — anndata profile (Python).** AnnData ⇄ L★; round-trip conformance on the shared-vocab core.
- **M3 — libstar + R; the cross-language demo.** C++ core (model + minimal zarr IO + primitives); R
  package via cpp11; `seurat` + `pagoda2` profiles. **Demo: `h5ad → L★ zarr` (Python) → `Seurat` (R)**,
  shared-vocab core lossless, off-vocab items reported. This demo is the proposal's core evidence.
- **M4 — conformance suite.** Formalize corpus + round-trip + cross-profile coverage; CI on Python & R.
- **M5 — collections.** `conos`, `cacoa` profiles; union axis, designs, gather; pagoda2.1 adoption.
- **M6 — performance & viewer.** Profile/optimize hot paths (share libstar into Python via pybind only
  if needed); wire the viewer to read the store.

---

## 10. Decisions

Settled (2026-06-11):

1. **Architecture:** option A — a shared **C++ core (libstar)** for the R and C++ bindings; **Python
   native** on zarr-python; the two IO paths kept consistent by the conformance suite.
2. **Core language:** **C++**, with a rolled minimal Zarr-v3 IO (blosc + zlib + nlohmann/json; no
   Bazel/TensorStore).
3. **v0 scope:** L★ zarr core + **AnnData + Seurat + pagoda2**, culminating in the cross-language
   `h5ad → L★ → Seurat` demo; conos/cacoa deferred to M5.
4. **License:** **MIT** (default; permissive, to maximize adoption as shared glue).
5. **Spec:** curate a tight `spec/lstar-spec.md` from the proposal alongside M1, so the schema is
   exercised by real code as it is written.

Dev-environment note: this workstation runs Python 3.8 (zarr-python v3 needs ≥3.11), so the Python
reference implementation targets **Zarr v2** here for now, written to be v3-ready; the C++ core and CI
target Zarr **v3 + sharding**.

---

## 11. Risks

- **Two Zarr IO implementations drift** (Python vs libstar) → mitigated by the conformance suite as CI.
- **C++/R packaging weight** undermines "lightweight" → mitigated by roll-minimal Zarr IO (no Bazel).
- **Vocabulary fragmentation** (proposal Part 5) is the standard-level risk, not a code risk, but
  lstar's job is to make the *reference* mappings concrete and tested so others converge on them.
- **Scope creep** into a compute framework → keep lstar to representation + conversion only.

---

## 12. Optimization & performance strategy

lstar is glue on hot paths (large sparse matrices, many fields, remote stores), so performance is a
design constraint, not an afterthought:

- **Lazy by default.** `read()` returns fields whose values are zarr-backed handles, materialized on
  access — converting or inspecting a dataset never loads the whole thing. (Python: keep the zarr
  array / optional dask; C++: chunk-on-demand.)
- **Streaming / chunked IO.** Read and write in chunks; never require the full matrix in RAM. Sparse
  matrices are chunked along the primary axis (CSC by gene-blocks, CSR by cell-blocks) so a single
  gene/cell or a block streams. CSR↔CSC operates block-wise where possible; a full transpose falls
  back to an out-of-core two-pass.
- **Multithreading.** The translation primitives (sparse transpose/convert, gather across collection
  samples, per-gene/per-cell summaries, codec encode/decode) parallelize over chunks/columns. The C++
  core uses OpenMP (matching pagoda2), thread count via a parameter/env with a cores−2 default cap;
  Python uses zarr's threaded chunk IO and drops heavy loops to the C++ core when present.
- **Zero-copy at the boundaries.** R↔C++ pass matrix buffers by pointer (cpp11/Rcpp) without copying;
  a future Python↔C++ path uses the buffer protocol / Arrow.
- **Codecs.** Production chunk compression via **blosc2** (fast, multithreaded); the minimal C++ IO
  supports **raw + gzip** (zlib) for portability and cross-impl interop, with blosc2 behind a flag.
  Embeddings may be int16-quantized for the viewer path.
- **No implicit densification.** Sparse stays sparse; gather/subset operate on sparse structure.
- **Collections gather, never concatenate.** `featureVector` over a collection streams per-sample
  slices and assembles them (partial-coverage aware); no global matrix is built.

Targets (revisit after benchmarks): convert a 1M-cell h5ad ⇄ L★ within ~native memory; per-gene fetch
= O(one chunk); CSR↔CSC of a 1M×30k sparse in seconds, multithreaded.

**Measured (2026-06-11), on real data.** Tabula Muris Senis Marrow — **40,220 cells × 20,138 genes,
77.6M nonzeros, float32** (`examples/real_perf.py`, `examples/roundtrip_chain.py`, `bench_colstats`,
`test_chunked`); 56-core box.
- **Lazy open.** `read(lazy=True)` on the 77.6M-nnz measure opens in **0.18 s / +9 MB** vs eager
  **4.6 s / +779 MB** (~87× less memory); on the 299 MB conos store, **+11 MB vs +368 MB**.
- **Streaming reduction.** Per-gene mean/var with on-the-fly `log1p` over all 20,138 genes runs in
  bounded **~165 MB** (block-streamed), matching a dense `np.var` ground truth to float32 precision;
  identical across thread counts.
- **Memory-lean precision.** Float32 measures stay float32 end to end (read, transpose, reduce) — no
  widening copy. Sums accumulate in float64 (numpy `bincount` / C++ `double` accumulators), so the
  moments are float64-accurate at the stored memory cost (low-precision storage, high-precision
  accumulation). The C++ kernels are templated on the value dtype to read float32 in place.
- **Multithreading, controllable from the call.** Both `lazy.stream_col_stats(n_threads=)` and
  `csc_col_mean_var(..., n_threads)` take the thread count (1 = serial, N, 0 = auto), results
  thread-invariant. On the real measure: C++ log1p reduction **2.46 s → 0.32 s (7.6×) at 16 threads**;
  the Python streamer ~**1.6–2× at 8 threads** (GIL-bound on the glue, numpy kernels release it).
  `bench_colstats` (16M-nnz synthetic): **5.2× plain / 6.6× log1p** at 16 threads.
- **Compression.** gzip-5 shrinks the conos raw-count store **314 → 66 MB (4.8×)**; the Marrow
  *normalized* float matrix only **1.6×** (high entropy) — ratio tracks the data, as expected.
- **Round-trip to the origin.** AnnData → L★ → AnnData over the real Marrow is a **fixed point across
  4 cycles** — `X`, `obs`, `obsm` (PCA/tSNE/UMAP), and `obsp` graphs return identical; `uns` is recorded
  in `dropped`, never silently lost. So a conversion chain of any length returns to the native format.
- **Portability.** The C++ core reads multi-chunk + gzip/zlib stores (any chunk grid, fill-padded edge
  chunks, float32/float64), writes a consolidated `.zmetadata`; Python/C++/R all read the same
  chunked+gzip store; `csc_to_csr` is an O(nnz) orientation flip that preserves the value dtype.

---

## 13. Progress log & reassessments

- **2026-06-11 — M1 done.** Python core (`Dataset`/`Axis`/`Field` + Zarr IO) round-trips; store layout
  matches Appendix A. Repo git-initialized (local only, no remote).
- **Reassessment M1 → M2.** Added the **SingleCellExperiment/MAE** profile (user request; Bioconductor
  adoption) and the optimization strategy (§12). Codec decision: minimal C++ Zarr IO uses raw+gzip for
  cross-impl interop; Python writer gains a `compressor` option; blosc2 deferred to production.
  Proceeding to M2 (AnnData profile — runnable now on the installed `anndata 0.8`).
- **2026-06-11 — M2 + M3 done.** AnnData profile round-trips (via the profile and through the zarr
  store); losses (uns) recorded, never silent. C++ `libstar` core (header-only model + minimal Zarr v2
  IO: single-chunk, uncompressed, utf8+offset strings; nlohmann/json) builds with g++9.4/cmake, and the
  **Python → C++ → Python cross-impl test is byte-faithful** — the format works across independent
  implementations. String encoding switched to utf8-bytes+offsets (cross-language-friendly); writer
  defaults to uncompressed chunks.
- **Validators & profile execution (answering: do we need validators/parsers?).** *Validators: yes.* A
  structural **store validator** (`lstar.validate`) is added (spans reference existing axes, shapes
  match axis lengths, relations span 2 axes, open-vocabulary roles/states warn). The profile
  round-trip, cross-profile, and cross-impl conformance tests are the *profile* validators (M5). *A
  profile-rule parser/DSL: no, by decision.* Mappings are logic-heavy (byDtype, varm↔obsm pairing,
  Seurat `DimReduc` exploding into 3 fields, conos gather, cacoa design), so profiles are code and the
  conformance suite — not a parser — guarantees they match the documented rule tables.
- **Reassessment M3 → next.** Next: the **R package** binding `libstar` (cpp11) + the
  **seurat/pagoda2/SCE** profiles, toward the headline cross-language `h5ad → L★ → Seurat` demo. C++
  multi-chunk reads, gzip/blosc codecs, `.zmetadata` consolidation, and the translation primitives
  (CSR↔CSC, gather, OpenMP) are deferred to the performance pass.
- **2026-06-11 — R package + 3 profiles + conformance done.** The R package binds `libstar` via cpp11
  (read+write); **Python ↔ C++ ↔ R** all read/write the same store faithfully. Profiles implemented and
  tested: **anndata** (Python), **seurat** (R, v5), **singlecellexperiment** (R). The headline
  **`AnnData (Python) ↔ L★ ↔ Seurat (R)`** demo runs (`examples/cross_language_demo.sh`), and a
  **cross-format conformance chain** `AnnData → L★ → Seurat → L★ → SCE → L★` preserves the
  shared-vocabulary core (`conformance/cross_format.sh`). A master runner (`conformance/run.sh`) builds
  the C++ core, installs the R package, and runs all Python + cross-impl + cross-format tests — all
  green. CI workflow added (pins `zarr<3` until the Python v3 port).
- **Reassessment → remaining.** Left: the **pagoda2** profile (pagoda2 is installed; needs a live
  Pagoda2 object to test), **lazy/streaming + perf** (Python lazy read; a C++ OpenMP translation
  primitive to exercise the threading model), and the deferred C++ items (multi-chunk, codecs, v3).
  Note: the Python reference impl uses the **Zarr v2 API**; a v3 port is the main portability follow-up.
- **2026-06-11 — multithreading demonstrated.** Added `csc_col_mean_var` (column-parallel OpenMP) to
  `libstar` with an on-the-fly `log1p`-normalized variant — the pagoda2 lazy-normalized-view pattern,
  which is compute-bound where the plain moments are bandwidth-bound. `bench_colstats` shows ~5.9× (plain)
  / ~7.2× (lognorm) at 16 threads; the earlier "1.00×" was a thread-state leak (`omp_set_num_threads(1)`
  not reset), fixed by passing explicit thread counts. Header re-vendored into the R package.
- **2026-06-11 — collection model (the L★ differentiator), driven by real data.** User flagged that the
  test datasets were toy-sized *and* that the realistic ones I reached for (Tabula Muris) are themselves
  *collections* flattened into one matrix — exactly the over-normalized shape L★ argues against. Built the
  **collection** representation and tested it on real collections: a `samples` axis + per-sample
  `cells.<s>`/`genes.<s>` axes & measures (heterogeneous: differing cells, possibly differing gene sets) +
  a union `cells` axis carrying the joint embedding/clusters/**graph (as a `relation`)**.
  - **Conos profile** (`R/R/profile_conos.R`, `write_conos`): ingests the real two-sample object
    (16,022 cells, 2 samples, 25.2M nnz, joint graph 312,912 nnz) → a **299 MB** L★ store that **Python
    reads & validates clean** (`examples/conos_collection_demo.R`). This is the headline collection result.
  - **Seurat v5 split-layer** collection: `read_seurat` recognizes a `split(assay, f=sample)` Assay5
    (`counts.<sample>` layers covering subsets of cells) as a collection (`kind="collection"`, per-sample
    axes + `samples` axis + `sample` design label). Demo + conformance added.
  - Portable **collection conformance** (`conformance/collection.sh`, no external data) + CI; master
    runner green end to end.
- **2026-06-11 — graceful version recognition (user request).** Profiles now *detect* the source format
  version and adapt instead of assuming one layout, recording `<format>@<version>` in `ds$profiles` and
  routing unrepresentable pieces to `dropped`: **Seurat** v3/v4 `Assay` vs v5 `Assay5` (`.seurat_versions`
  + a version-agnostic `.seurat_layer_access` with a `GetAssayData` fallback for SeuratObject < 5, both
  v3/v5 tested); **pagoda2** `getRawCounts()` accessor vs legacy `$counts` slot (recorded as
  `pagoda2-accessor`/`pagoda2-slot`); **AnnData** library version + the `.raw` slot kept on its own
  `genes_raw` axis when its gene set diverges (round-trips back into `adata.raw`); **SCE** version.
  Locked in by `python/tests/test_versions.py` (in the suite + CI).
- **2026-06-11 — lazy/streaming (Python) + C++ chunks/codecs/transpose.**
  - **Python lazy read.** `read(path, lazy=True)` leaves fields as `lazy.LazyDense`/`LazyCSX` proxies
    over the open zarr arrays (and prefers `open_consolidated`): the store opens without materializing
    the heavy arrays (+2 MB vs +133 MB on an 8M-nnz measure; +11 MB vs +368 MB on the conos store).
    `lazy.stream_col_stats` reduces a CSC measure by column block — per-gene mean/var with on-the-fly
    `log1p`, zero-aware (`ss = Σx² − (Σx)²/nrows`), bit-exact vs eager, in ~50 MB constant memory.
    `python/tests/test_lazy.py`, `examples/lazy_streaming_demo.py`.
  - **Python write** gained `chunk_elems` (chunk arrays along the first axis) and `compressor`
    (numcodecs, e.g. GZip) — chunking is what makes lazy reads touch only the needed blocks; gzip-5
    shrinks the conos store 314 → 66 MB (4.8×).
  - **C++ reader** now handles an arbitrary **chunk grid** (C-order, fill-padded edge chunks, missing
    chunk → fill 0) and **gzip/zlib** chunk decompression (via system zlib, `LSTAR_HAVE_ZLIB`; auto-detected
    header). Index arrays are dtype-normalized (`as_i64`) since scipy emits int32/int64 by size. The
    **writer** emits a consolidated `.zmetadata`. New primitive **`csc_to_csr`** (O(nnz) storage
    transpose) for the genes×cells orientation flips profiles need, without densifying.
  - Cross-impl proven on a **chunked + gzip** store: Python writes → C++ reads (and re-emits with
    `.zmetadata`) → Python `open_consolidated` re-reads identically; R reads it too.
    `conformance/chunked.sh` + `core/test/test_chunked.cpp`, in the suite + CI. R install switched to
    `--preclean` (R's make has no header-dep tracking, so a changed vendored header was silently stale).
- **Reassessment → next.** Lazy/streaming + the C++ chunk/codec/transpose work are in, tested, and
  benchmarked (see §12 Measured). Still open: a single-sample **pagoda2** `write_pagoda2` (the Conos path
  already reads per-sample Pagoda2 counts/PCA; factor it out); **blosc2** codec and **Zarr v3 + sharding**
  (the Python impl still uses the v2 API); a multithreaded/out-of-core full **transpose** and collection
  **gather**; a Python↔C++ zero-copy boundary.
