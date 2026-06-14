# Recipes (task → how)

Concrete, copy-pasteable solutions. Runnable end-to-end versions live in `examples/`.

## Convert h5ad → L★ → Seurat (cross-language)
```python
# Python
import anndata as ad, lstar
from lstar.profiles.anndata import read_anndata
lstar.write(read_anndata(ad.read_h5ad("a.h5ad")), "a.lstar.zarr")
```
```r
# R
ds <- lstar::lstar_read("a.lstar.zarr"); so <- lstar::write_seurat(ds)
```

## Round-trip back to the original format (variable length)
A chain of any length returns to the native format; what can't be carried is in `dropped`.
```python
from lstar.profiles.anndata import read_anndata, write_anndata
cur = adata
for _ in range(N):
    lstar.write(read_anndata(cur), "t.lstar.zarr")
    cur = write_anndata(lstar.read("t.lstar.zarr"))
# cur == adata on X/obs/obsm/obsp; uns recorded in ds.dropped (a fixed point)
```
Cross-format+language: `bash examples/roundtrip_xlang.sh <N> <file.h5ad>`
(AnnData → Seurat → SCE → … → AnnData).

## Ingest a collection of samples (don't flatten)
```r
ds <- lstar::write_conos(conos_obj)          # samples axis + per-sample axes/counts/pca +
                                             # union cells + sample label + joint graph(relation)
```
```r
# Seurat v5 split assay as a collection
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$sample)
ds <- lstar::read_seurat(obj)                # kind == "collection"
```

## Per-gene mean/variance at scale (lazy, streamed, threaded)
```python
ds = lstar.read("big.lstar.zarr", lazy=True)          # opens without reading the matrix
mean, var, nnz = lstar.stream_col_stats(ds.field("counts").values,
                                        lognorm=True, n_threads=0)   # all cores; C++ if available
hvg = np.argsort(var)[::-1][:2000]
```

## Write a chunked + compressed store (enables true streaming)
```python
import numcodecs
lstar.write(ds, "big.lstar.zarr", chunk_elems=1_000_000, compressor=numcodecs.GZip(5))
```
Chunking lets a lazy read touch only the blocks it needs; gzip shrinks raw-count stores ~5×
(normalized float matrices less — ratio tracks entropy).

## Check whether the fast path is active
```python
lstar.show_config()        # 'C++ accelerator ACTIVE' (+ OpenMP threads) or 'pure-Python fallback'
assert lstar.has_accel()   # in tests that must run the C++ engine
```
Force an engine for benchmarking: `engine="python"`/`"c++"` arg, or `LSTAR_ENGINE=python`.

## Validate a store / dataset
```python
errs = lstar.validate(lstar.read("s.lstar.zarr"))   # [] means valid
```

## Read an L★ store in C++ and reduce
```cpp
auto ds = lstar::read("s.lstar.zarr");
auto* f = ds.field("counts");                       // CSC measure
auto ip = lstar::as_i64(f->indptr);
auto s  = lstar::csc_col_mean_var(f->data.as<float>(), ip.data(),  // float32 in place
                                  f->shape[1], f->shape[0], 0, true);
```

## Write or extend an L★ store from JS/WASM (browser viewer)
```ts
import { writeStore, addToStore, type Compressor } from "lstar-js";
import createLstarKernels from "lstar-js/wasm";
const k = await createLstarKernels();
const gzip: Compressor = { id: "gzip", level: 1, compress: (b) => new Uint8Array(k.gzipCompress(b, 1)) };
await writeStore(store, dataset, { chunkElems: 1 << 18, compressor: gzip });   // chunked + gzip via WASM zlib
await addToStore(store, { fields: { od_score: { encoding: "dense", span: ["genes"], shape: [ng], data } } });
```
Every encoding writes (CSC/dense/categorical+factor/mask/partial/aux) and round-trips to Python/R/C++.
The compressor is **injected** (so the uncompressed path needs no WASM); `compressor: null`/omitted ->
one uncompressed chunk. Offsets/index are `<i8`; codes `<i4` (-1=missing). Verified by `conformance/js.sh`
(JS-write → Python/C++-read). Source: `js/core/writer.ts`, `js/wasm/lstar_wasm.cpp` (`gzipCompress`).

## Retest the profiles against the large LOCAL corpus (rare; keeps the synthetic CI tier honest)
```bash
bash conformance/sweep/retest_local.sh          # faithfulness guard + all cached-data sweeps + a TRIAGE
RETEST_HEAVY=1 bash conformance/sweep/retest_local.sh   # + GB-scale SeuratData sweeps
```
Two tiers: CI runs synthetic-only fixtures (`synth.py`, `LSTAR_SYNTHETIC_CORPUS=1`); the **local** corpus
(`testdata/`, gitignored) is what catches the long tail. The loop, run rarely: **(1)** the faithfulness
guard (`test_synth_faithful.py`) checks real-structure ⊆ synthetic — a gap means **update `synth.py`**;
**(2)** sweeps emit `/tmp/sweep_*.tsv` — a `FAIL`/`VALIDATE-ERR` is a **profile bug** (fix + add a fixture),
`LOADERR`/`SKIP` is just absent data; **(3)** refresh the coverage docs (`sweep/REPORT.md`, `sweep/CATALOG.md`,
`SUPPORT.md`) and re-grep for stale version/format enumerations. The orchestrator edits nothing — it names
what to update. Full directive: **`conformance/sweep/RETEST.md`**.

## Gotchas
- **Orientation:** L★ measures are cells×genes; Seurat/SCE are genes×cells (profiles transpose).
- **Big R in conformance → temp file, not `Rscript -e`:** a large `-e` program overflows R's ~8 KB
  command-line buffer and is silently ignored (`WARNING: '-e ...'`, runs nothing) → trips `pipefail` at a
  trailing `grep`. Quoted heredoc → temp file, `$RLIB` via `commandArgs`, `</dev/null` (see `seurat_v2.sh`).
- **dtype widths:** indptr/indices are int32 or int64 by size — use `as_i64` in C++; measures are
  often float32 — keep them float32, accumulate in float64.
- **Don't widen to gain precision** — it breaks the memory-lean contract; accumulate in float64 instead.
- **Single-chunk stores** can't stream lazily (a slice still reads the whole chunk) — write with
  `chunk_elems=` to stream.
- **R rebuilds:** after changing a vendored header, `R CMD INSTALL --preclean` (no header-dep tracking).
- **Relations** span two axes (e.g. `["cells","cells"]`); set `role="relation"` explicitly.

## Where things live
- Examples: `examples/*.py`, `examples/*.R`, `examples/*.sh`.
- Conformance: `conformance/run.sh` (master), plus `cross_format.sh`, `collection.sh`, `chunked.sh`.
- Design: `misc/Lstar_proposal.md`; plan + measured perf: `misc/plan1.md` §12.
- Docs: `docs/` (principles, model & format specs, worked examples).
