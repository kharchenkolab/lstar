# lstar — JavaScript / WebAssembly

The browser/Node side of lstar. The plan (see [`../docs`](../docs) and `pagoda2/misc/app_prop2.md`) is
a TypeScript data layer that reads **L★ Zarr stores** over HTTP or local disk, *sparingly* — fetching
the chunks/views a viewer needs — and runs the heavy compute in **WebAssembly**, reusing the same
`libstar` C++ core that backs the R (cpp11) and Python (pybind11) packages. So the numbers a browser
shows match R and Python exactly.

```
TypeScript (zarrita.js): fetch chunks, range reads, v3 sharding, store navigation, deck.gl, UI
        │  assembled typed arrays            ▲ results (typed arrays)
        ▼                                    │
WebAssembly (@lstar/wasm): decode/normalize/reduce/rank — the libstar kernels
```
WASM can't `fetch`; JS is slow at tight numeric loops — so I/O stays in TS, compute goes to WASM.

## Status

- **WASM kernels** (`dist/lstar_kernels.mjs`, Phase A): `libstar` compute primitives
  (`colMeanVar`, `cscToCsr`) compiled via Emscripten/embind. Verified in Node to match a dense
  reference *and* the Python kernel. Single-threaded for now (OpenMP/pthreads needs cross-origin
  isolation).
- **`core/` reader** (Phase B): a zarrita.js reader for L★ stores (any zarrita store — `FetchStore`
  over HTTP, a zip store, or `NodeFSStore` for Node). Reads the consolidated manifest and fetches
  fields *lazily*, including a single CSC **gene-column** as a slice (the viewer hot path). Works on
  chunked + gzip stores.
- **`core/` view** (Phase C): a framework-agnostic `LstarView` — `embedding`, `metadata` (categorical
  codes+categories / numeric), `geneExpression` (on-the-fly log1p of one gene's column), `colStats`
  (per-gene mean/var via WASM), `subsampleDE` (ranked genes A-vs-B), and a typed-array `Crossfilter`.

## API

```ts
import { openLstar, LstarView, Crossfilter, scalarToRGBA } from "lstar-js";
import { NodeFSStore } from "lstar-js/node-store";        // browser: new zarrita FetchStore(url)

const ds = await openLstar(new NodeFSStore("sample.lstar.zarr"));
ds.kind; ds.axisNames(); ds.fieldNames();                 // manifest
await ds.cscColumn("counts", geneIndex);                  // one gene's nonzeros, fetched as a slice

const view = new LstarView(ds);
const { data, n, dim } = await view.embedding("umap");    // positions for a point layer
const expr = await view.geneExpression("g3", { lognorm: true });   // per-cell scalar
const colors = scalarToRGBA(expr.values, expr.max);       // -> RGBA for the color attribute
const md = await view.metadata("leiden");                 // {kind:'categorical', codes, categories}
const { mean, var: variance } = await view.colStats({ lognorm: true });   // HVG (WASM)
const ranked = await view.subsampleDE(cellsA, cellsB);    // genes by |log fold change|
const sel = new Crossfilter(n).categorical(md.codes, [0]).selected();
```

Low-level kernels:
```js
import createLstarKernels from "lstar-js/wasm";
const M = await createLstarKernels();
M.colMeanVar(dataF64, indptrI32, nrows, /*n_threads*/ 1, /*lognorm*/ true);
M.cscToCsr(dataF64, indicesI32, indptrI32, nrows, ncols);
```

## Build & test

Needs [emsdk](https://emscripten.org/) (`EMSDK`, or `~/emsdk`), **Python ≥ 3.10** for `emcc`
(`LSTAR_EMCC_PYTHON=/path/to/python3.10` if the system one is older), a modern Node (emsdk bundles
one), and `npm install` (for zarrita).

```bash
LSTAR_EMCC_PYTHON=/path/to/python3.10 bash build.sh         # WASM -> dist/
PYTHONPATH=../python/src python3 test/make_store.py         # test store + references -> test/data/
node test/kernels.test.mjs
node --experimental-strip-types test/reader.test.ts
node --experimental-strip-types test/view.test.ts
```
Or all of it via `bash ../conformance/js.sh` (guarded; skips if emsdk/zarrita absent). `dist/` and
`test/data/` are git-ignored; `test/fixture_python.json` is a committed cross-language fixture.

The TypeScript is written in erasable syntax so Node runs it directly with `--experimental-strip-types`
(no build step for tests); a real consumer would type-check/bundle with `tsc`/esbuild.
