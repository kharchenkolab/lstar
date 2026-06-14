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
- **`core/` writer** (`writer.ts`): write a full L★ store from JS, or `addToStore` derived fields onto an
  existing (e.g. Python-written) one — **every encoding** both directions: CSC/dense, UTF-8, categorical
  (induces a factor axis), nullable mask, partial coverage, and the aux passthrough. Arrays are chunked
  (along axis 0) and optionally **gzip-compressed via the WASM zlib kernel** (`gzipCompress`); the
  compressor is dependency-injected so the uncompressed path needs no WASM. A consolidated `.zmetadata` is
  emitted. The bytes round-trip to the Python/R/C++ readers (`conformance/js.sh`: JS-write → Python/C++-read).

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

Write a store (chunked + gzip via the WASM zlib kernel), or extend an existing one:
```ts
import { writeStore, addToStore, type Compressor } from "lstar-js";
import createLstarKernels from "lstar-js/wasm";
const k = await createLstarKernels();
const gzip: Compressor = { id: "gzip", level: 1, compress: (b) => new Uint8Array(k.gzipCompress(b, 1)) };

await writeStore(new NodeFSStore("out.lstar.zarr"), {        // browser: any zarrita writable store
  axes: { cells: { labels: [...] }, genes: { labels: [...] } },
  fields: { counts: { encoding: "csc", span: ["cells","genes"], shape: [nc, ng], data, indices, indptr },
            leiden: { encoding: "categorical", span: ["cells"], codes, categories: ["A","B"], ordered: false } },
}, { chunkElems: 1 << 18, compressor: gzip });               // omit opts -> single uncompressed chunk

await addToStore(store, { fields: { od_score: { encoding: "dense", span: ["genes"], shape: [ng], data } },
                          profiles: ["viewer@0.1"] });        // append derived viewer fields
```

Low-level kernels:
```js
import createLstarKernels from "lstar-js/wasm";
const M = await createLstarKernels();
M.colMeanVar(dataF64, indptrI32, nrows, /*n_threads*/ 1, /*lognorm*/ true);
M.cscToCsr(dataF64, indicesI32, indptrI32, nrows, ncols);
M.gzipCompress(bytesU8, /*level*/ 1);   // gzip (RFC1952) for the write path
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
node --experimental-strip-types test/writer.test.ts            # write -> read-back (JS)
node --experimental-strip-types test/writer_make.ts /tmp/x.lstar.zarr   # write chunked+gzip, all encodings
PYTHONPATH=../python/src python3 test/writer_crossread.py /tmp/x.lstar.zarr   # Python reads it + validates
```

`bash ../conformance/js.sh` runs all of the above (kernels, reader, view, writer, and the JS-write →
Python-read cross-language gate) in one go.
Or all of it via `bash ../conformance/js.sh` (guarded; skips if emsdk/zarrita absent). `dist/` and
`test/data/` are git-ignored; `test/fixture_python.json` is a committed cross-language fixture.

The TypeScript is written in erasable syntax so Node runs it directly with `--experimental-strip-types`
(no build step for tests); a real consumer would type-check/bundle with `tsc`/esbuild.
