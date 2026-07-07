// Write a comprehensive L* store from JS -- chunked + gzip-compressed (via the WASM zlib kernel), with
// EVERY encoding: CSC measure, dense embedding, categorical (induced factor axis), nullable mask,
// partial coverage, and an aux passthrough. writer_crossread.py reads it back in Python (+ the C++
// reader) and asserts validate-clean + value equality. Usage: node --experimental-strip-types
// writer_make.ts <out-dir>. Pure-JS (no WASM) uncompressed mode if SkipWasm is set.
import { NodeFSStore } from "../core/node-store.ts";
import { writeStore, type Compressor } from "../core/writer.ts";

import createLstarKernels from "../dist/lstar_kernels.mjs";

const OUT = process.argv[2] || "/tmp/lstar-writer-cross.lstar.zarr";
const FORMAT = (process.argv[3] as "v2" | "v3") || "v2";   // on-disk Zarr format (default v2)

const k = await createLstarKernels();
const compressor: Compressor = { id: "gzip", level: 1, compress: (raw) => new Uint8Array(k.gzipCompress(raw, 1)) };

// CSC counts, 10 cells x 6 genes (column-major over genes)
const counts = {
  data: new Int32Array([1, 2, 3, 4, 5, 6, 7, 8]),
  indices: new Int32Array([0, 3, 1, 2, 5, 0, 7, 9]),  // cell rows
  indptr: new Int32Array([0, 2, 3, 5, 6, 8, 8]),      // per-gene (6 cols + 1)
};
const N = 10;

await writeStore(new NodeFSStore(OUT), {
  kind: "sample", profiles: ["jswriter@0.1"], dropped: [],
  axes: {
    cells: { labels: Array.from({ length: N }, (_, i) => `cell${i}`), role: "observation" },
    genes: { labels: ["g0", "g1", "g2", "g3", "g4", "g5"], role: "feature" },
    umap: { labels: ["umap0", "umap1"], origin: "derived", role: "coordinate" },
    celltype: { labels: ["T", "B", "NK"], origin: "derived", role: "factor", inducedBy: "celltype" },
  },
  fields: {
    counts: { role: "measure", span: ["cells", "genes"], encoding: "csc", state: "raw",
              shape: [N, 6], data: counts.data, indices: counts.indices, indptr: counts.indptr },
    umap: { role: "embedding", span: ["cells", "umap"], encoding: "dense", shape: [N, 2],
            data: Float32Array.from({ length: N * 2 }, (_, i) => i * 0.5) },
    celltype: { role: "label", span: ["cells"], encoding: "categorical", ordered: false,
                codes: new Int32Array([0, 1, 0, 2, 1, 0, -1, 2, 1, 0]), categories: ["T", "B", "NK"] },
    qc: { role: "measure", span: ["cells"], encoding: "dense", shape: [N],
          data: Float32Array.from({ length: N }, (_, i) => i + 0.25),
          mask: new Uint8Array([0, 0, 1, 0, 0, 0, 0, 1, 0, 0]) },   // cells 2 & 7 missing
    adt: { role: "measure", span: ["cells"], encoding: "dense", shape: [5],   // measured on 5 of 10 cells
           data: new Float32Array([10, 11, 12, 13, 14]),
           index: new BigInt64Array([0n, 2n, 4n, 6n, 8n]), indexAxis: "cells" },
  },
  aux: {
    "test.uns": {
      tree: { $obj: { n_pca: 50, method: "leiden", scores: { $array: "a0" }, names: { $strings: "a1", shape: [2] } } },
      arrays: [
        { id: "a0", kind: "dense", data: new Float64Array([1.1, 2.2, 3.3]) },
        { id: "a1", kind: "utf8", values: ["foo", "bar"] },
      ],
    },
  },
}, { chunkElems: 4, compressor }, FORMAT);   // chunkElems=4 forces multi-chunk; gzip via WASM

console.log(`wrote chunked+gzip ${FORMAT} store ->`, OUT);
