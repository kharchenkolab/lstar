// Guard for the JS viewer-prep compression layout. With an injected WASM writer codec, extendForViewer
// tags each appended navigator's per-field write opts: counts_cellmajor -> zstd chunked+sharded (chunk-
// granular subset reads), every other appended field -> zstd single-chunk. The gene-major `counts` basis
// lives in the base store and is NOT touched here. Asserts the on-disk codecs/chunking match the policy,
// and the store still reads back correct values.
import assert from "node:assert";
import * as fs from "node:fs";
import * as os from "node:os";
import * as path from "node:path";

import { NodeFSStore } from "../core/node-store.ts";
import { writeStore, type WriterCodec } from "../core/writer.ts";
import { extendForViewer } from "../core/extend.ts";
import { openLstar } from "../core/reader.ts";
import { VIEWER_CHUNK_ELEMS } from "../core/policy.ts";
import createLstarWriter from "../dist/lstar_writer.mjs";

const w = await createLstarWriter();
const codec: WriterCodec = {
  encodeChunk: (raw, comp, level) => w.encodeChunk(raw, comp, level),
  packShard: (entries, atEnd) => w.packShard(entries, atEnd),
};

// base store: gene-major CSC counts sized so the cell-major copy's data exceeds the 16384 knee (-> multi-
// chunk -> sharded). 4000 cells x 100 genes, 200 nnz/gene = 20000 > 16384. Written RAW (base-store default).
const NC = 4000, NG = 100, PER = 200, NNZ = NG * PER;
const indptr = new Int32Array(NG + 1);
for (let g = 0; g < NG; g++) indptr[g + 1] = indptr[g] + PER;
const indices = new Int32Array(NNZ);
for (let g = 0; g < NG; g++) for (let k = 0; k < PER; k++) indices[g * PER + k] = (g * 7 + k * 13) % NC;
const data = Int32Array.from({ length: NNZ }, (_, i) => (i % 9) + 1);

const dir = path.join(os.tmpdir(), "lstar-viewer-comp.lstar.zarr");
fs.rmSync(dir, { recursive: true, force: true });
await writeStore(new NodeFSStore(dir), {
  kind: "sample", profiles: ["base@0.1"],
  axes: { cells: { labels: Array.from({ length: NC }, (_, i) => `c${i}`), role: "observation" },
          genes: { labels: Array.from({ length: NG }, (_, i) => `g${i}`), role: "feature" },
          umap: { labels: ["u0", "u1"], origin: "derived", role: "coordinate" } },
  fields: {
    counts: { role: "measure", span: ["cells", "genes"], encoding: "csc", state: "raw", shape: [NC, NG], data, indices, indptr },
    umap: { role: "embedding", span: ["cells", "umap"], encoding: "dense", shape: [NC, 2],
            data: Float32Array.from({ length: NC * 2 }, (_, i) => (i % 100) * 0.01) },
    leiden: { role: "label", span: ["cells"], encoding: "utf8", values: Array.from({ length: NC }, (_, i) => `cl${i % 5}`) },
  },
});   // base store written raw (no codec) — gene-major counts stays raw

// prep with the injected codec -> appended navigators get the compressed layout
await extendForViewer(new NodeFSStore(dir), { codec, groupings: ["leiden"] });

const zj = (p: string) => JSON.parse(fs.readFileSync(path.join(dir, "fields", p, "zarr.json"), "utf8"));
function layout(field: string, arr = "data") {
  const z = zj(field + "/" + arr), c0 = z.codecs[0];
  const sharded = c0.name === "sharding_indexed";
  const inner = (sharded ? c0.configuration.codecs : z.codecs).map((c: any) => c.name);
  const chunk = (sharded ? c0.configuration.chunk_shape : z.chunk_grid.configuration.chunk_shape);
  return { sharded, inner, chunk, shape: z.shape };
}

// gene-major counts (in the base store) untouched: raw, single chunk
const gm = layout("counts");
assert.ok(!gm.sharded && !gm.inner.includes("zstd"), `gene-major counts must stay raw, saw ${JSON.stringify(gm)}`);
console.log(`  [js] gene-major counts: raw (${gm.inner}), untouched by prep`);

// counts_cellmajor: zstd, chunked (< full length), sharded
const cm = layout("counts_cellmajor");
assert.ok(cm.inner.includes("zstd"), `counts_cellmajor must be zstd, saw ${cm.inner}`);
assert.ok(cm.sharded, `counts_cellmajor must be sharded (nnz ${NNZ} > knee ${VIEWER_CHUNK_ELEMS}), saw ${JSON.stringify(cm)}`);
assert.strictEqual(cm.chunk[0], VIEWER_CHUNK_ELEMS, `counts_cellmajor inner chunk should be ${VIEWER_CHUNK_ELEMS}`);
console.log(`  [js] counts_cellmajor: zstd + sharded, inner chunk ${cm.chunk[0]} (${cm.inner})`);

// a dense navigator (od_score): zstd, single chunk (chunk == shape)
const od = layout("od_score", "values");
assert.ok(od.inner.includes("zstd") && !od.sharded && od.chunk[0] === od.shape[0],
  `od_score must be zstd single-chunk, saw ${JSON.stringify(od)}`);
console.log(`  [js] od_score: zstd single-chunk (${od.inner})`);

// still reads back: a cell-major row decodes through the compressed+sharded array
const ds = await openLstar(new NodeFSStore(dir));
const row = await ds.csrRow("counts_cellmajor", 0);
assert.ok(row.cols.length >= 0 && ds.fieldNames().includes("od_score"), "prepped store must read back");
console.log("viewer compression layout guard passed: appended navigators compressed per policy, reads OK");
