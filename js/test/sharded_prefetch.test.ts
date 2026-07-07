// Guard for the sharded whole-array PREFETCH (keysFor). Reading a whole sharded array must prefetch only
// the SHARD objects that exist (c/0..c/M-1), NOT the inner-chunk keys (c/0..c/N-1, N>=M). libzarr's
// ArrayMeta.chunk_shape is the INNER shape, so enumerating it would name c/M..c/N-1 that 404 —
// correctness-safe (decode reads only the shards) but it defeats sharding's fewer-HTTP-objects purpose.
// Here: write a store whose `counts/data` packs many inner chunks into ONE shard, read the whole field
// through a store that counts gets + misses, and assert NO miss on the data array's chunk objects.
import assert from "node:assert";
import * as fs from "node:fs";
import * as os from "node:os";
import * as path from "node:path";

import { NodeFSStore } from "../core/node-store.ts";
import { writeStore, type WriterCodec } from "../core/writer.ts";
import { openLstar } from "../core/reader.ts";
import type { LstarStore } from "../core/reader.ts";
import createLstarWriter from "../dist/lstar_writer.mjs";

// A store that records every get() and whether it was a miss (key absent) — the signal keysFor over-fetch
// produces (a fetch for an inner-chunk key that no shard object backs).
class CountingStore implements LstarStore {
  gets: string[] = [];
  misses: string[] = [];
  private inner: NodeFSStore;
  constructor(inner: NodeFSStore) { this.inner = inner; }
  async get(key: string): Promise<Uint8Array | undefined> {
    this.gets.push(key);
    const v = await this.inner.get(key);
    if (v === undefined) this.misses.push(key);
    return v;
  }
  async getRange(key: string, a: number, b: number) { return this.inner.getRange!(key, a, b); }
  async getSuffix(key: string, n: number) { return this.inner.getSuffix!(key, n); }
}

const out = path.join(os.tmpdir(), "lstar-shard-prefetch.lstar.zarr");
fs.rmSync(out, { recursive: true, force: true });

const w = await createLstarWriter();
const codec: WriterCodec = {
  encodeChunk: (raw, comp, level) => w.encodeChunk(raw, comp, level),
  packShard: (entries, atEnd) => w.packShard(entries, atEnd),
};

// 6 cells x 5 genes CSC, 8 nnz. chunkElems=2 -> data/indices are 4 inner chunks; shardElems=16 packs all
// of them into ONE shard object per array -> N=4 inner keys but M=1 store object.
const data = new Int32Array([1, 2, 3, 4, 5, 6, 7, 8]);
const indices = new Int32Array([0, 3, 1, 2, 5, 0, 4, 2]).map((x) => x % 6);
const indptr = new Int32Array([0, 2, 3, 5, 6, 8]);   // 5 genes
await writeStore(new NodeFSStore(out), {
  kind: "sample", profiles: ["shardpf@0.1"],
  axes: { cells: { labels: ["c0", "c1", "c2", "c3", "c4", "c5"], role: "observation" },
          genes: { labels: ["g0", "g1", "g2", "g3", "g4"], role: "feature" } },
  fields: { counts: { role: "measure", span: ["cells", "genes"], encoding: "csc", state: "raw",
                      shape: [6, 5], data, indices, indptr } },
}, { chunkElems: 2, compressor: { id: "gzip", level: 5 }, shardElems: 16, codec }, "v3");

// sanity: the data array is genuinely sharded, and has FEWER shard objects on disk than inner chunks
const dataDir = path.join(out, "fields/counts/data");
assert.strictEqual(JSON.parse(fs.readFileSync(path.join(dataDir, "zarr.json"), "utf8")).codecs[0].name,
  "sharding_indexed", "counts/data must be sharding_indexed");
const shardObjs = fs.existsSync(path.join(dataDir, "c")) ? fs.readdirSync(path.join(dataDir, "c")) : [];
console.log(`  [js] counts/data shard objects on disk: ${shardObjs.length} (${shardObjs.join(",")})`);

// whole-array read through the counting store
const store = new CountingStore(new NodeFSStore(out));
const ds = await openLstar(store);
const sp = await ds.fieldSparse("counts");
assert.deepStrictEqual([...sp.data], [1, 2, 3, 4, 5, 6, 7, 8], "values must round-trip");

const isDataChunk = (k: string) => /fields\/counts\/(data|indices)\/c\//.test(k);
const dataChunkGets = store.gets.filter(isDataChunk);
const dataChunkMisses = store.misses.filter(isDataChunk);
console.log(`  [js] counts data/indices chunk-object gets: ${dataChunkGets.length}, misses: ${dataChunkMisses.length}`);
assert.strictEqual(dataChunkMisses.length, 0,
  `keysFor over-fetch: reading a sharded array 404'd on ${dataChunkMisses.length} inner-chunk keys that no shard backs: ${[...new Set(dataChunkMisses)].join(", ")}`);
assert.ok(dataChunkGets.length > 0, "expected at least one shard-object read for the sharded field");

console.log("sharded prefetch guard passed: whole-array read fetched only shard objects (no inner-key 404s)");
