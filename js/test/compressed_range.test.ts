// Guard for the COMPRESSED byte-range path (_readRange -> _readCompressedRange -> readChunkDecoded).
// A compressed array must be sub-ranged at CHUNK granularity — decode only the covering chunk(s), not the
// whole array. The escape this catches (which a values-only test would miss): a compressed range path that
// silently reads the WHOLE array returns correct values, just slowly. So we assert BOTH (a) cscColumn ==
// the uncompressed twin, AND (b) reading one column fetches ~one data chunk object, not all of them.
// Cases: a column within one chunk, one in the PARTIAL last chunk, and one that STRADDLES a chunk boundary
// — for zstd unsharded AND zstd sharded.
import assert from "node:assert";
import * as fs from "node:fs";
import * as os from "node:os";
import * as path from "node:path";

import { NodeFSStore } from "../core/node-store.ts";
import { writeStore, type WriterCodec } from "../core/writer.ts";
import { openLstar } from "../core/reader.ts";
import type { LstarStore } from "../core/reader.ts";
import createLstarWriter from "../dist/lstar_writer.mjs";

// Counts a store's reads by key, so we can prove a column read touches ~one data chunk (not the whole array).
class CountingStore implements LstarStore {
  gets: string[] = [];
  private inner: NodeFSStore;
  constructor(inner: NodeFSStore) { this.inner = inner; }
  async get(key: string): Promise<Uint8Array | undefined> { this.gets.push(key); return this.inner.get(key); }
  async getRange(key: string, a: number, b: number): Promise<Uint8Array | undefined> {
    this.gets.push(key); return this.inner.getRange!(key, a, b);
  }
  async getSuffix(key: string, n: number): Promise<Uint8Array | undefined> {
    this.gets.push(key); return this.inner.getSuffix!(key, n);
  }
  reset() { this.gets = []; }
  // distinct chunk objects fetched under fields/<field>/data (v2 `N`, v3 `c/N`, or a shard object)
  dataChunkFetches(field: string): number {
    const re = new RegExp(`fields/${field}/data/(c/)?[0-9]`);
    return new Set(this.gets.filter((k) => re.test(k))).size;
  }
}

const w = await createLstarWriter();
const codec: WriterCodec = {
  encodeChunk: (raw, comp, level) => w.encodeChunk(raw, comp, level),
  packShard: (entries, atEnd) => w.packShard(entries, atEnd),
};

// 40 cells x 10 genes CSC. Per-gene nnz [8,3,7,5,9,2,6,4,8,3] (total 55). With chunkElems=8 the data/indices
// arrays are 7 chunks (8×6 + 7): gene 0 -> chunk 0 (within one), gene 2 -> data[11:18] straddles chunks
// 1|2, gene 9 -> data[52:55] in the partial last chunk 6.
const NNZ = [8, 3, 7, 5, 9, 2, 6, 4, 8, 3];
const indptr = new Int32Array(NNZ.length + 1);
for (let g = 0; g < NNZ.length; g++) indptr[g + 1] = indptr[g] + NNZ[g];
const total = indptr[NNZ.length];
const data = Int32Array.from({ length: total }, (_, i) => i + 1);
const indices = Int32Array.from({ length: total }, (_, i) => i % 40);
const CHUNK = 8;

function spec() {
  return {
    kind: "sample" as const, profiles: ["comprange@0.1"],
    axes: { cells: { labels: Array.from({ length: 40 }, (_, i) => `c${i}`), role: "observation" as const },
            genes: { labels: Array.from({ length: 10 }, (_, i) => `g${i}`), role: "feature" as const } },
    fields: { counts: { role: "measure" as const, span: ["cells", "genes"] as [string, string], encoding: "csc" as const,
                        state: "raw", shape: [40, 10] as number[], data, indices, indptr } },
  };
}

async function write(out: string, opts: { compressor: any; shardElems?: number }) {
  fs.rmSync(out, { recursive: true, force: true });
  await writeStore(new NodeFSStore(out), spec() as any,
    { chunkElems: CHUNK, compressor: opts.compressor, shardElems: opts.shardElems, codec }, "v3");
}

const base = path.join(os.tmpdir(), "lstar-comprange");
const REF = base + "-raw.lstar.zarr";       // uncompressed twin (reference values)
const ZS = base + "-zstd.lstar.zarr";        // zstd, chunked, unsharded
const ZSH = base + "-zstd-shard.lstar.zarr"; // zstd, chunked, sharded
await write(REF, { compressor: null });
await write(ZS, { compressor: { id: "zstd", level: 3 } });
await write(ZSH, { compressor: { id: "zstd", level: 3 }, shardElems: 32 });

// reference column values from the uncompressed store
const ref = await openLstar(new NodeFSStore(REF));
const refCols: Array<{ rows: number[]; vals: number[] }> = [];
for (let g = 0; g < 10; g++) {
  const c = await ref.cscColumn("counts", g);
  refCols.push({ rows: [...c.rows].map(Number), vals: [...c.vals].map(Number) });
}
assert.ok(refCols.some((c) => c.vals.length > 0), "reference must be non-empty");

// each compressed store: every column equals the uncompressed twin
for (const [label, dir] of [["zstd unsharded", ZS], ["zstd sharded", ZSH]] as const) {
  const cs = new CountingStore(new NodeFSStore(dir));
  const ds = await openLstar(cs);
  // genuinely compressed on disk?
  const zj = JSON.parse(fs.readFileSync(path.join(dir, "fields/counts/data/zarr.json"), "utf8"));
  const codecs = (zj.codecs[0].name === "sharding_indexed" ? zj.codecs[0].configuration.codecs : zj.codecs).map((c: any) => c.name);
  assert.ok(codecs.includes("zstd"), `${label}: counts/data should carry a zstd codec, saw ${codecs}`);

  // granularity FIRST, on a COLD cache: reading ONE within-a-chunk column (gene 0 -> data[0:8], chunk 0)
  // must fetch ~one data chunk object, NOT the whole array (7 chunks). This is the whole-array-cliff guard.
  // (Run before the value loop below, which would warm the decoded-chunk cache and hide the fetch.)
  cs.reset();
  await ds.cscColumn("counts", 0);
  const n = cs.dataChunkFetches("counts");
  console.log(`  [js] ${label}: cscColumn(gene 0) fetched ${n} data chunk object(s) [cold cache]`);
  assert.ok(n >= 1 && n <= 2, `${label}: expected ~1 data chunk fetch for a within-chunk column, got ${n} (whole-array read?)`);

  for (let g = 0; g < 10; g++) {
    const c = await ds.cscColumn("counts", g);
    assert.deepStrictEqual([...c.rows].map(Number), refCols[g].rows, `${label}: gene ${g} rows`);
    assert.deepStrictEqual([...c.vals].map(Number), refCols[g].vals, `${label}: gene ${g} vals`);
  }
  console.log(`  [js] ${label}: all 10 columns == uncompressed twin (incl. straddling g2 + partial-last g9)`);
}

console.log("compressed range guard passed: chunk-granular compressed byte-range == uncompressed, no whole-array read");
