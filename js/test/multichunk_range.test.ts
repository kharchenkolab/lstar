// Multi-chunk byte-range fast path — the landmine guard for `_readRange`'s chunk-offset math.
// `_readRange` now range-reads a MULTI-chunk uncompressed array when [lo,hi) lands within one chunk
// (chunk `floor(lo/chunkLen)`, offset relative to that chunk). An off-by-one at a chunk boundary would
// silently return WRONG values that still decode — so we build a CSC field whose `data`/`indices` are
// chunked with a PARTIAL LAST CHUNK, and assert `cscColumn` / `cscColumns` equal the whole-matrix ground
// truth (`fieldSparse`, read via zarrita, independent of the fast path) across FS-dir / FS-zip / HTTP-dir
// / HTTP-zip. Columns are laid out to exercise: within chunk 0 (== the old single-chunk path), within
// chunk 1 and the PARTIAL chunk 2 (the new `ci>0` offset math), a boundary-SPANNING column (falls back to
// the zarrita slice), and an EMPTY column. Self-contained — builds its own fixture via the JS writer.
import assert from "node:assert";
import * as fs from "node:fs";
import * as fsp from "node:fs/promises";
import * as http from "node:http";
import * as os from "node:os";
import * as path from "node:path";

import { openLstar } from "../core/reader.ts";
import { NodeFSStore } from "../core/node-store.ts";
import { HttpStore } from "../core/http-store.ts";
import { ZipStore, httpZipSource } from "../core/zip.ts";
import { writeStore, type DatasetSpec } from "../core/writer.ts";
import { nodeFileSource, writeStoreZip } from "../core/zip-node.ts";

// ---- a small CSC (gene-major) matrix: 60 cells × 8 genes, 37 nonzeros, chunked 16/chunk ----
// chunkElems=16 -> data/indices chunks [0,16) [16,32) [32,37)  (3 chunks; last is PARTIAL, 5 elems)
// column byte spans (from indptr):
//   g0 [0,3)    within chunk0        g4 [20,30)  within chunk1 (ci=1, NEW)
//   g1 [3,9)    within chunk0        g5 [30,32)  within chunk1 (ci=1, NEW)
//   g2 [9,11)   within chunk0        g6 [32,37)  within chunk2 (ci=2, PARTIAL, NEW)
//   g3 [11,20)  SPANS chunk0/1 -> fallback       g7 []         empty
const NCELLS = 60, NGENES = 8, CHUNK = 16;
const colRows: number[][] = [
  [1, 5, 10], [2, 3, 4, 6, 7, 8], [0, 9], [11, 12, 13, 14, 15, 16, 17, 18, 19],
  [20, 21, 22, 23, 24, 25, 26, 27, 28, 29], [30, 31], [32, 33, 34, 35, 36], [],
];
const indptr: number[] = [0];
const indices: number[] = [], data: number[] = [];
let v = 1;
for (const rows of colRows) { for (const r of rows) { indices.push(r); data.push(v++); } indptr.push(indices.length); }
assert.equal(indices.length, 37, "fixture nnz");
assert.equal(indptr.length, NGENES + 1);

const spec: DatasetSpec = {
  kind: "sample",
  axes: {
    cells: { labels: Array.from({ length: NCELLS }, (_, i) => "c" + i), origin: "observed" },
    genes: { labels: Array.from({ length: NGENES }, (_, i) => "g" + i), origin: "observed" },
  },
  fields: {
    counts: {
      role: "measure", span: ["cells", "genes"], encoding: "csc", state: "raw", shape: [NCELLS, NGENES],
      data: Float32Array.from(data), indices: Int32Array.from(indices), indptr: Int32Array.from(indptr),
      write: { chunkElems: CHUNK },   // force multi-chunk data/indices with a partial last chunk
    },
  },
};

// ground truth per column, straight from the arrays (what fieldSparse must reconstruct)
const truthCol = (c: number) => ({ rows: colRows[c].slice(), vals: colRows[c].map((_, k) => data[indptr[c] + k]) });

function serve(rootDir: string, zipPath: string): Promise<{ url: string; close: () => void }> {
  const server = http.createServer(async (req, res) => {
    try {
      const u = new URL(req.url!, "http://x");
      const p = u.pathname === "/store.zip" ? zipPath : path.join(rootDir, decodeURIComponent(u.pathname));
      let stat; try { stat = await fsp.stat(p); } catch { res.statusCode = 404; return res.end(); }
      if (!stat.isFile()) { res.statusCode = 404; return res.end(); }
      const size = stat.size;
      if (req.method === "HEAD") { res.setHeader("content-length", String(size)); res.statusCode = 200; return res.end(); }
      const range = req.headers.range;
      if (range) {
        const m = /bytes=(\d+)-(\d+)?/.exec(range)!;
        const start = parseInt(m[1], 10), end = m[2] ? parseInt(m[2], 10) : size - 1;
        res.statusCode = 206; res.setHeader("content-range", `bytes ${start}-${end}/${size}`);
        res.setHeader("content-length", String(end - start + 1));
        fs.createReadStream(p, { start, end }).pipe(res);
      } else { res.statusCode = 200; res.setHeader("content-length", String(size)); fs.createReadStream(p).pipe(res); }
    } catch (e) { res.statusCode = 500; res.end(String(e)); }
  });
  return new Promise((resolve) => server.listen(0, "127.0.0.1", () => {
    const a = server.address() as any; resolve({ url: `http://127.0.0.1:${a.port}`, close: () => server.close() });
  }));
}

const eq = (a: ArrayLike<number>, b: ArrayLike<number>, msg: string) =>
  assert.deepEqual(Array.from(a as any).map(Number), Array.from(b as any).map(Number), msg);

async function checkStore(label: string, store: any, rangeKeys?: string[]) {
  const ds = await openLstar(store);

  // ground truth via fieldSparse (whole-array reads through zarrita — independent of the fast path)
  const sp = await ds.fieldSparse("counts");
  const gp = Array.from(sp.indptr as any).map(Number);
  const colFromWhole = (c: number) => ({
    rows: Array.from(sp.indices as any).slice(gp[c], gp[c + 1]).map(Number),
    vals: Array.from(sp.data as any).slice(gp[c], gp[c + 1]).map(Number),
  });
  for (let c = 0; c < NGENES; c++) assert.deepEqual(colFromWhole(c), truthCol(c), `${label}: fieldSparse col ${c}`);

  // per-column cscColumn — exercises the multi-chunk fast path (within chunk0/1/2) + the spanning fallback (g3)
  for (let c = 0; c < NGENES; c++) {
    const col = await ds.cscColumn("counts", c);
    eq(col.rows, truthCol(c).rows, `${label}: cscColumn(${c}).rows`);
    eq(col.vals, truthCol(c).vals, `${label}: cscColumn(${c}).vals`);
  }

  // batched cscColumns — all columns in order
  const all = await ds.cscColumns("counts", [0, 1, 2, 3, 4, 5, 6, 7]);
  assert.equal(all.cols.length, NGENES);
  for (let c = 0; c < NGENES; c++) {
    eq(all.cols[c].rows, truthCol(c).rows, `${label}: cscColumns all[${c}].rows`);
    eq(all.cols[c].vals, truthCol(c).vals, `${label}: cscColumns all[${c}].vals`);
  }

  // batched cscColumns — SCATTERED subset, INPUT order, with a DUPLICATE and an EMPTY column
  const pick = [6, 0, 4, 6, 7];
  const got = await ds.cscColumns("counts", pick);
  assert.equal(got.cols.length, pick.length, `${label}: cscColumns returns one entry per input col (incl. dup)`);
  pick.forEach((c, i) => {
    eq(got.cols[i].rows, truthCol(c).rows, `${label}: cscColumns pick[${i}]=g${c}.rows`);
    eq(got.cols[i].vals, truthCol(c).vals, `${label}: cscColumns pick[${i}]=g${c}.vals`);
  });

  // confirm the MULTI-chunk fast path actually fired (a range read into a non-zero chunk index), not that
  // everything silently fell back. Chunk keys differ by format: v2 `.../data/1`, v3 `.../data/c/1` (default).
  if (rangeKeys) assert(rangeKeys.some((k) => /\/(data|indices)\/(?:c\/)?[12]$/.test(k)),
    `${label}: expected a byte-range read into chunk 1 or 2 (multi-chunk fast path) — saw: ${[...new Set(rangeKeys)].join(", ")}`);

  console.log(`  [js] multi-chunk range parity OK: ${label}`);
}

const tmp = await fsp.mkdtemp(path.join(os.tmpdir(), "lstar-mc-"));
try {
  const dir = path.join(tmp, "s.lstar.zarr");
  const zip = path.join(tmp, "s.lstar.zarr.zip");
  await writeStore(new NodeFSStore(dir), spec);
  await writeStoreZip(zip, spec);

  // FS dir — instrument getRange to prove the fast path fired on a non-zero chunk index
  const fsStore = new NodeFSStore(dir);
  const rangeKeys: string[] = [];
  await checkStore("FS dir  (NodeFSStore)", {
    get: fsStore.get.bind(fsStore),
    getRange: (k: string, a: number, b: number) => { rangeKeys.push(k); return fsStore.getRange!(k, a, b); },
  }, rangeKeys);

  await checkStore("FS zip  (nodeFileSource)", await ZipStore.open(nodeFileSource(zip), zip));

  const srv = await serve(dir, zip);
  try {
    await checkStore("HTTP dir(HttpStore)    ", new HttpStore(srv.url + "/"));
    await checkStore("HTTP zip(httpZipSource)", await ZipStore.open(httpZipSource(srv.url + "/store.zip")));
  } finally { srv.close(); }
} finally {
  await fsp.rm(tmp, { recursive: true, force: true });
}
console.log("multi-chunk range fast-path + cscColumns parity tests passed");
