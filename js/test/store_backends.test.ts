// Store-backend VALUE parity — the layer real consumers use. Every store backend, read through
// `openLstar`, must return byte-for-byte identical field *values* (not just open successfully): a
// categorical, a dense embedding, a full sparse matrix, AND a byte-range sparse column (the getRange
// fast path). Compared across FS-dir, FS-zip (nodeFileSource), HTTP-dir (HttpStore), and HTTP-zip
// (httpZipSource) — the four paths consumers actually use.
//
// This is the leg that would have caught the ZipStore leading-slash bug: zarrita forms chunk keys with
// a leading slash, and only value reads exercise them. "Opens" != "reads data correctly" — a store that
// returns undefined for every chunk still opens (metadata is slashless / from `.zmetadata`) and yields
// non-empty axes/fields; it just silently zeros all data. So we assert *decoded values*, not open.
//
//   usage: store_backends.test.ts <dir> <zip>
import assert from "node:assert";
import * as fs from "node:fs";
import * as fsp from "node:fs/promises";
import * as http from "node:http";
import * as path from "node:path";

import { openLstar } from "../core/reader.ts";
import { NodeFSStore } from "../core/node-store.ts";
import { HttpStore } from "../core/http-store.ts";
import { ZipStore, httpZipSource } from "../core/zip.ts";
import { nodeFileSource } from "../core/zip-node.ts";

const [, , dir, zip] = process.argv;

/** A tiny static file server with `Range` support: serves the store tree and the zip at /store.zip. */
function serve(rootDir: string, zipPath: string): Promise<{ url: string; close: () => void }> {
  const server = http.createServer(async (req, res) => {
    try {
      const u = new URL(req.url!, "http://x");
      const p = u.pathname === "/store.zip" ? zipPath : path.join(rootDir, decodeURIComponent(u.pathname));
      let stat;
      try { stat = await fsp.stat(p); } catch { res.statusCode = 404; return res.end(); }
      if (!stat.isFile()) { res.statusCode = 404; return res.end(); }
      const size = stat.size;
      if (req.method === "HEAD") { res.setHeader("content-length", String(size)); res.statusCode = 200; return res.end(); }
      const range = req.headers.range;
      if (range) {
        const m = /bytes=(\d+)-(\d+)?/.exec(range)!;
        const start = parseInt(m[1], 10);
        const end = m[2] ? parseInt(m[2], 10) : size - 1;
        res.statusCode = 206;
        res.setHeader("content-range", `bytes ${start}-${end}/${size}`);
        res.setHeader("content-length", String(end - start + 1));
        fs.createReadStream(p, { start, end }).pipe(res);
      } else {
        res.statusCode = 200;
        res.setHeader("content-length", String(size));
        fs.createReadStream(p).pipe(res);
      }
    } catch (e) { res.statusCode = 500; res.end(String(e)); }
  });
  return new Promise((resolve) => {
    server.listen(0, "127.0.0.1", () => {
      const addr = server.address() as any;
      resolve({ url: `http://127.0.0.1:${addr.port}`, close: () => server.close() });
    });
  });
}

/** The value surface a real consumer reads, through openLstar(store): axis labels, a categorical, a
 * dense embedding, a full sparse matrix, and a byte-range sparse column (exercises getRange). */
async function readValues(store: any) {
  const ds = await openLstar(store);
  const pca = await ds.fieldDense("pca");
  const sp = await ds.fieldSparse("data");
  const col0 = await ds.cscColumn("data", 0);
  // batched multi-column read: a scattered subset in INPUT order with a DUPLICATE — each entry must equal
  // the single-read cscColumn (asserted on the reference below), and be identical across backends.
  const ncols = (sp.shape as number[])[1] ?? 1;
  const pick = ncols >= 2 ? [0, 1, 0] : [0, 0];
  const batch = await ds.cscColumns("data", pick);
  return {
    cells: await ds.axisLabels("cells"),
    leiden: await ds.fieldStrings("leiden"),
    pca: Array.from(pca.data as any), pcaShape: pca.shape,
    spData: Array.from(sp.data as any), spIndices: Array.from(sp.indices as any),
    spIndptr: Array.from(sp.indptr as any), spShape: sp.shape,
    col0rows: Array.from(col0.rows as any), col0vals: Array.from(col0.vals as any),
    batchPick: pick,
    batchRows: batch.cols.map((c: any) => Array.from(c.rows as any)),
    batchVals: batch.cols.map((c: any) => Array.from(c.vals as any)),
  };
}

const ref = await readValues(new NodeFSStore(dir));
// the reference must be non-degenerate, or "equal to a collapsed store" would pass vacuously
assert(new Set(ref.leiden).size > 1, "reference leiden should have >1 distinct value (fixture too trivial?)");
assert(ref.pca.some((x: number) => x !== 0), "reference pca should not be all-zeros");
assert(ref.spData.some((x: number) => x !== 0), "reference sparse data should not be all-zeros");
// cscColumns must equal N× cscColumn (same fixture, single-read path) — the batched read is only a
// coalescing optimization, so every returned column is byte-identical to the per-column accessor.
{
  const ds = await openLstar(new NodeFSStore(dir));
  for (let i = 0; i < ref.batchPick.length; i++) {
    const one = await ds.cscColumn("data", ref.batchPick[i]);
    assert.deepEqual(ref.batchRows[i], Array.from(one.rows as any), `cscColumns[${i}].rows != cscColumn(${ref.batchPick[i]})`);
    assert.deepEqual(ref.batchVals[i], Array.from(one.vals as any), `cscColumns[${i}].vals != cscColumn(${ref.batchPick[i]})`);
  }
  console.log("  [js] cscColumns == N× cscColumn (batched read parity)");
}

const srv = await serve(dir, zip);
try {
  const backends: Array<[string, any]> = [
    ["FS zip   (nodeFileSource)", await ZipStore.open(nodeFileSource(zip), zip)],
    ["HTTP dir (HttpStore)     ", new HttpStore(srv.url + "/")],
    ["HTTP zip (httpZipSource) ", await ZipStore.open(httpZipSource(srv.url + "/store.zip"))],
  ];
  for (const [label, store] of backends) {
    const got = await readValues(store);
    assert.deepEqual(got, ref, `VALUE MISMATCH through openLstar: ${label} != FS dir (silent data collapse?)`);
    console.log(`  [js] values-through-reader equal: ${label} == FS dir`);
  }
  // belt-and-suspenders: ZipStore must tolerate a leading-slash key (the reader forms them). Chunk keys
  // differ by format: v3 (default) `.../values/c/0`, v2 `.../values/0` — probe whichever the store has.
  const zs = await ZipStore.open(nodeFileSource(zip), zip);
  const base = "fields/leiden/values";
  const key = (await zs.get(`${base}/c/0`)) ? `${base}/c/0` : `${base}/0`;
  const a = await zs.get(key), b = await zs.get("/" + key);
  assert(a && b && a.length === b.length && a.every((x, i) => x === b[i]),
    "ZipStore.get('/' + key) must equal ZipStore.get(key) — leading-slash tolerance");
  console.log("  [js] ZipStore tolerates a leading-slash key");
} finally {
  srv.close();
}
console.log("all store-backend value-parity tests passed");
