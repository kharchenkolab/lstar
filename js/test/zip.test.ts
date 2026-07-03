// JS single-file `.lstar.zarr.zip` conformance: ZipStore reads a Python-written STORED zip byte-for-byte
// like the directory store (seek-into-zip, incl. the getRange fast path), JS writes a STORED zip Python
// reads back, and the guardrails hold (DEFLATE rejected, ZIP64 read). Driven by conformance/zip_js.sh
// (or the Phase-6 matrix), which supplies the fixtures and cross-reads the JS-written zip in Python.
//
//   usage: zip.test.ts <dir> <pyzip> <outzip> <deflatezip> <zip64zip>
import assert from "node:assert";

import { openLstar } from "../core/reader.ts";
import { NodeFSStore } from "../core/node-store.ts";
import { ZipStore, readZipCentralDir } from "../core/zip.ts";
import { nodeFileSource, openLstarZip, writeStoreZip } from "../core/zip-node.ts";

const [, , dir, pyzip, outzip, deflatezip, zip64zip] = process.argv;

// (1) byte-level read parity: ZipStore.get(key) === NodeFSStore.get(key) for EVERY key, + getRange + open
{
  const idx = await readZipCentralDir(nodeFileSource(pyzip), pyzip);
  const zs = await ZipStore.open(nodeFileSource(pyzip), pyzip);
  const fss = new NodeFSStore(dir);
  let n = 0;
  for (const key of idx.keys()) {
    const a = await zs.get(key), b = await fss.get(key);
    assert(a && b, `absent: ${key}`);
    assert.equal(a.length, b.length, `length differs: ${key}`);
    for (let i = 0; i < a.length; i++) assert.equal(a[i], b[i], `byte ${i} differs: ${key}`);
    n++;
  }
  // getRange fast path: a chunk's head range == the head of its full bytes
  const chunkKey = [...idx.keys()].find((k) => /\/\d+(\.\d+)*$/.test(k));
  if (chunkKey) {
    const full = (await zs.get(chunkKey))!;
    const head = (await zs.getRange(chunkKey, 0, Math.min(8, full.length)))!;
    for (let i = 0; i < head.length; i++) assert.equal(head[i], full[i], `getRange head byte ${i}`);
  }
  const ds = await openLstar(zs);
  assert(ds.axes.size >= 1 && ds.fields.size >= 1, "openLstar(ZipStore) produced an empty dataset");
  console.log(`  [js] read parity: ${n} keys byte-identical (ZipStore == NodeFSStore); getRange + openLstar OK`);
}

// (2) write parity: JS writes a STORED zip that Python reads back (verified in Python by zip_js.sh)
{
  await writeStoreZip(outzip, {
    kind: "sample",
    axes: {
      cells: { labels: ["c0", "c1", "c2", "c3"], origin: "observed" },
      genes: { labels: ["g0", "g1", "g2"], origin: "observed" },
    },
    fields: {
      counts: {
        role: "measure", span: ["cells", "genes"], encoding: "csc", state: "raw",
        data: [1, 2, 3, 4], indices: [0, 1, 2, 3], indptr: [0, 1, 2, 4], shape: [4, 3],
      },
      umap: {
        role: "embedding", span: ["cells", "genes"], encoding: "dense", shape: [4, 3],
        data: [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11],
      },
      leiden: { role: "label", span: ["cells"], encoding: "utf8", values: ["a", "b", "a", "b"] },
    },
  });
  console.log("  [js] wrote", outzip, "(single-file STORED zip)");
}

// (3) guardrail: a DEFLATE-packed store is rejected with an actionable message
{
  let threw = false;
  try { await openLstarZip(deflatezip); }
  catch (e: any) { threw = /stored|deflate/i.test(String(e && e.message)); }
  assert(threw, "JS must reject a DEFLATE-packed .lstar.zarr.zip");
  console.log("  [js] rejects DEFLATE zip (guardrail)");
}

// (4) ZIP64: a ZIP64 STORED archive reads
{
  const ds = await openLstarZip(zip64zip);
  assert(ds.fields.size >= 1, "ZIP64 read produced no fields");
  console.log("  [js] reads a ZIP64 STORED archive");
}

console.log("all JS zip tests passed");
