// writer.test.ts — write an L* store from JS, read it back (JS), and exercise addToStore. The
// cross-language check (JS-write -> Python-read) runs in test/writer_crossread.py against the store
// this test leaves under test/.tmp-writer.
import assert from "node:assert";
import * as fs from "node:fs";
import * as os from "node:os";
import * as path from "node:path";

import { NodeFSStore } from "../core/node-store.ts";
import { writeStore, addToStore } from "../core/writer.ts";
import { openLstar } from "../core/reader.ts";

const OUT = process.env.WRITER_OUT || path.join(os.tmpdir(), "lstar-writer-test.lstar.zarr");

async function main() {
  fs.rmSync(OUT, { recursive: true, force: true });
  const store = new NodeFSStore(OUT);

  // a small CSC counts matrix (3 cells x 4 genes), an embedding, and a label
  const counts = {            // column-major (gene) CSC of a 3x4 matrix
    data: new Int32Array([1, 2, 3, 4, 5]),
    indices: new Int32Array([0, 2, 1, 0, 2]),  // cell rows
    indptr: new Int32Array([0, 2, 3, 4, 5]),   // per-gene column pointers (4 cols + 1)
  };
  await writeStore(store, {
    kind: "sample", profiles: ["test@0.1"],
    axes: {
      cells: { labels: ["c0", "c1", "c2"], role: "observation" },
      genes: { labels: ["g0", "g1", "g2", "g3"], role: "feature" },
      umap: { labels: ["umap0", "umap1"], origin: "derived", role: "coordinate" },
    },
    fields: {
      counts: { role: "measure", span: ["cells", "genes"], encoding: "csc", state: "raw",
                shape: [3, 4], data: counts.data, indices: counts.indices, indptr: counts.indptr },
      umap: { role: "embedding", span: ["cells", "umap"], encoding: "dense",
              shape: [3, 2], data: new Float32Array([0, 0, 1, 1, 2, 0.5]) },
      leiden: { role: "label", span: ["cells"], encoding: "utf8", values: ["a", "a", "b"] },
    },
  });

  // read it back with the JS reader
  const ds = await openLstar(store);
  assert.deepStrictEqual(ds.axisNames().sort(), ["cells", "genes", "umap"]);
  assert.deepStrictEqual(ds.fieldNames().sort(), ["counts", "leiden", "umap"]);
  assert.deepStrictEqual(await ds.axisLabels("genes"), ["g0", "g1", "g2", "g3"]);
  assert.deepStrictEqual(await ds.fieldStrings("leiden"), ["a", "a", "b"]);
  const emb = await ds.fieldDense("umap");
  assert.deepStrictEqual([...emb.data], [0, 0, 1, 1, 2, 0.5]);
  assert.deepStrictEqual(emb.shape, [3, 2]);
  const sp = await ds.fieldSparse("counts");
  assert.deepStrictEqual([...sp.data], [1, 2, 3, 4, 5]);
  assert.deepStrictEqual([...sp.indptr], [0, 2, 3, 4, 5]);
  assert.deepStrictEqual(sp.shape, [3, 4]);
  const col = await ds.cscColumn("counts", 0);             // gene 0 -> cells {0:1, 2:2}
  assert.deepStrictEqual([...col.rows], [0, 2]);
  assert.deepStrictEqual([...col.vals], [1, 2]);

  // addToStore: append a derived field + axis and update the manifest
  await addToStore(store, {
    axes: { groups_leiden: { labels: ["a", "b"], origin: "derived", role: "feature" } },
    fields: { od_score: { role: "measure", span: ["genes"], encoding: "dense", shape: [4],
                          data: new Float32Array([0.1, 0.2, 0.3, 0.4]) } },
    profiles: ["viewer@0.1"],
  });
  const ds2 = await openLstar(store);
  assert.ok(ds2.fieldNames().includes("od_score") && ds2.axisNames().includes("groups_leiden"));
  assert.ok(ds2.profiles.includes("viewer@0.1") && ds2.profiles.includes("test@0.1"));
  assert.deepStrictEqual([...(await ds2.fieldDense("od_score")).data], [0.1, 0.2, Math.fround(0.3), 0.4].map(Math.fround));

  console.log("writer.test: OK ->", OUT);
}

main().catch((e) => { console.error(e); process.exit(1); });
