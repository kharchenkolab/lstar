// reader_parallel.test.ts — a label/categorical field's independent component arrays are read CONCURRENTLY,
// not serially: fieldStrings reads values ‖ values_offsets, fieldCategorical reads codes ‖ categories ‖
// categories_offsets, axisLabels reads labels ‖ labels_offsets. On a hosted store each field is then ONE
// round-trip, not two/three — the fix for "clusters wait on the embedding" (the cluster colorBy paid a serial
// second round-trip on the viewer's first-paint path). Run: node --experimental-strip-types test/reader_parallel.test.ts
import assert from "node:assert";
import * as fs from "node:fs";
import * as os from "node:os";
import * as path from "node:path";

import { NodeFSStore } from "../core/node-store.ts";
import { writeStore } from "../core/writer.ts";
import { openLstar } from "../core/reader.ts";

// A NodeFSStore that records the MAX number of get()s in flight at once, with a small delay so overlapping
// reads actually register as concurrent. max === 1 ⇒ the reader serialized the component reads; ≥ 2 ⇒ parallel.
class CountingStore extends NodeFSStore {
  inFlight = 0; max = 0;
  reset() { this.max = 0; }
  private async track<T>(fn: () => Promise<T>): Promise<T> {
    this.inFlight++; this.max = Math.max(this.max, this.inFlight);
    try { await new Promise((r) => setTimeout(r, 5)); return await fn(); } finally { this.inFlight--; }
  }
  async get(key: string) { return this.track(() => super.get(key)); }
  async getRange(key: string, s: number, e: number) { return this.track(() => super.getRange(key, s, e)); }
}

async function main() {
  const base = fs.mkdtempSync(path.join(os.tmpdir(), "lstar-reader-parallel-"));
  const dir = path.join(base, "s.lstar.zarr");
  await writeStore(new NodeFSStore(dir), {
    kind: "sample", profiles: ["test@0.1"],
    axes: {
      cells: { labels: ["c0", "c1", "c2", "c3"], role: "observation" },
      genes: { labels: ["g0", "g1", "g2"], role: "feature" },
      umap: { labels: ["u0", "u1"], origin: "derived", role: "coordinate" },
    },
    fields: {
      counts: { role: "measure", span: ["cells", "genes"], encoding: "csc", state: "raw", shape: [4, 3],
                data: new Int32Array([1, 2, 3, 4]), indices: new Int32Array([0, 1, 2, 3]), indptr: new Int32Array([0, 2, 3, 4]) },
      umap: { role: "embedding", span: ["cells", "umap"], encoding: "dense", shape: [4, 2], data: new Float32Array([0, 0, 1, 1, 2, 0, 3, 1]) },
      leiden: { role: "label", span: ["cells"], encoding: "utf8", values: ["a", "a", "b", "b"] },
      ct: { role: "label", span: ["cells"], encoding: "categorical", codes: new Int32Array([0, 1, 0, 1]), categories: ["T", "B"] },
    },
  });

  const store = new CountingStore(dir);
  const ds = await openLstar(store);

  store.reset(); const s = await ds.fieldStrings("leiden");
  assert.deepStrictEqual(s, ["a", "a", "b", "b"], "values still decode correctly");
  assert.ok(store.max >= 2, `fieldStrings must read values ‖ values_offsets concurrently, got max=${store.max}`);

  store.reset(); const c = await ds.fieldCategorical("ct");
  assert.deepStrictEqual(c.categories, ["T", "B"], "categories still decode (stored order)");
  assert.deepStrictEqual([...c.codes], [0, 1, 0, 1], "codes still decode");
  assert.ok(store.max >= 2, `fieldCategorical must read codes ‖ categories ‖ offsets concurrently, got max=${store.max}`);

  store.reset(); const g = await ds.axisLabels("genes");
  assert.deepStrictEqual(g, ["g0", "g1", "g2"], "axis labels still decode");
  assert.ok(store.max >= 2, `axisLabels must read labels ‖ labels_offsets concurrently, got max=${store.max}`);

  fs.rmSync(base, { recursive: true, force: true });
  console.log("reader_parallel.test: OK (component reads concurrent + values correct)");
}

main().catch((e) => { console.error(e); process.exit(1); });
