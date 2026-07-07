// Guard for the reader-WASM instance lifecycle. The WASM module is a process-wide SINGLETON shared by
// every dataset (opening a fresh createLstarIO() per dataset accumulates WASM heaps and aborts a long
// session), while each dataset's native `Reader` is per-source and freed by dispose(). This asserts the
// contract: many datasets coexist + open cleanly; dispose() frees one without disturbing the others or
// the shared module; a dataset must not be read after dispose.
import assert from "node:assert";
import * as fs from "node:fs";
import * as os from "node:os";
import * as path from "node:path";

import { NodeFSStore } from "../core/node-store.ts";
import { writeStore } from "../core/writer.ts";
import { openLstar } from "../core/reader.ts";

const dir = path.join(os.tmpdir(), "lstar-dispose.lstar.zarr");
fs.rmSync(dir, { recursive: true, force: true });
await writeStore(new NodeFSStore(dir), {
  kind: "sample", profiles: ["dispose@0.1"],
  axes: { cells: { labels: ["c0", "c1", "c2"], role: "observation" },
          genes: { labels: ["g0", "g1", "g2", "g3"], role: "feature" } },
  fields: { counts: { role: "measure", span: ["cells", "genes"], encoding: "csc", state: "raw",
                      shape: [3, 4], data: new Int32Array([1, 2, 3, 4, 5]),
                      indices: new Int32Array([0, 2, 1, 0, 2]), indptr: new Int32Array([0, 2, 3, 4, 5]) } },
});   // v3 default

const store = new NodeFSStore(dir);

// two datasets over the same store coexist (shared module, isolated Readers) and both read correctly
const ds1 = await openLstar(store);
const ds2 = await openLstar(store);
assert.deepStrictEqual([...(await ds1.fieldSparse("counts")).data], [1, 2, 3, 4, 5]);
assert.deepStrictEqual([...(await ds2.fieldSparse("counts")).data], [1, 2, 3, 4, 5]);
console.log("  [js] two datasets on one shared WASM module read independently");

// dispose ds1: reading it must now fail (freed Reader); ds2 is unaffected and the module survives
ds1.dispose();
let threw = false;
try { await ds1.fieldSparse("counts"); } catch { threw = true; }
assert.ok(threw, "reading a disposed dataset must throw (its Reader was freed)");
assert.deepStrictEqual([...(await ds2.fieldSparse("counts")).data], [1, 2, 3, 4, 5],
  "disposing one dataset must not disturb another");
ds1.dispose();   // idempotent
console.log("  [js] dispose() frees one dataset; others + the shared module survive; idempotent");

// a NEW dataset opens fine after disposing an earlier one (module was not torn down)
const ds3 = await openLstar(store);
assert.deepStrictEqual([...(await ds3.fieldSparse("counts")).data], [1, 2, 3, 4, 5]);
ds2.dispose(); ds3.dispose();

// churn many datasets sequentially — the singleton bounds this to one WASM heap (per-open modules would
// accumulate and, per the pagoda3 report, abort after several).
for (let i = 0; i < 24; i++) {
  const d = await openLstar(store);
  assert.strictEqual((await d.fieldSparse("counts")).data.length, 5);
  d.dispose();
}
console.log("  [js] opened+disposed 24 datasets on one shared module without abort");
console.log("reader dispose/lifecycle guard passed");
