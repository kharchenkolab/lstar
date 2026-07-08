// Regression: a DENSE int64/uint64 field (a BigInt typed array in JS) must read through fieldAsCsc. The
// dense->CSC path assigned a BigInt into a Float64Array (throws "Cannot convert a BigInt value to a
// number") and compared it to the number 0 (`0n !== 0` is true, so zeros miscount). Caught on real Visium
// spatial data (int64 `spatial` coords / `in_tissue`) by conformance/v3_corpus.sh; guarded here in CI
// (which has no corpus). Fixed by coercing with Number() in denseToCscArrays.
import assert from "node:assert";
import * as fs from "node:fs";
import * as os from "node:os";
import * as path from "node:path";

import { NodeFSStore } from "../core/node-store.ts";
import { writeStore } from "../core/writer.ts";
import { openLstar } from "../core/reader.ts";

const dir = fs.mkdtempSync(path.join(os.tmpdir(), "lstar-bigint-"));
// 3 cells x 4 genes, DENSE, int64 (-> BigInt64Array on the JS read). Column-major nonzeros give a known CSC.
const data = new BigInt64Array([
  0n, 2n, 0n, 5n,   // cell 0
  1n, 0n, 0n, 0n,   // cell 1
  0n, 0n, 3n, 0n,   // cell 2
]);
await writeStore(new NodeFSStore(dir), {
  kind: "sample",
  axes: { cells: { labels: ["c0", "c1", "c2"], role: "observation" },
          genes: { labels: ["g0", "g1", "g2", "g3"], role: "feature" } },
  fields: { m: { role: "measure", span: ["cells", "genes"], encoding: "dense", state: "raw", shape: [3, 4], data } },
} as any);

const ds = await openLstar(new NodeFSStore(dir));
const csc = await ds.fieldAsCsc("m");   // dense -> CSC: exercises denseToCscArrays on a BigInt source
// gene 0 -> {cell1:1}; gene 1 -> {cell0:2}; gene 2 -> {cell2:3}; gene 3 -> {cell0:5}. nnz=4.
assert.strictEqual(csc.indptr.length, 5, "indptr length");
assert.deepStrictEqual([...csc.indptr].map(Number), [0, 1, 2, 3, 4], "per-gene column pointers");
assert.deepStrictEqual([...csc.data].map(Number), [1, 2, 3, 5], "nonzero values (BigInt coerced to number)");
assert.deepStrictEqual([...csc.indices].map(Number), [1, 0, 2, 0], "row indices");
console.log("dense int64 (BigInt) -> CSC: OK (values coerced, zeros not miscounted)");
