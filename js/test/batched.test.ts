// Stage 5: csrRows batched-row reads. Returns a re-based CSR submatrix for a cell selection, coalescing
// rows into a few merged byte-range requests (not one per row). Verifies correctness for contiguous,
// scattered, duplicated, and empty selections, and that coalescing reduces the request count.
import * as fs from "node:fs";
import * as os from "node:os";
import * as path from "node:path";

import { NodeFSStore } from "../core/node-store.ts";
import { writeStore } from "../core/writer.ts";
import { openLstar } from "../core/reader.ts";

let fail = 0;
const check = (name: string, ok: boolean) => { console.log(`  ${ok ? "OK" : "FAIL"}  ${name}`); if (!ok) fail++; };

class RangeCountingStore {
  ranges = 0; inner: any;
  constructor(inner: any) { this.inner = inner; }
  get(key: string, opts?: unknown) { return this.inner.get(key, opts); }
  async getRange(key: string, s: number, e: number) { this.ranges++; return this.inner.getRange(key, s, e); }
  set(k: string, v: Uint8Array) { return this.inner.set(k, v); }
  delete(k: string) { return this.inner.delete(k); }
}

const NR = 8, NC = 6;
const dense = [
  [0, 2, 0, 0, 1, 0],
  [3, 0, 0, 4, 0, 0],
  [0, 0, 5, 0, 0, 6],
  [0, 0, 0, 0, 0, 0],   // empty row
  [7, 8, 0, 0, 9, 0],
  [0, 0, 1, 2, 0, 3],
  [4, 0, 0, 0, 0, 0],
  [0, 5, 0, 6, 0, 7],
];
function toCsr() {
  const data: number[] = [], indices: number[] = [], indptr: number[] = [0];
  for (let r = 0; r < NR; r++) { for (let c = 0; c < NC; c++) if (dense[r][c]) { data.push(dense[r][c]); indices.push(c); } indptr.push(data.length); }
  return { data: Float64Array.from(data), indices: Int32Array.from(indices), indptr: Int32Array.from(indptr) };
}
const csr = toCsr();

// reconstruct a dense row from a csrRows result (i-th requested row)
function rowFromSub(sub: any, i: number): number[] {
  const out = new Array(NC).fill(0);
  for (let k = sub.indptr[i]; k < sub.indptr[i + 1]; k++) out[Number(sub.indices[k])] = Number(sub.data[k]);
  return out;
}
const matches = (sub: any, rows: number[]) =>
  rows.every((r, i) => rowFromSub(sub, i).join() === dense[r].join());

const dir = fs.mkdtempSync(path.join(os.tmpdir(), "lstar-batched-"));
await writeStore(new NodeFSStore(dir), {
  kind: "sample",
  axes: { cells: { labels: Array.from({ length: NR }, (_, i) => "c" + i), role: "observation" },
          genes: { labels: Array.from({ length: NC }, (_, i) => "g" + i), role: "feature" } },
  fields: { counts_cm: { role: "measure", span: ["cells", "genes"], encoding: "csr", state: "raw",
                         shape: [NR, NC], data: csr.data, indices: csr.indices, indptr: csr.indptr } },
} as any);

const ds = await openLstar(new NodeFSStore(dir));

// contiguous, scattered, including-empty-row, duplicates, single, empty
const sels: number[][] = [[1, 2, 3], [0, 4, 7], [3], [5, 5, 6], [0], []];
for (const sel of sels) {
  const sub = await ds.csrRows("counts_cm", sel);
  check(`csrRows [${sel.join(",")}] correct`, matches(sub, sel) && sub.rows.join() === sel.join() && sub.indptr.length === sel.length + 1);
}

// equivalence with per-row csrRow for a scattered selection
{
  const sel = [7, 2, 5];
  const sub = await ds.csrRows("counts_cm", sel);
  let ok = true;
  for (let i = 0; i < sel.length; i++) {
    const single = await ds.csrRow("counts_cm", sel[i]);
    const a = sub.indptr[i], b = sub.indptr[i + 1];
    ok = ok && Array.from(sub.indices.slice(a, b)).map(Number).join() === Array.from(single.cols).map(Number).join() &&
               Array.from(sub.data.slice(a, b)).map(Number).join() === Array.from(single.vals).map(Number).join();
  }
  check("csrRows == per-row csrRow (scattered)", ok);
}

// coalescing: a contiguous run of rows uses far fewer getRange calls than one-per-row
{
  const rc = new RangeCountingStore(new NodeFSStore(dir));
  const ds2 = await openLstar(rc);
  await ds2.csrRows("counts_cm", [1, 2, 3, 4, 5]);  // contiguous run -> 1 indptr + 1 data + 1 indices
  check("coalesced contiguous run uses <= 3 getRange (not 1/row)", rc.ranges <= 3);
}

console.log(fail === 0 ? "\nbatched OK" : `\nbatched FAIL: ${fail}`);
process.exit(fail === 0 ? 0 : 1);
