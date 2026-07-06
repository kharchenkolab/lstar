// Stage 3: byte-range sparse reads. cscColumn/csrRow must (a) return the correct nonzeros, (b) give
// identical results on the byte-range fast path and the zarrita slice fallback, (c) actually take the
// fast path when getRange + consolidated metadata are present, and (d) fall back on a multi-chunk store.
import * as fs from "node:fs";
import * as os from "node:os";
import * as path from "node:path";

import { NodeFSStore } from "../core/node-store.ts";
import { writeStore } from "../core/writer.ts";
import { openLstar } from "../core/reader.ts";

let fail = 0;
const check = (name: string, ok: boolean) => { console.log(`  ${ok ? "OK" : "FAIL"}  ${name}`); if (!ok) fail++; };
const arr = (a: any) => Array.from(a as any).map(Number);

// Counts getRange calls + chunk gets, to prove which path was taken. A chunk key ends with "/<int>".
class RangeCountingStore {
  ranges = 0; chunkGets = 0; inner: any;
  constructor(inner: any) { this.inner = inner; }
  async get(key: string, opts?: unknown) {
    const norm = key[0] === "/" ? key.slice(1) : key;
    if (/\/\d+(?:\.\d+)*$/.test(norm)) this.chunkGets++;
    return this.inner.get(key, opts);
  }
  async getRange(key: string, s: number, e: number) { this.ranges++; return this.inner.getRange(key, s, e); }
  set(k: string, v: Uint8Array) { return this.inner.set(k, v); }
  delete(k: string) { return this.inner.delete(k); }
}
// A store WITHOUT getRange -> forces the slice fallback.
class NoRangeStore {
  inner: any;
  constructor(inner: any) { this.inner = inner; }
  get(key: string, opts?: unknown) { return this.inner.get(key, opts); }
  set(k: string, v: Uint8Array) { return this.inner.set(k, v); }
  delete(k: string) { return this.inner.delete(k); }
}

// dense (row-major, nrows x ncols) -> CSC (column-major) arrays
function denseToCsc(dense: number[][], nrows: number, ncols: number) {
  const data: number[] = [], indices: number[] = [], indptr: number[] = [0];
  for (let c = 0; c < ncols; c++) {
    for (let r = 0; r < nrows; r++) if (dense[r][c] !== 0) { data.push(dense[r][c]); indices.push(r); }
    indptr.push(data.length);
  }
  return { data: Float64Array.from(data), indices: Int32Array.from(indices), indptr: Int32Array.from(indptr) };
}
function denseToCsr(dense: number[][], nrows: number, ncols: number) {
  const data: number[] = [], indices: number[] = [], indptr: number[] = [0];
  for (let r = 0; r < nrows; r++) {
    for (let c = 0; c < ncols; c++) if (dense[r][c] !== 0) { data.push(dense[r][c]); indices.push(c); }
    indptr.push(data.length);
  }
  return { data: Float64Array.from(data), indices: Int32Array.from(indices), indptr: Int32Array.from(indptr) };
}

const NR = 4, NC = 5;
const dense = [
  [0, 3, 0, 0, 7],
  [1, 0, 0, 5, 0],
  [0, 0, 0, 0, 9],   // col 2 is all-zero (empty column/row edge case)
  [2, 4, 0, 6, 0],
];
const csc = denseToCsc(dense, NR, NC);
const csr = denseToCsr(dense, NR, NC);

async function buildStore(dir: string, opts?: any) {
  await writeStore(new NodeFSStore(dir), {
    kind: "sample",
    axes: { cells: { labels: ["c0", "c1", "c2", "c3"], role: "observation" },
            genes: { labels: ["g0", "g1", "g2", "g3", "g4"], role: "feature" } },
    fields: {
      counts: { role: "measure", span: ["cells", "genes"], encoding: "csc", state: "raw",
                shape: [NR, NC], data: csc.data, indices: csc.indices, indptr: csc.indptr },
      counts_cm: { role: "measure", span: ["cells", "genes"], encoding: "csr", state: "raw",
                   shape: [NR, NC], data: csr.data, indices: csr.indices, indptr: csr.indptr },
    },
  } as any, opts);
}

// reference: nonzeros of a dense column / row
const refCol = (c: number) => dense.map((row, r) => [r, row[c]]).filter(([, v]) => v !== 0);
const refRow = (r: number) => dense[r].map((v, c) => [c, v]).filter(([, v]) => v !== 0);

// ---- fast path (getRange + consolidated) ----
const dir = fs.mkdtempSync(path.join(os.tmpdir(), "lstar-sparse-"));
await buildStore(dir);
{
  const rc = new RangeCountingStore(new NodeFSStore(dir));
  const ds = await openLstar(rc);
  let ok = true;
  for (let c = 0; c < NC; c++) {
    const { rows, vals } = await ds.cscColumn("counts", c);
    const ref = refCol(c);
    ok = ok && arr(rows).join() === ref.map(([r]) => r).join() && arr(vals).join() === ref.map(([, v]) => v).join();
  }
  check("cscColumn fast-path correct (all columns incl. empty)", ok);
  check("cscColumn took the byte-range path (getRange called, no chunk gets)", rc.ranges > 0 && rc.chunkGets === 0);

  let rok = true;
  for (let r = 0; r < NR; r++) {
    const { cols, vals } = await ds.csrRow("counts_cm", r);
    const ref = refRow(r);
    rok = rok && arr(cols).join() === ref.map(([c]) => c).join() && arr(vals).join() === ref.map(([, v]) => v).join();
  }
  check("csrRow fast-path correct (all rows)", rok);
}

// ---- equivalence: fast path == slice fallback (no-getRange store) ----
{
  const fast = await openLstar(new NodeFSStore(dir));
  const slow = await openLstar(new NoRangeStore(new NodeFSStore(dir)));
  let same = true;
  for (let c = 0; c < NC; c++) {
    const a = await fast.cscColumn("counts", c), b = await slow.cscColumn("counts", c);
    same = same && arr(a.rows).join() === arr(b.rows).join() && arr(a.vals).join() === arr(b.vals).join();
  }
  check("cscColumn fast == slow (equivalence)", same);
}

// ---- multi-chunk store: the byte-range fast path still applies per-chunk (chunk index = floor(lo/
// chunkLen)); a column whose span CROSSES a chunk boundary falls back to a whole-array read. Correct
// either way. ----
{
  const mdir = fs.mkdtempSync(path.join(os.tmpdir(), "lstar-sparse-mc-"));
  await buildStore(mdir, { chunkElems: 2 });   // uncompressed but multi-chunk along axis 0
  const rc = new RangeCountingStore(new NodeFSStore(mdir));
  const ds = await openLstar(rc);
  let ok = true;
  for (let c = 0; c < NC; c++) {
    const { rows, vals } = await ds.cscColumn("counts", c);
    const ref = refCol(c);
    ok = ok && arr(rows).join() === ref.map(([r]) => r).join() && arr(vals).join() === ref.map(([, v]) => v).join();
  }
  check("cscColumn correct on multi-chunk store", ok);
  check("multi-chunk still takes the byte-range fast path where a column fits a chunk", rc.ranges > 0);
}

console.log(fail === 0 ? "\nsparse OK" : `\nsparse FAIL: ${fail}`);
process.exit(fail === 0 ? 0 : 1);
