// Stage 4: per-field compression. One writeStore call can keep the gene-major copy raw single-chunk
// (byte-range fast path) while gzip-compressing the cell-major copy (sequential bulk reads) -- via
// FieldSpec.write. Verifies the on-disk .zarray compressors and that both copies read back correctly
// (the compressed one exercises zarrita's gzip decode on the fallback path).
import * as fs from "node:fs";
import * as os from "node:os";
import * as path from "node:path";
import * as zlib from "node:zlib";

import { NodeFSStore } from "../core/node-store.ts";
import { writeStore } from "../core/writer.ts";
import { openLstar } from "../core/reader.ts";

let fail = 0;
const check = (name: string, ok: boolean) => { console.log(`  ${ok ? "OK" : "FAIL"}  ${name}`); if (!ok) fail++; };
const arr = (a: any) => Array.from(a as any).map(Number);

const gzip = { id: "gzip" as const, level: 6, compress: (raw: Uint8Array) => new Uint8Array(zlib.gzipSync(raw, { level: 6 })) };

const NR = 4, NC = 5;
const dense = [
  [0, 3, 0, 0, 7],
  [1, 0, 0, 5, 0],
  [0, 0, 0, 0, 9],
  [2, 4, 0, 6, 0],
];
function toCsc() {
  const data: number[] = [], indices: number[] = [], indptr: number[] = [0];
  for (let c = 0; c < NC; c++) { for (let r = 0; r < NR; r++) if (dense[r][c]) { data.push(dense[r][c]); indices.push(r); } indptr.push(data.length); }
  return { data: Float64Array.from(data), indices: Int32Array.from(indices), indptr: Int32Array.from(indptr) };
}
function toCsr() {
  const data: number[] = [], indices: number[] = [], indptr: number[] = [0];
  for (let r = 0; r < NR; r++) { for (let c = 0; c < NC; c++) if (dense[r][c]) { data.push(dense[r][c]); indices.push(c); } indptr.push(data.length); }
  return { data: Float64Array.from(data), indices: Int32Array.from(indices), indptr: Int32Array.from(indptr) };
}
const csc = toCsc(), csr = toCsr();

const dir = fs.mkdtempSync(path.join(os.tmpdir(), "lstar-comp-"));
await writeStore(new NodeFSStore(dir), {
  kind: "sample",
  axes: { cells: { labels: ["c0", "c1", "c2", "c3"], role: "observation" },
          genes: { labels: ["g0", "g1", "g2", "g3", "g4"], role: "feature" } },
  fields: {
    // gene-major: raw single chunk (no per-field write) -> byte-range fast path
    counts: { role: "measure", span: ["cells", "genes"], encoding: "csc", state: "raw",
              shape: [NR, NC], data: csc.data, indices: csc.indices, indptr: csc.indptr },
    // cell-major: gzip-compressed (per-field override)
    counts_cm: { role: "measure", span: ["cells", "genes"], encoding: "csr", state: "raw",
                 shape: [NR, NC], data: csr.data, indices: csr.indices, indptr: csr.indptr,
                 write: { compressor: gzip } },
  },
} as any, undefined, "v2");   // this test inspects v2 .zarray per-field compressors

// ---- on-disk compressors are asymmetric ----
const za = (p: string) => JSON.parse(fs.readFileSync(path.join(dir, p, ".zarray"), "utf8"));
check("gene-major counts/data is uncompressed", za("fields/counts/data").compressor === null);
check("cell-major counts_cm/data is gzip", za("fields/counts_cm/data")?.compressor?.id === "gzip");

// ---- both copies read back correctly ----
const ds = await openLstar(new NodeFSStore(dir));
// gene-major via fast path
let ok = true;
for (let c = 0; c < NC; c++) {
  const { rows, vals } = await ds.cscColumn("counts", c);
  const ref = dense.map((row, r) => [r, row[c]]).filter(([, v]) => v);
  ok = ok && arr(rows).join() === ref.map(([r]) => r).join() && arr(vals).join() === ref.map(([, v]) => v).join();
}
check("raw gene-major reads correctly (fast path)", ok);

// cell-major via fallback -> exercises zarrita gzip decode
let rok = true;
for (let r = 0; r < NR; r++) {
  const { cols, vals } = await ds.csrRow("counts_cm", r);
  const ref = dense[r].map((v, c) => [c, v]).filter(([, v]) => v);
  rok = rok && arr(cols).join() === ref.map(([c]) => c).join() && arr(vals).join() === ref.map(([, v]) => v).join();
}
check("gzip cell-major reads correctly (zarrita decode)", rok);

// fieldSparse on the compressed field decodes the whole thing
const sp = await ds.fieldSparse("counts_cm");
check("fieldSparse on gzip field decodes", arr(sp.data).join() === arr(csr.data).join() && arr(sp.indptr).join() === arr(csr.indptr).join());

console.log(fail === 0 ? "\ncompression OK" : `\ncompression FAIL: ${fail}`);
process.exit(fail === 0 ? 0 : 1);
