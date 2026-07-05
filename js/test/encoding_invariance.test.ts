// Encoding-invariance for the LIVE-viewer measure compute: the SAME matrix, written as dense / csc / csr,
// must give identical colStats (HVG ranking) and subsampleDE (interactive DE). This is the coverage that
// was missing when JS read the measure via a sparse-hardcoded path — a dense primary measure (SCE
// logcounts / scaled AnnData X) threw NotFoundError and viewer compute silently broke. Metamorphic: no
// Python reference needed — encoding must not change the answer, and a direct dense reference pins correctness.
import * as fs from "node:fs";
import * as os from "node:os";
import * as path from "node:path";

import { openLstar } from "../core/reader.ts";
import { NodeFSStore } from "../core/node-store.ts";
import { writeStore } from "../core/writer.ts";
import { LstarView } from "../core/view.ts";

let fail = 0;
const check = (name: string, ok: boolean) => { console.log(`  ${ok ? "OK" : "FAIL"}  ${name}`); if (!ok) fail++; };
const approx = (a: ArrayLike<number>, b: ArrayLike<number>, rtol = 1e-9, atol = 1e-9): boolean => {
  if (a.length !== b.length) return false;
  for (let i = 0; i < a.length; i++)
    if (Math.abs(Number(a[i]) - Number(b[i])) > atol + rtol * Math.abs(Number(b[i]))) return false;
  return true;
};

const NR = 40, NC = 12;
// deterministic (no Math.random): ~40% zeros, cell-parity-dependent so A/B DE has signal
const val = (i: number, j: number): number => {
  const h = ((i * 73856093) ^ (j * 19349663)) >>> 0;
  return h % 5 < 2 ? 0 : (h % 13) + (i % 2 ? 1 : 0);
};
const dense = new Float64Array(NR * NC);
for (let i = 0; i < NR; i++) for (let j = 0; j < NC; j++) dense[i * NC + j] = val(i, j);

function denseToCsc(d: ArrayLike<number>, nr: number, nc: number) {
  const indptr = new Int32Array(nc + 1);
  for (let c = 0; c < nc; c++) { let k = 0; for (let r = 0; r < nr; r++) if (d[r * nc + c] !== 0) k++; indptr[c + 1] = indptr[c] + k; }
  const data = new Float64Array(indptr[nc]), indices = new Int32Array(indptr[nc]); let w = 0;
  for (let c = 0; c < nc; c++) for (let r = 0; r < nr; r++) { const v = d[r * nc + c]; if (v !== 0) { data[w] = v; indices[w] = r; w++; } }
  return { data, indices, indptr };
}
function denseToCsr(d: ArrayLike<number>, nr: number, nc: number) {
  const indptr = new Int32Array(nr + 1);
  for (let r = 0; r < nr; r++) { let k = 0; for (let c = 0; c < nc; c++) if (d[r * nc + c] !== 0) k++; indptr[r + 1] = indptr[r] + k; }
  const data = new Float64Array(indptr[nr]), indices = new Int32Array(indptr[nr]); let w = 0;
  for (let r = 0; r < nr; r++) for (let c = 0; c < nc; c++) { const v = d[r * nc + c]; if (v !== 0) { data[w] = v; indices[w] = c; w++; } }
  return { data, indices, indptr };
}

const csc = denseToCsc(dense, NR, NC);
const csr = denseToCsr(dense, NR, NC);
const fieldFor: Record<string, any> = {
  dense: { role: "measure", span: ["cells", "genes"], encoding: "dense", state: "raw", shape: [NR, NC], data: dense },
  csc:   { role: "measure", span: ["cells", "genes"], encoding: "csc",   state: "raw", shape: [NR, NC], ...csc },
  csr:   { role: "measure", span: ["cells", "genes"], encoding: "csr",   state: "raw", shape: [NR, NC], ...csr },
};

const tmp = fs.mkdtempSync(path.join(os.tmpdir(), "lstar-enc-"));
async function build(enc: string): Promise<LstarView> {
  const out = path.join(tmp, enc + ".lstar.zarr");
  await writeStore(new NodeFSStore(out), {
    kind: "sample", profiles: [], dropped: [],
    axes: {
      cells: { labels: Array.from({ length: NR }, (_, i) => `c${i}`), role: "observation" },
      genes: { labels: Array.from({ length: NC }, (_, j) => `g${j}`), role: "feature" },
    },
    fields: { X: fieldFor[enc] },
  });
  return new LstarView(await openLstar(new NodeFSStore(out)));
}

const A = Array.from({ length: NR / 2 }, (_, i) => i * 2);        // even cells
const B = Array.from({ length: NR / 2 }, (_, i) => i * 2 + 1);    // odd cells
type Res = { mean: Float64Array; var: Float64Array; nnz: Float64Array; meanA: Float64Array; meanB: Float64Array };
async function run(v: LstarView): Promise<Res> {
  const cs = await v.colStats({ field: "X", lognorm: true });
  const de = await v.subsampleDE(A, B, { field: "X", lognorm: true });
  const meanA = new Float64Array(NC), meanB = new Float64Array(NC);
  for (const r of de) { meanA[r.gene] = r.meanA; meanB[r.gene] = r.meanB; }
  return { mean: cs.mean as any, var: cs.var as any, nnz: cs.nnz as any, meanA, meanB };
}

const R: Record<string, Res> = {};
for (const enc of ["dense", "csc", "csr"]) R[enc] = await run(await build(enc));

// (1) encoding-invariance: dense == csc == csr on every statistic (the guard #101 was missing)
for (const enc of ["csc", "csr"]) {
  check(`colStats mean  ${enc}==dense`, approx(R[enc].mean, R.dense.mean));
  check(`colStats var   ${enc}==dense`, approx(R[enc].var, R.dense.var));
  check(`colStats nnz   ${enc}==dense`, approx(R[enc].nnz, R.dense.nnz));
  check(`subsampleDE meanA ${enc}==dense`, approx(R[enc].meanA, R.dense.meanA));
  check(`subsampleDE meanB ${enc}==dense`, approx(R[enc].meanB, R.dense.meanB));
}

// (2) direct correctness reference (guards against all-three-consistently-wrong): per-gene zero-aware
// log1p mean + nnz computed straight from the dense matrix must match colStats.
const refMean = new Float64Array(NC), refNnz = new Float64Array(NC);
for (let j = 0; j < NC; j++) { let s = 0, n = 0; for (let i = 0; i < NR; i++) { const v = dense[i * NC + j]; s += Math.log1p(v); if (v !== 0) n++; } refMean[j] = s / NR; refNnz[j] = n; }
check("colStats mean == direct dense reference", approx(R.csc.mean, refMean, 1e-9, 1e-9));
check("colStats nnz  == direct dense reference", approx(R.csc.nnz, refNnz));

// (3) positive assertions: the compute actually produced non-trivial results (not silently empty)
check("nnz has nonzeros (compute ran)", Array.from(R.dense.nnz).some((x) => x > 0));
check("subsampleDE produced per-gene means", Array.from(R.dense.meanA).some((x) => x > 0));

fs.rmSync(tmp, { recursive: true, force: true });
console.log(fail === 0 ? "\nencoding-invariance OK" : `\nencoding-invariance FAIL: ${fail}`);
process.exit(fail === 0 ? 0 : 1);
