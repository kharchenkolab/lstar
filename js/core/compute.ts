// The viewer-compute RECIPE, single-sourced. Pure array-in / results-out reductions over the libstar
// WASM kernels — the per-(group,gene) sufficient stats, 1-vs-rest markers, pagoda2 overdispersion, and
// A-vs-B DE that the viewer is built from. Shared by lstar's `extend.ts` (precompute-once prep), lstar's
// `LstarView` (live queries), and pagoda3's browser view, so a fix lands once and the conformance suite
// covers all three. Mirrors what Python/R single-source in `viewer.py`.
//
// Deliberately NO store / IO / render coupling: the caller reads the measure (the whole matrix via
// `ds.fieldAsCsc`, or a subset's rows via `reader.csrRows`) and hands CSC arrays + codes + cell-sets in.
// Orchestration — workers, SharedArrayBuffer, caches, sampling caps, subset *selection*, annotation
// overlays — stays with the caller. The subset fast path is which rows you read, not a concern here.
import createLstarKernels from "../dist/lstar_kernels.mjs";

/** A CSC measure as the kernels expect it: cells × genes, column-major over genes. `data`/`indices`/
 * `indptr` are the CSC triple; `ncells`/`ngenes` are the shape. (`ds.fieldAsCsc(name)` returns the triple
 * + `shape=[ncells,ngenes]` for any on-disk encoding.) */
export interface CscMeasure {
  data: ArrayLike<number>;
  indptr: ArrayLike<number>;
  indices: ArrayLike<number>;
  ncells: number;
  ngenes: number;
}

let _kernels: Promise<any> | null = null;
/** The libstar WASM kernels, memoized per module instance. Pass your own handle to any primitive (e.g. a
 * worker-local instance you already loaded) to skip this. */
export function kernels(): Promise<any> { return (_kernels ??= createLstarKernels()); }

// `Number` as the map fn keeps these safe for BigInt64Array/BigUint64Array (int64 zarr arrays decode to
// BigInt typed arrays, which don't implicitly convert) as well as ordinary numeric arrays. Widening
// float32→float64 is lossless and the kernels compute in double, so wrapping never changes a result.
const toF64 = (a: any): Float64Array => (a instanceof Float64Array ? a : Float64Array.from(a, Number));
const toI32 = (a: any): Int32Array => (a instanceof Int32Array ? a : Int32Array.from(a, Number));

/** Per-gene zero-aware mean / variance / nnz over a CSC measure (the HVG-ranking inputs), via colMeanVar.
 * Computed over `log1p(x)` when `lognorm` (default true), else raw. */
export async function colStats(
  m: { data: ArrayLike<number>; indptr: ArrayLike<number>; ncells: number },
  opts: { lognorm?: boolean } = {}, M?: any,
): Promise<{ mean: Float64Array; var: Float64Array; nnz: Int32Array }> {
  M ??= await kernels();
  return M.colMeanVar(toF64(m.data), toI32(m.indptr), m.ncells, 1, opts.lognorm ?? true);
}

/** Per-gene pagoda2 overdispersion (adjustVariance: LOWESS residual of log(var)~log(mean) + F-test) from
 * PRE-REDUCED per-gene stats — the downstream-math half, layout-independent. A caller that reduced the
 * measure itself (e.g. per-gene mean/var/n_expr over a cell-major subset, so it never formed a gene-major
 * CSC) calls this directly; the reduction stays theirs, the algorithm stays here. `nnz` is per-gene n_expr
 * (the F-test's degrees of freedom). Thin over the WASM kernel, which coerces the arrays. */
export async function overdispersionFromStats(
  mean: ArrayLike<number>, variance: ArrayLike<number>, nnz: ArrayLike<number>, M?: any,
): Promise<Float64Array> {
  M ??= await kernels();
  return M.overdispersion(mean, variance, nnz);
}

/** Per-gene pagoda2 overdispersion over a CSC measure — the reduce (colMeanVar) + downstream
 * ({@link overdispersionFromStats}) together. Global over the whole measure, or over a subset if the caller
 * read only a subset's rows. Callers who reduce cell-major should reduce themselves and call
 * {@link overdispersionFromStats} instead of transposing to gene-major CSC. */
export async function overdispersionScore(
  m: { data: ArrayLike<number>; indptr: ArrayLike<number>; ncells: number },
  opts: { lognorm?: boolean } = {}, M?: any,
): Promise<Float64Array> {
  M ??= await kernels();
  const cmv = await colStats(m, opts, M);
  return overdispersionFromStats(cmv.mean, cmv.var, cmv.nnz, M);
}

/** Per-(group, gene) sufficient stats (group-major K × ngenes): `{sum, sumsq, n_expr}` over `log1p(x)`
 * when `lognorm`, else raw. `codes` maps each cell to a group in [0, K) (-1 = ungrouped). */
export async function groupSufficientStats(
  m: CscMeasure, codes: Int32Array, K: number, opts: { lognorm?: boolean } = {}, M?: any,
): Promise<{ sum: Float64Array; sumsq: Float64Array; n_expr: Float64Array }> {
  M ??= await kernels();
  return M.colSumByGroup(toF64(m.data), toI32(m.indptr), toI32(m.indices), m.ncells, m.ngenes, codes, K, opts.lognorm ?? true);
}

/** Group sizes (cells per code, ignoring -1) — the `nper` that {@link markers} needs. */
export function groupSizes(codes: ArrayLike<number>, K: number): Int32Array {
  const nper = new Int32Array(K);
  for (let i = 0; i < codes.length; i++) { const c = Number(codes[i]); if (c >= 0) nper[c]++; }
  return nper;
}

/** 1-vs-rest marker table (gene-major ngenes × K): `{lfc, padj}` from per-(group,gene) sufficient stats
 * (`sum`, `n_expr` from {@link groupSufficientStats}) and group sizes `nper` ({@link groupSizes}). */
export async function markers(
  sum: ArrayLike<number>, n_expr: ArrayLike<number>, nper: Int32Array,
  K: number, ngenes: number, ncells: number, M?: any,
): Promise<{ lfc: Float64Array; padj: Float64Array }> {
  M ??= await kernels();
  return M.markersOneVsRest(sum, n_expr, nper, K, ngenes, ncells);
}

/** Rank genes distinguishing cell set A from B: per-group mean over `log1p(x)` (when `lognorm`), sorted by
 * |log fold change|, with per-group means. A pure JS reduction (no kernel) — the caller passes the measure
 * whose rows it wants compared (the whole matrix, or a subset it read via `csrRows`) and the already-capped
 * A/B cell-id sets (sampling caps are the caller's policy). */
export function deAvsB(
  m: CscMeasure, A: number[], B: number[], opts: { lognorm?: boolean } = {},
): Array<{ gene: number; meanA: number; meanB: number; lfc: number }> {
  const lognorm = opts.lognorm ?? true;
  const { data, indices, indptr, ncells, ngenes } = m;
  const inA = new Uint8Array(ncells), inB = new Uint8Array(ncells);
  for (const c of A) inA[c] = 1;
  for (const c of B) inB[c] = 1;
  const sumA = new Float64Array(ngenes), sumB = new Float64Array(ngenes);
  for (let j = 0; j < ngenes; j++) {
    for (let k = Number(indptr[j]); k < Number(indptr[j + 1]); k++) {
      const r = Number(indices[k]);
      const v = lognorm ? Math.log1p(Number(data[k])) : Number(data[k]);
      if (inA[r]) sumA[j] += v; else if (inB[r]) sumB[j] += v;
    }
  }
  const out: Array<{ gene: number; meanA: number; meanB: number; lfc: number }> = [];
  for (let j = 0; j < ngenes; j++) {
    const meanA = sumA[j] / A.length, meanB = sumB[j] / B.length;
    out.push({ gene: j, meanA, meanB, lfc: meanA - meanB });
  }
  out.sort((x, y) => Math.abs(y.lfc) - Math.abs(x.lfc));
  return out;
}
