// Phase C — a framework-agnostic viewer query API over an L* store.
//
// These are the operations the next-gen pagoda2 viewer needs (app_prop2.md): navigate an embedding,
// color by gene or metadata, rank genes, and cross-filter selections. Numeric reductions go to the
// WASM kernels (the same libstar core as R/Python); I/O comes from the reader. This is the generic
// L*-store implementation — a future "viewer profile" (precomputed summaries / cell-major DE panel)
// will make groupStats/DE instant; here they read the measure, which is fine at moderate scale.
import createLstarKernels from "../dist/lstar_kernels.mjs";

import type { LstarDataset } from "./reader.ts";

export interface ColStats { mean: Float64Array; var: Float64Array; nnz: Int32Array; }
export type Metadata =
  | { kind: "categorical"; codes: Int32Array; categories: string[] }
  | { kind: "numeric"; values: Float32Array };

function toF64(a: any): Float64Array { return a instanceof Float64Array ? a : Float64Array.from(a); }
function toI32(a: any): Int32Array { return a instanceof Int32Array ? a : Int32Array.from(a); }
function toF32(a: any): Float32Array { return a instanceof Float32Array ? a : Float32Array.from(a); }

export class LstarView {
  ds: LstarDataset;
  private kernels: any = null;

  constructor(ds: LstarDataset) { this.ds = ds; }
  private async M(): Promise<any> {
    if (!this.kernels) this.kernels = await createLstarKernels();
    return this.kernels;
  }

  /** Embedding coordinates as a flat Float32Array (n_cells × d, C-order) for a point layer. */
  async embedding(name = "umap"): Promise<{ data: Float32Array; n: number; dim: number }> {
    const { data, shape } = await this.ds.fieldDense(name);
    return { data: toF32(data), n: shape[0], dim: shape[1] ?? 1 };
  }

  /** Metadata for coloring/filtering: a label field → codes+categories, a numeric field → values. */
  async metadata(name: string): Promise<Metadata> {
    const meta = this.ds.fields.get(name);
    if (!meta) throw new Error("no field " + name);
    if (meta.encoding === "categorical") {        // codes + categories are stored directly
      const { codes, categories } = await this.ds.fieldCategorical(name);
      return { kind: "categorical", codes, categories };
    }
    if (meta.encoding === "utf8") {
      const strings = await this.ds.fieldStrings(name);
      const categories: string[] = [];
      const idx = new Map<string, number>();
      const codes = new Int32Array(strings.length);
      for (let i = 0; i < strings.length; i++) {
        let c = idx.get(strings[i]);
        if (c === undefined) { c = categories.length; categories.push(strings[i]); idx.set(strings[i], c); }
        codes[i] = c;
      }
      return { kind: "categorical", codes, categories };
    }
    const { data } = await this.ds.fieldDense(name);
    return { kind: "numeric", values: toF32(data) };
  }

  /**
   * Per-cell expression of one gene (0 for non-expressing), normalized on the fly — the hot path:
   * fetch only that gene's CSC column, transform its nonzeros. Returns the scalar field; map it to
   * colors with `scalarToRGBA`.
   */
  async geneExpression(gene: string, opts: { lognorm?: boolean; field?: string } = {}):
      Promise<{ values: Float32Array; max: number; col: number }> {
    const lognorm = opts.lognorm ?? true;
    const field = opts.field ?? "counts";
    const genes = await this.ds.axisLabels("genes");
    const col = genes.indexOf(gene);
    if (col < 0) throw new Error("no gene " + gene);
    const ncells = this.ds.axes.get("cells")!.length;
    const { rows, vals } = await this.ds.cscColumn(field, col);
    const out = new Float32Array(ncells);
    let max = 0;
    for (let k = 0; k < rows.length; k++) {
      const v = lognorm ? Math.log1p(vals[k]) : vals[k];
      out[rows[k]] = v;
      if (v > max) max = v;
    }
    return { values: out, max, col };
  }

  /** Per-gene zero-aware mean/variance (HVG ranking) over the whole measure, via the WASM kernel. */
  async colStats(opts: { lognorm?: boolean; field?: string } = {}): Promise<ColStats> {
    const lognorm = opts.lognorm ?? true;
    const sp = await this.ds.fieldSparse(opts.field ?? "counts");
    if (sp.fmt !== "csc") throw new Error("colStats needs a csc measure, got " + sp.fmt);
    const M = await this.M();
    return M.colMeanVar(toF64(sp.data), toI32(sp.indptr), sp.shape[0], 1, lognorm);
  }

  /**
   * Rank genes distinguishing cell set A from B. Generic-store implementation: a single pass over the
   * measure accumulating per-group sums (subsample each group to bound cost). Returns genes sorted by
   * |log fold change| with per-group means. (A viewer profile's cell-major DE panel makes this a few
   * row reads — a later phase.)
   */
  async subsampleDE(cellsA: number[], cellsB: number[],
                    opts: { lognorm?: boolean; field?: string; maxPerGroup?: number } = {}):
      Promise<Array<{ gene: number; meanA: number; meanB: number; lfc: number }>> {
    const lognorm = opts.lognorm ?? true;
    const cap = opts.maxPerGroup ?? Infinity;
    const A = cap === Infinity ? cellsA : cellsA.slice(0, cap);
    const B = cap === Infinity ? cellsB : cellsB.slice(0, cap);
    const sp = await this.ds.fieldSparse(opts.field ?? "counts");
    const [nrows, ncols] = sp.shape;
    const inA = new Uint8Array(nrows), inB = new Uint8Array(nrows);
    for (const c of A) inA[c] = 1;
    for (const c of B) inB[c] = 1;
    const { data, indices, indptr } = sp;
    const sumA = new Float64Array(ncols), sumB = new Float64Array(ncols);
    for (let j = 0; j < ncols; j++) {
      for (let k = Number(indptr[j]); k < Number(indptr[j + 1]); k++) {
        const r = indices[k];
        const v = lognorm ? Math.log1p(data[k]) : data[k];
        if (inA[r]) sumA[j] += v; else if (inB[r]) sumB[j] += v;
      }
    }
    const out = [];
    for (let j = 0; j < ncols; j++) {
      const meanA = sumA[j] / A.length, meanB = sumB[j] / B.length;
      out.push({ gene: j, meanA, meanB, lfc: meanA - meanB });
    }
    out.sort((x, y) => Math.abs(y.lfc) - Math.abs(x.lfc));
    return out;
  }
}

/** A typed-array bitmap crossfilter over cells: facets combine by AND; the mask drives selection. */
export class Crossfilter {
  n: number;
  mask: Uint8Array;
  constructor(n: number) { this.n = n; this.mask = new Uint8Array(n).fill(1); }
  reset(): this { this.mask.fill(1); return this; }
  categorical(codes: Int32Array, keep: number[]): this {
    const set = new Set(keep);
    for (let i = 0; i < this.n; i++) if (!set.has(codes[i])) this.mask[i] = 0;
    return this;
  }
  range(values: ArrayLike<number>, lo: number, hi: number): this {
    for (let i = 0; i < this.n; i++) if (values[i] < lo || values[i] > hi) this.mask[i] = 0;
    return this;
  }
  selected(): number[] {
    const out: number[] = [];
    for (let i = 0; i < this.n; i++) if (this.mask[i]) out.push(i);
    return out;
  }
  count(): number { let c = 0; for (let i = 0; i < this.n; i++) c += this.mask[i]; return c; }
}

/** Map a scalar field to RGBA bytes with a simple perceptual ramp (viewer color attribute). */
export function scalarToRGBA(values: ArrayLike<number>, max: number): Uint8Array {
  const n = values.length;
  const rgba = new Uint8Array(n * 4);
  const m = max > 0 ? max : 1;
  for (let i = 0; i < n; i++) {
    const t = Math.max(0, Math.min(1, values[i] / m));
    rgba[i * 4] = Math.round(255 * t);                 // R rises
    rgba[i * 4 + 1] = Math.round(64 + 120 * t);        // G
    rgba[i * 4 + 2] = Math.round(255 * (1 - t));       // B falls
    rgba[i * 4 + 3] = 255;
  }
  return rgba;
}
