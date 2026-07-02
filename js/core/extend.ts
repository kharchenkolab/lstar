// viewer@0.1 store optimization — the JS twin of Python's `lstar.extend_for_viewer` / `pagoda3.write_viewer`.
// A base L* store (just `counts` + labels + embedding) is fully *viewable*, but DE / variable-genes recompute from
// the whole matrix each session. This computes the precomputed "navigators" a viewer-optimized store carries —
//   • counts_cellmajor (CSR copy)    + counts_cellmajor_order   — substrate for scope compute / per-cell range reads
//   • od_score                       — global per-gene overdispersion residual (log1p)
//   • stats_<g>_{sum,sumsq,nexpr}    — per-group sufficient stats over an induced groups_<g> axis
//   • markers_<g>_{lfc,padj}         — 1-vs-rest marker tables
// using the SAME libstar WASM kernels the live viewer runs (so prepped == live), then `addToStore`s them under the
// `viewer@0.1` profile. Pure Node — no Python / numpy. Kernels: ../wasm/lstar_wasm.cpp.
import { openLstar, type LstarDataset } from "./reader.ts";
import { addToStore, type FieldSpec, type AxisSpec, type LstarWritableStore } from "./writer.ts";
import { selectCountsBasis } from "./basis.ts";
import createLstarKernels from "../dist/lstar_kernels.mjs";

const VIEWER_PROFILE = "viewer@0.1";

export interface ExtendOptions {
  groupings?: string[];   // categorical label fields to build stats/markers for; default = auto-detect
  markers?: boolean;      // also compute 1-vs-rest marker tables (default true)
  counts?: string;        // force the count measure (else auto: a raw-state measure, name "counts" as fallback)
  basis?: string;         // "lognorm" to prep (approximately) from an already log-normalized measure
}

// Per-cell label codes + category names, from EITHER a categorical-encoded field (codes stored) or a utf8 label
// field (strings → derived codes) — both are valid L* label encodings; the viewer treats them the same.
async function labelCodes(ds: LstarDataset, name: string): Promise<{ codes: Int32Array; categories: string[] }> {
  if (ds.field(name)?.encoding === "categorical") { const c = await ds.fieldCategorical(name); return { codes: c.codes, categories: c.categories }; }
  const strings = await ds.fieldStrings(name);
  const idx = new Map<string, number>(), categories: string[] = [], codes = new Int32Array(strings.length);
  for (let i = 0; i < strings.length; i++) { const s = strings[i]; let c = idx.get(s); if (c === undefined) { c = categories.length; categories.push(s); idx.set(s, c); } codes[i] = c; }
  return { codes, categories };
}

// Per-cell label fields (categorical or utf8) with 2..60 distinct values; clustering / cell-type names sort first.
async function detectGroupings(ds: LstarDataset): Promise<string[]> {
  const out: string[] = [];
  for (const name of ds.fieldNames()) {
    const f = ds.field(name);
    if (!f || (f.encoding !== "categorical" && f.encoding !== "utf8") || (f.span?.length ?? 0) !== 1) continue;
    try { const { categories } = await labelCodes(ds, name); if (categories.length >= 2 && categories.length <= 60) out.push(name); } catch { /* not a clean label */ }
  }
  const rank = (n: string) => (/leiden|cluster|cell.?type|louvain/i.test(n) ? 0 : 1);
  return out.sort((a, b) => rank(a) - rank(b));
}

// The primary embedding that keys the within-cluster (Hilbert) locality order: an `embedding`-role field
// over the cell axis with >=2 dims, preferring umap. Mirrors Python/R detect_embedding so the shared core
// reorder gets the same secondary key on every surface. null when no embedding is present.
function detectEmbedding(ds: LstarDataset, cellAxis: string): string | null {
  const cands: string[] = [];
  for (const name of ds.fieldNames()) {
    const f = ds.field(name);
    if (!f || f.role !== "embedding" || !f.span || f.span[0] !== cellAxis) continue;
    const shp = f.shape;
    if (!shp || shp.length < 2 || shp[1] < 2) continue;
    cands.push(name);
  }
  if (!cands.length) return null;
  cands.sort((a, b) => (/umap/i.test(a) ? 0 : 1) - (/umap/i.test(b) ? 0 : 1) || (a < b ? -1 : a > b ? 1 : 0));
  return cands[0];
}

// First 2 embedding dims as a row-major Float64Array (ncells x 2) for viewerCellOrder.
async function readEmbedding2(ds: LstarDataset, name: string, ncells: number): Promise<Float64Array> {
  const { data, shape } = await ds.fieldDense(name);
  const cols = shape[1] ?? 1;
  const out = new Float64Array(ncells * 2);
  for (let i = 0; i < ncells; i++) { out[2 * i] = data[i * cols]; out[2 * i + 1] = data[i * cols + 1]; }
  return out;
}

// Physically reorder CSR rows so physical row p holds cell perm[p]. Deterministic gather -- matches the
// scipy/Matrix row reorders on the other surfaces given the same perm.
function reorderCsrRows(csr: { data: any; indices: any; indptr: any }, perm: Int32Array): { data: any; indices: any; indptr: Int32Array } {
  const n = perm.length, indptr = csr.indptr, indices = csr.indices, data = csr.data;
  const newIndptr = new Int32Array(n + 1);
  for (let p = 0; p < n; p++) { const c = perm[p]; newIndptr[p + 1] = newIndptr[p] + (indptr[c + 1] - indptr[c]); }
  const nnz = newIndptr[n];
  const newData = new (data.constructor as any)(nnz);
  const newIndices = new (indices.constructor as any)(nnz);
  let w = 0;
  for (let p = 0; p < n; p++) { const c = perm[p]; for (let k = indptr[c]; k < indptr[c + 1]; k++) { newData[w] = data[k]; newIndices[w] = indices[k]; w++; } }
  return { data: newData, indices: newIndices, indptr: newIndptr };
}

/** Add the `viewer@0.1` navigator fields to an existing store, in place. `store` must be read+write (e.g. NodeFSStore). */
export async function extendForViewer(store: LstarWritableStore, opts: ExtendOptions = {}): Promise<void> {
  const ds = await openLstar(store as any);
  // Select the count basis by content/state (raw preferred, log1p'd), not the literal name "counts";
  // `counts=`/`basis=` override, and a clear error lists present measures when there's no raw basis.
  const { field: counts, log1p } = selectCountsBasis(ds, { counts: opts.counts, basis: opts.basis });
  const [ncells, ngenes] = ds.field(counts)!.shape as number[];
  const [cellAxis, geneAxis] = ds.field(counts)!.span;
  const basisState = ds.field(counts)!.state ?? "raw";

  const M: any = await createLstarKernels();

  // Normalize counts to CSC -- the layout every kernel expects. A converted store can carry AnnData's
  // native CSR; Python/R normalize via scipy/Matrix, the browser via the shared csrToCsc kernel, so
  // extendForViewer accepts either encoding on every surface (was: JS threw on CSR -- the reported bug).
  const raw = await ds.fieldSparse(counts);
  let cscData = raw.data, cscIndices = raw.indices, cscIndptr = raw.indptr;
  if (raw.fmt === "csr") {
    const c = M.csrToCsc(raw.data, raw.indices, raw.indptr, ncells, ngenes);
    cscData = c.data; cscIndices = c.indices; cscIndptr = c.indptr;
  } else if (raw.fmt !== "csc") {
    throw new Error("extendForViewer: `" + counts + "` must be CSC or CSR (cells x genes), got " + raw.fmt);
  }

  let groupings = (opts.groupings ?? await detectGroupings(ds)).filter((g) => ds.hasField(g));
  if (!groupings.length) throw new Error("extendForViewer: no categorical grouping found (pass {groupings:[...]})");
  const markers = opts.markers ?? true;
  const axes: Record<string, AxisSpec> = {}, fields: Record<string, FieldSpec> = {};
  const prov = { cache: VIEWER_PROFILE };

  // 1) global od_score — per-gene mean/var of log1p over all cells -> pagoda2 overdispersion residual.
  const cmv = M.colMeanVar(cscData, cscIndptr, ncells, 1, log1p);                 // {mean, var, nnz} per gene
  fields["od_score"] = { role: "measure", span: [geneAxis], encoding: "dense", shape: [ngenes], data: M.overdispersion(cmv.mean, cmv.var, cmv.nnz), provenance: { ...prov, method: "viewer.od", basis: log1p ? "log1p" : "lognorm-input" } };

  // 2) per grouping — sufficient stats (group-major K x ngenes) + 1-vs-rest markers (gene-major ngenes x K).
  for (const g of groupings) {
    const { codes, categories } = await labelCodes(ds, g);
    const K = categories.length, gaxis = "groups_" + g;
    axes[gaxis] = { labels: categories, origin: "derived", role: "feature" };
    const s = M.colSumByGroup(cscData, cscIndptr, cscIndices, ncells, ngenes, codes, K, log1p);
    fields["stats_" + g + "_sum"]   = { role: "measure", span: [gaxis, geneAxis], encoding: "dense", shape: [K, ngenes], data: s.sum,    provenance: prov };
    fields["stats_" + g + "_sumsq"] = { role: "measure", span: [gaxis, geneAxis], encoding: "dense", shape: [K, ngenes], data: s.sumsq,  provenance: prov };
    fields["stats_" + g + "_nexpr"] = { role: "measure", span: [gaxis, geneAxis], encoding: "dense", shape: [K, ngenes], data: s.n_expr, provenance: prov };
    if (markers) {
      const nper = new Int32Array(K); for (let i = 0; i < ncells; i++) { const c = codes[i]; if (c >= 0) nper[c]++; }
      const mk = M.markersOneVsRest(s.sum, s.n_expr, nper, K, ngenes, ncells);
      const mp = { ...prov, method: "viewer.markers", test: "1-vs-rest" };
      fields["markers_" + g + "_lfc"]  = { role: "measure", span: [geneAxis, gaxis], encoding: "dense", shape: [ngenes, K], data: mk.lfc,  provenance: mp };
      fields["markers_" + g + "_padj"] = { role: "measure", span: [geneAxis, gaxis], encoding: "dense", shape: [ngenes, K], data: mk.padj, provenance: mp };
    }
  }

  // 3) cell-major (CSR) counts, physically reordered via the SHARED core reorder (cluster code, then
  //    Hilbert of the embedding when present) -- byte-identical to the Python/R surfaces. pos_of[cell] =
  //    physical row; perm is its inverse (physical row -> cell) for the CSR row gather.
  const primaryCodes = (await labelCodes(ds, groupings[0])).codes;
  const embName = detectEmbedding(ds, cellAxis);
  const emb = embName ? await readEmbedding2(ds, embName, ncells) : null;
  const posOf: Int32Array = M.viewerCellOrder(primaryCodes, emb, ncells, 1024);
  const perm = new Int32Array(ncells); for (let i = 0; i < ncells; i++) perm[posOf[i]] = i;
  const csr = reorderCsrRows(M.cscToCsr(cscData, cscIndices, cscIndptr, ncells, ngenes), perm);
  const order = new Float64Array(ncells); for (let i = 0; i < ncells; i++) order[i] = posOf[i];
  fields["counts_cellmajor_order"] = { role: "measure", span: [cellAxis], encoding: "dense", state: "permutation", shape: [ncells], data: order, provenance: { ...prov, method: "viewer.reorder", curve: emb ? "hilbert" : "cluster" } };
  fields["counts_cellmajor"] = { role: "measure", span: [cellAxis, geneAxis], encoding: "csr", state: basisState, shape: [ncells, ngenes], data: csr.data, indices: csr.indices, indptr: csr.indptr, provenance: prov };

  await addToStore(store, { axes, fields, profiles: [VIEWER_PROFILE] });
}
