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
import createLstarKernels from "../dist/lstar_kernels.mjs";

const VIEWER_PROFILE = "viewer@0.1";

export interface ExtendOptions {
  groupings?: string[];   // categorical label fields to build stats/markers for; default = auto-detect
  markers?: boolean;      // also compute 1-vs-rest marker tables (default true)
  counts?: string;        // raw-counts measure to summarize (default "counts"; e.g. "X" for an AnnData .X)
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

/** Add the `viewer@0.1` navigator fields to an existing store, in place. `store` must be read+write (e.g. NodeFSStore). */
export async function extendForViewer(store: LstarWritableStore, opts: ExtendOptions = {}): Promise<void> {
  const ds = await openLstar(store as any);
  const counts = opts.counts ?? "counts";
  if (!ds.hasField(counts)) throw new Error("extendForViewer: store has no `" + counts + "` measure");
  const sp = await ds.fieldSparse(counts);
  if (sp.fmt !== "csc") throw new Error("extendForViewer: `" + counts + "` must be CSC (cells x genes), got " + sp.fmt);
  const [ncells, ngenes] = sp.shape;
  const [cellAxis, geneAxis] = ds.field(counts)!.span;

  let groupings = (opts.groupings ?? await detectGroupings(ds)).filter((g) => ds.hasField(g));
  if (!groupings.length) throw new Error("extendForViewer: no categorical grouping found (pass {groupings:[...]})");
  const markers = opts.markers ?? true;

  const M: any = await createLstarKernels();
  const axes: Record<string, AxisSpec> = {}, fields: Record<string, FieldSpec> = {};
  const prov = { cache: VIEWER_PROFILE };

  // 1) global od_score — per-gene mean/var of log1p over all cells -> pagoda2 overdispersion residual.
  const cmv = M.colMeanVar(sp.data, sp.indptr, ncells, 1, true);                 // {mean, var, nnz} per gene
  fields["od_score"] = { role: "measure", span: [geneAxis], encoding: "dense", shape: [ngenes], data: M.overdispersion(cmv.mean, cmv.var, cmv.nnz), provenance: { ...prov, method: "viewer.od", basis: "log1p" } };

  // 2) per grouping — sufficient stats (group-major K x ngenes) + 1-vs-rest markers (gene-major ngenes x K).
  for (const g of groupings) {
    const { codes, categories } = await labelCodes(ds, g);
    const K = categories.length, gaxis = "groups_" + g;
    axes[gaxis] = { labels: categories, origin: "derived", role: "feature" };
    const s = M.colSumByGroup(sp.data, sp.indptr, sp.indices, ncells, ngenes, codes, K, true);
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

  // 3) cell-major (CSR) counts copy + order. v1 keeps cells in cell order with an identity order so the store reads
  //    as fully optimized; the hybrid locality reorder (Hilbert over the embedding — a remote-host read-coalescing
  //    win) is a follow-up.
  const csr = M.cscToCsr(sp.data, sp.indices, sp.indptr, ncells, ngenes);
  const order = new Float64Array(ncells); for (let i = 0; i < ncells; i++) order[i] = i;
  fields["counts_cellmajor_order"] = { role: "measure", span: [cellAxis], encoding: "dense", state: "permutation", shape: [ncells], data: order, provenance: { ...prov, method: "viewer.reorder", curve: "identity" } };
  fields["counts_cellmajor"] = { role: "measure", span: [cellAxis, geneAxis], encoding: "csr", state: "raw", shape: [ncells, ngenes], data: csr.data, indices: csr.indices, indptr: csr.indptr, provenance: prov };

  await addToStore(store, { axes, fields, profiles: [VIEWER_PROFILE] });
}
