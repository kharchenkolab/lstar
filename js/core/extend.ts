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
import { MIN_GROUPS, MAX_GROUPS, groupingRank, embeddingRank, HILBERT_GRID } from "./policy.ts";
import createLstarKernels from "../dist/lstar_kernels.mjs";

const VIEWER_PROFILE = "viewer@0.1";

export interface ExtendOptions {
  groupings?: string[];   // categorical label fields to build stats/markers for; default = auto-detect
  primary?: string;       // the grouping the VIEWER OPENS ON: hoisted to the front (the counts_cellmajor reorder
                          // key + summarized first), and COMPOSES with auto-detect (the rest are still prepped) —
                          // which `groupings` alone can't express. Default: the first detected grouping is primary.
  markers?: boolean;      // also compute 1-vs-rest marker tables (default true)
  counts?: string;        // force the count measure (else auto: a raw-state measure, name "counts" as fallback)
  basis?: string;         // "lognorm" to prep (approximately) from an already log-normalized measure
  order?: string;         // "hybrid" (default) locality reorder + _order; "none" keeps rows in cell order
}

// Per-cell label codes + category names, from EITHER a categorical-encoded field (codes stored) or a utf8 label
// field (strings → derived codes) — both are valid L* label encodings; the viewer treats them the same.
async function labelCodes(ds: LstarDataset, name: string): Promise<{ codes: Int32Array; categories: string[] }> {
  // Decode to per-cell strings from EITHER a categorical (stored codes+categories) or a utf8 label field.
  let strings: string[];
  if (ds.field(name)?.encoding === "categorical") {
    const c = await ds.fieldCategorical(name);
    strings = Array.from(c.codes, (k) => (k >= 0 ? c.categories[k] : ""));
  } else {
    strings = await ds.fieldStrings(name);
  }
  // Category order MUST be sorted-unique to match Python (np.unique) and R (sort(unique)) so the induced
  // groups_<g> axis + the group codes align field-for-field across surfaces. Was: first-seen (utf8) /
  // stored (categorical) order -- a silent divergence that permuted stats/markers rows and the reorder key.
  const categories = Array.from(new Set(strings)).sort();
  const idx = new Map(categories.map((s, i) => [s, i] as [string, number]));
  const codes = new Int32Array(strings.length);
  for (let i = 0; i < strings.length; i++) codes[i] = idx.get(strings[i])!;
  return { codes, categories };
}

// Per-cell label fields (categorical or utf8) with 2..60 distinct values; clustering / cell-type names sort first.
async function detectGroupings(ds: LstarDataset, cellAxis: string): Promise<string[]> {
  const out: string[] = [];
  for (const name of ds.fieldNames()) {
    const f = ds.field(name);
    // Must be a label over the CELL axis: `span[0] === cellAxis`, not merely 1-D. Without this, a 1-D
    // label over the GENE axis (e.g. a `highly_variable` gene flag with 2..60 levels) was picked as a
    // grouping in JS but rejected by Py/R -> ngenes codes fed to a cells kernel -> out-of-bounds/garbage.
    // skip Seurat's active-idents mirror (subtype "active_ident" — a UI-state copy of the current identity,
    // usually == a clustering already present); it's not a separate grouping (matches Python/R).
    if (!f || f.subtype === "active_ident" || (f.encoding !== "categorical" && f.encoding !== "utf8") || (f.span?.length ?? 0) !== 1 || f.span![0] !== cellAxis) continue;
    try { const { categories } = await labelCodes(ds, name); if (categories.length >= MIN_GROUPS && categories.length <= MAX_GROUPS) out.push(name); } catch { /* not a clean label */ }
  }
  // preferred names (clustering / cell-type) first by list position, then alphabetical -- identical to
  // Python _detect_groupings / R .detect_groupings (single source: policy.ts <-> viewer_policy.json).
  return out.sort((a, b) => groupingRank(a) - groupingRank(b) || (a < b ? -1 : a > b ? 1 : 0));
}

// The primary embedding that keys the within-cluster (Hilbert) locality order: an `embedding`-role field
// over the cell axis with >=2 dims, preferring umap. Mirrors Python/R detect_embedding so the shared core
// reorder gets the same secondary key on every surface. null when no embedding is present.
function detectEmbedding(ds: LstarDataset, cellAxis: string): string | null {
  const cands: string[] = [];
  for (const name of ds.fieldNames()) {
    const f = ds.field(name);
    // Dense fields don't carry `shape` in the field manifest (the zarr array holds it), so key off the
    // span: an embedding over the cell axis whose 2nd axis has >=2 dims. (Was: `f.shape[1]` -> undefined
    // for a dense embedding -> every embedding skipped -> cluster-only order, diverging from Python/R.)
    if (!f || f.role !== "embedding" || !f.span || f.span[0] !== cellAxis || f.span.length < 2) continue;
    if (ds.axisLength(f.span[1]) < 2) continue;
    cands.push(name);
  }
  if (!cands.length) return null;
  cands.sort((a, b) => embeddingRank(a) - embeddingRank(b) || (a < b ? -1 : a > b ? 1 : 0));
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

/** Dense (flat, C-order, cells x genes) -> CSC arrays, dropping zeros — matches scipy's `csc_matrix(dense)`
 * so the viewer kernels (which expect CSC and treat absent entries as zero) yield identical stats for a
 * dense primary measure as for a native-sparse one. Column-major: one CSC column per gene. */
function denseToCsc(dense: ArrayLike<number>, ncells: number, ngenes: number): { data: Float64Array; indices: Int32Array; indptr: Int32Array } {
  const indptr = new Int32Array(ngenes + 1);
  for (let g = 0; g < ngenes; g++) {
    let c = 0;
    for (let i = 0; i < ncells; i++) if (dense[i * ngenes + g] !== 0) c++;
    indptr[g + 1] = indptr[g] + c;
  }
  const nnz = indptr[ngenes];
  const data = new Float64Array(nnz), indices = new Int32Array(nnz);
  let w = 0;
  for (let g = 0; g < ngenes; g++)
    for (let i = 0; i < ncells; i++) { const v = dense[i * ngenes + g]; if (v !== 0) { data[w] = v; indices[w] = i; w++; } }
  return { data, indices, indptr };
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
  // Read the counts basis as CSC regardless of its on-disk encoding — dense (`/values`), CSR, or CSC —
  // mirroring Python (`X.tocsc() if issparse else csc_matrix(X)`) and R (C++ core densify->CSC). A DENSE
  // primary measure (an SCE logcounts assay, or a scaled/dense AnnData X) lives at `/values`, not the sparse
  // `/data`; reading it via fieldSparse threw NotFoundError and silently skipped viewer-opt (the reported bug).
  let cscData: any, cscIndices: any, cscIndptr: any;
  if (ds.field(counts)!.encoding === "dense") {
    const dv = await ds.fieldDense(counts);                          // flat C-order (cells x genes)
    const csc = denseToCsc(dv.data, ncells, ngenes);                 // column-major, zeros dropped (== scipy csc)
    cscData = csc.data; cscIndices = csc.indices; cscIndptr = csc.indptr;
  } else {
    const raw = await ds.fieldSparse(counts);
    cscData = raw.data; cscIndices = raw.indices; cscIndptr = raw.indptr;
    if (raw.fmt === "csr") {
      const c = M.csrToCsc(raw.data, raw.indices, raw.indptr, ncells, ngenes);
      cscData = c.data; cscIndices = c.indices; cscIndptr = c.indptr;
    } else if (raw.fmt !== "csc") {
      throw new Error("extendForViewer: `" + counts + "` must be dense, CSC, or CSR (cells x genes), got " + raw.fmt);
    }
  }

  let groupings = (opts.groupings ?? await detectGroupings(ds, cellAxis)).filter((g) => ds.hasField(g));
  // Hoist the viewer's primary grouping to the front (guaranteed present): it keys the counts_cellmajor
  // reorder + is summarized first. Composes with auto-detect — the other groupings are still prepped.
  if (opts.primary != null) {
    if (!ds.hasField(opts.primary)) throw new Error("extendForViewer: primary `" + opts.primary + "` is not a field");
    // must be a 1-D grouping over the CELL axis (else a cryptic reorder crash); span==[cellAxis] is the check
    // identical across Py/R/JS (their detection predicates differ, but this structural one does not).
    const psp = ds.fields.get(opts.primary)?.span;
    if (!psp || psp.length !== 1 || psp[0] !== cellAxis)
      throw new Error("extendForViewer: primary `" + opts.primary + "` must be a grouping over the cell axis `" + cellAxis + "` (a 1-D label)");
    groupings = [opts.primary, ...groupings.filter((g) => g !== opts.primary)];
  }
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

  // 3) cell-major (CSR) counts. order="hybrid" (default) physically reorders rows via the SHARED core
  //    reorder (cluster code, then Hilbert of the embedding) + writes counts_cellmajor_order -- byte-
  //    identical to Python/R. order="none" keeps rows in cell order and omits _order (parity with Python).
  let csr: { data: any; indices: any; indptr: any };
  if ((opts.order ?? "hybrid") !== "none") {
    const primaryCodes = (await labelCodes(ds, groupings[0])).codes;
    const embName = detectEmbedding(ds, cellAxis);
    const emb = embName ? await readEmbedding2(ds, embName, ncells) : null;
    const posOf: Int32Array = M.viewerCellOrder(primaryCodes, emb, ncells, HILBERT_GRID);
    const perm = new Int32Array(ncells); for (let i = 0; i < ncells; i++) perm[posOf[i]] = i;   // physical row -> cell
    csr = reorderCsrRows(M.cscToCsr(cscData, cscIndices, cscIndptr, ncells, ngenes), perm);
    const order = new Float64Array(ncells); for (let i = 0; i < ncells; i++) order[i] = posOf[i];
    fields["counts_cellmajor_order"] = { role: "measure", span: [cellAxis], encoding: "dense", state: "permutation", shape: [ncells], data: order, provenance: { ...prov, method: "viewer.reorder", curve: emb ? "hilbert" : "cluster", group: groupings[0] } };
  } else {
    csr = M.cscToCsr(cscData, cscIndices, cscIndptr, ncells, ngenes);   // no reorder, no _order field
  }
  fields["counts_cellmajor"] = { role: "measure", span: [cellAxis, geneAxis], encoding: "csr", state: basisState, shape: [ncells, ngenes], data: csr.data, indices: csr.indices, indptr: csr.indptr, provenance: prov };

  await addToStore(store, { axes, fields, profiles: [VIEWER_PROFILE] });
}
