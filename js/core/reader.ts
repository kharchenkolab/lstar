// @lstar/core — a lazy reader for L* Zarr stores, over zarrita.js.
//
// Opens a store (HTTP/local/zip), reads the consolidated L* metadata, and exposes axes/fields with
// values fetched only when asked — including a single CSC gene-column (the viewer's hot path). The
// heavy numeric work belongs in the WASM kernels (see view.ts); this module is I/O + assembly.
import * as zarr from "zarrita";

/** Minimal store contract (zarrita-compatible): fetch one object by key, or undefined if absent.
 * `getRange` is optional: a store that can serve a byte sub-range of an object (HTTP `Range`, a file
 * `pread`) enables the reader's sub-chunk fast path (one gene/cell = a few KB instead of the whole
 * uncompressed chunk). Stores without it still work — the reader falls back to whole-chunk reads. */
export interface LstarStore {
  get(key: string): Promise<Uint8Array | undefined>;
  getRange?(key: string, start: number, end: number): Promise<Uint8Array | undefined>;
}

export interface AxisMeta {
  name: string;
  origin: string;
  role?: string;            // observation | feature | coordinate | factor | ...
  induced_by?: string;      // the field this axis was induced from (a `factor` axis's categorical label)
  length: number;
}

export interface FieldMeta {
  name: string;
  role?: string;
  span: string[];
  encoding: string; // dense | csc | csr | coo | utf8
  state?: string;
  subtype?: string;
  shape?: number[];
  nullable?: boolean;   // carries a `mask` array (1 == missing): nullable Int/boolean/string
}

const TD = new TextDecoder();
const TE = new TextEncoder();

// Keys whose bytes are Zarr metadata (served from consolidated metadata when available). Everything
// else (chunk keys ".../0", label/value byte arrays) is data and goes to the underlying store.
const META_RE = /(?:\.zgroup|\.zattrs|\.zarray|zarr\.json)$/;

/**
 * Wraps a store so every metadata read is served from a parsed consolidated `.zmetadata` map instead
 * of the network: a present key returns its bytes, an ABSENT metadata key returns undefined *without*
 * a request (this is what suppresses zarrita's per-node v3 `zarr.json` probes — the difference between
 * one open request and ~80). Data chunk reads (and byte-range reads) pass straight through.
 */
class ConsolidatedStore {
  inner: LstarStore;
  meta: Record<string, unknown>;
  constructor(inner: LstarStore, meta: Record<string, unknown>) { this.inner = inner; this.meta = meta; }
  async get(key: string, _opts?: unknown): Promise<Uint8Array | undefined> {
    const norm = key[0] === "/" ? key.slice(1) : key;
    if (META_RE.test(norm)) {
      const v = this.meta[norm];
      return v === undefined ? undefined : TE.encode(JSON.stringify(v));
    }
    return this.inner.get(key);
  }
  async getRange(key: string, start: number, end: number): Promise<Uint8Array | undefined> {
    return this.inner.getRange?.(key, start, end);
  }
}

// Zarr v2 little-endian dtype -> [typed-array constructor, itemsize]. The writer only emits these LE
// dtypes; the byte-range fast path decodes raw bytes through this table (we assume LE — true on every
// platform Node/browsers run, and the writer always tags "<").
const DTYPE: Record<string, [any, number]> = {
  "<f8": [Float64Array, 8], "<f4": [Float32Array, 4],
  "<i4": [Int32Array, 4], "<i8": [BigInt64Array, 8], "|i1": [Int8Array, 1],
  "<i2": [Int16Array, 2], "|u1": [Uint8Array, 1], "<u4": [Uint32Array, 4], "<u8": [BigUint64Array, 8],
};

// Copy raw little-endian bytes into a freshly-allocated typed array (the copy guarantees alignment —
// a Range body may start at an arbitrary offset — and length is always an itemsize multiple).
function decodeTyped(bytes: Uint8Array, ctor: any, isize: number): any {
  const out = new ctor(Math.floor(bytes.byteLength / isize));
  new Uint8Array(out.buffer).set(bytes.subarray(0, out.byteLength));
  return out;
}

function decodeStrings(bytes: Uint8Array, offsets: ArrayLike<number | bigint>): string[] {
  const n = offsets.length - 1;
  const out: string[] = new Array(n);
  for (let i = 0; i < n; i++) {
    out[i] = TD.decode(bytes.subarray(Number(offsets[i]), Number(offsets[i + 1])));
  }
  return out;
}

export class LstarDataset {
  rawStore: any;                          // the underlying store (data chunks, byte-range reads)
  store: any;                             // == rawStore, or a ConsolidatedStore wrapping it
  root: any;
  meta: Record<string, any> = {};         // parsed `.zmetadata` map (empty if the store has none)
  kind = "sample";
  specVersion = "0.1";
  profiles: string[] = [];
  dropped: string[] = [];
  auxNames: string[] = [];   // namespaces of the lossless passthrough subtree (uns/@misc)
  axes = new Map<string, AxisMeta>();
  fields = new Map<string, FieldMeta>();

  constructor(store: any) {
    this.rawStore = store;
    this.store = store;
    this.root = zarr.root(store);
  }

  private _open(p: string) {
    return zarr.open(this.root.resolve(p), { kind: "array" });
  }
  private async _get(p: string, sel?: any) {
    return zarr.get(await this._open(p), sel);
  }

  async init(): Promise<this> {
    // Consolidated open: one `.zmetadata` read up front serves every group/array metadata object, so
    // opening the manifest + axes + fields is a single request instead of ~80 (and gives cscColumn the
    // per-array dtype/chunks/compressor it needs for the byte-range fast path). Absent or malformed ->
    // per-object reads (still correct), e.g. an older store or one extended without refreshing it.
    const zmeta = await this.rawStore.get(".zmetadata");
    if (zmeta) {
      try {
        const parsed = JSON.parse(TD.decode(zmeta));
        this.meta = parsed?.metadata ?? {};
        if (Object.keys(this.meta).length) {
          this.store = new ConsolidatedStore(this.rawStore, this.meta);
          this.root = zarr.root(this.store);
        }
      } catch { /* malformed consolidated metadata -> fall back to per-object reads */ }
    }
    const grp = await zarr.open(this.root, { kind: "group" });
    const m = (grp.attrs as any).lstar;
    if (!m) throw new Error("not an L* store (no 'lstar' root attribute)");
    this.kind = m.kind ?? "sample";
    this.specVersion = m.spec_version ?? "0.1";
    this.profiles = m.profiles ?? [];
    this.dropped = m.dropped ?? [];
    this.auxNames = m.passthrough ?? [];
    for (const name of m.axes as string[]) {
      const ax = await zarr.open(this.root.resolve("axes/" + name), { kind: "group" });
      const offs = await this._open("axes/" + name + "/labels_offsets");
      const lm = (ax.attrs as any).lstar ?? {};
      this.axes.set(name, { name, origin: lm.origin ?? "observed", role: lm.role,
                            induced_by: lm.induced_by ?? undefined,
                            length: (offs.shape[0] as number) - 1 });
    }
    for (const name of m.fields as string[]) {
      const f = await zarr.open(this.root.resolve("fields/" + name), { kind: "group" });
      const lm = (f.attrs as any).lstar ?? {};
      this.fields.set(name, { name, role: lm.role, span: lm.span ?? [], encoding: lm.encoding,
                              state: lm.state ?? undefined, subtype: lm.subtype ?? undefined,
                              shape: lm.shape, ordered: lm.ordered,
                              nullable: lm.nullable ?? undefined } as any);
    }
    return this;
  }

  axisNames(): string[] { return [...this.axes.keys()]; }
  fieldNames(): string[] { return [...this.fields.keys()]; }
  hasField(name: string): boolean { return this.fields.has(name); }
  field(name: string): FieldMeta | undefined { return this.fields.get(name); }
  axisLength(name: string): number { return this.axes.get(name)?.length ?? 0; }

  async axisLabels(name: string): Promise<string[]> {
    const bytes = (await this._get("axes/" + name + "/labels")).data as Uint8Array;
    const offs = (await this._get("axes/" + name + "/labels_offsets")).data as any;
    return decodeStrings(bytes, offs);
  }

  /** A dense field's values (flat, C-order) + shape. */
  async fieldDense(name: string): Promise<{ data: any; shape: number[] }> {
    const c = await this._get("fields/" + name + "/values");
    return { data: c.data, shape: c.shape as number[] };
  }

  /** A utf8 (label) field's string values. */
  async fieldStrings(name: string): Promise<string[]> {
    const bytes = (await this._get("fields/" + name + "/values")).data as Uint8Array;
    const offs = (await this._get("fields/" + name + "/values_offsets")).data as any;
    return decodeStrings(bytes, offs);
  }

  /**
   * A lossless-passthrough namespace (e.g. "anndata.uns"), reconstructed into a plain JS object: the
   * stored JSON `tree` with its array-leaf references resolved (a numpy structured array becomes a
   * `{field: array}` object). Read-only — for inspection / promoting a recognized structure to a field.
   */
  async aux(ns: string): Promise<any> {
    const g = await zarr.open(this.root.resolve("passthrough/" + ns), { kind: "group" });
    const lm = (g.attrs as any).lstar ?? {};
    const tree = typeof lm.tree === "string" ? JSON.parse(lm.tree) : lm.tree;
    const byId: Record<string, any> = {};
    for (const a of (lm.arrays ?? []) as Array<{ id: string; kind: string }>) {
      if (a.kind === "utf8") {
        const bytes = (await this._get("passthrough/" + ns + "/" + a.id)).data as Uint8Array;
        const offs = (await this._get("passthrough/" + ns + "/" + a.id + "_offsets")).data as any;
        byId[a.id] = decodeStrings(bytes, offs);
      } else {
        byId[a.id] = (await this._get("passthrough/" + ns + "/" + a.id)).data;
      }
    }
    const build = (node: any): any => {
      if (node === null || typeof node !== "object") return node;
      if ("$obj" in node) { const o: any = {}; for (const k in node.$obj) o[k] = build(node.$obj[k]); return o; }
      if ("$list" in node) return (node.$list as any[]).map(build);
      if ("$array" in node) return byId[node.$array];
      if ("$strings" in node) return byId[node.$strings];
      if ("$record" in node) { const o: any = {}; for (const k in node.$record.fields) o[k] = build(node.$record.fields[k]); return o; }
      if ("$bytes" in node) return node.$bytes;
      if ("$dropped" in node) return null;
      return node;
    };
    return build(tree);
  }

  /** A nullable field's validity mask (`1 == missing`), or null when the field carries no nulls. */
  async fieldMask(name: string): Promise<Uint8Array | null> {
    if (!(this.fields.get(name) as any)?.nullable) return null;
    return Uint8Array.from((await this._get("fields/" + name + "/mask")).data as any);
  }

  /** A categorical (factor) field: integer codes (-1 = missing) + the category labels + ordered. */
  async fieldCategorical(name: string): Promise<{ codes: Int32Array; categories: string[]; ordered: boolean }> {
    const codes = Int32Array.from((await this._get("fields/" + name + "/codes")).data as any);
    const bytes = (await this._get("fields/" + name + "/categories")).data as Uint8Array;
    const offs = (await this._get("fields/" + name + "/categories_offsets")).data as any;
    return { codes, categories: decodeStrings(bytes, offs), ordered: !!(this.fields.get(name) as any)?.ordered };
  }

  /** A sparse (csc/csr) field's raw arrays + shape. Reads the whole field. */
  async fieldSparse(name: string): Promise<{ data: any; indices: any; indptr: any; shape: number[]; fmt: string }> {
    const meta = this.fields.get(name)!;
    const data = (await this._get("fields/" + name + "/data")).data;
    const indices = (await this._get("fields/" + name + "/indices")).data;
    const indptr = (await this._get("fields/" + name + "/indptr")).data;
    return { data, indices, indptr, shape: meta.shape as number[], fmt: meta.encoding };
  }

  /**
   * Read elements [lo, hi) of a 1-D array as a typed array. Fast path: when the store supports byte
   * ranges AND the array is a single uncompressed, unfiltered chunk (known from consolidated metadata),
   * issue ONE `getRange` over the exact bytes [lo·itemsize, hi·itemsize) of chunk "0" — a gene/cell is
   * a few KB instead of the whole chunk. Otherwise fall back to a zarrita slice (which decompresses /
   * stitches chunks / works without consolidated metadata). Equivalent results either way.
   */
  private async _rangeOrSlice(arrPath: string, lo: number, hi: number): Promise<any> {
    const za: any = this.meta[arrPath + "/.zarray"];
    const dt = za ? DTYPE[za.dtype] : undefined;
    const singleChunk = za && Array.isArray(za.chunks) && za.chunks.length === 1 &&
      za.chunks[0] >= (za.shape?.[0] ?? 0);
    if (typeof this.store.getRange === "function" && dt && singleChunk &&
        za.compressor == null && (za.filters == null || (Array.isArray(za.filters) && za.filters.length === 0))) {
      const [ctor, isize] = dt;
      const bytes = await this.store.getRange(arrPath + "/0", lo * isize, hi * isize);
      if (bytes) return decodeTyped(bytes, ctor, isize);
    }
    return (await this._get(arrPath, [zarr.slice(lo, hi)])).data;
  }

  /**
   * One CSC column (e.g. a gene's expression across cells) — the viewer's hot path. Returns the
   * nonzero `rows` (cell indices) and `vals`, reading only indptr[col..col+1] and the corresponding
   * indices/data ranges (a single byte-range each on a consolidated, uncompressed store).
   */
  async cscColumn(name: string, col: number): Promise<{ rows: any; vals: any }> {
    const base = "fields/" + name;
    const ip = await this._rangeOrSlice(base + "/indptr", col, col + 2);
    const a = Number(ip[0]), b = Number(ip[1]);
    if (b <= a) return { rows: new Int32Array(0), vals: new Float64Array(0) };
    const rows = await this._rangeOrSlice(base + "/indices", a, b);
    const vals = await this._rangeOrSlice(base + "/data", a, b);
    return { rows, vals };
  }

  /**
   * One CSR row (e.g. a cell's expression across genes), symmetric to {@link cscColumn}. Returns the
   * nonzero `cols` (feature indices) and `vals` for the row, via a single byte-range each on the fast
   * path. (The cell-major counts copy stores cells as rows — this is its per-cell accessor.)
   */
  async csrRow(name: string, row: number): Promise<{ cols: any; vals: any }> {
    const base = "fields/" + name;
    const ip = await this._rangeOrSlice(base + "/indptr", row, row + 2);
    const a = Number(ip[0]), b = Number(ip[1]);
    if (b <= a) return { cols: new Int32Array(0), vals: new Float64Array(0) };
    const cols = await this._rangeOrSlice(base + "/indices", a, b);
    const vals = await this._rangeOrSlice(base + "/data", a, b);
    return { cols, vals };
  }

  // Cell -> physical-row permutation for a row-reordered cell-major field (a locality order), or null when the field is
  // stored in canonical cell order. Convention: a sibling field `<name>_order` holds each cell's physical row index; its
  // presence flags `<name>` as reordered. Loaded once per field. (A later format rev could carry this in the field's own
  // metadata rather than by name.)
  private _rowMapCache = new Map<string, Int32Array | null>();
  private async _rowMap(name: string): Promise<Int32Array | null> {
    if (this._rowMapCache.has(name)) return this._rowMapCache.get(name)!;
    let map: Int32Array | null = null;
    const orderField = name + "_order";
    if (this.fields.has(orderField)) {
      const d = (await this.fieldDense(orderField)).data as ArrayLike<number>;
      map = d instanceof Int32Array ? d : Int32Array.from(d as any);
    }
    this._rowMapCache.set(name, map);
    return map;
  }

  /**
   * Rows of a CSR field for a cell selection, as a re-based CSR submatrix `{ data, indices, indptr,
   * rows }` (indptr over the requested rows, in input order). Selected rows are coalesced into runs so
   * the data/indices are read in a few merged byte-range requests instead of one per row — DE /
   * overdispersion on a small selection then touches a few MB of the cell-major copy, not the whole
   * matrix. Contiguous rows always merge (no waste); a `gap` of intervening elements is bridged to
   * amortize request latency (those extra bytes are fetched but not assembled into the output).
   */
  async csrRows(name: string, rowIndices: number[], gap = 4096, onProgress?: (done: number, total: number) => void): Promise<{ data: any; indices: any; indptr: Int32Array; rows: number[] }> {
    const base = "fields/" + name;
    if (rowIndices.length === 0) return { data: new Float64Array(0), indices: new Int32Array(0), indptr: Int32Array.from([0]), rows: [] };
    const nrows = (this.fields.get(name)!.shape as number[])[0];
    const indptr = await this._rangeOrSlice(base + "/indptr", 0, nrows + 1);   // small; read once per batch
    const at = (r: number) => Number(indptr[r]);

    // LOCALITY reorder: if this field's rows are stored in a permuted (e.g. Hilbert/cluster) order, map each requested
    // cell id -> its PHYSICAL row. The physical rows are what we coalesce + read, so a topologically-coherent selection
    // (a cluster, a lasso) becomes a few contiguous runs instead of thousands. Output stays in REQUESTED (cell) order,
    // so callers are unaffected. Identity when the field isn't reordered.
    const rowMap = await this._rowMap(name);
    const phys = rowMap ? rowIndices.map((c) => rowMap[c]) : rowIndices;

    // coalesce sorted-unique (physical) rows into runs (inclusive row ranges) to merge byte-range requests
    const uniq = [...new Set(phys)].sort((x, y) => x - y);
    const runs: { r0: number; r1: number }[] = [];
    for (const r of uniq) {
      const last = runs[runs.length - 1];
      if (last && at(r) - at(last.r1 + 1) <= gap) last.r1 = r;
      else runs.push({ r0: r, r1: r });
    }
    // one merged data/indices read per run; report progress as runs complete (done/total runs) for a fetch UI
    let done = 0;
    const blocks = await Promise.all(runs.map(async (run) => {
      const a = at(run.r0), b = at(run.r1 + 1);
      const data = b > a ? await this._rangeOrSlice(base + "/data", a, b) : new Float64Array(0);
      const indices = b > a ? await this._rangeOrSlice(base + "/indices", a, b) : new Int32Array(0);
      onProgress?.(++done, runs.length);
      return { r0: run.r0, r1: run.r1, a, data, indices };
    }));
    const blockOf = (r: number) => blocks.find((bl) => r >= bl.r0 && r <= bl.r1)!;

    // assemble the requested rows (input/cell order) into a re-based CSR submatrix — `phys` is each requested cell's
    // physical row (== the cell id when the field isn't reordered), so the output row order matches `rowIndices`.
    const outData: number[] = [], outIdx: number[] = [], outPtr: number[] = [0];
    for (const r of phys) {
      const a = at(r), b = at(r + 1), bl = blockOf(r), off = a - bl.a;
      for (let k = 0; k < b - a; k++) { outData.push(Number(bl.data[off + k])); outIdx.push(Number(bl.indices[off + k])); }
      outPtr.push(outData.length);
    }
    return { data: Float64Array.from(outData), indices: Int32Array.from(outIdx), indptr: Int32Array.from(outPtr), rows: rowIndices.slice() };
  }
}

/** Open an L* store and read its manifest. `store` is any zarrita-compatible store. */
export async function openLstar(store: any): Promise<LstarDataset> {
  return new LstarDataset(store).init();
}
