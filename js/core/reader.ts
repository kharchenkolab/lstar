// @lstar/core — a lazy reader for L* Zarr stores, over the libzarr WASM core (see wasm-source.ts).
//
// Opens a store (HTTP/local/zip), reads the consolidated L* metadata through the SAME libzarr core that
// backs R/Python (so v2 and v3 are read by one recipe, not a JS reimplementation), and exposes axes/
// fields with values fetched only when asked — including a single CSC gene-column (the viewer's hot
// path, a byte-range read straight off the store). The heavy numeric work belongs in the WASM kernels
// (see view.ts); this module is I/O + assembly.

/** Minimal store contract: fetch one object by key, or undefined if absent.
 * `getRange` is optional: a store that can serve a byte sub-range of an object (HTTP `Range`, a file
 * `pread`) enables the reader's sub-chunk fast path (one gene/cell = a few KB instead of the whole
 * uncompressed chunk). Stores without it still work — the reader falls back to whole-chunk reads. */
export interface LstarStore {
  get(key: string): Promise<Uint8Array | undefined>;
  getRange?(key: string, start: number, end: number): Promise<Uint8Array | undefined>;
  /** Read the LAST `n` bytes of an object (a suffix range), or undefined if absent. Enables the
   * byte-range fast path on SHARDED arrays: a v3 shard's index sits at the end of the shard object,
   * so the reader suffix-reads the index, then range-reads just the wanted chunk's bytes. Optional —
   * a store without it falls back to a whole-array read on sharded arrays (correct, not streamed). */
  getSuffix?(key: string, n: number): Promise<Uint8Array | undefined>;
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
  provenance?: any;     // free-form origin metadata (written by Python/C++/R; surfaced for inspection)
  coverage?: string;    // "full" | "partial" — partial fields carry an `index` array into `index_axis`
  index_axis?: string;  // for partial coverage: the axis the `index` positions refer to
}

const TD = new TextDecoder();

/** CSR arrays (nrows x ncols) -> CSC arrays. A one-time layout normalization on read (the heavy compute
 * kernels stay in WASM); lets every measure consumer work in CSC regardless of on-disk orientation. */
function csrToCscArrays(data: ArrayLike<number>, indices: ArrayLike<number>, indptr: ArrayLike<number>,
                        nrows: number, ncols: number): { data: Float64Array; indices: Int32Array; indptr: Int32Array } {
  const nnz = indptr[nrows] as number;
  const colptr = new Int32Array(ncols + 1);
  for (let k = 0; k < nnz; k++) colptr[(indices[k] as number) + 1]++;
  for (let c = 0; c < ncols; c++) colptr[c + 1] += colptr[c];
  const outData = new Float64Array(nnz), outIdx = new Int32Array(nnz);
  const next = Int32Array.from(colptr.subarray(0, ncols));            // running write cursor per column
  for (let r = 0; r < nrows; r++)
    for (let k = indptr[r] as number; k < (indptr[r + 1] as number); k++) {
      const col = indices[k] as number, dst = next[col]++;
      outData[dst] = data[k] as number; outIdx[dst] = r;
    }
  return { data: outData, indices: outIdx, indptr: colptr };
}

/** Dense (flat, C-order, nrows x ncols) -> CSC arrays, dropping zeros — matches scipy's `csc_matrix(dense)`,
 * so a dense primary measure yields identical stats to a native-sparse one. Column-major: one CSC column per gene. */
function denseToCscArrays(dense: ArrayLike<number | bigint>, nrows: number, ncols: number): { data: Float64Array; indices: Int32Array; indptr: Int32Array } {
  // `Number(...)` coerces an int64/uint64 dense field (a BigInt typed array) to number: assigning a BigInt
  // into the Float64 `data`, or comparing it to the number `0`, would otherwise throw / miscount zeros
  // (`0n !== 0` is true). No-op for a numeric source; matches the Float64 target either way.
  const indptr = new Int32Array(ncols + 1);
  for (let c = 0; c < ncols; c++) { let cnt = 0; for (let r = 0; r < nrows; r++) if (Number(dense[r * ncols + c]) !== 0) cnt++; indptr[c + 1] = indptr[c] + cnt; }
  const nnz = indptr[ncols];
  const data = new Float64Array(nnz), indices = new Int32Array(nnz);
  let w = 0;
  for (let c = 0; c < ncols; c++)
    for (let r = 0; r < nrows; r++) { const v = Number(dense[r * ncols + c]); if (v !== 0) { data[w] = v; indices[w] = r; w++; } }
  return { data, indices, indptr };
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
  store: any;                             // == rawStore (kept for the byte-range fast path's getRange)
  meta: Record<string, any> = {};         // per-array metadata (from the source's consolidated md; byte-range fast path)
  kind = "sample";
  specVersion = "0.1";
  profiles: string[] = [];
  dropped: string[] = [];
  auxNames: string[] = [];   // namespaces of the lossless passthrough subtree (uns/@misc)
  axes = new Map<string, AxisMeta>();
  fields = new Map<string, FieldMeta>();
  src: any;   // the libzarr/WASM read source — whole-array + group-metadata reads go through it (the one reader)

  constructor(store: any, src: any) {
    this.rawStore = store;
    this.store = store;          // == the underlying store; the byte-range fast path reads chunks off it
    this.src = src;
  }

  private async _get(p: string): Promise<any> {   // a whole array via libzarr -> { data, shape }
    return this.src.array(p);
  }
  // A group's L* attributes ("" = root); libzarr reads both v2 (.zattrs) and v3 (zarr.json).
  private async _groupLstar(path: string): Promise<any> {
    return this.src.groupLstar(path);
  }
  // A 1-D array's length — used for axis lengths from labels_offsets. Metadata only (no chunk read).
  private async _arrayShape(p: string): Promise<number[]> {
    return this.src.arrayShape(p);
  }

  async init(): Promise<this> {
    // The libzarr source loads + expands the store's consolidated metadata once (into `src.meta`, which
    // the byte-range fast path reads for per-array dtype/chunks), then serves every group/array metadata
    // read through the WASM core — one read up front instead of ~80, and reads both v2 and v3.
    await this.src.init();
    this.meta = this.src.meta ?? {};
    const m = await this._groupLstar("");
    if (!m) throw new Error("not an L* store at this location — no L* manifest ('lstar' root attribute; check the store URL / ?store= path).");
    this.kind = m.kind ?? "sample";
    this.specVersion = m.spec_version ?? "0.1";
    this.profiles = m.profiles ?? [];
    this.dropped = m.dropped ?? [];
    this.auxNames = m.passthrough ?? [];
    for (const name of m.axes as string[]) {
      const lm = (await this._groupLstar("axes/" + name)) ?? {};
      const shape = await this._arrayShape("axes/" + name + "/labels_offsets");
      this.axes.set(name, { name, origin: lm.origin ?? "observed", role: lm.role,
                            induced_by: lm.induced_by ?? undefined,
                            length: (shape[0] as number) - 1 });
    }
    for (const name of m.fields as string[]) {
      const lm = (await this._groupLstar("fields/" + name)) ?? {};
      this.fields.set(name, { name, role: lm.role, span: lm.span ?? [], encoding: lm.encoding,
                              state: lm.state ?? undefined, subtype: lm.subtype ?? undefined,
                              shape: lm.shape, ordered: lm.ordered,
                              nullable: lm.nullable ?? undefined,
                              provenance: lm.provenance, coverage: lm.coverage,
                              index_axis: lm.index_axis } as any);
    }
    return this;
  }

  /** Release the dataset's native reader state (the WASM `Reader` + its byte cache). The shared WASM
   * module stays alive for other datasets. Call this when done with a dataset in a long-lived session
   * that opens many — otherwise each dataset's C++ `Reader` accumulates in the module heap. Idempotent;
   * the dataset must not be read after disposing. A no-op for a source without `dispose` (e.g. a mock). */
  dispose(): void { this.src?.dispose?.(); }

  axisNames(): string[] { return [...this.axes.keys()]; }
  fieldNames(): string[] { return [...this.fields.keys()]; }
  hasField(name: string): boolean { return this.fields.has(name); }
  field(name: string): FieldMeta | undefined { return this.fields.get(name); }
  axisLength(name: string): number { return this.axes.get(name)?.length ?? 0; }

  async axisLabels(name: string): Promise<string[]> {
    // labels + labels_offsets are independent (both just read a chunk; the array metadata is already loaded) —
    // read CONCURRENTLY so a labelled axis costs ONE round-trip, not two serial ones.
    const [bytes, offs] = await Promise.all([
      this._get("axes/" + name + "/labels"), this._get("axes/" + name + "/labels_offsets"),
    ]);
    return decodeStrings(bytes.data as Uint8Array, offs.data as any);
  }

  /** A dense field's values (flat, C-order) + shape. */
  async fieldDense(name: string): Promise<{ data: any; shape: number[] }> {
    const c = await this._get("fields/" + name + "/values");
    return { data: c.data, shape: c.shape as number[] };
  }

  /** A partial-coverage field's `index` array (0-based positions into `index_axis`), or undefined if the
   * field has no `index` (full coverage). Same read path as fieldDense — surfaces coverage for a consumer. */
  async fieldIndex(name: string): Promise<any | undefined> {
    const c = await this._get("fields/" + name + "/index").catch(() => undefined);
    return c?.data;
  }

  /** A utf8 (label) field's string values. */
  async fieldStrings(name: string): Promise<string[]> {
    // values + values_offsets are independent — read CONCURRENTLY so a label field costs ONE round-trip, not
    // two serial ones. On a hosted store this is on the viewer's first-paint critical path (the cluster colour
    // + every facet column); the serial reads made "clusters wait" for a second round-trip that had no dependency.
    const [bytes, offs] = await Promise.all([
      this._get("fields/" + name + "/values"), this._get("fields/" + name + "/values_offsets"),
    ]);
    return decodeStrings(bytes.data as Uint8Array, offs.data as any);
  }

  /**
   * A lossless-passthrough namespace (e.g. "anndata.uns"), reconstructed into a plain JS object: the
   * stored JSON `tree` with its array-leaf references resolved (a numpy structured array becomes a
   * `{field: array}` object). Read-only — for inspection / promoting a recognized structure to a field.
   */
  async aux(ns: string): Promise<any> {
    const lm = (await this._groupLstar("passthrough/" + ns)) ?? {};
    const tree = typeof lm.tree === "string" ? JSON.parse(lm.tree) : lm.tree;
    const byId: Record<string, any> = {};
    for (const a of (lm.arrays ?? []) as Array<{ id: string; kind: string }>) {
      if (a.kind === "utf8") {
        const [bytes, offs] = await Promise.all([
          this._get("passthrough/" + ns + "/" + a.id), this._get("passthrough/" + ns + "/" + a.id + "_offsets"),
        ]);
        byId[a.id] = decodeStrings(bytes.data as Uint8Array, offs.data as any);
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
    // codes + categories + categories_offsets are independent — read CONCURRENTLY (one round-trip, not three serial).
    const [codesC, bytesC, offsC] = await Promise.all([
      this._get("fields/" + name + "/codes"),
      this._get("fields/" + name + "/categories"),
      this._get("fields/" + name + "/categories_offsets"),
    ]);
    return { codes: Int32Array.from(codesC.data as any), categories: decodeStrings(bytesC.data as Uint8Array, offsC.data as any),
             ordered: !!(this.fields.get(name) as any)?.ordered };
  }

  /** A sparse (csc/csr) field's raw arrays + shape. Reads the whole field. */
  async fieldSparse(name: string): Promise<{ data: any; indices: any; indptr: any; shape: number[]; fmt: string }> {
    const meta = this.fields.get(name)!;
    const data = (await this._get("fields/" + name + "/data")).data;
    const indices = (await this._get("fields/" + name + "/indices")).data;
    const indptr = (await this._get("fields/" + name + "/indptr")).data;
    return { data, indices, indptr, shape: meta.shape as number[], fmt: meta.encoding };
  }

  /** A 2-D measure as CSC arrays, whatever its on-disk encoding — dense (`/values`), CSR, or CSC. The
   * SINGLE point where a measure is read for compute, so no consumer (extendForViewer, view.colStats,
   * view.subsampleDE) special-cases encoding and can diverge — a DENSE primary measure (SCE logcounts /
   * scaled AnnData X) lives at `/values`, not the sparse `/data`, and reading it as sparse threw
   * NotFoundError. Mirrors Python `X.tocsc() if issparse else csc_matrix(X)` and R's Matrix coercion. */
  async fieldAsCsc(name: string): Promise<{ data: any; indices: any; indptr: any; shape: number[] }> {
    if (this.fields.get(name)?.encoding === "dense") {
      const dv = await this.fieldDense(name);
      const [nr, nc] = dv.shape as number[];
      return { ...denseToCscArrays(dv.data, nr, nc), shape: dv.shape };
    }
    const sp = await this.fieldSparse(name);
    if (sp.fmt === "csc") return { data: sp.data, indices: sp.indices, indptr: sp.indptr, shape: sp.shape };
    if (sp.fmt === "csr") {
      const [nr, nc] = sp.shape;
      return { ...csrToCscArrays(sp.data, sp.indices, sp.indptr, nr, nc), shape: sp.shape };
    }
    throw new Error("fieldAsCsc: `" + name + "` must be dense, CSC, or CSR (2-D cells x genes), got " + sp.fmt);
  }

  // Bounded LRU read cache. Re-running compute over the SAME cells re-issues the SAME byte-range reads, so cache the
  // decoded slices keyed by (array, lo, hi). Holds only RANGE reads (csrRows blocks, gene columns) — the whole-matrix
  // fieldSparse path uses _get, not this — so it stays bounded and never pins the full matrix. LRU-evicted to a budget.
  private _rcache = new Map<string, { v: any; bytes: number }>();
  private _rcacheBytes = 0;
  private _rcacheBudget = 256 * 1024 * 1024;   // 256MB default
  /** Cap (bytes) for the in-memory read cache; LRU-evicts down to it. 0 disables + clears it. */
  setReadCacheBudget(bytes: number): void { this._rcacheBudget = Math.max(0, Math.floor(bytes)); if (this._rcacheBudget === 0) { this._rcache.clear(); this._rcacheBytes = 0; } else this._evictReadCache(); }
  private _csrConcurrency = 12;   // max simultaneous range reads in csrRows — bounds connection use on a scattered (non-coalescing) selection
  /** Max in-flight range reads per csrRows batch (default 12). Lower for a flaky/limited server, higher for HTTP/2. */
  setCsrConcurrency(n: number): void { this._csrConcurrency = Math.max(1, Math.floor(n)); }
  readCacheStats(): { entries: number; bytes: number; budget: number } { return { entries: this._rcache.size, bytes: this._rcacheBytes, budget: this._rcacheBudget }; }
  private _evictReadCache(): void {
    while (this._rcacheBytes > this._rcacheBudget && this._rcache.size) {
      const k = this._rcache.keys().next().value as string; const e = this._rcache.get(k)!; this._rcache.delete(k); this._rcacheBytes -= e.bytes;
    }
  }
  // cache wrapper around the range reader: LRU-bump on hit, store + evict on miss (skip a slice too big to ever fit).
  private async _rangeOrSlice(arrPath: string, lo: number, hi: number): Promise<any> {
    const ck = arrPath + ":" + lo + ":" + hi;
    const hit = this._rcache.get(ck);
    if (hit) { this._rcache.delete(ck); this._rcache.set(ck, hit); return hit.v; }
    const v = await this._readRange(arrPath, lo, hi);
    const bytes = v?.byteLength || 0;
    if (this._rcacheBudget > 0 && bytes > 0 && bytes <= this._rcacheBudget) { this._rcache.set(ck, { v, bytes }); this._rcacheBytes += bytes; this._evictReadCache(); }
    return v;
  }

  /**
   * Read elements [lo, hi) of a 1-D array as a typed array. Two chunk-granular fast paths, both on v2+v3
   * (libzarr `arrayInfo` reports dtype/chunk-shape/raw-ness/sharded from either format; the fetch loop
   * stays in JS, the Zarr interpretation in libzarr):
   *   • UNCOMPRESSED + `getRange` + single covering chunk → ONE exact byte-range read of that chunk (no
   *     decode) — a gene/cell is a few KB. Sharded: resolved through the shard index (`_readShardedChunkRange`).
   *   • COMPRESSED → fetch + decode only the chunk(s) covering [lo,hi) (`_readCompressedRange`) — O(chunk),
   *     not O(array). Removes the old "compressed ⇒ whole-array" cliff, so a compressed (+sharded) store
   *     stays byte-range-readable; decoded chunks are cached (viewer chunk-locality).
   * Anything without a fast path (0-d, a fill/missing chunk, an unrangeable sharded store) falls through to
   * a whole-array read + slice via libzarr. Edge (partial last) chunks are stored full-size + fill-padded,
   * so a sub-range is always in-bounds within the chunk. Equivalent results either way.
   */
  private async _readRange(arrPath: string, lo: number, hi: number): Promise<any> {
    if (hi > lo) {
      const info = await this.src.arrayInfo(arrPath);      // {dtype, itemsize, chunkShape, uncompressed, sharded} — v2 or v3
      const dt = DTYPE[info.dtype];
      const oneD = info.chunkShape.length === 1 && info.chunkShape[0] > 0;
      if (dt && oneD) {
        const [ctor, isize] = dt;
        const chunkLen = info.chunkShape[0];
        const ci = Math.floor(lo / chunkLen);
        if (info.uncompressed) {
          // exact byte range of the covering chunk (no decode) — needs getRange + a single covering chunk
          if (typeof this.store.getRange === "function" && hi <= (ci + 1) * chunkLen) {
            const off = lo - ci * chunkLen;
            if (!info.sharded) {
              const key = this.src.chunkKey(arrPath, [ci]);  // libzarr gives the leaf key; join under arrPath
              const bytes = await this.store.getRange(arrPath + "/" + key, off * isize, (off + (hi - lo)) * isize);
              if (bytes) return decodeTyped(bytes, ctor, isize);
            } else {
              const got = await this._readShardedChunkRange(arrPath, ci, off, hi - lo, ctor, isize);
              if (got) return got;
            }
          }
        } else {
          const got = await this._readCompressedRange(arrPath, lo, hi, ci, chunkLen, ctor, info.sharded);
          if (got) return got;
        }
      }
    }
    const whole = (await this.src.array(arrPath)).data;    // no fast path (fill chunk, 0-d, decode miss) -> whole-array + slice
    return (whole as any).slice(lo, hi);
  }

  /**
   * Compressed byte-range: decode the chunk(s) covering [lo,hi) via `readChunkDecoded` (which fetches only
   * the covering chunk object / the chunk's bytes inside its shard, and caches the decoded chunk) and
   * return exactly [lo,hi). One covering chunk is the hot case; a span crossing chunk boundaries decodes
   * each covered chunk and concatenates — still O(covered chunks), never the whole array. Returns undefined
   * (→ whole-array fallback) if any covering chunk is a fill/missing chunk or the store can't range-read it.
   */
  private async _readCompressedRange(arrPath: string, lo: number, hi: number, ci: number, chunkLen: number,
                                     ctor: any, sharded: boolean): Promise<any | undefined> {
    const cj = Math.floor((hi - 1) / chunkLen);
    if (ci === cj) {
      const chunk = await this.src.readChunkDecoded(arrPath, ci, sharded);
      return chunk ? chunk.slice(lo - ci * chunkLen, hi - ci * chunkLen) : undefined;
    }
    const out = new ctor(hi - lo);
    for (let c = ci; c <= cj; c++) {
      const chunk = await this.src.readChunkDecoded(arrPath, c, sharded);
      if (!chunk) return undefined;                        // fill/missing or unrangeable -> whole-array fallback
      const cStart = c * chunkLen;
      const a = Math.max(lo, cStart), b = Math.min(hi, cStart + chunkLen);
      out.set(chunk.subarray(a - cStart, b - cStart), a - lo);
    }
    return out;
  }

  /**
   * Stream `n` elements at offset `off` within inner chunk `ci` of a SHARDED (v3) array, reading only
   * the shard's index + the chunk's own bytes — not the whole shard. Resolves through libzarr's shard
   * math: `shardLocate` names the shard object + index layout; a suffix (or leading) read fetches the
   * small index; `shardEntry` decodes it to the chunk's [offset, nbytes) within the shard; one final
   * `getRange` reads the wanted elements. Returns undefined (→ caller falls back to a whole-array read,
   * always correct) when the store can't suffix-read an end-located index, a fetch misses, or the chunk
   * is a fill (missing) chunk.
   */
  private async _readShardedChunkRange(arrPath: string, ci: number, off: number, n: number, ctor: any, isize: number): Promise<any | undefined> {
    const loc = this.src.shardLocate(arrPath, [ci]);        // { shardKey (leaf), intra, indexSize, indexAtEnd }
    const shardKey = arrPath + "/" + loc.shardKey;
    let idx: Uint8Array | undefined;
    if (loc.indexAtEnd) {
      if (typeof this.store.getSuffix !== "function") return undefined;   // no suffix read -> whole-array fallback
      idx = await this.store.getSuffix(shardKey, loc.indexSize);
    } else {
      idx = await this.store.getRange!(shardKey, 0, loc.indexSize);
    }
    if (!idx || idx.byteLength < loc.indexSize) return undefined;
    const e = this.src.shardEntry(arrPath, idx, loc.intra);  // { offset, nbytes, missing }
    if (e.missing) return undefined;
    const start = e.offset + off * isize;
    const bytes = await this.store.getRange!(shardKey, start, start + n * isize);
    return bytes ? decodeTyped(bytes, ctor, isize) : undefined;
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
   * Many CSC columns in ONE batched read — the column-major twin of {@link csrRows}, for a gene panel
   * (a dotplot/heatmap's marker columns, or a subset recompute). The requested columns' byte spans are
   * COALESCED into a few merged range requests (bounded by `_csrConcurrency`) instead of one request per
   * column — so N gene columns cost a handful of parallel reads rather than N reads into the SAME chunk
   * object, which a browser serializes on its per-URL cache lock. Each returned column is byte-identical
   * to `cscColumn(name, col)`; results are in INPUT order (a repeated column re-slices the same block).
   * `gap` bridges small inter-column gaps to amortize latency (extra bytes fetched, not assembled).
   */
  async cscColumns(name: string, cols: number[], gap = 4096): Promise<{ cols: { rows: any; vals: any }[] }> {
    const base = "fields/" + name;
    if (cols.length === 0) return { cols: [] };
    const ncols = (this.fields.get(name)!.shape as number[])[1];
    const indptr = await this._rangeOrSlice(base + "/indptr", 0, ncols + 1);   // small; read once per batch
    const at = (c: number) => Number(indptr[c]);
    // coalesce the sorted-unique columns into runs (contiguous column ranges) so nearby columns' data/indices
    // are read in a few merged byte-range requests; a large gap between columns starts a new run.
    const uniq = [...new Set(cols)].sort((x, y) => x - y);
    const runs: { c0: number; c1: number }[] = [];
    for (const c of uniq) {
      const last = runs[runs.length - 1];
      if (last && at(c) - at(last.c1 + 1) <= gap) last.c1 = c;
      else runs.push({ c0: c, c1: c });
    }
    // one merged data/indices read per run, bounded concurrency (a scattered gene set doesn't coalesce, so
    // `runs` can be large — cap in-flight reads instead of exhausting the connection pool, exactly like csrRows).
    const blocks: { c0: number; c1: number; a: number; indices: any; data: any }[] = new Array(runs.length);
    let next = 0;
    const worker = async () => {
      while (next < runs.length) {
        const i = next++; const run = runs[i];
        const a = at(run.c0), b = at(run.c1 + 1);
        blocks[i] = { c0: run.c0, c1: run.c1, a,
          indices: b > a ? await this._rangeOrSlice(base + "/indices", a, b) : new Int32Array(0),
          data:    b > a ? await this._rangeOrSlice(base + "/data", a, b)    : new Float64Array(0) };
      }
    };
    await Promise.all(Array.from({ length: Math.min(this._csrConcurrency, runs.length) }, worker));
    const blockOf = (c: number) => blocks.find((bl) => c >= bl.c0 && c <= bl.c1)!;
    // slice each requested column back out (INPUT order) — identical bytes to cscColumn(name, col).
    const out = cols.map((c) => {
      const a = at(c), b = at(c + 1);
      if (b <= a) return { rows: new Int32Array(0), vals: new Float64Array(0) };
      const bl = blockOf(c), off = a - bl.a, n = b - a;
      return { rows: (bl.indices as any).slice(off, off + n), vals: (bl.data as any).slice(off, off + n) };
    });
    return { cols: out };
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
   * Locality-aware subsample of `cellIds` down to ~`maxRows` cell ids — for an approximate (ranking-grade) compute that
   * doesn't need every cell. On a REORDERED field, pick `windows` evenly-spaced CONTIGUOUS runs of physical rows, so the
   * subsequent csrRows read is ~`windows` coalesced reads (a few MB) rather than the whole selection; the windows are
   * spread across the selection's physical extent (≈ its spatial extent under a Hilbert order) to stay representative.
   * On a canonical field (no locality) returns a uniform stride sample. Returns all ids unchanged when ≤ maxRows.
   */
  async sampleRows(name: string, cellIds: number[], maxRows: number, windows = 16): Promise<number[]> {
    if (!(maxRows > 0) || cellIds.length <= maxRows) return cellIds.slice();
    const map = await this._rowMap(name);
    if (!map) { const out: number[] = [], step = cellIds.length / maxRows; for (let i = 0; i < maxRows; i++) out.push(cellIds[Math.floor(i * step)]); return out; }
    const pairs = cellIds.map((c) => [map[c], c] as [number, number]).sort((a, b) => a[0] - b[0]);
    const W = Math.max(1, Math.min(windows, maxRows)), perWin = Math.ceil(maxRows / W), seg = pairs.length / W;
    const out: number[] = [];
    for (let w = 0; w < W; w++) { const start = Math.floor(w * seg); for (let i = 0; i < perWin && start + i < pairs.length; i++) out.push(pairs[start + i][1]); }
    return out;
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
    // one merged data/indices read per run; report progress as runs complete (done/total runs) for a fetch UI.
    // BOUNDED concurrency: a selection that's ORTHOGONAL to the row order (a sample/condition across a cluster reorder,
    // or a sparse scatter) doesn't coalesce, so `runs` can be in the thousands. Firing them all at once exhausts the
    // browser/server connection pool and the whole batch fails fast ("Failed to fetch"). Cap in-flight reads instead;
    // a worst-case scattered read then just streams in waves rather than crashing.
    let done = 0;
    const blocks: { r0: number; r1: number; a: number; data: any; indices: any }[] = new Array(runs.length);
    let next = 0;
    const worker = async () => {
      while (next < runs.length) {
        const i = next++; const run = runs[i];
        const a = at(run.r0), b = at(run.r1 + 1);
        const data = b > a ? await this._rangeOrSlice(base + "/data", a, b) : new Float64Array(0);
        const indices = b > a ? await this._rangeOrSlice(base + "/indices", a, b) : new Int32Array(0);
        blocks[i] = { r0: run.r0, r1: run.r1, a, data, indices };
        onProgress?.(++done, runs.length);
      }
    };
    await Promise.all(Array.from({ length: Math.min(this._csrConcurrency, runs.length) }, worker));
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

/**
 * Open an L* store and read its manifest, through the libzarr WASM core — the SAME reader R and Python
 * use, so v2 and v3 stores are read by one recipe (no JS zarr reimplementation). `store` is any object
 * with the async `get(key)` / optional `getRange(key,start,end)` contract (see {@link LstarStore},
 * {@link NodeFSStore}, {@link HttpStore}). Whole-array + metadata reads go through libzarr; the byte-range
 * streaming hot path (cscColumn/csrRow) reads chunk sub-ranges straight off the store's `getRange`.
 */
export async function openLstar(store: any): Promise<LstarDataset> {
  const { WasmSource } = await import("./wasm-source.ts");
  return new LstarDataset(store, new WasmSource(store)).init();
}

/** @deprecated libzarr is now the default reader — use {@link openLstar}. Kept as an alias. */
export const openLstarWasm = openLstar;
