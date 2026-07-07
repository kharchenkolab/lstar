// WasmSource — the libzarr-backed read primitive for LstarDataset, replacing the zarrita reimplementation
// for whole-array + group-metadata reads. The SAME C++ core that R/Python use, compiled to WASM: it reads
// both v2 and v3 stores, so the viewer stops carrying a second, JS-only zarr decoder. The libzarr Store is
// synchronous, so this prefetches (async) the exact keys the reader will touch into a cache the WASM
// callback reads from (prefetch-then-sync-decode; no Asyncify). The byte-range streaming hot path stays in
// reader.ts on the raw store — this covers manifest, group attrs, and whole-array reads.
import createLstarIO from "../dist/lstar_io.mjs";
import type { LstarStore } from "./reader.ts";

// ONE shared Emscripten instance across every dataset. The module is stateless — pure chunk-key/shard/
// codec functions plus a per-source `Reader` bound to that source's store callback — so there is nothing
// to isolate per dataset. Instantiating `createLstarIO()` per open would allocate a fresh WASM heap each
// time and a long session that loads many datasets would accumulate them until it aborts; a singleton
// bounds it to one heap. Per-dataset native state (the `Reader`) is freed by `WasmSource.dispose()`.
let _modPromise: Promise<any> | null = null;
function ioModule(): Promise<any> { return (_modPromise ??= createLstarIO()); }

const TD = new TextDecoder();
const TE = new TextEncoder();

// Zarr LE dtype -> [typed-array constructor, itemsize] (the writer only emits these; assume LE).
const DTYPE: Record<string, [any, number]> = {
  "<f8": [Float64Array, 8], "<f4": [Float32Array, 4],
  "<i4": [Int32Array, 4], "<i8": [BigInt64Array, 8], "|i1": [Int8Array, 1],
  "<i2": [Int16Array, 2], "|u1": [Uint8Array, 1], "<u4": [Uint32Array, 4], "<u8": [BigUint64Array, 8],
};
function decodeTyped(bytes: Uint8Array, ctor: any, isize: number): any {
  const out = new ctor(Math.floor(bytes.byteLength / isize));
  new Uint8Array(out.buffer).set(bytes.subarray(0, out.byteLength));
  return out;
}

export class WasmSource {
  store: LstarStore;
  meta: Record<string, any> = {};        // consolidated array metadata (byte-range fast path reads this)
  private M: any;
  private reader: any;
  private cache = new Map<string, Uint8Array | null>();
  private haveConsolidated = false;      // true once .zmetadata / inline-v3 md is expanded into the cache

  constructor(store: LstarStore) { this.store = store; }

  async init(): Promise<void> {
    this.M = await ioModule();
    this.reader = new this.M.Reader((key: string) => this.cache.get(key) ?? null);
    await this._loadConsolidated();
  }

  // Free this dataset's native state: the embind `Reader` object (a C++ allocation inside the shared
  // module — NOT garbage-collected until deleted) and the prefetch byte cache. The shared WASM module
  // itself stays alive for other datasets. Idempotent; the source must not be used after disposing.
  dispose(): void {
    try { this.reader?.delete?.(); } catch { /* already deleted */ }
    this.reader = null;
    this.cache.clear();
  }

  // Fetch one object into the cache (once). A miss is cached as null so libzarr's "missing chunk = fill"
  // path is honoured without re-fetching.
  private async _put(key: string): Promise<void> {
    if (this.cache.has(key)) return;
    this.cache.set(key, (await this.store.get(key)) ?? null);
  }

  // Load the store's consolidated metadata and expand it into per-node cache entries, so a fresh per-node
  // open (Group::open / Array::open) reads its metadata straight from cache. v3 first (root zarr.json with
  // the inline convention), then v2 (.zmetadata). Also records array metadata in `meta` for the byte-range
  // fast path. Stores without consolidation still work — `_ensure` fetches per-node metadata lazily.
  private async _loadConsolidated(): Promise<void> {
    // v2 (.zmetadata) first — the current default format, so the common path is ONE store read and no
    // per-node metadata probes. Only if it's absent do we look for a v3 root zarr.json (whose inline
    // convention carries the whole tree). The two are mutually exclusive on disk.
    const zm = await this.store.get(".zmetadata");
    if (zm) {
      this.cache.set(".zmetadata", zm);
      const md = JSON.parse(TD.decode(zm)).metadata ?? {};
      for (const [key, val] of Object.entries(md as Record<string, any>)) {
        this.cache.set(key, TE.encode(JSON.stringify(val)));
        if (key.endsWith("/.zarray")) this.meta[key] = val;
      }
      this.haveConsolidated = true;
      return;
    }
    const zj = await this.store.get("zarr.json");
    if (zj) {
      this.cache.set("zarr.json", zj);
      const cm = JSON.parse(TD.decode(zj))?.consolidated_metadata?.metadata;
      if (cm) {
        for (const [node, m] of Object.entries(cm as Record<string, any>)) {
          this.cache.set(node + "/zarr.json", TE.encode(JSON.stringify(m)));
          if (m.node_type === "array") this.meta[node + "/zarr.json"] = m;
        }
        this.haveConsolidated = true;
      }
    }
  }

  // Make sure a node's metadata documents are in cache. When the store is consolidated they already are
  // (from _loadConsolidated), so this is a no-op — crucially, it must NOT probe absent keys (e.g. zarr.json
  // on a v2 group), which would fire a store read per node and defeat the one-read consolidated open. Only
  // a non-consolidated store falls through to per-node metadata fetches.
  private async _ensure(path: string): Promise<void> {
    if (this.haveConsolidated) return;
    const p = path ? path + "/" : "";
    await Promise.all([p + "zarr.json", p + ".zgroup", p + ".zattrs", p + ".zarray"].map((k) => this._put(k)));
  }

  // A group's L* attributes ("" = root; the manifest lives at root under "lstar").
  async groupLstar(path: string): Promise<any> {
    await this._ensure(path);
    return JSON.parse(this.reader.groupAttrs(path)).lstar;
  }

  // An array's shape from metadata ONLY — no chunk reads (used for axis lengths; a large offsets array
  // must not be materialized just to read its length).
  async arrayShape(path: string): Promise<number[]> {
    await this._ensure(path);
    return this.reader.shape(path) as number[];
  }

  // The byte-range fast path's descriptor for an array (dtype, itemsize, chunk shape, uncompressed) —
  // format-agnostic (v2/v3), computed by libzarr. Cached per path (a gene panel reads many columns of
  // the same array). Metadata only; no chunk reads.
  private _info = new Map<string, { dtype: string; itemsize: number; chunkShape: number[]; uncompressed: boolean; sharded: boolean }>();
  async arrayInfo(path: string): Promise<{ dtype: string; itemsize: number; chunkShape: number[]; uncompressed: boolean; sharded: boolean }> {
    let v = this._info.get(path);
    if (!v) { await this._ensure(path); v = this.reader.arrayInfo(path); this._info.set(path, v!); }
    return v!;
  }
  // The store key of one chunk in the array's own encoding (v2 "0" / v3 "c/0"). Sync — the metadata is
  // already cached (call after arrayInfo/_ensure). libzarr owns the chunk-key math.
  chunkKey(path: string, idx: number[]): string {
    return this.reader.chunkKey(path, idx);
  }

  // Sharded byte-range resolve (v3 sharding). Both sync (metadata already cached via arrayInfo/_ensure)
  // — libzarr owns the shard-index math; JS owns the two fetches (index, then chunk). shardLocate: which
  // shard object holds inner chunk `idx` (leaf key — JS prepends the array path, like chunkKey), the
  // chunk's slot, and the index layout (byte size + at-end). shardEntry: decode the fetched index bytes
  // -> the chunk's [offset, nbytes) within the shard (missing == a fill chunk).
  shardLocate(path: string, idx: number[]): { shardKey: string; intra: number; indexSize: number; indexAtEnd: boolean } {
    return this.reader.shardLocate(path, idx);
  }
  shardEntry(path: string, indexBytes: Uint8Array, intra: number): { offset: number; nbytes: number; missing: boolean } {
    return this.reader.shardEntry(path, indexBytes, intra);
  }

  // A whole array -> { data: TypedArray, shape } (matches zarrita's chunk shape/data contract).
  async array(path: string): Promise<{ data: any; shape: number[] }> {
    await this._ensure(path);
    const keys: string[] = this.reader.keysFor(path);       // metadata + every chunk key, libzarr's own scheme
    await Promise.all(keys.map((k) => this._put(k)));
    const a = this.reader.array(path);
    const [ctor, isize] = DTYPE[a.dtype] ?? [Uint8Array, 1];
    return { data: decodeTyped(a.bytes, ctor, isize), shape: a.shape as number[] };
  }
}
