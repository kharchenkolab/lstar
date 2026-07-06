// WasmSource — the libzarr-backed read primitive for LstarDataset, replacing the zarrita reimplementation
// for whole-array + group-metadata reads. The SAME C++ core that R/Python use, compiled to WASM: it reads
// both v2 and v3 stores, so the viewer stops carrying a second, JS-only zarr decoder. The libzarr Store is
// synchronous, so this prefetches (async) the exact keys the reader will touch into a cache the WASM
// callback reads from (prefetch-then-sync-decode; no Asyncify). The byte-range streaming hot path stays in
// reader.ts on the raw store — this covers manifest, group attrs, and whole-array reads.
import createLstarIO from "../dist/lstar_io.mjs";
import type { LstarStore } from "./reader.ts";

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

  constructor(store: LstarStore) { this.store = store; }

  async init(): Promise<void> {
    this.M = await createLstarIO();
    this.reader = new this.M.Reader((key: string) => this.cache.get(key) ?? null);
    await this._loadConsolidated();
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
    const zj = await this.store.get("zarr.json");
    if (zj) {
      this.cache.set("zarr.json", zj);
      const cm = JSON.parse(TD.decode(zj))?.consolidated_metadata?.metadata;
      if (cm) for (const [node, m] of Object.entries(cm as Record<string, any>)) {
        this.cache.set(node + "/zarr.json", TE.encode(JSON.stringify(m)));
        if (m.node_type === "array") this.meta[node + "/zarr.json"] = m;
      }
      return;
    }
    const zm = await this.store.get(".zmetadata");
    if (zm) {
      this.cache.set(".zmetadata", zm);
      const md = JSON.parse(TD.decode(zm)).metadata ?? {};
      for (const [key, val] of Object.entries(md as Record<string, any>)) {
        this.cache.set(key, TE.encode(JSON.stringify(val)));
        if (key.endsWith("/.zarray")) this.meta[key] = val;
      }
    }
  }

  // Make sure a node's metadata documents are in cache (no-op when consolidated md already populated them).
  private async _ensure(path: string): Promise<void> {
    const p = path ? path + "/" : "";
    await Promise.all([p + "zarr.json", p + ".zgroup", p + ".zattrs", p + ".zarray"].map((k) => this._put(k)));
  }

  // A group's L* attributes ("" = root; the manifest lives at root under "lstar").
  async groupLstar(path: string): Promise<any> {
    await this._ensure(path);
    return JSON.parse(this.reader.groupAttrs(path)).lstar;
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
