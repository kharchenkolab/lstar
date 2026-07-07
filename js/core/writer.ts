// @lstar/core — write L* Zarr stores (the write side of reader.ts). Emits Zarr v2 (single chunk
// per array, uncompressed — readable by zarrita and by the R/Python lstar readers). Use `writeStore`
// for a full dataset, or `addToStore` to append derived fields/axes to an existing store and update
// the root manifest (what the pagoda3 prep does — extend a base store with viewer navigators).
//
// Chunking/compression are intentionally simple here (single uncompressed chunk); the base store's
// own fields keep their original chunking when `addToStore` leaves them untouched.

/** Write side of the store contract: put one object by key (mkdir-p handled by the store).
 * `delete` is optional and used only to drop stale consolidated metadata after `addToStore`. */
export interface LstarWritableStore {
  get(key: string): Promise<Uint8Array | undefined>;
  set(key: string, value: Uint8Array): Promise<void>;
  delete?(key: string): Promise<void>;
}

export interface AxisSpec {
  labels: (string | number)[]; origin?: string; role?: string;
  inducedBy?: string;   // for a derived `factor` axis: the categorical field whose categories it mirrors
}
export interface FieldSpec {
  role?: string; span: string[]; encoding: "dense" | "csc" | "csr" | "utf8" | "categorical";
  state?: string | null; subtype?: string | null; shape?: number[] | null;
  // dense: `data` (flat C-order) + `shape`; utf8: `values` (strings); csc/csr: data+indices+indptr + `shape`;
  // categorical: `codes` (i4, -1=missing) + `categories` (strings) + `ordered`.
  data?: ArrayLike<number>; values?: (string | number)[];
  indices?: ArrayLike<number>; indptr?: ArrayLike<number>;
  codes?: ArrayLike<number>; categories?: (string | number)[]; ordered?: boolean;
  mask?: Uint8Array;                          // nullable: 1 == missing (over span[0])
  index?: BigInt64Array | ArrayLike<number>;  // partial coverage: 0-based positions into `indexAxis`
  indexAxis?: string;
  write?: WriteOptions;                       // per-field chunk/compressor override (else the call-level opts)
  provenance?: Record<string, unknown>;       // free-form metadata persisted to the field's lstar attrs
}
// Aux passthrough (the lossless uns/@misc subtree): an opaque JSON `tree` plus the array leaves it
// references by id. Leaf grammar in `tree`: {"$array":id} for a dense leaf, {"$strings":id,...} for a
// utf8 leaf (the rest -- $obj/$list/scalars -- are inline). The writer stores `tree` verbatim as a string.
export interface AuxArraySpec { id: string; kind: "dense" | "utf8"; data?: ArrayLike<number>; values?: (string | number)[]; }
export interface AuxSpec { tree: unknown; arrays: AuxArraySpec[]; }
export interface DatasetSpec {
  kind?: string; specVersion?: string; profiles?: string[]; dropped?: string[];
  axes: Record<string, AxisSpec>; fields: Record<string, FieldSpec>;
  aux?: Record<string, AuxSpec>;
}

/** A chunk compressor (dependency-injected so writer.ts never imports a WASM module — the caller wires
 * `compress`, e.g. the WASM `gzipCompress` or a pure-JS gzip). `id`/`level` go into the codec metadata.
 * With a `codec` (the WASM writer, below) `compress` is unused — `codec.encodeChunk` handles gzip AND
 * zstd by `id`; without one, `compress` is applied per chunk (gzip only, the pre-WASM path). */
export interface Compressor { id: "gzip" | "zstd"; level: number; compress?: (raw: Uint8Array) => Uint8Array; }
/** The libzarr write-side pure functions (from the lstar_writer WASM module), dependency-injected so the
 * drift-prone bytes stay in libzarr: `encodeChunk` runs a chunk through [bytes, gzip|zstd]; `packShard`
 * assembles a shard object from slot-ordered encoded chunks (index + crc32c). Same core the reader
 * decodes with, so a JS-written store is byte-consistent with the C++/Python writers. */
export interface WriterCodec {
  encodeChunk: (raw: Uint8Array, compressor: "none" | "gzip" | "zstd", level: number) => Uint8Array;
  packShard: (entries: (Uint8Array | null)[], indexAtEnd: boolean) => Uint8Array;
}
/** Write knobs: `chunkElems` splits arrays multi-chunk along axis 0; `compressor` compresses each chunk;
 * `shardElems` (v3 only, needs `codec`) packs ~that many elements' chunks into each shard object; `codec`
 * is the injected WASM writer (enables zstd + sharding). All omitted -> a single uncompressed chunk per
 * array (the default, byte-identical to before). */
export interface WriteOptions { chunkElems?: number; compressor?: Compressor | null; shardElems?: number; codec?: WriterCodec; }

const ENC = new TextEncoder();
const TD = new TextDecoder();
const j = (o: unknown) => ENC.encode(JSON.stringify(o));

/** On-disk Zarr format. v2 (default, byte-identical to before) writes .zarray/.zgroup/.zattrs + a
 * consolidated .zmetadata; v3 writes per-node zarr.json + a root zarr.json carrying the manifest AND
 * an inline consolidated_metadata map, with chunk keys under `c/`. Both are read by the C++/R/JS/Python
 * libzarr cores. `addToStore` writes whichever format the base store already uses (auto-detected), so
 * appending viewer navigators never mixes formats. */
export type Fmt = "v2" | "v3";
// v2 little-endian dtype string -> Zarr v3 `data_type` name (endian moves into the `bytes` codec).
const V3_DTYPE: Record<string, string> = {
  "<f8": "float64", "<f4": "float32", "<i4": "int32", "<i8": "int64",
  "|i1": "int8", "|u1": "uint8", "<i2": "int16", "<u4": "uint32", "<u8": "uint64",
};

/** Order consolidated-metadata keys parent-before-child (by path depth, then lexicographically) so a
 * strict consolidated reader can nest them. zarr-python 3's parser walks the flat map top-down and
 * assumes every parent group appears before its children (and that siblings group consecutively); an
 * `addToStore` that simply appended new keys broke that. zarr v2 readers are order-insensitive, so
 * this is purely additive robustness that keeps the one-request consolidated open working everywhere. */
function sortedMeta(meta: Record<string, unknown>): Record<string, unknown> {
  const depth = (k: string) => (k.match(/\//g) || []).length;
  const out: Record<string, unknown> = {};
  for (const k of Object.keys(meta).sort((a, b) => depth(a) - depth(b) || (a < b ? -1 : a > b ? 1 : 0))) out[k] = meta[k];
  return out;
}

function dtypeOf(a: ArrayLike<number>): [string, number, Uint8Array] {
  if (a instanceof Float64Array) return ["<f8", 8, new Uint8Array(a.buffer, a.byteOffset, a.byteLength)];
  if (a instanceof Float32Array) return ["<f4", 4, new Uint8Array(a.buffer, a.byteOffset, a.byteLength)];
  if (a instanceof Int32Array) return ["<i4", 4, new Uint8Array(a.buffer, a.byteOffset, a.byteLength)];
  if (a instanceof Int8Array) return ["|i1", 1, new Uint8Array(a.buffer, a.byteOffset, a.byteLength)];
  if (a instanceof Uint8Array) return ["|u1", 1, a];
  if (a instanceof BigInt64Array) return ["<i8", 8, new Uint8Array(a.buffer, a.byteOffset, a.byteLength)];
  // a plain number[] or other -> materialize as f8
  const f = Float64Array.from(a as ArrayLike<number>);
  return ["<f8", 8, new Uint8Array(f.buffer)];
}

// Chunk shape: split ONLY along axis 0 (so each chunk is a contiguous C-order byte slice -- mirrors the
// C++ write_array / Python _chunks_for). `chunkElems` is a target element budget per chunk; the inner
// dims stay whole. Omitted -> the whole array is one chunk.
function chunkShapeFor(shape: number[], chunkElems?: number): number[] {
  if (!shape.length) return [1];
  if (!chunkElems || chunkElems <= 0) return shape.slice();
  const inner = shape.slice(1).reduce((p, c) => p * c, 1);
  const rows = Math.max(1, Math.floor(chunkElems / Math.max(1, inner)));
  return rows >= shape[0] ? shape.slice() : [rows, ...shape.slice(1)];
}

// Shard (outer) shape: pack k WHOLE inner chunks along axis 0 into one shard object (~shardElems elements).
// Mirrors the C++ shard_shape_for / Python _shards_for. null (unsharded) for a 0-d/unchunked array or a
// single-chunk one (nothing to pack).
function shardShapeFor(shape: number[], chunks: number[], shardElems?: number): number[] | null {
  if (!shardElems || shardElems <= 0 || !shape.length) return null;
  const chunkRows = chunks[0];
  if (!(chunkRows > 0) || chunkRows >= shape[0]) return null;    // single chunk -> nothing to shard
  const inner = shape.slice(1).reduce((p, c) => p * c, 1);
  const numChunks = Math.ceil(shape[0] / chunkRows);
  const targetRows = Math.floor(shardElems / Math.max(1, inner));
  const k = Math.max(1, Math.min(numChunks, Math.floor(targetRows / chunkRows)));
  return [k * chunkRows, ...shape.slice(1)];
}

// One Zarr array, chunked along axis 0 and optionally compressed. Each chunk is padded to the full chunk
// size with fill_value 0 BEFORE compression (edge chunks are full-size; the reader trims via `shape`). With
// no opts this writes a single uncompressed chunk -- byte-identical to the old writer. v3 + `shardElems` +
// an injected `codec` packs the inner chunks into shard objects (sharding_indexed). The drift-prone bytes
// (compression, shard index + crc32c) go through the injected libzarr codec; JS owns layout + keys + meta.
async function writeArray(store: LstarWritableStore, base: string, a: ArrayLike<number>, shape: number[],
                          opts?: WriteOptions, fmt: Fmt = "v3"): Promise<void> {
  const [dtype, isize, raw] = dtypeOf(a);
  const sh = shape.length ? shape : [a.length];
  const chunks = chunkShapeFor(sh, opts?.chunkElems).map((s) => Math.max(s, 1));
  const compressor = opts?.compressor ?? null;
  const codec = opts?.codec;
  // encode one (padded) chunk: via the injected WASM codec (gzip/zstd, same core the reader decodes with)
  // if present, else the pre-WASM injected gzip `compress`, else raw bytes.
  const enc = (block: Uint8Array): Uint8Array =>
    codec ? codec.encodeChunk(block, compressor ? compressor.id : "none", compressor?.level ?? 5)
    : compressor?.compress ? compressor.compress(block) : block;
  const inner = chunks.slice(1).reduce((p, c) => p * c, 1);     // elements per axis-0 row (inner dims)
  const chunkBytes = chunks[0] * inner * isize;
  const nRows = sh[0] ?? 0;
  const nChunks = nRows <= 0 ? 1 : Math.ceil(nRows / chunks[0]);
  const chunkBlock = (ci: number): Uint8Array => {              // padded raw bytes of inner chunk `ci`
    const block = new Uint8Array(chunkBytes);
    block.set(raw.subarray(ci * chunkBytes, Math.min((ci + 1) * chunkBytes, raw.byteLength)));
    return block;
  };

  // v3 SHARDED path: pack inner chunks into shard objects (fewer files, still range-readable via the index).
  const shardShape = fmt === "v3" && codec ? shardShapeFor(sh, chunks, opts?.shardElems) : null;
  if (shardShape) {
    const innerCodecs: unknown[] = [{ name: "bytes", configuration: { endian: "little" } }];
    if (compressor) innerCodecs.push({ name: compressor.id, configuration: { level: compressor.level } });
    await store.set(base + "/zarr.json", j({
      zarr_format: 3, node_type: "array", shape: sh, data_type: V3_DTYPE[dtype] ?? "uint8",
      chunk_grid: { name: "regular", configuration: { chunk_shape: shardShape } },   // outer grid = shards
      chunk_key_encoding: { name: "default", configuration: { separator: "/" } },
      codecs: [{ name: "sharding_indexed", configuration: {
        chunk_shape: chunks,                                    // inner chunk shape
        codecs: innerCodecs,
        index_codecs: [{ name: "bytes", configuration: { endian: "little" } }, { name: "crc32c" }],
        index_location: "end",
      } }],
      fill_value: 0,
    }));
    const perShard = shardShape[0] / chunks[0];                 // inner chunks per shard (integer)
    const nShards = Math.ceil(nChunks / perShard);
    const stail = sh.slice(1).map(() => "0").join("/");         // outer key inner dims (only axis 0 sharded)
    for (let s = 0; s < nShards; s++) {
      const entries: (Uint8Array | null)[] = [];
      for (let slot = 0; slot < perShard; slot++) {
        const ci = s * perShard + slot;
        entries.push(ci < nChunks ? enc(chunkBlock(ci)) : null);   // slots beyond the array are fill (missing)
      }
      const shardBytes = codec!.packShard(entries, /*indexAtEnd=*/true);
      if (shardBytes.length) await store.set(base + "/c/" + (stail ? s + "/" + stail : String(s)), shardBytes);
    }
    return;
  }

  // UNSHARDED path (v2 or v3): one object per chunk. Chunk BYTES are the same as v2 for both formats
  // (LE contiguous, optionally compressed) — only the metadata + chunk-key scheme differ.
  if (fmt === "v3") {
    const codecs: unknown[] = [{ name: "bytes", configuration: { endian: "little" } }];
    if (compressor) codecs.push({ name: compressor.id, configuration: { level: compressor.level } });
    await store.set(base + "/zarr.json", j({
      zarr_format: 3, node_type: "array", shape: sh, data_type: V3_DTYPE[dtype] ?? "uint8",
      chunk_grid: { name: "regular", configuration: { chunk_shape: chunks } },
      chunk_key_encoding: { name: "default", configuration: { separator: "/" } },
      codecs, fill_value: 0,
    }));
  } else {
    await store.set(base + "/.zarray", j({
      zarr_format: 2, shape: sh, chunks, dtype,
      compressor: compressor ? { id: compressor.id, level: compressor.level } : null,
      fill_value: 0, order: "C", filters: null,
    }));
  }
  const sep = fmt === "v3" ? "/" : ".";                         // v3 default key encoding uses "/" + a "c" prefix
  const tail = sh.slice(1).map(() => "0").join(sep);            // only axis 0 is chunked -> inner dims are 0
  for (let ci = 0; ci < nChunks; ci++) {
    const payload = enc(chunkBlock(ci));
    const leaf = tail ? ci + sep + tail : String(ci);
    await store.set(base + "/" + (fmt === "v3" ? "c/" + leaf : leaf), payload);
  }
}

// A group node: v2 writes a bare `.zgroup` (+ `.zattrs` carrying the L* attrs when `lstar` is given);
// v3 writes a single `zarr.json` (node_type group, L* attrs inline). `lstar === undefined` is a bare
// container group (axes/, fields/, models/, passthrough/) with no L* attributes.
async function writeGroup(store: LstarWritableStore, base: string, lstar: unknown, fmt: Fmt = "v3"): Promise<void> {
  const p = base ? base + "/" : "";
  if (fmt === "v3") {
    await store.set(p + "zarr.json", j({ zarr_format: 3, node_type: "group",
      attributes: lstar !== undefined ? { lstar } : {} }));
  } else {
    await store.set(p + ".zgroup", j({ zarr_format: 2 }));
    if (lstar !== undefined) await store.set(p + ".zattrs", j({ lstar }));
  }
}

// A utf8 string array -> a concatenated `|u1` byte array + an `<i8` offsets array (length n+1).
async function writeStrings(store: LstarWritableStore, base: string, valuesKey: string, offsetsKey: string,
                            strings: (string | number)[], opts?: WriteOptions, fmt: Fmt = "v3"): Promise<void> {
  const parts = strings.map((s) => ENC.encode(String(s)));
  const total = parts.reduce((n, p) => n + p.length, 0);
  // offsets are `<i8` (int64) -- the cross-language contract (Python/C++ readers expect int64); the JS
  // reader reads them width-agnostically. A length-(n+1) prefix-sum into the concatenated UTF-8 bytes.
  const bytes = new Uint8Array(total), offs = new BigInt64Array(strings.length + 1);
  let o = 0;
  for (let i = 0; i < parts.length; i++) { bytes.set(parts[i], o); o += parts[i].length; offs[i + 1] = BigInt(o); }
  await writeArray(store, base + "/" + valuesKey, bytes, [total], opts, fmt);
  await writeArray(store, base + "/" + offsetsKey, offs, [strings.length + 1], opts, fmt);
}

async function writeAxis(store: LstarWritableStore, name: string, ax: AxisSpec, opts?: WriteOptions, fmt: Fmt = "v3"): Promise<void> {
  await writeGroup(store, "axes/" + name, { kind: "axis", origin: ax.origin ?? "observed", role: ax.role ?? "", induced_by: ax.inducedBy ?? null, provenance: {} }, fmt);
  await writeStrings(store, "axes/" + name, "labels", "labels_offsets", ax.labels, opts, fmt);
}

async function writeField(store: LstarWritableStore, name: string, f: FieldSpec, opts?: WriteOptions, fmt: Fmt = "v3"): Promise<void> {
  const enc = f.encoding ?? "dense";
  const base = "fields/" + name;
  const partial = f.index != null;
  // Per-field write options override the call-level default. This is what makes the asymmetric layout
  // possible: a gene-major copy stays raw single-chunk (random per-gene byte-range access), while a
  // cell-major copy of the same counts is gzip-chunked (sequential bulk reads) -- set `write` on each.
  const fopts = f.write ?? opts;
  await writeGroup(store, base, {
    kind: "field", role: f.role ?? "", span: f.span ?? [], encoding: enc,
    state: f.state ?? null, subtype: f.subtype ?? null, shape: f.shape ?? null,
    coverage: partial ? "partial" : "full", index_axis: partial ? (f.indexAxis ?? null) : undefined,
    nullable: f.mask != null ? true : undefined,
    ordered: enc === "categorical" ? (f.ordered ?? false) : undefined,
    provenance: f.provenance ?? {},
  }, fmt);
  if (enc === "categorical") {
    await writeArray(store, base + "/codes", f.codes!, [f.codes!.length], fopts, fmt);
    await writeStrings(store, base, "categories", "categories_offsets", f.categories!, fopts, fmt);
  } else if (enc === "csc" || enc === "csr") {
    await writeArray(store, base + "/data", f.data!, [f.data!.length], fopts, fmt);
    await writeArray(store, base + "/indices", f.indices!, [f.indices!.length], fopts, fmt);
    await writeArray(store, base + "/indptr", f.indptr!, [f.indptr!.length], fopts, fmt);
  } else if (enc === "utf8") {
    await writeStrings(store, base, "values", "values_offsets", f.values!, fopts, fmt);
  } else {
    await writeArray(store, base + "/values", f.data!, f.shape ?? [f.data!.length], fopts, fmt);
  }
  if (f.mask != null) await writeArray(store, base + "/mask", f.mask, [f.mask.length], fopts, fmt);
  if (f.index != null) await writeArray(store, base + "/index", f.index, [f.index.length], fopts, fmt);
}

// The `passthrough/<ns>` lossless passthrough: the opaque `tree` (stored as a string) + the array leaves.
async function writeAux(store: LstarWritableStore, ns: string, aux: AuxSpec, opts?: WriteOptions, fmt: Fmt = "v3"): Promise<void> {
  const base = "passthrough/" + ns;
  await writeGroup(store, base, { kind: "passthrough", tree: JSON.stringify(aux.tree), arrays: aux.arrays.map((a) => ({ id: a.id, kind: a.kind })) }, fmt);
  for (const a of aux.arrays) {
    if (a.kind === "utf8") await writeStrings(store, base, a.id, a.id + "_offsets", a.values!, opts, fmt);
    else await writeArray(store, base + "/" + a.id, a.data!, [a.data!.length], opts, fmt);
  }
}

// Record a node's metadata as it's written, keyed for the consolidated map: v2 keeps the store key
// (".../.zarray" etc.); v3 keys by NODE PATH (strip the trailing "zarr.json"), and the ROOT ("") is
// excluded — it is the carrier of the inline consolidated map, not a member of it.
function recordMeta(meta: Record<string, unknown>, key: string, value: Uint8Array, fmt: Fmt): void {
  if (fmt === "v3") {
    if (key.endsWith("zarr.json")) {
      const node = key.slice(0, -"zarr.json".length).replace(/\/$/, "");
      if (node) meta[node] = JSON.parse(TD.decode(value));
    }
  } else if (key.endsWith(".zarray") || key.endsWith(".zattrs") || key.endsWith(".zgroup")) {
    meta[key] = JSON.parse(TD.decode(value));
  }
}

/** Write a complete L* store (root manifest + axes + fields + aux + consolidated metadata). `format`
 * is the on-disk Zarr format: "v3" (default; per-node zarr.json + a root zarr.json carrying the manifest
 * and an inline consolidated_metadata map) or "v2" (legacy .zarray/.zgroup/.zattrs + a consolidated
 * .zmetadata). */
export async function writeStore(store: LstarWritableStore, ds: DatasetSpec, opts?: WriteOptions, format: Fmt = "v3"): Promise<void> {
  // Transparently record every node's metadata as it's written, so we can emit the consolidated map at
  // the end (one-read opens; strict readers). Readers also fall back to per-node metadata.
  const meta: Record<string, unknown> = {};
  const rec: LstarWritableStore = {
    get: (k) => store.get(k),
    set: async (k, v) => { recordMeta(meta, k, v, format); return store.set(k, v); },
    delete: store.delete ? (k) => store.delete!(k) : undefined,
  };
  const manifest = {
    kind: ds.kind ?? "sample", spec_version: ds.specVersion ?? "0.1",
    profiles: ds.profiles ?? [], dropped: ds.dropped ?? [],
    axes: Object.keys(ds.axes), fields: Object.keys(ds.fields), passthrough: Object.keys(ds.aux ?? {}),
  };
  await writeGroup(rec, "", manifest, format);
  // the `axes/`/`fields/`/`models/`/`passthrough/` container groups must carry their own group marker, or
  // strict readers (Python `zarr`) won't navigate into them (zarrita is lenient and resolves by path).
  // `models/` is always present (empty) so a JS-written top-level tree matches the Python/C++/R one.
  await writeGroup(rec, "axes", undefined, format);
  await writeGroup(rec, "fields", undefined, format);
  await writeGroup(rec, "models", undefined, format);
  if (ds.aux && Object.keys(ds.aux).length) await writeGroup(rec, "passthrough", undefined, format);
  for (const [name, ax] of Object.entries(ds.axes)) await writeAxis(rec, name, ax, opts, format);
  for (const [name, f] of Object.entries(ds.fields)) await writeField(rec, name, f, opts, format);
  for (const [ns, a] of Object.entries(ds.aux ?? {})) await writeAux(rec, ns, a, opts, format);
  if (format === "v3") {
    // v3 consolidation is INLINE in the root zarr.json (the map excludes the root itself); overwrite the
    // root node written above with the same manifest plus the consolidated_metadata map.
    await store.set("zarr.json", j({ zarr_format: 3, node_type: "group", attributes: { lstar: manifest },
      consolidated_metadata: { kind: "inline", must_understand: false, metadata: sortedMeta(meta) } }));
  } else {
    await store.set(".zmetadata", j({ zarr_consolidated_format: 1, metadata: sortedMeta(meta) }));
  }
}

/** Append derived axes/fields to an existing store and update the root manifest (+ profiles). AUTO-DETECTS
 * the base store's on-disk format — v3 if a root `zarr.json` group node is present, else v2 — and writes
 * matching nodes + refreshes that format's consolidated metadata, so appending viewer navigators never
 * mixes formats (a v2 base gets v2 nodes, a v3 base gets v3 nodes). */
export async function addToStore(store: LstarWritableStore, add: { axes?: Record<string, AxisSpec>; fields?: Record<string, FieldSpec>; profiles?: string[] }, opts?: WriteOptions): Promise<void> {
  // detect the base format from its root node.
  const rootV3 = await store.get("zarr.json");
  let fmt: Fmt, root: any, m: any;
  if (rootV3) {
    root = JSON.parse(TD.decode(rootV3));
    if (root.node_type !== "group") throw new Error("addToStore: root zarr.json is not a group (not an L* store)");
    fmt = "v3"; m = (root.attributes ??= {}).lstar;
    if (!m) throw new Error("addToStore: no L* manifest in root zarr.json (not an L* store)");
  } else {
    const raw = await store.get(".zattrs");
    if (!raw) throw new Error("addToStore: no root zarr.json / .zattrs (not an L* store)");
    fmt = "v2"; root = JSON.parse(TD.decode(raw)); m = root.lstar;
  }
  // Record every node we (re)write so we can refresh the store's consolidated metadata (instead of
  // dropping it) — keeping the one-request consolidated open valid after an extend.
  const newMeta: Record<string, unknown> = {};
  const rec: LstarWritableStore = {
    get: (k) => store.get(k),
    set: async (k, v) => { recordMeta(newMeta, k, v, fmt); return store.set(k, v); },
    delete: store.delete ? (k) => store.delete!(k) : undefined,
  };
  for (const [name, ax] of Object.entries(add.axes ?? {})) {
    if (!m.axes.includes(name)) m.axes.push(name);
    await writeAxis(rec, name, ax, opts, fmt);
  }
  for (const [name, f] of Object.entries(add.fields ?? {})) {
    if (!m.fields.includes(name)) m.fields.push(name);
    await writeField(rec, name, f, opts, fmt);
  }
  for (const p of add.profiles ?? []) if (!(m.profiles ?? (m.profiles = [])).includes(p)) m.profiles.push(p);
  if (fmt === "v3") {
    // v3: manifest + inline consolidated map both live in the root zarr.json. `m` was mutated in place,
    // so `root` already carries the updated manifest; merge the new nodes into the consolidated map.
    const cm = root.consolidated_metadata?.metadata ?? {};
    root.consolidated_metadata = { kind: "inline", must_understand: false, metadata: sortedMeta({ ...cm, ...newMeta }) };
    await store.set("zarr.json", j(root));
  } else {
    const rootAttrs = j(root);
    await store.set(".zattrs", rootAttrs);
    newMeta[".zattrs"] = JSON.parse(TD.decode(rootAttrs));
    // merge the new keys (+ updated root manifest) into the base `.zmetadata` if present; else leave it
    // absent — readers fall back to per-object metadata.
    const existing = await store.get(".zmetadata");
    if (existing) {
      const cm = JSON.parse(TD.decode(existing));
      cm.metadata = sortedMeta({ ...(cm.metadata ?? {}), ...newMeta });
      await store.set(".zmetadata", j(cm));
    }
  }
}
