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

/** A chunk compressor (dependency-injected so writer.ts never imports the WASM/zlib module — the caller
 * wires `compress`, e.g. the WASM `gzipCompress` or a pure-JS gzip). `id`/`level` go into `.zarray`'s
 * `compressor`; `compress` is applied per (padded) chunk. Only "gzip" (RFC1952) is decodable by all
 * readers unchanged. */
export interface Compressor { id: "gzip"; level: number; compress: (raw: Uint8Array) => Uint8Array; }
/** Write knobs: `chunkElems` splits arrays into multi-chunk along axis 0 (target elements/chunk);
 * `compressor` gzips each chunk. Both omitted -> a single uncompressed chunk per array (the default,
 * byte-identical to before). */
export interface WriteOptions { chunkElems?: number; compressor?: Compressor | null; }

const ENC = new TextEncoder();
const j = (o: unknown) => ENC.encode(JSON.stringify(o));

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

// One Zarr v2 array, chunked along axis 0 and optionally gzip-compressed. Each chunk is padded to the
// full chunk size with fill_value 0 BEFORE compression (v2 edge chunks are full-size); the reader trims
// via `shape`. With no opts this writes a single uncompressed chunk -- byte-identical to the old writer.
async function writeArray(store: LstarWritableStore, base: string, a: ArrayLike<number>, shape: number[],
                          opts?: WriteOptions): Promise<void> {
  const [dtype, isize, raw] = dtypeOf(a);
  const sh = shape.length ? shape : [a.length];
  const chunks = chunkShapeFor(sh, opts?.chunkElems).map((s) => Math.max(s, 1));
  const compressor = opts?.compressor ?? null;
  await store.set(base + "/.zarray", j({
    zarr_format: 2, shape: sh, chunks, dtype,
    compressor: compressor ? { id: compressor.id, level: compressor.level } : null,
    fill_value: 0, order: "C", filters: null,
  }));
  const inner = chunks.slice(1).reduce((p, c) => p * c, 1);     // elements per axis-0 row (inner dims)
  const chunkBytes = chunks[0] * inner * isize;
  const nRows = sh[0] ?? 0;
  const nChunks = nRows <= 0 ? 1 : Math.ceil(nRows / chunks[0]);
  const tail = sh.slice(1).map(() => "0").join(".");            // only axis 0 is chunked -> ".0..." inner
  for (let ci = 0; ci < nChunks; ci++) {
    const block = new Uint8Array(chunkBytes);                   // zero-padded to full chunk size
    block.set(raw.subarray(ci * chunkBytes, Math.min((ci + 1) * chunkBytes, raw.byteLength)));
    const payload = compressor ? compressor.compress(block) : block;
    await store.set(base + "/" + (tail ? ci + "." + tail : String(ci)), payload);
  }
}

async function writeGroup(store: LstarWritableStore, base: string, lstar: unknown): Promise<void> {
  const p = base ? base + "/" : "";
  await store.set(p + ".zgroup", j({ zarr_format: 2 }));
  await store.set(p + ".zattrs", j({ lstar }));
}

// A utf8 string array -> a concatenated `|u1` byte array + an `<i8` offsets array (length n+1).
async function writeStrings(store: LstarWritableStore, base: string, valuesKey: string, offsetsKey: string,
                            strings: (string | number)[], opts?: WriteOptions): Promise<void> {
  const parts = strings.map((s) => ENC.encode(String(s)));
  const total = parts.reduce((n, p) => n + p.length, 0);
  // offsets are `<i8` (int64) -- the cross-language contract (Python/C++ readers expect int64); the JS
  // reader reads them width-agnostically. A length-(n+1) prefix-sum into the concatenated UTF-8 bytes.
  const bytes = new Uint8Array(total), offs = new BigInt64Array(strings.length + 1);
  let o = 0;
  for (let i = 0; i < parts.length; i++) { bytes.set(parts[i], o); o += parts[i].length; offs[i + 1] = BigInt(o); }
  await writeArray(store, base + "/" + valuesKey, bytes, [total], opts);
  await writeArray(store, base + "/" + offsetsKey, offs, [strings.length + 1], opts);
}

async function writeAxis(store: LstarWritableStore, name: string, ax: AxisSpec, opts?: WriteOptions): Promise<void> {
  await writeGroup(store, "axes/" + name, { kind: "axis", origin: ax.origin ?? "observed", role: ax.role ?? "", induced_by: ax.inducedBy ?? null, provenance: {} });
  await writeStrings(store, "axes/" + name, "labels", "labels_offsets", ax.labels, opts);
}

async function writeField(store: LstarWritableStore, name: string, f: FieldSpec, opts?: WriteOptions): Promise<void> {
  const enc = f.encoding ?? "dense";
  const base = "fields/" + name;
  const partial = f.index != null;
  await writeGroup(store, base, {
    kind: "field", role: f.role ?? "", span: f.span ?? [], encoding: enc,
    state: f.state ?? null, subtype: f.subtype ?? null, shape: f.shape ?? null,
    coverage: partial ? "partial" : "full", index_axis: partial ? (f.indexAxis ?? null) : undefined,
    nullable: f.mask != null ? true : undefined,
    ordered: enc === "categorical" ? (f.ordered ?? false) : undefined,
    provenance: {},
  });
  if (enc === "categorical") {
    await writeArray(store, base + "/codes", f.codes!, [f.codes!.length], opts);
    await writeStrings(store, base, "categories", "categories_offsets", f.categories!, opts);
  } else if (enc === "csc" || enc === "csr") {
    await writeArray(store, base + "/data", f.data!, [f.data!.length], opts);
    await writeArray(store, base + "/indices", f.indices!, [f.indices!.length], opts);
    await writeArray(store, base + "/indptr", f.indptr!, [f.indptr!.length], opts);
  } else if (enc === "utf8") {
    await writeStrings(store, base, "values", "values_offsets", f.values!, opts);
  } else {
    await writeArray(store, base + "/values", f.data!, f.shape ?? [f.data!.length], opts);
  }
  if (f.mask != null) await writeArray(store, base + "/mask", f.mask, [f.mask.length], opts);
  if (f.index != null) await writeArray(store, base + "/index", f.index, [f.index.length], opts);
}

// The `passthrough/<ns>` lossless passthrough: the opaque `tree` (stored as a string) + the array leaves.
async function writeAux(store: LstarWritableStore, ns: string, aux: AuxSpec, opts?: WriteOptions): Promise<void> {
  const base = "passthrough/" + ns;
  await writeGroup(store, base, { kind: "passthrough", tree: JSON.stringify(aux.tree), arrays: aux.arrays.map((a) => ({ id: a.id, kind: a.kind })) });
  for (const a of aux.arrays) {
    if (a.kind === "utf8") await writeStrings(store, base, a.id, a.id + "_offsets", a.values!, opts);
    else await writeArray(store, base + "/" + a.id, a.data!, [a.data!.length], opts);
  }
}

/** Write a complete L* store (root manifest + axes + fields + aux + consolidated `.zmetadata`). */
export async function writeStore(store: LstarWritableStore, ds: DatasetSpec, opts?: WriteOptions): Promise<void> {
  // Transparently record every .zgroup/.zattrs/.zarray as it's written, so we can emit a consolidated
  // `.zmetadata` at the end (one-read opens; strict readers). Readers also fall back to per-object meta.
  const meta: Record<string, unknown> = {};
  const rec: LstarWritableStore = {
    get: (k) => store.get(k),
    set: async (k, v) => {
      if (k.endsWith(".zarray") || k.endsWith(".zattrs") || k.endsWith(".zgroup")) meta[k] = JSON.parse(new TextDecoder().decode(v));
      return store.set(k, v);
    },
    delete: store.delete ? (k) => store.delete!(k) : undefined,
  };
  await writeGroup(rec, "", {
    kind: ds.kind ?? "sample", spec_version: ds.specVersion ?? "0.1",
    profiles: ds.profiles ?? [], dropped: ds.dropped ?? [],
    axes: Object.keys(ds.axes), fields: Object.keys(ds.fields), passthrough: Object.keys(ds.aux ?? {}),
  });
  // the `axes/`/`fields/`/`passthrough/` container groups must carry their own .zgroup marker, or strict
  // readers (Python `zarr`) won't navigate into them (zarrita is lenient and resolves by path).
  await rec.set("axes/.zgroup", j({ zarr_format: 2 }));
  await rec.set("fields/.zgroup", j({ zarr_format: 2 }));
  if (ds.aux && Object.keys(ds.aux).length) await rec.set("passthrough/.zgroup", j({ zarr_format: 2 }));
  for (const [name, ax] of Object.entries(ds.axes)) await writeAxis(rec, name, ax, opts);
  for (const [name, f] of Object.entries(ds.fields)) await writeField(rec, name, f, opts);
  for (const [ns, a] of Object.entries(ds.aux ?? {})) await writeAux(rec, ns, a, opts);
  await store.set(".zmetadata", j({ zarr_consolidated_format: 1, metadata: meta }));
}

/** Append derived axes/fields to an existing store and update the root manifest (+ profiles). */
export async function addToStore(store: LstarWritableStore, add: { axes?: Record<string, AxisSpec>; fields?: Record<string, FieldSpec>; profiles?: string[] }, opts?: WriteOptions): Promise<void> {
  const raw = await store.get(".zattrs");
  if (!raw) throw new Error("addToStore: no root .zattrs (not an L* store)");
  const root = JSON.parse(new TextDecoder().decode(raw));
  const m = root.lstar;
  for (const [name, ax] of Object.entries(add.axes ?? {})) {
    if (!m.axes.includes(name)) m.axes.push(name);
    await writeAxis(store, name, ax, opts);
  }
  for (const [name, f] of Object.entries(add.fields ?? {})) {
    if (!m.fields.includes(name)) m.fields.push(name);
    await writeField(store, name, f, opts);
  }
  for (const p of add.profiles ?? []) if (!(m.profiles ?? (m.profiles = [])).includes(p)) m.profiles.push(p);
  await store.set(".zattrs", j(root));
  // The base store may carry consolidated metadata (.zmetadata) that now lists a stale field set;
  // drop it so readers fall back to the live per-object metadata we just wrote (both zarrita's
  // withConsolidated and Python's _open_root fall back gracefully when it's absent).
  if (store.delete) await store.delete(".zmetadata");
}
