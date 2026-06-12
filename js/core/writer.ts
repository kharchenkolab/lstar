// @lstar/core — write L* Zarr stores (the write side of reader.ts). Emits Zarr v2 (single chunk
// per array, uncompressed — readable by zarrita and by the R/Python lstar readers). Use `writeStore`
// for a full dataset, or `addToStore` to append derived fields/axes to an existing store and update
// the root manifest (what the pagoda3 prep does — extend a base store with viewer navigators).
//
// Chunking/compression are intentionally simple here (single uncompressed chunk); the base store's
// own fields keep their original chunking when `addToStore` leaves them untouched.

/** Write side of the store contract: put one object by key (mkdir-p handled by the store). */
export interface LstarWritableStore {
  get(key: string): Promise<Uint8Array | undefined>;
  set(key: string, value: Uint8Array): Promise<void>;
}

export interface AxisSpec { labels: (string | number)[]; origin?: string; role?: string; }
export interface FieldSpec {
  role?: string; span: string[]; encoding: "dense" | "csc" | "csr" | "utf8";
  state?: string | null; subtype?: string | null; shape?: number[] | null;
  // dense: `data` (flat C-order) + `shape`; utf8: `values` (strings); csc/csr: data+indices+indptr + `shape`.
  data?: ArrayLike<number>; values?: (string | number)[];
  indices?: ArrayLike<number>; indptr?: ArrayLike<number>;
}
export interface DatasetSpec {
  kind?: string; specVersion?: string; profiles?: string[]; dropped?: string[];
  axes: Record<string, AxisSpec>; fields: Record<string, FieldSpec>;
}

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

// One Zarr v2 array as a single chunk (chunks == shape, no compressor). Last/only chunk is padded
// to the chunk size with fill_value 0 (here chunk == shape, so padding only kicks in for empty dims).
async function writeArray(store: LstarWritableStore, base: string, a: ArrayLike<number>, shape: number[]): Promise<void> {
  const [dtype, isize, raw] = dtypeOf(a);
  const chunks = shape.length ? shape.map((s) => Math.max(s, 1)) : [1];
  const sh = shape.length ? shape : [a.length];
  await store.set(base + "/.zarray", j({
    zarr_format: 2, shape: sh, chunks, dtype, compressor: null, fill_value: 0, order: "C", filters: null,
  }));
  const nbytes = chunks.reduce((p, c) => p * c, 1) * isize;
  let buf = raw;
  if (buf.byteLength < nbytes) { const padded = new Uint8Array(nbytes); padded.set(buf); buf = padded; }
  const key = sh.map(() => "0").join(".");
  await store.set(base + "/" + key, buf);
}

async function writeGroup(store: LstarWritableStore, base: string, lstar: unknown): Promise<void> {
  const p = base ? base + "/" : "";
  await store.set(p + ".zgroup", j({ zarr_format: 2 }));
  await store.set(p + ".zattrs", j({ lstar }));
}

// A utf8 string array -> a concatenated `|u1` byte array + an `<i4` offsets array (length n+1).
async function writeStrings(store: LstarWritableStore, base: string, valuesKey: string, offsetsKey: string,
                            strings: (string | number)[]): Promise<void> {
  const parts = strings.map((s) => ENC.encode(String(s)));
  const total = parts.reduce((n, p) => n + p.length, 0);
  const bytes = new Uint8Array(total), offs = new Int32Array(strings.length + 1);
  let o = 0;
  for (let i = 0; i < parts.length; i++) { bytes.set(parts[i], o); o += parts[i].length; offs[i + 1] = o; }
  await writeArray(store, base + "/" + valuesKey, bytes, [total]);
  await writeArray(store, base + "/" + offsetsKey, offs, [strings.length + 1]);
}

async function writeAxis(store: LstarWritableStore, name: string, ax: AxisSpec): Promise<void> {
  await writeGroup(store, "axes/" + name, { kind: "axis", origin: ax.origin ?? "observed", role: ax.role ?? "", induced_by: null, provenance: {} });
  await writeStrings(store, "axes/" + name, "labels", "labels_offsets", ax.labels);
}

async function writeField(store: LstarWritableStore, name: string, f: FieldSpec): Promise<void> {
  const enc = f.encoding ?? "dense";
  const base = "fields/" + name;
  await writeGroup(store, base, {
    kind: "field", role: f.role ?? "", span: f.span ?? [], encoding: enc,
    state: f.state ?? null, subtype: f.subtype ?? null, shape: f.shape ?? null, coverage: "full", provenance: {},
  });
  if (enc === "csc" || enc === "csr") {
    await writeArray(store, base + "/data", f.data!, [f.data!.length]);
    await writeArray(store, base + "/indices", f.indices!, [f.indices!.length]);
    await writeArray(store, base + "/indptr", f.indptr!, [f.indptr!.length]);
  } else if (enc === "utf8") {
    await writeStrings(store, base, "values", "values_offsets", f.values!);
  } else {
    await writeArray(store, base + "/values", f.data!, f.shape ?? [f.data!.length]);
  }
}

/** Write a complete L* store (root manifest + axes + fields). */
export async function writeStore(store: LstarWritableStore, ds: DatasetSpec): Promise<void> {
  await writeGroup(store, "", {
    kind: ds.kind ?? "sample", spec_version: ds.specVersion ?? "0.1",
    profiles: ds.profiles ?? [], dropped: ds.dropped ?? [],
    axes: Object.keys(ds.axes), fields: Object.keys(ds.fields),
  });
  // the `axes/` and `fields/` container groups must carry their own .zgroup marker, or strict
  // readers (Python `zarr`) won't navigate into them (zarrita is lenient and resolves by path).
  await store.set("axes/.zgroup", j({ zarr_format: 2 }));
  await store.set("fields/.zgroup", j({ zarr_format: 2 }));
  for (const [name, ax] of Object.entries(ds.axes)) await writeAxis(store, name, ax);
  for (const [name, f] of Object.entries(ds.fields)) await writeField(store, name, f);
}

/** Append derived axes/fields to an existing store and update the root manifest (+ profiles). */
export async function addToStore(store: LstarWritableStore, add: { axes?: Record<string, AxisSpec>; fields?: Record<string, FieldSpec>; profiles?: string[] }): Promise<void> {
  const raw = await store.get(".zattrs");
  if (!raw) throw new Error("addToStore: no root .zattrs (not an L* store)");
  const root = JSON.parse(new TextDecoder().decode(raw));
  const m = root.lstar;
  for (const [name, ax] of Object.entries(add.axes ?? {})) {
    if (!m.axes.includes(name)) m.axes.push(name);
    await writeAxis(store, name, ax);
  }
  for (const [name, f] of Object.entries(add.fields ?? {})) {
    if (!m.fields.includes(name)) m.fields.push(name);
    await writeField(store, name, f);
  }
  for (const p of add.profiles ?? []) if (!(m.profiles ?? (m.profiles = [])).includes(p)) m.profiles.push(p);
  await store.set(".zattrs", j(root));
}
