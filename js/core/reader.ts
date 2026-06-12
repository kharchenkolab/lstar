// @lstar/core — a lazy reader for L* Zarr stores, over zarrita.js.
//
// Opens a store (HTTP/local/zip), reads the consolidated L* metadata, and exposes axes/fields with
// values fetched only when asked — including a single CSC gene-column (the viewer's hot path). The
// heavy numeric work belongs in the WASM kernels (see view.ts); this module is I/O + assembly.
import * as zarr from "zarrita";

/** Minimal store contract (zarrita-compatible): fetch one object by key, or undefined if absent. */
export interface LstarStore {
  get(key: string): Promise<Uint8Array | undefined>;
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
}

const TD = new TextDecoder();

function decodeStrings(bytes: Uint8Array, offsets: ArrayLike<number | bigint>): string[] {
  const n = offsets.length - 1;
  const out: string[] = new Array(n);
  for (let i = 0; i < n; i++) {
    out[i] = TD.decode(bytes.subarray(Number(offsets[i]), Number(offsets[i + 1])));
  }
  return out;
}

export class LstarDataset {
  store: any;
  root: any;
  kind = "sample";
  specVersion = "0.1";
  profiles: string[] = [];
  dropped: string[] = [];
  axes = new Map<string, AxisMeta>();
  fields = new Map<string, FieldMeta>();

  constructor(store: any) {
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
    const grp = await zarr.open(this.root, { kind: "group" });
    const m = (grp.attrs as any).lstar;
    if (!m) throw new Error("not an L* store (no 'lstar' root attribute)");
    this.kind = m.kind ?? "sample";
    this.specVersion = m.spec_version ?? "0.1";
    this.profiles = m.profiles ?? [];
    this.dropped = m.dropped ?? [];
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
                              shape: lm.shape, ordered: lm.ordered } as any);
    }
    return this;
  }

  axisNames(): string[] { return [...this.axes.keys()]; }
  fieldNames(): string[] { return [...this.fields.keys()]; }

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
   * One CSC column (e.g. a gene's expression across cells), fetched as a slice — the hot path.
   * Returns the nonzero `rows` (cell indices) and `vals`. Reads only indptr[col..col+1] and the
   * corresponding data/indices ranges (a few chunks on a chunked store).
   */
  async cscColumn(name: string, col: number): Promise<{ rows: any; vals: any }> {
    const ip = await this._open("fields/" + name + "/indptr");
    const slice = await zarr.get(ip, [zarr.slice(col, col + 2)]);
    const a = Number(slice.data[0]), b = Number(slice.data[1]);
    if (b <= a) return { rows: new Int32Array(0), vals: new Float64Array(0) };
    const rows = (await this._get("fields/" + name + "/indices", [zarr.slice(a, b)])).data;
    const vals = (await this._get("fields/" + name + "/data", [zarr.slice(a, b)])).data;
    return { rows, vals };
  }
}

/** Open an L* store and read its manifest. `store` is any zarrita-compatible store. */
export async function openLstar(store: any): Promise<LstarDataset> {
  return new LstarDataset(store).init();
}
