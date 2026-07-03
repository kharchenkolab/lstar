// STORED single-file `.lstar.zarr.zip` codec — read (seek-into-zip) + write (pack). No deflate: zarr
// chunks are already codec-compressed, and only a STORED entry stays byte-range-readable inside the
// archive — the point of a hosted single file. A `ZipStore` reads a chunk by issuing ONE range read
// (an HTTP `Range` into the hosted zip, or a file `pread`) at the entry's offset, no decompression.
// Mirrors the Python/C++ codec; ZIP64-aware. Browser-safe: no node imports here (node file/dir helpers
// live in ./zip-node.ts).
import type { LstarStore } from "./reader.ts";

/** A single-file byte source the ZipStore reads ranges of: a local file (`pread`) or a URL (`Range`). */
export interface ByteSource {
  size(): Promise<number>;
  range(start: number, end: number): Promise<Uint8Array>; // [start, end)
}

const encoder = new TextEncoder();
const decoder = new TextDecoder();

const u16 = (b: Uint8Array, o: number) => b[o] | (b[o + 1] << 8);
const u32 = (b: Uint8Array, o: number) =>
  (b[o] | (b[o + 1] << 8) | (b[o + 2] << 16) | (b[o + 3] << 24)) >>> 0;
const u64 = (b: Uint8Array, o: number) => u32(b, o + 4) * 0x100000000 + u32(b, o); // Number, exact < 2^53

const w16 = (v: number) => Uint8Array.of(v & 0xff, (v >>> 8) & 0xff);
const w32 = (v: number) => Uint8Array.of(v & 0xff, (v >>> 8) & 0xff, (v >>> 16) & 0xff, (v >>> 24) & 0xff);
const w64 = (v: number) => {
  const lo = v >>> 0, hi = Math.floor(v / 0x100000000) >>> 0;
  return Uint8Array.of(lo & 0xff, (lo >>> 8) & 0xff, (lo >>> 16) & 0xff, (lo >>> 24) & 0xff,
                       hi & 0xff, (hi >>> 8) & 0xff, (hi >>> 16) & 0xff, (hi >>> 24) & 0xff);
};

function concat(arrs: Uint8Array[]): Uint8Array {
  let n = 0;
  for (const a of arrs) n += a.length;
  const out = new Uint8Array(n);
  let o = 0;
  for (const a of arrs) { out.set(a, o); o += a.length; }
  return out;
}

const CRC_TABLE = (() => {
  const t = new Uint32Array(256);
  for (let i = 0; i < 256; i++) {
    let c = i;
    for (let k = 0; k < 8; k++) c = (c & 1) ? (0xEDB88320 ^ (c >>> 1)) : (c >>> 1);
    t[i] = c >>> 0;
  }
  return t;
})();

/** Standard CRC-32 (poly 0xEDB88320) over `bytes`, as an unsigned 32-bit number. */
export function crc32(bytes: Uint8Array): number {
  let c = 0xffffffff;
  for (let i = 0; i < bytes.length; i++) c = CRC_TABLE[(c ^ bytes[i]) & 0xff] ^ (c >>> 8);
  return (c ^ 0xffffffff) >>> 0;
}

interface ZEntry { lho: number; size: number; dataOff?: number; }

/** Parse a zip's central directory (ZIP64-aware) into name -> {local-header offset, size}. STORED-only:
 * a DEFLATE-compressed entry throws a clear, actionable error (a hosted single-file store must be STORED
 * so its chunks stay byte-range-readable). */
export async function readZipCentralDir(src: ByteSource, forwhat = "zip"): Promise<Map<string, ZEntry>> {
  const fsize = await src.size();
  const taillen = Math.min(fsize, 65557);
  const tail = await src.range(fsize - taillen, fsize);
  let e = -1;
  for (let i = tail.length - 22; i >= 0; i--) if (u32(tail, i) === 0x06054b50) { e = i; break; }
  if (e < 0) throw new Error(`${forwhat}: not a zip (no end-of-central-directory record)`);
  let nEntries = u16(tail, e + 10), cdSize = u32(tail, e + 12), cdOff = u32(tail, e + 16);
  if (nEntries === 0xffff || cdSize === 0xffffffff || cdOff === 0xffffffff) { // ZIP64
    const loc = e - 20;
    if (loc < 0 || u32(tail, loc) !== 0x07064b50)
      throw new Error(`${forwhat}: ZIP64 end-of-central-directory locator missing`);
    const z64off = u64(tail, loc + 8);
    const z = await src.range(z64off, z64off + 56);
    if (u32(z, 0) !== 0x06064b50) throw new Error(`${forwhat}: bad ZIP64 EOCD record`);
    nEntries = u64(z, 32); cdSize = u64(z, 40); cdOff = u64(z, 48);
  }
  const cd = await src.range(cdOff, cdOff + cdSize);
  const idx = new Map<string, ZEntry>();
  const deflated: string[] = [];
  let p = 0;
  for (let i = 0; i < nEntries; i++) {
    if (p + 46 > cd.length || u32(cd, p) !== 0x02014b50)
      throw new Error(`${forwhat}: corrupt central directory`);
    const method = u16(cd, p + 10);
    let usize = u32(cd, p + 24);
    const csizeRaw = u32(cd, p + 20);
    const nlen = u16(cd, p + 28), elen = u16(cd, p + 30), clen = u16(cd, p + 32);
    let lho = u32(cd, p + 42);
    const name = decoder.decode(cd.subarray(p + 46, p + 46 + nlen));
    let ep = p + 46 + nlen; const eend = ep + elen;                  // ZIP64 extra fills 0xFFFFFFFF fields
    while (ep + 4 <= eend) {
      const hid = u16(cd, ep), hsz = u16(cd, ep + 2); let q = ep + 4;
      if (hid === 0x0001) {
        if (usize === 0xffffffff && q + 8 <= eend) { usize = u64(cd, q); q += 8; }
        if (csizeRaw === 0xffffffff && q + 8 <= eend) { q += 8; }   // compressed (== usize for STORED)
        if (lho === 0xffffffff && q + 8 <= eend) { lho = u64(cd, q); q += 8; }
      }
      ep += 4 + hsz;
    }
    if (method !== 0) deflated.push(name);
    idx.set(name, { lho, size: usize });
    p += 46 + nlen + elen + clen;
  }
  if (deflated.length)
    throw new Error(
      `${forwhat}: this .lstar.zarr.zip is DEFLATE-compressed (${deflated.length} entries, e.g. ` +
      `'${deflated[0]}') — a hosted single-file store must be written STORED so its chunks stay ` +
      "byte-range-readable. Repack it STORED (lstar convert, or `zip -0 -r`).");
  return idx;
}

/** A read-only L* store backed by a single STORED zip: `get`/`getRange` resolve a key to its byte range
 * in the archive and issue ONE range read against the underlying source (the seek-into-zip fast path). */
export class ZipStore implements LstarStore {
  private src: ByteSource;
  private idx: Map<string, ZEntry>;
  private constructor(src: ByteSource, idx: Map<string, ZEntry>) { this.src = src; this.idx = idx; }

  static async open(src: ByteSource, forwhat = "zip"): Promise<ZipStore> {
    return new ZipStore(src, await readZipCentralDir(src, forwhat));
  }

  /** Resolve an entry's data offset lazily (the LOCAL header's name/extra lengths can differ from the
   * central record's), caching it — so open() costs only the central-directory reads, not one per entry. */
  private async dataOffset(e: ZEntry): Promise<number> {
    if (e.dataOff !== undefined) return e.dataOff;
    const h = await this.src.range(e.lho, e.lho + 30);
    e.dataOff = e.lho + 30 + u16(h, 26) + u16(h, 28);
    return e.dataOff;
  }

  async get(key: string): Promise<Uint8Array | undefined> {
    const e = this.idx.get(key);
    if (!e) return undefined;
    const off = await this.dataOffset(e);
    return this.src.range(off, off + e.size);
  }

  async getRange(key: string, start: number, end: number): Promise<Uint8Array | undefined> {
    const e = this.idx.get(key);
    if (!e) return undefined;
    const off = await this.dataOffset(e);
    const s = off + Math.max(0, start), en = off + Math.min(e.size, Math.max(0, end));
    return this.src.range(s, Math.max(s, en));
  }
}

/** A `ByteSource` over an HTTP(S) URL: `Range` for sub-ranges (falls back to slicing a 200 body), and
 * `content-length` (or a `content-range` probe) for the size. This is what makes a hosted single-file
 * `.lstar.zarr.zip` range-readable — the reason the archive must be STORED. */
export function httpZipSource(url: string, fetchImpl?: typeof fetch): ByteSource {
  const doFetch = fetchImpl ?? globalThis.fetch.bind(globalThis);
  let cachedSize: number | undefined;
  return {
    async size() {
      if (cachedSize !== undefined) return cachedSize;
      const head = await doFetch(url, { method: "HEAD" });
      const len = head.headers.get("content-length");
      if (head.ok && len) return (cachedSize = parseInt(len, 10));
      const probe = await doFetch(url, { headers: { Range: "bytes=0-0" } });
      const cr = probe.headers.get("content-range"); // "bytes 0-0/12345"
      if (cr && cr.includes("/")) return (cachedSize = parseInt(cr.split("/")[1], 10));
      throw new Error("httpZipSource: cannot determine size of " + url);
    },
    async range(start, end) {
      const res = await doFetch(url, { headers: { Range: `bytes=${start}-${end - 1}` } });
      if (!res.ok && res.status !== 206) throw new Error("httpZipSource: range fetch failed " + res.status);
      const buf = new Uint8Array(await res.arrayBuffer());
      return res.status === 200 ? buf.subarray(start, end) : buf; // 200 = server ignored Range
    },
  };
}

/** Pack `entries` (key -> bytes) into ONE STORED zip (ZIP64 when needed). Metadata (`.z*`) is placed
 * first so a reader hits the manifest early. In-memory — for node, ./zip-node.ts packs a directory. */
export function packStoredZip(entries: Array<[string, Uint8Array]>): Uint8Array {
  const isMeta = (n: string) => n.slice(n.lastIndexOf("/") + 1).startsWith(".z");
  const sorted = [...entries].sort((a, b) => {
    const am = isMeta(a[0]), bm = isMeta(b[0]);
    if (am !== bm) return am ? -1 : 1;
    return a[0] < b[0] ? -1 : a[0] > b[0] ? 1 : 0;
  });
  const pieces: Uint8Array[] = [];
  const cds: { name: Uint8Array; crc: number; size: number; off: number }[] = [];
  let offset = 0;
  for (const [nameStr, data] of sorted) {
    const name = encoder.encode(nameStr);
    const crc = crc32(data), size = data.length, z64 = size >= 0xffffffff;
    const extra = z64 ? concat([w16(0x0001), w16(16), w64(size), w64(size)]) : new Uint8Array(0);
    const lh = concat([w32(0x04034b50), w16(z64 ? 45 : 20), w16(0), w16(0), w16(0), w16(0x21), w32(crc),
      w32(z64 ? 0xffffffff : size), w32(z64 ? 0xffffffff : size),
      w16(name.length), w16(extra.length), name, extra]);
    pieces.push(lh, data);
    cds.push({ name, crc, size, off: offset });
    offset += lh.length + data.length;
  }
  const cdStart = offset;
  const cdPieces: Uint8Array[] = [];
  for (const c of cds) {
    const zsz = c.size >= 0xffffffff, zoff = c.off >= 0xffffffff, z64 = zsz || zoff;
    const body: Uint8Array[] = [];
    if (zsz) { body.push(w64(c.size), w64(c.size)); }
    if (zoff) body.push(w64(c.off));
    const bodyLen = body.reduce((n, a) => n + a.length, 0);
    const extra = z64 ? concat([w16(0x0001), w16(bodyLen), ...body]) : new Uint8Array(0);
    cdPieces.push(concat([w32(0x02014b50), w16(z64 ? 45 : 20), w16(z64 ? 45 : 20), w16(0), w16(0), w16(0),
      w16(0x21), w32(c.crc), w32(zsz ? 0xffffffff : c.size), w32(zsz ? 0xffffffff : c.size),
      w16(c.name.length), w16(extra.length), w16(0), w16(0), w16(0), w32(0),
      w32(zoff ? 0xffffffff : c.off), c.name, extra]));
  }
  const cd = concat(cdPieces);
  const cdSize = cd.length, nrec = cds.length;
  const needZ64 = nrec >= 0xffff || cdStart >= 0xffffffff || cdSize >= 0xffffffff;
  const tailPieces: Uint8Array[] = [];
  if (needZ64) {
    const z64eocd = cdStart + cdSize;
    tailPieces.push(concat([w32(0x06064b50), w64(44), w16(45), w16(45), w32(0), w32(0),
      w64(nrec), w64(nrec), w64(cdSize), w64(cdStart),
      w32(0x07064b50), w32(0), w64(z64eocd), w32(1)]));
  }
  tailPieces.push(concat([w32(0x06054b50), w16(0), w16(0),
    w16(nrec >= 0xffff ? 0xffff : nrec), w16(nrec >= 0xffff ? 0xffff : nrec),
    w32(cdSize >= 0xffffffff ? 0xffffffff : cdSize), w32(cdStart >= 0xffffffff ? 0xffffffff : cdStart),
    w16(0)]));
  return concat([...pieces, cd, ...tailPieces]);
}
