// A minimal HTTP(S) store for L* Zarr stores, with byte-range support so the reader can fetch one
// gene/cell as a few KB (a `Range` request) instead of the whole uncompressed chunk. zarrita-compatible:
// `get(key)` resolves a key against a base URL and returns its bytes (or undefined on 404). `getRange`
// issues a `Range` request and gracefully falls back to slicing the full body when the server ignores
// Range (responds 200 instead of 206) -- correct either way, only the fast path saves bandwidth.
import type { LstarStore } from "./reader.ts";

export class HttpStore implements LstarStore {
  base: string;
  private _fetch: typeof fetch;
  constructor(base: string, opts?: { fetch?: typeof fetch }) {
    this.base = base.endsWith("/") ? base : base + "/";          // predictable URL join
    // bind to the global: an unbound `globalThis.fetch` throws "Illegal invocation" in browsers (fetch
    // requires `this === window`). Node tolerates the unbound call, so this only bit under a real browser.
    this._fetch = opts?.fetch ?? globalThis.fetch.bind(globalThis);
  }
  private url(key: string): string {
    return this.base + (key.startsWith("/") ? key.slice(1) : key);
  }
  async get(key: string, _opts?: unknown): Promise<Uint8Array | undefined> {
    const res = await this._fetch(this.url(key));
    if (res.status === 404) return undefined;
    if (!res.ok) throw new Error(`HTTP ${res.status} for ${this.url(key)}`);
    return new Uint8Array(await res.arrayBuffer());
  }
  async getRange(key: string, start: number, end: number): Promise<Uint8Array | undefined> {
    if (end <= start) return new Uint8Array(0);
    // HTTP byte ranges are inclusive on both ends: [start, end) -> bytes=start-(end-1).
    const res = await this._fetch(this.url(key), { headers: { Range: `bytes=${start}-${end - 1}` } });
    if (res.status === 404) return undefined;
    if (!res.ok && res.status !== 206) throw new Error(`HTTP ${res.status} for ${this.url(key)}`);
    const buf = new Uint8Array(await res.arrayBuffer());
    // 206 Partial Content: the body is exactly the requested slice. 200 (server/CDN ignored Range):
    // the body is the whole object -> slice it to honor the [start, end) contract.
    if (res.status === 206) return buf;
    return buf.subarray(start, Math.min(end, buf.length));
  }
  /** Read the LAST `n` bytes of an object (a suffix `Range: bytes=-n`), or undefined on 404. Backs the
   * byte-range fast path on sharded arrays (the v3 shard index lives at the object's end). Falls back to
   * slicing the tail of the full body when a server ignores Range (200) — correct, just no bandwidth win. */
  async getSuffix(key: string, n: number): Promise<Uint8Array | undefined> {
    if (n <= 0) return new Uint8Array(0);
    const res = await this._fetch(this.url(key), { headers: { Range: `bytes=-${n}` } });
    if (res.status === 404) return undefined;
    if (!res.ok && res.status !== 206) throw new Error(`HTTP ${res.status} for ${this.url(key)}`);
    const buf = new Uint8Array(await res.arrayBuffer());
    if (res.status === 206) return buf;                     // exactly the last n bytes
    return buf.subarray(Math.max(0, buf.length - n));       // 200: whole body -> take the tail
  }
}
