// Guard: httpZipSource must issue EVERY fetch with `cache: "no-store"`. All reads target the one archive
// URL, and a browser serializes concurrent same-URL fetches under the default cache mode (a per-cache-entry
// lock) — so a hosted zip's cold open would be strictly serial. `no-store` opts out of that lock. The lock
// is browser-only (node's fetch has no HTTP cache), so the real effect can't be asserted here; instead we
// assert the INTENT — the option is passed on range + size fetches — which is the thing a regression drops.
import assert from "node:assert";

import { httpZipSource } from "../core/zip.ts";

const calls: Array<{ url: any; opts: any }> = [];
const mockFetch = (async (url: any, opts: any) => {
  calls.push({ url, opts });
  const h = new Map<string, string>([["content-length", "1024"]]);
  return {
    ok: true, status: 206,
    headers: { get: (k: string) => h.get(k.toLowerCase()) ?? null },
    arrayBuffer: async () => new Uint8Array(16).buffer,
  } as any;
}) as unknown as typeof fetch;

const src = httpZipSource("https://example.test/store.lstar.zarr.zip", mockFetch);
await src.range(0, 16);
await src.size();

assert(calls.length >= 2, `expected range + size fetches, got ${calls.length}`);
for (const c of calls) {
  assert(c.opts && c.opts.cache === "no-store",
    `httpZipSource must fetch with cache:"no-store" (else the browser same-URL cache lock serializes ` +
    `concurrent zip reads); got cache=${JSON.stringify(c.opts?.cache)} on ${JSON.stringify(c.opts)}`);
}
console.log(`  [js] httpZipSource: all ${calls.length} fetch(es) use cache:"no-store"`);
console.log("httpZipSource no-store guard passed");
