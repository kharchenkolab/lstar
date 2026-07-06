// Stage 1: the store byte-range contract. NodeFSStore.getRange does a positional file read; HttpStore
// .getRange issues a `Range` request and falls back to slicing the full body when a server ignores it.
// getSuffix reads the LAST n bytes (a `bytes=-n` suffix range) — the sharded byte-range path uses it to
// fetch a shard's trailing index; same positional-read / 200-fallback contract as getRange.
import * as fs from "node:fs";
import * as os from "node:os";
import * as path from "node:path";
import * as http from "node:http";

import { NodeFSStore } from "../core/node-store.ts";
import { HttpStore } from "../core/http-store.ts";

let fail = 0;
const check = (name: string, ok: boolean) => { console.log(`  ${ok ? "OK" : "FAIL"}  ${name}`); if (!ok) fail++; };
const eqBytes = (a?: Uint8Array, b?: Uint8Array) =>
  !!a && !!b && a.length === b.length && a.every((x, i) => x === b[i]);

// A deterministic 1 KiB payload.
const payload = new Uint8Array(1024);
for (let i = 0; i < payload.length; i++) payload[i] = (i * 7 + 3) & 0xff;

// ---- NodeFSStore.getRange ----
const dir = fs.mkdtempSync(path.join(os.tmpdir(), "lstar-range-"));
fs.writeFileSync(path.join(dir, "blob"), payload);
const nfs = new NodeFSStore(dir);

check("nfs get whole", eqBytes(await nfs.get("blob"), payload));
check("nfs getRange slice", eqBytes(await nfs.getRange("blob", 100, 200), payload.subarray(100, 200)));
check("nfs getRange to-EOF", eqBytes(await nfs.getRange("blob", 1000, 1024), payload.subarray(1000, 1024)));
check("nfs getRange empty", eqBytes(await nfs.getRange("blob", 50, 50), new Uint8Array(0)));
check("nfs getRange missing -> undefined", (await nfs.getRange("nope", 0, 4)) === undefined);
check("nfs getSuffix last 16", eqBytes(await nfs.getSuffix("blob", 16), payload.subarray(1024 - 16)));
check("nfs getSuffix n > size -> whole", eqBytes(await nfs.getSuffix("blob", 5000), payload));
check("nfs getSuffix 0 -> empty", eqBytes(await nfs.getSuffix("blob", 0), new Uint8Array(0)));
check("nfs getSuffix missing -> undefined", (await nfs.getSuffix("nope", 8)) === undefined);

// ---- HttpStore.getRange against a Range-honoring server (206) ----
const rangeServer = http.createServer((req, res) => {
  if (req.url !== "/blob") { res.statusCode = 404; res.end(); return; }
  const range = req.headers["range"];
  const suffix = range && /bytes=-(\d+)/.exec(String(range));
  if (suffix) {
    const n = Math.min(Number(suffix[1]), payload.length);
    const s = payload.length - n;
    res.statusCode = 206;
    res.setHeader("Content-Range", `bytes ${s}-${payload.length - 1}/${payload.length}`);
    res.end(Buffer.from(payload.subarray(s)));
  } else if (range) {
    const m = /bytes=(\d+)-(\d+)/.exec(String(range))!;
    const s = Number(m[1]), e = Number(m[2]); // inclusive
    res.statusCode = 206;
    res.setHeader("Content-Range", `bytes ${s}-${e}/${payload.length}`);
    res.end(Buffer.from(payload.subarray(s, e + 1)));
  } else {
    res.statusCode = 200;
    res.end(Buffer.from(payload));
  }
});
await new Promise<void>((r) => rangeServer.listen(0, "127.0.0.1", r));
const rport = (rangeServer.address() as any).port;
const hs = new HttpStore(`http://127.0.0.1:${rport}/`);

check("http get whole", eqBytes(await hs.get("blob"), payload));
check("http get missing -> undefined", (await hs.get("nope")) === undefined);
check("http getRange (206) slice", eqBytes(await hs.getRange("blob", 100, 200), payload.subarray(100, 200)));
check("http getRange (206) to-EOF", eqBytes(await hs.getRange("blob", 1000, 1024), payload.subarray(1000, 1024)));
check("http getRange empty", eqBytes(await hs.getRange("blob", 7, 7), new Uint8Array(0)));
check("http getSuffix (206) last 16", eqBytes(await hs.getSuffix("blob", 16), payload.subarray(1024 - 16)));
check("http getSuffix (206) missing -> undefined", (await hs.getSuffix("nope", 8)) === undefined);
rangeServer.close();

// ---- HttpStore.getRange against a server that IGNORES Range (always 200 full body) ----
const fullServer = http.createServer((req, res) => {
  if (req.url !== "/blob") { res.statusCode = 404; res.end(); return; }
  res.statusCode = 200; // ignore Range entirely
  res.end(Buffer.from(payload));
});
await new Promise<void>((r) => fullServer.listen(0, "127.0.0.1", r));
const fport = (fullServer.address() as any).port;
const hs2 = new HttpStore(`http://127.0.0.1:${fport}`); // note: no trailing slash -> normalized

check("http getRange (200 fallback) slices correctly", eqBytes(await hs2.getRange("blob", 100, 200), payload.subarray(100, 200)));
check("http getSuffix (200 fallback) takes the tail", eqBytes(await hs2.getSuffix("blob", 16), payload.subarray(1024 - 16)));
fullServer.close();

console.log(fail === 0 ? "\nrange OK" : `\nrange FAIL: ${fail}`);
process.exit(fail === 0 ? 0 : 1);
