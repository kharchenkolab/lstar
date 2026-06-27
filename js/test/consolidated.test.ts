// Stage 2: consolidated-metadata open. With a `.zmetadata`, opening the manifest + axes + fields must
// hit the underlying store for ZERO metadata objects (all served from the parsed map) — one read total.
// Without it, the reader falls back to per-object metadata reads (still correct). addToStore refreshes
// `.zmetadata` so the consolidated open stays valid after an extend.
import * as fs from "node:fs";
import * as os from "node:os";
import * as path from "node:path";

import { NodeFSStore } from "../core/node-store.ts";
import { writeStore, addToStore } from "../core/writer.ts";
import { openLstar } from "../core/reader.ts";

let fail = 0;
const check = (name: string, ok: boolean) => { console.log(`  ${ok ? "OK" : "FAIL"}  ${name}`); if (!ok) fail++; };
const META_RE = /(?:\.zgroup|\.zattrs|\.zarray|zarr\.json)$/;

/** Wraps a store and counts how many METADATA objects (and `.zmetadata`) reach the underlying store. */
class CountingStore {
  meta = 0; zmeta = 0; data = 0;
  inner: any;
  constructor(inner: any) { this.inner = inner; }
  async get(key: string, opts?: unknown) {
    const norm = key[0] === "/" ? key.slice(1) : key;
    if (norm.endsWith(".zmetadata")) this.zmeta++;
    else if (META_RE.test(norm)) this.meta++;
    else this.data++;
    return this.inner.get(key, opts);
  }
  set(key: string, v: Uint8Array) { return this.inner.set(key, v); }
  delete(key: string) { return this.inner.delete(key); }
  getRange(key: string, s: number, e: number) { return this.inner.getRange(key, s, e); }
}

const spec = {
  kind: "sample", profiles: ["test@0.1"],
  axes: {
    cells: { labels: ["c0", "c1", "c2"], role: "observation" },
    genes: { labels: ["g0", "g1", "g2", "g3"], role: "feature" },
  },
  fields: {
    counts: { role: "measure", span: ["cells", "genes"], encoding: "csc" as const, state: "raw",
              shape: [3, 4], data: new Int32Array([1, 2, 3, 4, 5]),
              indices: new Int32Array([0, 2, 1, 0, 2]), indptr: new Int32Array([0, 2, 3, 4, 5]) },
    leiden: { role: "label", span: ["cells"], encoding: "utf8" as const, values: ["a", "a", "b"] },
  },
};

const dir = fs.mkdtempSync(path.join(os.tmpdir(), "lstar-consol-"));
await writeStore(new NodeFSStore(dir), spec as any);

// ---- consolidated open: zero metadata objects fetched from the underlying store ----
{
  const cs = new CountingStore(new NodeFSStore(dir));
  const ds = await openLstar(cs);
  check("consolidated: 1 .zmetadata read", cs.zmeta === 1);
  check("consolidated: 0 per-object metadata reads", cs.meta === 0);
  check("consolidated: manifest correct", ds.fieldNames().sort().join() === "counts,leiden" &&
        ds.axes.get("genes")!.length === 4);
  // reads still work (data passes through the wrapper)
  check("consolidated: axisLabels", (await ds.axisLabels("genes")).join() === "g0,g1,g2,g3");
  check("consolidated: fieldStrings", (await ds.fieldStrings("leiden")).join() === "a,a,b");
  const col = await ds.cscColumn("counts", 0);
  check("consolidated: cscColumn", Array.from(col.rows).map(Number).join() === "0,2" &&
        Array.from(col.vals).map(Number).join() === "1,2");
  check("consolidated: meta map populated", Object.keys(ds.meta).length > 0);
}

// ---- fallback: no .zmetadata -> per-object reads, still correct ----
{
  const noconsol = path.join(dir, "..", path.basename(dir) + "-noconsol");
  fs.cpSync(dir, noconsol, { recursive: true });
  fs.rmSync(path.join(noconsol, ".zmetadata"));
  const cs = new CountingStore(new NodeFSStore(noconsol));
  const ds = await openLstar(cs);
  check("fallback: .zmetadata absent -> per-object metadata reads happen", cs.meta > 0);
  check("fallback: manifest still correct", ds.fieldNames().sort().join() === "counts,leiden");
  check("fallback: fieldStrings still correct", (await ds.fieldStrings("leiden")).join() === "a,a,b");
  check("fallback: meta map empty", Object.keys(ds.meta).length === 0);
}

// ---- addToStore refreshes .zmetadata: consolidated open stays valid + sees the new field ----
{
  const store = new NodeFSStore(dir);
  await addToStore(store, {
    fields: { od_score: { role: "measure", span: ["genes"], encoding: "dense", shape: [4],
                          data: new Float32Array([0.1, 0.2, 0.3, 0.4]) } },
    profiles: ["viewer@0.1"],
  });
  check("addToStore: .zmetadata still present", fs.existsSync(path.join(dir, ".zmetadata")));
  const cs = new CountingStore(new NodeFSStore(dir));
  const ds = await openLstar(cs);
  check("addToStore: still consolidated (0 per-object metadata reads)", cs.meta === 0);
  check("addToStore: new field visible via consolidated open", ds.fieldNames().includes("od_score"));
  check("addToStore: new field reads correctly",
        [...(await ds.fieldDense("od_score")).data].map(Math.fround).join() ===
        [0.1, 0.2, 0.3, 0.4].map(Math.fround).join());
  check("addToStore: profiles merged", ds.profiles.includes("viewer@0.1") && ds.profiles.includes("test@0.1"));
}

console.log(fail === 0 ? "\nconsolidated OK" : `\nconsolidated FAIL: ${fail}`);
process.exit(fail === 0 ? 0 : 1);
