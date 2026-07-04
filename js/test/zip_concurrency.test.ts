// Store-backend throughput parity: a single-file zip must read a first-screen of fields with the same
// concurrency and ~the same number of round-trips as the directory store — NOT an extra serial
// local-header hop (2× requests) per read. Regression guard for the ZipStore read amplification bug.
//
// Builds one multi-field store (embeddings + labels) as a directory and as a zip, fires ALL field-value
// reads concurrently through `openLstar` with an instrumented IO layer, and asserts: (a) the zip issues
// no more than ~1.5× the directory's reads (a per-read local-header hop would double them), and (b) the
// zip preserves read concurrency (max-in-flight ≈ #fields, not 1). Self-contained — no external fixture.
import assert from "node:assert";
import * as fsp from "node:fs/promises";
import * as os from "node:os";
import * as path from "node:path";

import { openLstar } from "../core/reader.ts";
import { NodeFSStore } from "../core/node-store.ts";
import { ZipStore } from "../core/zip.ts";
import { writeStore, type DatasetSpec } from "../core/writer.ts";
import { nodeFileSource, writeStoreZip } from "../core/zip-node.ts";

const D = 20;
const sleep = (ms: number) => new Promise((r) => setTimeout(r, ms));

function instrument() {
  let inflight = 0, maxInflight = 0, reads = 0;
  const wrap = <F extends (...a: any[]) => Promise<any>>(fn: F): F =>
    (async (...args: any[]) => {
      reads++; inflight++; maxInflight = Math.max(maxInflight, inflight);
      try { await sleep(D); return await fn(...args); } finally { inflight--; }
    }) as F;
  return { wrap, stats: () => ({ reads, maxInflight }) };
}

async function readConcurrently(store: any) {
  const ds = await openLstar(store);
  const jobs: Promise<any>[] = [];
  for (const [name, f] of ds.fields as Map<string, any>) {
    if (f.encoding === "dense") jobs.push(ds.fieldDense(name));
    else if (f.encoding === "utf8") jobs.push(ds.fieldStrings(name));
  }
  const t0 = performance.now();
  await Promise.all(jobs);
  return { ms: performance.now() - t0, nFields: jobs.length };
}

// ---- build a multi-field store (6 embeddings + 4 labels) as a directory AND a zip ----
const n = 1000;
const emb = ["umap", "pca", "tsne", "diffmap", "phate", "fa2"];
const labs = ["leiden", "louvain", "celltype", "batch"];
const spec: DatasetSpec = {
  kind: "sample",
  axes: { cells: { labels: Array.from({ length: n }, (_, i) => "c" + i), origin: "observed" },
          emb2: { labels: ["d0", "d1"], origin: "derived" } },
  fields: {},
};
for (const e of emb)
  spec.fields[e] = { role: "embedding", span: ["cells", "emb2"], encoding: "dense", shape: [n, 2],
                     data: Float32Array.from({ length: n * 2 }, (_, i) => (i * 7 % 100) / 10) };
for (const l of labs)
  spec.fields[l] = { role: "label", span: ["cells"], encoding: "utf8",
                     values: Array.from({ length: n }, (_, i) => l + (i % 7)) };

const tmp = await fsp.mkdtemp(path.join(os.tmpdir(), "lstar-zipperf-"));
try {
  const dir = path.join(tmp, "s.lstar.zarr");
  const zip = path.join(tmp, "s.lstar.zarr.zip");
  await writeStore(new NodeFSStore(dir), spec);
  await writeStoreZip(zip, spec);

  const di = instrument();
  const dir0 = new NodeFSStore(dir);
  const dres = await readConcurrently({ get: di.wrap(dir0.get.bind(dir0)), getRange: di.wrap(dir0.getRange!.bind(dir0)) });
  const ds = di.stats();

  const zi = instrument();
  const zsrc = nodeFileSource(zip);
  const zstore = await ZipStore.open({ size: () => zsrc.size(), range: zi.wrap(zsrc.range.bind(zsrc)) }, zip);
  const zres = await readConcurrently(zstore);
  const zs = zi.stats();

  console.log(`  dir: ${dres.nFields} fields | ${ds.reads} reads | max-in-flight ${ds.maxInflight} | ${dres.ms.toFixed(0)}ms`);
  console.log(`  zip: ${zres.nFields} fields | ${zs.reads} reads | max-in-flight ${zs.maxInflight} | ${zres.ms.toFixed(0)}ms`);

  // (a) no read amplification: a per-read local-header hop would ~double the zip's reads
  assert(zs.reads <= Math.ceil(ds.reads * 1.5) + 2,
    `zip issues too many reads (${zs.reads}) vs dir (${ds.reads}) — a per-read local-header hop is back?`);
  // (b) concurrency preserved: reads overlap, not serialized one-round-trip-at-a-time
  assert(zs.maxInflight >= Math.ceil(zres.nFields / 2),
    `zip serialized reads (max-in-flight ${zs.maxInflight}, expected ~${zres.nFields}) — reads not overlapping`);
  console.log("  zip throughput parity OK (no amplification, concurrency preserved)");
} finally {
  await fsp.rm(tmp, { recursive: true, force: true });
}
console.log("zip concurrency/throughput test passed");
