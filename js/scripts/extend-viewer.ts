// Optimize an L* store for the viewer (add the viewer@0.1 navigators) — the Node CLI twin of
// `pagoda3.write_viewer` / `lstar convert … --viewer`, no Python/numpy.
//   node --experimental-strip-types extend-viewer.ts <store> [--basis lognorm] [--out DST] [grouping ...]
// <store> and --out DST may each be a directory OR a single-file `.lstar.zarr.zip`. Without --out the
// store is extended IN PLACE (directory only) — operate on a COPY, not a fixture you want pristine.
// With --out, the input is left pristine and the extended store is written to DST (packed STORED when
// DST ends in `.zip`, e.g. `extend-viewer.ts sample.lstar.zarr --out sample.viewer.lstar.zarr.zip`).
import * as fsp from "node:fs/promises";
import * as os from "node:os";
import * as path from "node:path";

import { NodeFSStore } from "../core/node-store.ts";
import { extendForViewer } from "../core/extend.ts";
import { packStoredZipDir, extractZipToDir } from "../core/zip-node.ts";

const input = process.argv[2];
if (!input) {
  console.error("usage: extend-viewer.ts <store> [--basis lognorm] [--primary G] [--out DST] [grouping ...]");
  process.exit(2);
}
const rest = process.argv.slice(3);
let basis: string | undefined;                          // "lognorm" to prep from a log measure (no raw counts)
const bi = rest.indexOf("--basis");
if (bi >= 0) { basis = rest[bi + 1]; rest.splice(bi, 2); }
let out: string | undefined;                            // write here instead of extending in place
const oi = rest.indexOf("--out");
if (oi >= 0) { out = rest[oi + 1]; rest.splice(oi, 2); }
let primary: string | undefined;                        // the grouping the viewer opens on (hoisted to front)
const pi = rest.indexOf("--primary");
if (pi >= 0) { primary = rest[pi + 1]; rest.splice(pi, 2); }
const groupings = rest;

const opts: { groupings?: string[]; basis?: string; primary?: string } = {};
if (groupings.length) opts.groupings = groupings;
if (basis) opts.basis = basis;
if (primary) opts.primary = primary;

const t = Date.now();
const note = `${groupings.length ? " [" + groupings.join(", ") + "]" : ""}${basis ? " basis=" + basis : ""}${primary ? " primary=" + primary : ""}`;

if (!out) {
  if (input.endsWith(".zip"))
    throw new Error("in-place extend needs a directory store; pass --out DST to read a .zip and write elsewhere");
  await extendForViewer(new NodeFSStore(input), opts);   // extend the directory in place
  console.error(`extended ${input}${note} in ${Date.now() - t}ms`);
} else {
  // extend a working copy (input stays pristine), then emit DST as a directory or a STORED zip
  const work = await fsp.mkdtemp(path.join(os.tmpdir(), "lstar-ext-"));
  try {
    if (input.endsWith(".zip")) await extractZipToDir(input, work);
    else await fsp.cp(input, work, { recursive: true });
    await extendForViewer(new NodeFSStore(work), opts);
    if (out.endsWith(".zip")) {
      await packStoredZipDir(work, out);
    } else {
      await fsp.rm(out, { recursive: true, force: true });
      await fsp.cp(work, out, { recursive: true });
    }
  } finally {
    await fsp.rm(work, { recursive: true, force: true });
  }
  console.error(`extended ${input} -> ${out}${note} in ${Date.now() - t}ms`);
}
