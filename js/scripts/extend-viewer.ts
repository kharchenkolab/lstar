// Optimize an L* store for the viewer (add the viewer@0.1 navigators) IN PLACE — the Node CLI twin of
// `pagoda3.write_viewer` / `lstar convert … --viewer`, no Python/numpy.
//   node --experimental-strip-types extend-viewer.ts <store-dir> [--basis lognorm] [grouping ...]
// Operate on a COPY, not a fixture you want to keep pristine.
import { NodeFSStore } from "../core/node-store.ts";
import { extendForViewer } from "../core/extend.ts";

const dir = process.argv[2];
if (!dir) { console.error("usage: extend-viewer.ts <store-dir> [--basis lognorm] [grouping ...]"); process.exit(2); }
const rest = process.argv.slice(3);
let basis: string | undefined;                          // "lognorm" to prep from a log measure (no raw counts)
const bi = rest.indexOf("--basis");
if (bi >= 0) { basis = rest[bi + 1]; rest.splice(bi, 2); }
const groupings = rest;

const opts: { groupings?: string[]; basis?: string } = {};
if (groupings.length) opts.groupings = groupings;
if (basis) opts.basis = basis;

const store = new NodeFSStore(dir);
const t = Date.now();
await extendForViewer(store, opts);
console.error(`extended ${dir}${groupings.length ? " [" + groupings.join(", ") + "]" : ""}${basis ? " basis=" + basis : ""} in ${Date.now() - t}ms`);
