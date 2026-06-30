// Optimize an L* store for the viewer (add the viewer@0.1 navigators) IN PLACE — the Node CLI twin of
// `pagoda3.write_viewer` / `lstar convert … --viewer`, no Python/numpy.
//   node --experimental-strip-types extend-viewer.ts <store-dir> [grouping ...]
// Operate on a COPY, not a fixture you want to keep pristine.
import { NodeFSStore } from "../core/node-store.ts";
import { extendForViewer } from "../core/extend.ts";

const dir = process.argv[2];
if (!dir) { console.error("usage: extend-viewer.ts <store-dir> [grouping ...]"); process.exit(2); }
const groupings = process.argv.slice(3);

const store = new NodeFSStore(dir);
const t = Date.now();
await extendForViewer(store, groupings.length ? { groupings } : {});
console.error(`extended ${dir}${groupings.length ? " [" + groupings.join(", ") + "]" : ""} in ${Date.now() - t}ms`);
