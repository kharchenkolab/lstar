// Reverse cross-language leg (the pagoda3-prep scenario): JS `addToStore` appends a derived field +
// factor axis to a store the PYTHON writer produced, then Python reads it back. Usage:
// node --experimental-strip-types writer_extend.ts <python-store-dir>  (operate on a COPY, not the fixture).
import { NodeFSStore } from "../core/node-store.ts";
import { addToStore } from "../core/writer.ts";
import { openLstar } from "../core/reader.ts";

const STORE = process.argv[2];
if (!STORE) { console.error("usage: writer_extend.ts <store-dir>"); process.exit(2); }
const store = new NodeFSStore(STORE);

// read an existing axis to size the derived field correctly (so validate stays clean)
const ds = await openLstar(store);
const genes = await ds.axisLabels("genes");
const od = Float32Array.from({ length: genes.length }, (_, i) => i * 0.1);

await addToStore(store, {
  // a derived navigator axis (plain labelled axis -- not factor-induced, so no inducing field needed)
  axes: { od_groups: { labels: ["lo", "hi"], origin: "derived", role: "feature" } },
  fields: { od_score: { role: "measure", span: ["genes"], encoding: "dense", shape: [genes.length], data: od } },
  profiles: ["viewer@0.1"],
});
console.log("extended (addToStore) ->", STORE, `od_score over ${genes.length} genes`);
