// Real-data format-invariance for the libzarr WASM reader: on a real dataset, openLstarWasm must read
// the v2 and v3 copies to IDENTICAL values across the full L* API (dense/sparse/utf8/categorical), field
// for field. Absolute correctness (vs the zarr-python reference) is checked separately, at the array
// level, by io_parity.mjs; this checks the higher-level L* assembly is format-invariant. Zarrita is NOT
// used as the oracle here: it mis-reads boolean (`|b1`) dense fields as NaN (a real bug the libzarr reader
// does not share), so comparing against it would spuriously fail on any dataset carrying a bool field.
// Args: <v2dir> <v3dir>.
import { openLstarWasm } from "../core/reader.ts";
import { NodeFSStore } from "../core/node-store.ts";

const [v2dir, v3dir] = process.argv.slice(2);
let fail = 0;
const check = (n, ok) => { if (!ok) { console.log(`  FAIL  ${n}`); fail++; } };
const numEq = (a, b) => {   // NaN==NaN (real DE tables carry NaN placeholders)
  if (a.length !== b.length) return false;
  for (let i = 0; i < a.length; i++) { const x = Number(a[i]), y = Number(b[i]); if (x !== y && !(Number.isNaN(x) && Number.isNaN(y))) return false; }
  return true;
};
const strEq = (a, b) => JSON.stringify(a) === JSON.stringify(b);

async function fieldVals(ds, name) {                     // generic read dispatched on encoding
  const enc = ds.field(name).encoding;
  if (enc === "csc" || enc === "csr" || enc === "dense") { const c = await ds.fieldAsCsc(name); return { kind: "num", v: c.data }; }
  if (enc === "utf8") return { kind: "str", v: await ds.fieldStrings(name) };
  if (enc === "categorical") { const c = await ds.fieldCategorical(name); return { kind: "cat", codes: c.codes, cats: c.categories }; }
  return { kind: "skip" };
}

const w2 = await openLstarWasm(new NodeFSStore(v2dir));   // libzarr (v2)
const w3 = await openLstarWasm(new NodeFSStore(v3dir));   // libzarr (v3), same data

check("fieldNames", strEq(w3.fieldNames(), w2.fieldNames()));
check("axisNames", strEq(w3.axisNames(), w2.axisNames()));

let nfields = 0;
for (const name of w2.fieldNames()) {
  const [r2, r3] = await Promise.all([fieldVals(w2, name), fieldVals(w3, name)]);
  if (r2.kind === "skip") continue;
  nfields++;
  if (r2.kind === "num") check(`${name}`, numEq(r3.v, r2.v));
  else if (r2.kind === "str") check(`${name}`, strEq(r3.v, r2.v));
  else if (r2.kind === "cat") check(`${name}`, numEq(r3.codes, r2.codes) && strEq(r3.cats, r2.cats));
}
const ax = w2.axisNames()[0];
check(`axisLabels(${ax})`, strEq(await w3.axisLabels(ax), await w2.axisLabels(ax)));

console.log(fail === 0
  ? `openLstarWasm reads v2 and v3 to identical values across ${nfields} fields`
  : `  ${fail} format-invariance failures over ${nfields} fields`);
process.exit(fail === 0 ? 0 : 1);
