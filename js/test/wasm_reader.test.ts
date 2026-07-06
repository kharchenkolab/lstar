// The libzarr WASM reader (openLstarWasm) must drive the full L* API identically to the zarrita reader
// (openLstar) on a v2 store, AND read a v3 store of the SAME data to the same values -- so the viewer can
// retire the zarrita reimplementation and gain v3 reads. Compares manifest + every encoding (dense/csc/csr
// /utf8/categorical/nullable/partial/aux) + the CSC byte-range hot path. Stores come from v3_gen.py (v2
// seed) + the C++ test_v3 writer (v3 copy); paths are argv[2] (v2) and argv[3] (v3).
import { openLstar, openLstarWasm } from "../core/reader.ts";
import { NodeFSStore } from "../core/node-store.ts";

const [v2dir, v3dir] = process.argv.slice(2);
let fail = 0;
const check = (name: string, ok: boolean) => { if (!ok) { console.log(`  FAIL  ${name}`); fail++; } };

function arrEq(a: ArrayLike<any>, b: ArrayLike<any>): boolean {
  if (a.length !== b.length) return false;
  for (let i = 0; i < a.length; i++) if (Number(a[i]) !== Number(b[i])) return false;
  return true;
}
const strEq = (a: any[], b: any[]) => JSON.stringify(a) === JSON.stringify(b);   // for string arrays (labels/names)

const zar = await openLstar(new NodeFSStore(v2dir));       // zarrita reference (v2)
const w2 = await openLstarWasm(new NodeFSStore(v2dir));    // libzarr (v2)
const w3 = await openLstarWasm(new NodeFSStore(v3dir));    // libzarr (v3), same data

// manifest
for (const [tag, ds] of [["v2", w2], ["v3", w3]] as const) {
  check(`${tag} kind`, ds.kind === zar.kind);
  check(`${tag} axisNames`, strEq(ds.axisNames(), zar.axisNames()));
  check(`${tag} fieldNames`, strEq(ds.fieldNames(), zar.fieldNames()));
  check(`${tag} axisLength(cells)`, ds.axisLength("cells") === zar.axisLength("cells"));
}

// per-encoding reads, each compared to the zarrita reference
for (const [tag, ds] of [["v2", w2], ["v3", w3]] as const) {
  check(`${tag} axisLabels(cells)`, strEq(await ds.axisLabels("cells"), await zar.axisLabels("cells")));
  check(`${tag} axisLabels(genes)`, strEq(await ds.axisLabels("genes"), await zar.axisLabels("genes")));

  const dPca = await ds.fieldDense("pca"), zPca = await zar.fieldDense("pca");
  check(`${tag} dense pca`, arrEq(dPca.data, zPca.data) && JSON.stringify(dPca.shape) === JSON.stringify(zPca.shape));

  const dC = await ds.fieldAsCsc("counts"), zC = await zar.fieldAsCsc("counts");   // csc
  check(`${tag} csc counts`, arrEq(dC.data, zC.data) && arrEq(dC.indices, zC.indices) && arrEq(dC.indptr, zC.indptr));

  const dD = await ds.fieldAsCsc("data"), zD = await zar.fieldAsCsc("data");        // csr -> csc
  check(`${tag} csr data`, arrEq(dD.data, zD.data) && arrEq(dD.indices, zD.indices));

  check(`${tag} utf8 barcode`, strEq(await ds.fieldStrings("barcode"), await zar.fieldStrings("barcode")));

  const dCat = await ds.fieldCategorical("leiden"), zCat = await zar.fieldCategorical("leiden");
  check(`${tag} categorical leiden`, arrEq(dCat.codes, zCat.codes) && JSON.stringify(dCat.categories) === JSON.stringify(zCat.categories) && dCat.ordered === zCat.ordered);

  const dM = await ds.fieldMask("n_counts"), zM = await zar.fieldMask("n_counts");  // nullable
  check(`${tag} nullable mask n_counts`, (dM == null) === (zM == null) && (dM == null || arrEq(dM, zM!)));

  const dI = await ds.fieldIndex("adt"), zI = await zar.fieldIndex("adt");           // partial
  check(`${tag} partial index adt`, (dI == null) === (zI == null) && (dI == null || arrEq(dI, zI)));

  const dCol = await ds.cscColumn("counts", 5), zCol = await zar.cscColumn("counts", 5);  // byte-range hot path
  check(`${tag} cscColumn(counts,5)`, arrEq(dCol.rows, zCol.rows) && arrEq(dCol.vals, zCol.vals));

  const dAux = await ds.aux("anndata.uns"), zAux = await zar.aux("anndata.uns");
  check(`${tag} aux anndata.uns`, JSON.stringify(dAux.params) === JSON.stringify(zAux.params) && arrEq(dAux.scores, zAux.scores));
}

console.log(fail === 0
  ? "  openLstarWasm == openLstar across all encodings (v2) + reads v3 to the same values"
  : `  ${fail} parity failures`);
process.exit(fail === 0 ? 0 : 1);
