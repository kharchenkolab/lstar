// Parity gate for the WASM/libzarr reader: read every array of an L* store through the WASM core (a
// sync-fs callback stands in for the browser's prefetch cache) and check dtype + shape + raw bytes against
// the zarr-python reference dump. Also checks the root manifest. Runs for both a v2 and a v3 store, so the
// one reader is proven to read exactly what Python does in either format -- the point of retiring zarrita.
import createLstarIO from "../dist/lstar_io.mjs";
import fs from "node:fs";

const [storeDir, dumpPath] = process.argv.slice(2);
const dump = JSON.parse(fs.readFileSync(dumpPath, "utf8"));
const M = await createLstarIO();
const rd = new M.Reader((key) => {
  const p = `${storeDir}/${key}`;
  return fs.existsSync(p) ? new Uint8Array(fs.readFileSync(p)) : null;
});

let fail = 0;
const check = (name, ok) => { if (!ok) { console.log(`  FAIL  ${name}`); fail++; } };

// manifest (root attrs) matches
const man = JSON.parse(rd.rootAttrs()).lstar;
check("manifest fields", JSON.stringify(man.fields) === JSON.stringify(dump.manifest.fields));
check("manifest axes", JSON.stringify(man.axes) === JSON.stringify(dump.manifest.axes));

// every array: dtype + shape + raw bytes identical to zarr-python
let n = 0;
for (const [path, exp] of Object.entries(dump.arrays)) {
  const a = rd.array(path);
  const shapeOk = JSON.stringify(a.shape) === JSON.stringify(exp.shape);
  const dtypeOk = a.dtype === exp.dtype;
  const b64 = Buffer.from(a.bytes).toString("base64");
  const bytesOk = b64 === exp.b64;
  check(`${path} dtype`, dtypeOk);
  check(`${path} shape`, shapeOk);
  check(`${path} bytes`, bytesOk);
  n++;
}
console.log(fail === 0
  ? `  WASM/libzarr reader == zarr-python across ${n} arrays + manifest`
  : `  ${fail} mismatches over ${n} arrays`);
process.exit(fail === 0 ? 0 : 1);
