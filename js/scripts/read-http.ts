// Authenticity check for pagoda3's standalone server: open a store over HTTP using the SAME reader the
// browser viewer runs (@lstar/core openLstar + HttpStore byte-range reads) — pointed at pagoda3.serve.
//   node --experimental-strip-types read-http.ts <base-url>
import { openLstar } from "../core/index.ts";
import { HttpStore } from "../core/http-store.ts";

const base = process.argv[2] || "http://127.0.0.1:9305/";
const ds = await openLstar(new HttpStore(base));

const names = ds.fieldNames();
console.error("opened", base, "·", names.length, "fields");
console.error("  has counts_cellmajor:", ds.hasField("counts_cellmajor"),
              "| counts:", ds.hasField("counts"));

// embedding — a dense [cells, 2] field; read it via the real reader (range chunks under the hood)
const embName = names.find((n) => ds.field(n)?.role === "embedding");
if (embName) {
  const e = await ds.fieldDense(embName);
  let finite = 0;
  for (let i = 0; i < Math.min(e.data.length, 2000); i++) if (Number.isFinite(e.data[i])) finite++;
  console.error(`  embedding '${embName}' shape=${e.shape} first=(${e.data[0]?.toFixed?.(3)}, ${e.data[1]?.toFixed?.(3)}) finite/2000=${finite}`);
}

// a gene column from counts via byte-range (the gene-coloring fast path)
const sp = await ds.fieldSparse("counts");
console.error(`  counts ${sp.fmt} shape=${sp.shape} nnz=${sp.data.length}`);

console.error("READ OK — the viewer's reader reads pagoda3.serve over HTTP range");
