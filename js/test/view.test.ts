// Phase C: the viewer query API must produce the values a viewer renders, matching Python references.
import * as fs from "node:fs";
import * as path from "node:path";
import { fileURLToPath } from "node:url";

import { openLstar } from "../core/reader.ts";
import { NodeFSStore } from "../core/node-store.ts";
import { LstarView, Crossfilter, scalarToRGBA } from "../core/view.ts";

const HERE = path.dirname(fileURLToPath(import.meta.url));
const expected = JSON.parse(fs.readFileSync(path.join(HERE, "data", "expected.json"), "utf8"));

function approx(a: ArrayLike<number>, b: ArrayLike<number>, rtol = 1e-5, atol = 1e-6): boolean {
  if (a.length !== b.length) return false;
  for (let i = 0; i < a.length; i++)
    if (Math.abs(Number(a[i]) - Number(b[i])) > atol + rtol * Math.abs(Number(b[i]))) return false;
  return true;
}
let fail = 0;
const check = (name: string, ok: boolean) => { console.log(`  ${ok ? "OK" : "FAIL"}  ${name}`); if (!ok) fail++; };

const ds = await openLstar(new NodeFSStore(path.join(HERE, "data", "sample.lstar.zarr")));
const view = new LstarView(ds);

// embedding
const emb = await view.embedding("umap");
check("embedding(umap)", emb.n === expected.n_cells && emb.dim === 2 && approx(emb.data, expected.umap));

// metadata: categorical (leiden) + numeric (n_umi)
const md = await view.metadata("leiden");
check("metadata(leiden) categorical", md.kind === "categorical" &&
      (md as any).categories.length === 3 &&
      Array.from((md as any).codes).every((c: number, i: number) =>
        (md as any).categories[c] === expected.leiden[i]));
const num = await view.metadata("n_umi");
check("metadata(n_umi) numeric", num.kind === "numeric" && approx((num as any).values, expected.n_umi));

// gene coloring: scatter the gene's CSC column, log1p; compare to a dense reference built from expected
{
  const gc = expected.gene_col;
  const ref = new Float32Array(expected.n_cells);
  for (let k = 0; k < gc.rows.length; k++) ref[gc.rows[k]] = Math.log1p(gc.vals[k]);
  const genes = await ds.axisLabels("genes");
  const ge = await view.geneExpression(genes[gc.index], { lognorm: true });
  check("geneExpression(g) values", approx(ge.values, ref));
  const rgba = scalarToRGBA(ge.values, ge.max);
  check("scalarToRGBA shape", rgba.length === expected.n_cells * 4);
}

// per-gene stats via WASM, vs Python stream_col_stats
{
  const cs = await view.colStats({ lognorm: true });
  const ref = expected.colstats_lognorm;
  check("colStats mean/var/nnz (WASM)",
        approx(cs.mean, ref.mean, 1e-9) && approx(cs.var, ref.var, 1e-7) && approx(cs.nnz, ref.nnz));
}

// DE: per-group means for leiden A vs B vs Python reference
{
  const de = expected.de_ref;
  const ranked = await view.subsampleDE(de.cellsA, de.cellsB, { lognorm: true });
  const meanA = new Float64Array(expected.n_genes), meanB = new Float64Array(expected.n_genes);
  for (const r of ranked) { meanA[r.gene] = r.meanA; meanB[r.gene] = r.meanB; }
  check("subsampleDE group means", approx(meanA, de.meanA, 1e-9) && approx(meanB, de.meanB, 1e-9));
  check("subsampleDE ranked by |lfc|", Math.abs(ranked[0].lfc) >= Math.abs(ranked[ranked.length - 1].lfc));
}

// crossfilter: leiden == "A" selects exactly the A cells
{
  const m = await view.metadata("leiden");
  const aCode = (m as any).categories.indexOf("A");
  const cf = new Crossfilter(expected.n_cells).categorical((m as any).codes, [aCode]);
  const nA = expected.leiden.filter((x: string) => x === "A").length;
  check("crossfilter categorical", cf.count() === nA && cf.selected().every((i: number) => expected.leiden[i] === "A"));
}

console.log(fail === 0 ? "\nview OK" : `\nview FAIL: ${fail}`);
process.exit(fail === 0 ? 0 : 1);
