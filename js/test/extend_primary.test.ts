// extend_primary.test.ts — `extendForViewer`'s `primary` option: the grouping the viewer opens on is hoisted
// to the front (the counts_cellmajor locality-reorder key + summarized first) and COMPOSES with auto-detect
// (the other groupings are still prepped) — which `groupings` alone can't express. Mirrors the Python test
// (python/tests/test_viewer.py::test_primary_*). Run: node --experimental-strip-types test/extend_primary.test.ts
import assert from "node:assert";
import * as fs from "node:fs";
import * as os from "node:os";
import * as path from "node:path";

import { NodeFSStore } from "../core/node-store.ts";
import { writeStore } from "../core/writer.ts";
import { openLstar } from "../core/reader.ts";
import { extendForViewer } from "../core/extend.ts";

// A 6-cell x 4-gene store with TWO groupings (leiden + cell_type) so the reorder-key choice is observable.
async function freshStore(dir: string): Promise<NodeFSStore> {
  fs.rmSync(dir, { recursive: true, force: true });
  const store = new NodeFSStore(dir);
  await writeStore(store, {
    kind: "sample", profiles: ["test@0.1"],
    axes: {
      cells: { labels: ["c0", "c1", "c2", "c3", "c4", "c5"], role: "observation" },
      genes: { labels: ["g0", "g1", "g2", "g3"], role: "feature" },
      umap: { labels: ["umap0", "umap1"], origin: "derived", role: "coordinate" },
    },
    fields: {
      counts: { role: "measure", span: ["cells", "genes"], encoding: "csc", state: "raw", shape: [6, 4],
                data: new Int32Array([1, 2, 3, 4, 5, 6, 7, 8]),   // per-gene columns
                indices: new Int32Array([0, 1, 2, 3, 4, 5, 0, 5]),
                indptr: new Int32Array([0, 2, 4, 6, 8]) },
      umap: { role: "embedding", span: ["cells", "umap"], encoding: "dense", shape: [6, 2],
              data: new Float32Array([0, 0, 1, 0, 0, 1, 1, 1, 2, 2, 2, 0]) },
      // detection prefers "leiden" over "cell_type" (viewer_policy preferred_groupings order)
      leiden:    { role: "label", span: ["cells"], encoding: "utf8", values: ["a", "a", "b", "b", "c", "c"] },
      cell_type: { role: "label", span: ["cells"], encoding: "utf8", values: ["T", "B", "T", "B", "NK", "NK"] },
    },
  });
  return store;
}

async function main() {
  const base = fs.mkdtempSync(path.join(os.tmpdir(), "lstar-extend-primary-"));

  // 1) primary="cell_type": both groupings prepped (compose w/ auto-detect), reorder keyed on cell_type.
  {
    const store = await freshStore(path.join(base, "primary"));
    await extendForViewer(store, { primary: "cell_type" });
    const ds = await openLstar(store);
    for (const g of ["cell_type", "leiden"]) {
      assert.ok(ds.fieldNames().includes("stats_" + g + "_sum"), `stats_${g}_sum prepped`);
      assert.ok(ds.fieldNames().includes("markers_" + g + "_lfc"), `markers_${g}_lfc prepped`);
    }
    assert.strictEqual((ds.field("counts_cellmajor_order") as any).provenance?.group, "cell_type",
      "reorder keyed on the primary grouping");
  }

  // 2) default (no primary): reorder keys on the auto-detected first grouping — leiden is preferred over
  //    cell_type, so the default primary is leiden. This is exactly what `primary=` lets the viewer override.
  {
    const store = await freshStore(path.join(base, "default"));
    await extendForViewer(store, {});
    const ds = await openLstar(store);
    assert.strictEqual((ds.field("counts_cellmajor_order") as any).provenance?.group, "leiden",
      "default reorder keys on the first detected grouping (leiden)");
  }

  // 3) an unknown primary is a clear error (not a silent no-op); and a field that isn't a cell grouping
  //    (umap: a 2-D embedding over cells) is rejected too — a clear error, not a cryptic reorder crash.
  {
    const store = await freshStore(path.join(base, "bad"));
    await assert.rejects(() => extendForViewer(store, { primary: "not_a_field" }), /primary/,
      "unknown primary rejects");
    const store2 = await freshStore(path.join(base, "bad2"));
    await assert.rejects(() => extendForViewer(store2, { primary: "umap" }), /cell axis/,
      "a non-grouping primary (2-D embedding) rejects");
  }

  // 4) Seurat's active-idents mirror (subtype "active_ident") is NOT a viewer grouping — even though it's a
  //    valid 1-D cell label duplicating the clustering (read_seurat's `ident` field). Matches Python/R.
  {
    const dir = path.join(base, "ident");
    fs.rmSync(dir, { recursive: true, force: true });
    const store = new NodeFSStore(dir);
    await writeStore(store, {
      kind: "sample", profiles: ["test@0.1"],
      axes: {
        cells: { labels: ["c0", "c1", "c2", "c3", "c4", "c5"], role: "observation" },
        genes: { labels: ["g0", "g1", "g2", "g3"], role: "feature" },
        umap: { labels: ["umap0", "umap1"], origin: "derived", role: "coordinate" },
      },
      fields: {
        counts: { role: "measure", span: ["cells", "genes"], encoding: "csc", state: "raw", shape: [6, 4],
                  data: new Int32Array([1, 2, 3, 4, 5, 6, 7, 8]), indices: new Int32Array([0, 1, 2, 3, 4, 5, 0, 5]), indptr: new Int32Array([0, 2, 4, 6, 8]) },
        umap: { role: "embedding", span: ["cells", "umap"], encoding: "dense", shape: [6, 2], data: new Float32Array([0, 0, 1, 0, 0, 1, 1, 1, 2, 2, 2, 0]) },
        leiden: { role: "label", span: ["cells"], encoding: "utf8", values: ["a", "a", "b", "b", "c", "c"] },
        ident: { role: "label", span: ["cells"], encoding: "utf8", subtype: "active_ident", values: ["a", "a", "b", "b", "c", "c"] },
      },
    });
    await extendForViewer(store, {});
    const ds = await openLstar(store);
    assert.ok(ds.fieldNames().includes("stats_leiden_sum"), "leiden is a grouping");
    assert.ok(!ds.fieldNames().includes("stats_ident_sum"), "the active_ident mirror is NOT a viewer grouping (matches Py/R)");
  }

  fs.rmSync(base, { recursive: true, force: true });
  console.log("extend_primary.test: OK");
}

main().catch((e) => { console.error(e); process.exit(1); });
