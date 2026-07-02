// selectCountsBasis — pick the viewer-prep count basis by STATE, not the literal name "counts".
// Mirrors Python lstar.viewer._select_counts_basis (test_convert_viewer_basis.py).
import { test } from "node:test";
import assert from "node:assert/strict";
import { selectCountsBasis } from "../core/basis.ts";

const M = (span: string[] = ["cells", "genes"], role = "measure", state?: string) => ({ role, span, state } as any);
const ds = (f: Record<string, any>) => ({ fieldNames: () => Object.keys(f), field: (n: string) => f[n] });

test("default: prefers a measure named 'counts' (raw → log1p)", () => {
  const r = selectCountsBasis(ds({ counts: M(["cells", "genes"], "measure", "raw"), X: M(["cells", "genes"], "measure", "scaled") }));
  assert.deepEqual(r, { field: "counts", log1p: true });
});

test("default: no 'counts' name → picks a raw-state measure (e.g. an AnnData .X)", () => {
  const r = selectCountsBasis(ds({ X: M(["cells", "genes"], "measure", "raw") }));
  assert.deepEqual(r, { field: "X", log1p: true });
});

test("no raw basis (scaled .X + lognorm .raw, the pbmc3k tutorial) → clear error listing measures", () => {
  const call = () => selectCountsBasis(ds({ X: M(["cells", "genes"], "measure", "scaled"), rawX: M(["cells", "genes"], "measure", "lognorm") }));
  assert.throws(call, /no raw counts/);
  assert.throws(call, /X\[scaled\]/);          // present measures + state are listed
  assert.throws(call, /basis="lognorm"/);      // and it offers the way forward
});

test("counts= forces a measure (and rejects an unknown one)", () => {
  const r = selectCountsBasis(ds({ X: M(["cells", "genes"], "measure", "scaled") }), { counts: "X" });
  assert.deepEqual(r, { field: "X", log1p: true });     // a forced non-lognorm measure is still log1p'd
  assert.throws(() => selectCountsBasis(ds({ X: M() }), { counts: "nope" }), /counts="nope" is not a measure/);
});

test("basis='lognorm' picks a log-normalized measure, used as-is (log1p false)", () => {
  const r = selectCountsBasis(ds({ X: M(["cells", "genes"], "measure", "scaled"), data: M(["cells", "genes"], "measure", "lognorm") }), { basis: "lognorm" });
  assert.deepEqual(r, { field: "data", log1p: false });
  // errors only when there's neither a lognorm-state measure NOR a conventional lognorm name (X/data/logcounts)
  assert.throws(() => selectCountsBasis(ds({ counts: M(["cells", "genes"], "measure", "raw") }), { basis: "lognorm" }), /no log-normalized measure/);
});

test("ignores non-measures and 1-D fields (only cells×genes measures are candidates)", () => {
  const r = selectCountsBasis(ds({ leiden: M(["cells"], "label"), od_score: M(["genes"], "measure", "raw"), counts: M(["cells", "genes"], "measure", "raw") }));
  assert.equal(r.field, "counts");
});

test("multi-sample cell axis (cells.s1) still qualifies as a cells×genes measure", () => {
  const r = selectCountsBasis(ds({ X: M(["cells.s1", "genes"], "measure", "raw") }));
  assert.deepEqual(r, { field: "X", log1p: true });
});
