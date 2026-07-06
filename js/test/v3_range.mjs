// v3-aware byte-range fast path (the pagoda3 ship-blocker): on an UNCOMPRESSED store, cscColumn/csrRow
// must issue exact byte-range reads (not whole-array reads) on BOTH v2 and v3 — the v3 case forms the
// chunk key in the v3 default encoding (c/0), computed by libzarr, with the fetch loop in JS — and return
// identical values across formats. Args: <v2dir> <v3dir> (same data, uncompressed, one v2 + one v3).
import { openLstar } from "../core/reader.ts";
import { NodeFSStore } from "../core/node-store.ts";

const [v2dir, v3dir] = process.argv.slice(2);
let fail = 0;
const check = (n, ok) => { if (!ok) { console.log(`  FAIL  ${n}`); fail++; } };

// counts byte-range (getRange) vs whole-chunk (get on a chunk key) so we can prove which path ran.
class Counting {
  constructor(inner) { this.inner = inner; this.ranges = 0; this.chunkGets = 0; }
  get(k) { if (/\/(c\/)?\d+(?:[./]\d+)*$/.test(k)) this.chunkGets++; return this.inner.get(k); }
  getRange(k, s, e) { this.ranges++; return this.inner.getRange(k, s, e); }
}
const eq = (a, b) => { if (a.length !== b.length) return false; for (let i = 0; i < a.length; i++) if (Number(a[i]) !== Number(b[i])) return false; return true; };

const cols = {};
for (const [tag, dir] of [["v2", v2dir], ["v3", v3dir]]) {
  const st = new Counting(new NodeFSStore(dir));
  const ds = await openLstar(st);
  const col = await ds.cscColumn("counts", 7);
  const row = await ds.csrRow("counts_cellmajor", 11);
  // byte-range fast path took: getRange used, no whole-chunk gets
  check(`${tag} cscColumn+csrRow via byte-range (getRange used, no chunk gets)`, st.ranges > 0 && st.chunkGets === 0);
  check(`${tag} cscColumn(7) has nonzeros`, col.rows.length > 0);
  cols[tag] = { col, row };
}
// identical values across formats
check("cscColumn v2 == v3 (rows)", eq(cols.v2.col.rows, cols.v3.col.rows));
check("cscColumn v2 == v3 (vals)", eq(cols.v2.col.vals, cols.v3.col.vals));
check("csrRow v2 == v3 (cols)", eq(cols.v2.row.cols, cols.v3.row.cols));
check("csrRow v2 == v3 (vals)", eq(cols.v2.row.vals, cols.v3.row.vals));

console.log(fail === 0
  ? "  v3-aware byte-range: cscColumn/csrRow take the fast path on v2 AND v3; values identical"
  : `  ${fail} v3 byte-range failures`);
process.exit(fail === 0 ? 0 : 1);
