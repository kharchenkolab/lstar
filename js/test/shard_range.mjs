// Sharding byte-range hot path (B1). On an UNCOMPRESSED store a sharded array's inner chunks live INSIDE
// shard objects, so the plain chunk-key range read does not map 1:1. The reader resolves the chunk through
// the shard index (libzarr shardLocate/shardEntry): suffix-read the small index, then range-read exactly the
// chunk's bytes — so cscColumn on a SHARDED store (a) returns the same nonzeros as the UNSHARDED store, and
// (b) STREAMS: it uses getSuffix + getRange and never whole-object-reads a shard. Args: <unsharded> <sharded>.
import { openLstar } from "../core/reader.ts";
import { NodeFSStore } from "../core/node-store.ts";

const [udir, sdir] = process.argv.slice(2);
let fail = 0;
const check = (n, ok) => { if (!ok) { console.log(`  FAIL  ${n}`); fail++; } };
const eq = (a, b) => { if (a.length !== b.length) return false; for (let i = 0; i < a.length; i++) if (Number(a[i]) !== Number(b[i])) return false; return true; };
const isChunk = (k) => /\/c\/\d+$/.test(k);   // a data-chunk / shard-object key ".../c/<n>"

// Counting wrapper around the sharded store, to prove which store ops the read path used.
class Counting {
  constructor(inner) { this.inner = inner; this.gets = []; this.ranges = 0; this.suffixes = 0; }
  async get(k) { this.gets.push(k); return this.inner.get(k); }
  async getRange(k, s, e) { this.ranges++; return this.inner.getRange(k, s, e); }
  async getSuffix(k, n) { this.suffixes++; return this.inner.getSuffix(k, n); }
}

const u = await openLstar(new NodeFSStore(udir));
const cs = new Counting(new NodeFSStore(sdir));
const s = await openLstar(cs);

// (a) VALUE parity across columns spanning several chunks (so we hit chunks at various intra-shard slots).
const ncols = (u.field("counts").shape)[1];
let ncok = 0;
for (const c of [0, 1, 3, 7, Math.floor(ncols / 2), ncols - 1]) {
  if (c >= ncols) continue;
  const cu = await u.cscColumn("counts", c), cx = await s.cscColumn("counts", c);
  check(`cscColumn(${c}) sharded == unsharded`, eq(cu.rows, cx.rows) && eq(cu.vals, cx.vals));
  ncok++;
}

// (b) STREAMING: disable the LRU read cache (else a re-read serves from memory with no store ops), read a
// FRESH deep column, and assert the shard-resolve path fired — index via getSuffix, chunk bytes via
// getRange, and NO whole-shard object read.
s.setReadCacheBudget(0);
cs.gets.length = 0; cs.ranges = 0; cs.suffixes = 0;
const deep = Math.min(ncols - 1, Math.max(2, Math.floor(ncols * 0.84)));   // a column deep inside a shard (intra > 0)
await s.cscColumn("counts", deep);
const wholeShardGets = cs.gets.filter(isChunk).length;
check("sharded read STREAMS: suffix-reads the shard index (getSuffix used)", cs.suffixes >= 1);
check("sharded read STREAMS: range-reads chunk bytes (getRange used)", cs.ranges >= 1);
check("sharded read STREAMS: never whole-object-reads a shard", wholeShardGets === 0);

console.log(fail === 0
  ? `  sharded byte-range reads == unsharded across ${ncok} columns, and STREAM (index via suffix, chunk via range)`
  : `  ${fail} sharded byte-range failures`);
process.exit(fail === 0 ? 0 : 1);
