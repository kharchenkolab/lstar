"""Lossless passthrough of a format's untyped long-tail (AnnData `uns`, Seurat `@misc`/`@commands`).

Upgrades `dropped` from *name-only* to **verbatim preserve + reproduce**. The store carries the payload
as a *self-describing* subtree -- a JSON skeleton (`tree`) whose array leaves are references into a flat
`arrays` manifest -- so it stays **inspectable** (a recognized structure can be promoted to a typed
field later) and round-trips through every reader using L*'s own primitives (JSON attrs + typed arrays).
The core (C++/R) round-trips `tree` + the manifest *verbatim*, never interpreting it; only the
originating profile walks the tree to rebuild the live object.

Tree grammar (each node):
  - a JSON scalar (None / bool / int / float / str)                        -- inline
  - {"$obj":   {key: node, ...}}                                           -- a dict
  - {"$list":  [node, ...]}                                                -- a list/tuple
  - {"$array": id}                                                         -- a numeric/bool ndarray
  - {"$strings": id, "shape": [...]}                                       -- a string ndarray
  - {"$record": {"fields": {name: node}, "length": n}}                     -- a numpy structured array
  - {"$bytes": "<latin1>"} / {"$dropped": "<typename>"}                    -- raw bytes / unrepresentable

`to_store(obj) -> (tree, arrays)` where `arrays` is a list of `{"id", "kind": "dense"|"utf8", "data"}`;
`from_store(tree, arrays) -> obj` is the exact inverse. Genuinely unrepresentable leaves (sparse,
opaque objects) become `$dropped` -- recorded, never silently lost.
"""
import numpy as np


def to_store(obj):
    arrays = []

    def emit_array(a):
        a = np.asarray(a)
        i = "a%d" % len(arrays)
        if a.dtype.kind in ("U", "S", "O"):
            arrays.append({"id": i, "kind": "utf8", "data": a})
            return {"$strings": i, "shape": [int(x) for x in a.shape]}
        arrays.append({"id": i, "kind": "dense", "data": a})
        return {"$array": i}

    def emit(v):
        if v is None or isinstance(v, (bool, str)):
            return v
        if isinstance(v, np.bool_):
            return bool(v)
        if isinstance(v, (int, np.integer)):
            return int(v)
        if isinstance(v, (float, np.floating)):
            return float(v)
        if isinstance(v, bytes):
            return {"$bytes": v.decode("latin1")}
        if isinstance(v, np.ndarray):
            if v.ndim == 0:
                return emit(v.item())
            if v.dtype.names:                              # structured / record array
                fields = {n: emit_array(np.asarray(v[n])) for n in v.dtype.names}
                return {"$record": {"fields": fields, "length": int(v.shape[0])}}
            return emit_array(v)
        if isinstance(v, dict):
            return {"$obj": {str(k): emit(val) for k, val in v.items()}}
        if isinstance(v, (list, tuple)):
            return {"$list": [emit(x) for x in v]}
        try:                                               # a numpy scalar of some other kind
            return emit(np.asarray(v).item()) if np.ndim(v) == 0 else {"$dropped": type(v).__name__}
        except Exception:
            return {"$dropped": type(v).__name__}

    return emit(obj), arrays


def from_store(tree, arrays):
    by_id = {a["id"]: a["data"] for a in arrays}

    def build(node):
        if not isinstance(node, dict):
            return node                                    # JSON scalar
        if "$obj" in node:
            return {k: build(v) for k, v in node["$obj"].items()}
        if "$list" in node:
            return [build(x) for x in node["$list"]]
        if "$array" in node:
            return np.asarray(by_id[node["$array"]])
        if "$strings" in node:
            arr = np.asarray(by_id[node["$strings"]], dtype=str)
            shp = node.get("shape")
            return arr.reshape(shp) if shp and len(shp) > 1 else arr
        if "$record" in node:
            rec = node["$record"]
            cols = {n: np.asarray(build(fn)) for n, fn in rec["fields"].items()}
            dt = np.dtype([(n, cols[n].dtype) for n in cols])
            out = np.empty(int(rec["length"]), dtype=dt)
            for n in cols:
                out[n] = cols[n]
            return out
        if "$bytes" in node:
            return node["$bytes"].encode("latin1")
        if "$dropped" in node:
            return None
        return {k: build(v) for k, v in node.items()}      # defensive: a bare dict

    return build(tree)
