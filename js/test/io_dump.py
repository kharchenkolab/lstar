#!/usr/bin/env python3
# Dump every array in an L* store (v2 or v3) as {path: {dtype, shape, b64}} plus the root manifest, so a
# Node harness can check the WASM/libzarr reader against the zarr-python reference byte-for-byte.
import os, sys, json, base64, zarr, numpy as np

store = sys.argv[1]

def array_keys(root):
    keys = []
    for dp, _, fns in os.walk(root):
        marker = "zarr.json" if "zarr.json" in fns else (".zarray" if ".zarray" in fns else None)
        if not marker:
            continue
        rel = os.path.relpath(dp, root)
        if rel == ".":
            continue
        if marker == ".zarray" or json.load(open(os.path.join(dp, "zarr.json"))).get("node_type") == "array":
            keys.append(rel.replace(os.sep, "/"))
    return keys

g = zarr.open_group(store, mode="r")
out = {"manifest": g.attrs["lstar"], "arrays": {}}
for k in array_keys(store):
    a = np.ascontiguousarray(np.asarray(g[k][...]))
    out["arrays"][k] = {"dtype": a.dtype.str, "shape": list(a.shape),
                        "b64": base64.b64encode(a.tobytes()).decode()}
json.dump(out, open(sys.argv[2], "w"))
print(f"  [py] dumped {len(out['arrays'])} arrays from {os.path.basename(store)}")
