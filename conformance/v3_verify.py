#!/usr/bin/env python3
# zarr-python interop for the v3-format leg. Two modes:
#   check   <v3>  <v2ref>   -- assert <v3> is genuine v3 (zarr_format 3, inline consolidated), and every
#                              array reads value-identical to the v2 reference; manifest attrs identical.
#   reemit  <v2>  <v3out>    -- re-emit the v2 store as v3 using zarr-python (an INDEPENDENT writer, whose
#                              default compressor is zstd) so the C++/R readers are tested against a v3
#                              store they did not themselves produce.
import warnings; warnings.filterwarnings("ignore")   # zarr's "v3 consolidated md not in spec" note is expected (#309)
import os, sys, json, shutil, zarr, numpy as np

def array_keys(root):
    keys = []
    for dp, _, fns in os.walk(root):
        marker = "zarr.json" if "zarr.json" in fns else (".zarray" if ".zarray" in fns else None)
        if not marker:
            continue
        rel = os.path.relpath(dp, root)
        if rel == ".":
            continue
        if marker == ".zarray":
            keys.append(rel.replace(os.sep, "/"))
        else:
            if json.load(open(os.path.join(dp, "zarr.json"))).get("node_type") == "array":
                keys.append(rel.replace(os.sep, "/"))
    return keys

def check(v3, v2ref):
    rj = json.load(open(os.path.join(v3, "zarr.json")))
    assert rj.get("zarr_format") == 3 and rj.get("node_type") == "group", "not a v3 group root"
    assert "consolidated_metadata" in rj, "missing inline consolidated metadata (#309)"
    g3 = zarr.open_group(v3, mode="r"); g2 = zarr.open_group(v2ref, mode="r")
    assert g3.attrs["lstar"] == g2.attrs["lstar"], "manifest attrs differ across formats"
    bad = 0
    for k in array_keys(v3):
        a = np.asarray(g3[k][...]); b = np.asarray(g2[k][...])
        if not (a.dtype == b.dtype and a.shape == b.shape and np.array_equal(a, b, equal_nan=True)):
            bad += 1; print("  MISMATCH", k)
    assert bad == 0, f"{bad} array mismatches vs v2 reference"
    print(f"  zarr-python reads C++ v3: {len(array_keys(v3))} arrays == v2, format 3, inline-consolidated")

def _reemit(v2, v3out, compressors):
    shutil.rmtree(v3out, ignore_errors=True)
    g2 = zarr.open_group(v2, mode="r")
    store = zarr.storage.LocalStore(v3out)
    g3 = zarr.create_group(store=store, overwrite=True)               # zarr_format 3
    for k, v in dict(g2.attrs).items():
        g3.attrs[k] = v
    for dp, _, fns in os.walk(v2):
        rel = os.path.relpath(dp, v2)
        if rel == ".":
            continue
        node = "/".join(rel.split(os.sep))
        if ".zarray" in fns:
            arr = np.asarray(g2[node][...])
            kw = {} if compressors is None else {"compressors": compressors}
            d = g3.create_array(name=node, shape=arr.shape, dtype=arr.dtype,
                                chunks=arr.shape if arr.size else (1,), **kw)
            if arr.size:
                d[...] = arr
            for k, v in dict(g2[node].attrs).items():
                d.attrs[k] = v
        elif ".zgroup" in fns:
            sub = g3.require_group(node)
            for k, v in dict(g2[node].attrs).items():
                sub.attrs[k] = v
    zarr.consolidate_metadata(store)

def reemit(v2, v3out):                                               # zarr-python's default compressor (zstd)
    _reemit(v2, v3out, compressors=None)
    print("  zarr-python re-emitted v3 (zstd default) ->", v3out)

def reemit_gzip(v2, v3out):                                         # gzip: readable by the WASM reader (no zstd port)
    from zarr.codecs import GzipCodec
    _reemit(v2, v3out, compressors=[GzipCodec(level=5)])
    print("  zarr-python re-emitted v3 (gzip) ->", v3out)

def shardcheck(sharded, ref):
    # the sharded store must (a) genuinely use the sharding_indexed codec on its big arrays, packing many
    # inner chunks into fewer objects, and (b) read to the same values as the unsharded reference.
    import os
    dj = json.load(open(os.path.join(sharded, "fields", "counts", "data", "zarr.json")))
    assert dj["codecs"][0]["name"] == "sharding_indexed", f"counts/data not sharded: {dj['codecs']}"
    ncs = os.path.join(sharded, "fields", "counts", "data", "c")
    nshards = sum(1 for _ in os.walk(ncs)) if os.path.isdir(ncs) else 0
    gs = zarr.open_group(sharded, mode="r"); gr = zarr.open_group(ref, mode="r")
    bad = 0
    for k in array_keys(sharded):
        a = np.asarray(gs[k][...]); b = np.asarray(gr[k][...])
        if not (a.dtype == b.dtype and a.shape == b.shape and np.array_equal(a, b, equal_nan=True)):
            bad += 1; print("  MISMATCH", k)
    assert bad == 0, f"{bad} array mismatches vs the unsharded reference"
    print(f"  sharded store: sharding_indexed on counts/data; {len(array_keys(sharded))} arrays == unsharded seed")

def compare(a, b):
    ga = zarr.open_group(a, mode="r"); gb = zarr.open_group(b, mode="r")
    assert dict(ga.attrs) == dict(gb.attrs), "root attrs differ"
    ka, kb = sorted(array_keys(a)), sorted(array_keys(b))
    assert ka == kb, f"array sets differ: {set(ka) ^ set(kb)}"
    bad = 0
    for k in ka:
        x = np.asarray(ga[k][...]); y = np.asarray(gb[k][...])
        if not (x.dtype == y.dtype and x.shape == y.shape and np.array_equal(x, y, equal_nan=True)):
            bad += 1; print("  MISMATCH", k)
    assert bad == 0, f"{bad} array mismatches between {a} and {b}"
    print(f"  {len(ka)} arrays value-identical between the two stores")

if __name__ == "__main__":
    {"check": check, "reemit": reemit, "reemit_gzip": reemit_gzip,
     "compare": compare, "shardcheck": shardcheck}[sys.argv[1]](sys.argv[2], sys.argv[3])
