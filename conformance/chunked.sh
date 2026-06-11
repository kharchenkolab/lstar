#!/usr/bin/env bash
# Chunked + compressed cross-impl conformance: the C++ core must read a multi-chunk,
# gzip-compressed store written by Python, and its output (with a consolidated .zmetadata) must
# re-open in Python. Also exercises the csc<->csr transpose primitive. This is the portability
# path that the single-chunk-uncompressed MVP could not handle.
set -e
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
IN=/tmp/conf_chunked.lstar.zarr
OUT=/tmp/conf_chunked_cpp.lstar.zarr

read NNZ SUM < <(PYTHONPATH="$ROOT/python/src" python3 - "$IN" <<'PY'
import sys, warnings; warnings.filterwarnings("ignore")
import numpy as np, scipy.sparse as sp, numcodecs, lstar
rng = np.random.default_rng(7)
X = sp.csc_matrix(sp.random(400, 250, density=0.08, format="csc", random_state=rng))
X.data = X.data * 10 + 0.3
ds = lstar.Dataset(kind="sample")
ds.add_axis("cells", [f"cell{i}" for i in range(400)])
ds.add_axis("genes", [f"g{i}" for i in range(250)])
ds.add_axis("pca", [f"PC{i}" for i in range(10)], role="coordinate")
ds.add_field("counts", X, role="measure", span=["cells", "genes"], state="raw")
ds.add_field("pca", rng.standard_normal((400, 10)).astype("f4"), role="embedding", span=["cells", "pca"])
ds.add_field("leiden", np.array([f"c{i%5}" for i in range(400)]), role="label", span=["cells"])
lstar.write(ds, sys.argv[1], compressor=numcodecs.GZip(5), chunk_elems=3000)
print(int(X.nnz), float(X.data.sum()))
PY
)
echo "  [py] wrote chunked+gzip store (nnz=$NNZ)"

"$ROOT/core/build/test_chunked" "$IN" "$OUT" "$NNZ" "$SUM"

PYTHONPATH="$ROOT/python/src" python3 - "$IN" "$OUT" <<'PY'
import sys, os, warnings; warnings.filterwarnings("ignore")
import numpy as np, zarr, lstar
assert os.path.exists(sys.argv[2] + "/.zmetadata"), "C++ did not write .zmetadata"
zarr.open_consolidated(sys.argv[2], mode="r")
a, b = lstar.read(sys.argv[1]), lstar.read(sys.argv[2])
assert (a.fields["counts"].values != b.fields["counts"].values).nnz == 0
assert np.allclose(np.asarray(a.fields["pca"].values), np.asarray(b.fields["pca"].values))
assert (np.asarray(a.fields["leiden"].values) == np.asarray(b.fields["leiden"].values)).all()
print("  [py] re-opened C++ output (consolidated); counts/pca/leiden identical")
PY
echo "chunked conformance PASSED."
