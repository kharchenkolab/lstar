"""Cross-implementation conformance: Python -> C++ -> Python must be identical.

  Python writes a store -> the C++ test_crossimpl reads it and writes a copy ->
  Python reads the C++-written copy and asserts equality with the original.

Requires the C++ binary built at core/build/test_crossimpl (cmake --build).
Run: PYTHONPATH=python/src python3 python/tests/test_crossimpl.py
"""
import os
import subprocess
import sys
import tempfile

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from test_roundtrip import make_ds  # noqa: E402

from lstar.zarr_io import write, read  # noqa: E402

_HERE = os.path.dirname(__file__)
CPP_BIN = os.environ.get(
    "LSTAR_CPP_BIN", os.path.abspath(os.path.join(_HERE, "..", "..", "core", "build", "test_crossimpl")))


def run():
    if not os.path.exists(CPP_BIN):
        raise SystemExit("C++ binary not found at %s (build core first)" % CPP_BIN)
    d = tempfile.mkdtemp()
    pin = os.path.join(d, "py.lstar.zarr")
    pout = os.path.join(d, "cpp.lstar.zarr")
    ds = make_ds()
    write(ds, pin)

    r = subprocess.run([CPP_BIN, pin, pout], capture_output=True, text=True)
    print(r.stdout.strip())
    if r.returncode != 0:
        raise SystemExit("C++ test_crossimpl failed:\n" + r.stderr)

    ds2 = read(pout)  # the C++-written store, read back in Python
    assert set(ds2.axes) == set(ds.axes)
    assert set(ds2.fields) == set(ds.fields)
    assert (np.asarray(ds2.axis("cells").labels) == np.asarray(ds.axis("cells").labels)).all()
    assert ds2.axis("pca").origin == "derived"
    assert ds2.field("counts").encoding == "csc" and ds2.field("counts").state == "raw"
    assert (ds2.field("counts").values.toarray() == ds.field("counts").values.toarray()).all()
    assert np.allclose(ds2.field("pca").values, ds.field("pca").values)
    assert (ds2.field("leiden").values == ds.field("leiden").values).all()
    print("cross-impl OK: Python -> C++ -> Python round-trip is byte-faithful")


def test_crossimpl():
    run()


if __name__ == "__main__":
    run()
