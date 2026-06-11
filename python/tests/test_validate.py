"""Validator tests: well-formed datasets pass; malformed ones are caught."""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from test_roundtrip import make_ds  # noqa: E402

import lstar  # noqa: E402


def _errors(ds):
    return [i for i in lstar.validate(ds) if i.startswith("ERROR")]


def run():
    # the standard sample validates clean
    ds = make_ds()
    assert _errors(ds) == [], _errors(ds)
    lstar.validate(ds, strict=True)  # does not raise

    # an AnnData-derived dataset validates clean too
    try:
        import anndata  # noqa: F401
        sys.path.insert(0, os.path.dirname(__file__))
        from test_anndata_profile import make_adata
        ds_a = lstar.read_anndata(make_adata())
        assert _errors(ds_a) == [], _errors(ds_a)
    except ImportError:
        pass

    # span referencing a missing axis -> ERROR
    bad = make_ds()
    bad.fields["pca"].span = ["cells", "no_such_axis"]
    assert any("unknown axes" in e for e in _errors(bad)), _errors(bad)

    # shape mismatch -> ERROR
    bad2 = make_ds()
    bad2.add_field("wrong", np.zeros((50, 5)), role="measure", span=["genes", "pca"])
    # genes=50 ok but pca=10 != 5 -> error
    assert any("wrong" in e for e in _errors(bad2)), _errors(bad2)

    # strict raises
    try:
        lstar.validate(bad, strict=True)
        raise AssertionError("strict validate should have raised")
    except ValueError:
        pass

    print("validate OK: clean datasets pass; missing-axis and shape-mismatch caught")


def test_validate():
    run()


if __name__ == "__main__":
    run()
