"""Induction (model.md "Three induction rules"): a categorical `label` induces a **factor axis** whose
labels ARE the label's categories, so per-group results become ordinary fields over it. Covers eager
auto-induce, canonical identity (reuse on identical labels / error on a clash), the round-trip of the
induced axis + `induced_by`, and `validate()`'s drift check.

Run: PYTHONPATH=python/src python3 python/tests/test_induce.py
"""
import os
import tempfile

import numpy as np

import lstar
from lstar import Categorical


def _store():
    return os.path.join(tempfile.mkdtemp(), "induce.lstar.zarr")


def test_auto_induce_creates_factor_axis():
    ds = lstar.Dataset(kind="sample")
    ds.add_axis("cells", [f"c{i}" for i in range(6)])
    ds.add_field("leiden", Categorical(np.array([0, 1, 2, -1, 1, 0]),
                                       np.array(["A", "B", "C"]), ordered=True), span=["cells"])
    # a categorical label induced its factor axis (bare name, derived, role=factor, labels=categories)
    assert "leiden" in ds.axes
    ax = ds.axis("leiden")
    assert ax.role == "factor" and ax.origin == "derived" and ax.induced_by == "leiden"
    assert list(ax.labels) == ["A", "B", "C"]
    assert not lstar.validate(ds)
    # axis and field are separate namespaces -- the cells->category map and the category set coexist
    assert "leiden" in ds.fields and ds.field("leiden").role == "label"
    print("auto-induce: categorical label -> factor axis (labels = categories, induced_by set)")


def test_plain_string_label_does_not_induce():
    ds = lstar.Dataset(kind="sample")
    ds.add_axis("cells", [f"c{i}" for i in range(3)])
    ds.add_field("name", np.array(["x", "y", "z"]), span=["cells"])      # plain string, not categorical
    assert ds.field("name").role == "label" and ds.field("name").encoding != "categorical"
    assert "name" not in ds.axes                                        # no factor axis for free-text
    print("plain string label stays utf8 -> no factor axis induced")


def test_induced_axis_roundtrips():
    ds = lstar.Dataset(kind="sample")
    ds.add_axis("cells", [f"c{i}" for i in range(5)])
    ds.add_field("ct", Categorical(np.array([0, 1, 0, 2, 1]),
                                   np.array(["T", "B", "NK"])), span=["cells"])
    p = _store()
    lstar.write(ds, p)
    ds2 = lstar.read(p)
    assert "ct" in ds2.axes and ds2.axis("ct").role == "factor"
    assert ds2.axis("ct").induced_by == "ct"
    assert list(ds2.axis("ct").labels) == ["T", "B", "NK"]
    assert not lstar.validate(ds2)                                      # consistency survives the store
    print("induced factor axis round-trips through the store (labels + induced_by + factor role)")


def test_canonical_identity_reuse_and_collision():
    # idempotent: re-inducing the same field returns the same axis, no duplicate
    ds = lstar.Dataset(kind="sample")
    ds.add_axis("cells", [f"c{i}" for i in range(4)])
    ds.add_field("g", Categorical(np.array([0, 1, 0, 1]), np.array(["lo", "hi"])), span=["cells"])
    assert ds.induce("g") == "g" and ds.induce("g") == "g"             # reuse, not re-create
    n_axes = len(ds.axes)
    assert ds.induce("g") == "g" and len(ds.axes) == n_axes

    # a hand-declared axis with the SAME labels is adopted (induced_by linked) -- independent results align
    ds2 = lstar.Dataset(kind="sample")
    ds2.add_axis("cells", [f"c{i}" for i in range(4)])
    ds2.add_axis("grp", ["lo", "hi"], origin="derived", role="factor")  # pre-existing, no induced_by
    ds2.add_field("grp", Categorical(np.array([0, 1, 0, 1]), np.array(["lo", "hi"])), span=["cells"])
    assert ds2.axis("grp").induced_by == "grp"                         # adopted on add

    # a name clash with DIFFERENT labels: explicit induce raises, auto-induce skips + validate warns
    ds3 = lstar.Dataset(kind="sample")
    ds3.add_axis("cells", [f"c{i}" for i in range(4)])
    ds3.add_axis("clash", ["p", "q", "r"], role="factor")             # different labels
    ds3.add_field("clash", Categorical(np.array([0, 1, 0, 1]), np.array(["lo", "hi"])), span=["cells"])
    assert ds3.axis("clash").labels.tolist() == ["p", "q", "r"]        # auto-induce did NOT overwrite
    try:
        ds3.induce("clash"); assert False, "expected a collision error"
    except ValueError as e:
        assert "different labels" in str(e)
    assert any("name clash" in i for i in lstar.validate(ds3))         # surfaced, never silent
    print("canonical identity: idempotent reuse, hand-axis adoption, collision raises + validate warns")


def test_validate_catches_induced_axis_drift():
    ds = lstar.Dataset(kind="sample")
    ds.add_axis("cells", [f"c{i}" for i in range(4)])
    ds.add_field("k", Categorical(np.array([0, 1, 0, 1]), np.array(["a", "b"])), span=["cells"])
    assert not lstar.validate(ds)
    ds.axis("k").labels = np.array(["a", "ZZZ"])                        # drift the induced axis
    assert any("induced-axis drift" in i for i in lstar.validate(ds))
    print("validate catches drift between a factor axis and its inducing field's categories")


if __name__ == "__main__":
    test_auto_induce_creates_factor_axis()
    test_plain_string_label_does_not_induce()
    test_induced_axis_roundtrips()
    test_canonical_identity_reuse_and_collision()
    test_validate_catches_induced_axis_drift()
