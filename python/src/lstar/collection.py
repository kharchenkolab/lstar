"""Assemble a *collection* of heterogeneous samples (not an aligned tensor) from per-sample sources.

A collection keeps each sample's own data --- per-sample ``cells.<s>``/``genes.<s>`` axes and
``<field>.<s>`` measures, possibly over **divergent or even disjoint** gene sets --- alongside a
``samples`` axis and a *union* ``cells`` axis carrying the joint analysis (a shared embedding, clusters,
a graph). ``collection_from`` builds exactly the structure the Conos / Seurat-v5-split profiles produce,
but from any list of per-sample objects, so the common "I have N separately-processed samples" workflow
is first-class rather than something you hand-assemble.
"""
from __future__ import annotations

import numpy as np

from .model import DERIVED, Categorical, Dataset


def _as_dataset(obj):
    """Read one per-sample source into an L* :class:`Dataset` (passthrough if already one)."""
    if isinstance(obj, Dataset):
        return obj
    try:
        import anndata as ad
        if isinstance(obj, ad.AnnData):
            from .profiles.anndata import read_anndata
            return read_anndata(obj)
    except ImportError:
        pass
    try:
        import mudata as md
        if isinstance(obj, md.MuData):
            from .profiles.mudata import read_mudata
            return read_mudata(obj)
    except ImportError:
        pass
    raise TypeError(
        f"collection_from: each sample must be an lstar.Dataset or AnnData/MuData, got {type(obj).__name__}")


def collection_from(samples, joint=None, sample_field="sample", prefix_cells=True):
    """Build a ``kind='collection'`` :class:`Dataset` from per-sample sources.

    Parameters
    ----------
    samples : dict ``{name: Dataset|AnnData|MuData}`` or a list of those (auto-named ``s0, s1, ...``).
        Each sample keeps its own ``cells.<name>``/``genes.<name>`` axes and ``<field>.<name>`` fields ---
        gene sets may overlap, differ, or be entirely disjoint across samples.
    joint : optional dict ``{name: value}`` over the **union** cells --- the integration outputs. A 2-D
        array becomes a joint embedding; a categorical/label becomes a clustering (inducing a factor axis);
        an ``(n_union x n_union)`` sparse matrix becomes a cell-cell relation (graph).
    sample_field : name of the design label recording each union cell's sample (default ``"sample"``).
    prefix_cells : prefix each cell label with its sample name so union labels are unique (default True).
    """
    items = list(samples.items()) if isinstance(samples, dict) else [(f"s{i}", o) for i, o in enumerate(samples)]
    if not items:
        raise ValueError("collection_from: no samples given")

    ds = Dataset(kind="collection")
    union_cells, sample_of, n_per = [], [], []
    for s, src in items:
        dss = _as_dataset(src)
        if "cells" not in dss.axes or "genes" not in dss.axes:
            raise ValueError(f"collection_from: sample {s!r} has no 'cells'/'genes' axes")
        cell_labels = list(np.asarray(dss.axis("cells").labels, dtype=str))
        if prefix_cells:
            cell_labels = [f"{s}_{c}" for c in cell_labels]
        rename = {}
        for ax in dss.axes:                                   # namespace every axis as "<axis>.<s>"
            a = dss.axis(ax)
            rename[ax] = f"{ax}.{s}"
            if a.role == "factor" and a.induced_by:           # re-induced below by its (renamed) label field
                continue
            labs = cell_labels if ax == "cells" else list(np.asarray(a.labels, dtype=str))
            ds.add_axis(rename[ax], labs, origin=a.origin, role=a.role)
        for fn, fl in dss.fields.items():                     # namespace every field as "<field>.<s>"
            ds.add_field(f"{fn}.{s}", fl.values, role=fl.role, span=[rename[a] for a in (fl.span or [])],
                         state=fl.state, encoding=fl.encoding, subtype=fl.subtype, mask=fl.mask,
                         coverage=fl.coverage, index=fl.index,
                         index_axis=rename.get(fl.index_axis) if fl.index_axis else None,
                         provenance=dict(fl.provenance) if fl.provenance else None)
        union_cells.extend(cell_labels)
        sample_of.extend([s] * len(cell_labels))
        n_per.append(len(cell_labels))

    ds.add_axis("samples", [s for s, _ in items], role="sample")
    ds.add_field("n_cells", np.asarray(n_per, dtype="int64"), role="measure", span=["samples"])

    if len(set(union_cells)) != len(union_cells):
        raise ValueError("collection_from: union cell labels are not unique (use prefix_cells=True)")
    ds.add_axis("cells", union_cells, origin=DERIVED, role="observation")
    ds.add_field(sample_field, np.asarray(sample_of, dtype=str), role="label", span=["cells"], subtype="design")

    for name, val in (joint or {}).items():
        _add_joint(ds, name, val)
    return ds


def _add_joint(ds, name, val):
    """Type one joint (union-cells) field: sparse -> relation, 2-D -> embedding, else clustering label."""
    import scipy.sparse as sp
    if sp.issparse(val):
        ds.add_field(name, val, role="relation", span=["cells", "cells"], subtype="knn")
    elif isinstance(val, Categorical) or (not np.isscalar(val) and np.asarray(val).ndim == 2):
        if isinstance(val, Categorical):
            ds.add_field(name, val, role="label", span=["cells"])   # induces a factor axis
        else:
            v = np.asarray(val)
            ds.add_axis(name, [f"{name}{i}" for i in range(v.shape[1])], origin=DERIVED, role="coordinate")
            ds.add_field(name, v, role="embedding", span=["cells", name])
    else:
        ds.add_field(name, val, role="label", span=["cells"])       # 1-D clustering -> induces a factor axis
