"""lstar core in-memory model: Dataset, Axis, Field.

A dataset is a set of Axes (entities you index by) and Fields (typed data over axes).
This is the Python reference implementation of the L* model; see misc/Lstar_proposal.md.

Only `values` is required when adding a field; role / span / encoding are inferred
(the resolved view), and may be overridden to unlock role-specific behavior.
"""
from dataclasses import dataclass, field as _dcfield
from typing import Any, List, Optional

import numpy as np

OBSERVED = "observed"
DERIVED = "derived"


@dataclass
class Axis:
    name: str
    labels: Any
    origin: str = OBSERVED
    role: Optional[str] = None          # observation | feature | coordinate | None
    induced_by: Optional[str] = None
    provenance: dict = _dcfield(default_factory=dict)

    def __post_init__(self):
        self.labels = np.asarray(self.labels)

    def __len__(self) -> int:
        return int(self.labels.shape[0])


@dataclass
class Field:
    name: str
    values: Any                          # np.ndarray (dense) or scipy.sparse (csr/csc/coo)
    role: Optional[str] = None
    span: Optional[List[str]] = None
    state: Optional[str] = None
    encoding: Optional[str] = None       # dense | csr | csc | coo
    coverage: str = "full"
    directed: Optional[bool] = None
    weighted: Optional[bool] = None
    subtype: Optional[str] = None
    uncertainty: Optional[str] = None
    provenance: dict = _dcfield(default_factory=dict)


class Dataset:
    """A set of axes and fields. The unit of L* interchange."""

    def __init__(self, kind="sample", spec_version="0.1"):
        self.kind = kind
        self.spec_version = spec_version
        self.axes = {}        # name -> Axis
        self.fields = {}      # name -> Field
        self.models = {}      # name -> (later milestones)
        self.profiles = []    # which profiles wrote this dataset
        self.dropped = []     # native locations a profile could not represent (loss is recorded)

    # ---- axes ----
    def add_axis(self, name, labels, origin=OBSERVED, role=None,
                 induced_by=None, provenance=None):
        self.axes[name] = Axis(name, labels, origin, role, induced_by, provenance or {})
        return self.axes[name]

    def axis(self, name):
        return self.axes[name]

    # ---- fields ----
    def add_field(self, name, values, role=None, span=None, state=None, encoding=None,
                  coverage="full", directed=None, weighted=None, subtype=None,
                  uncertainty=None, provenance=None):
        span = self._infer_span(values, span)
        role = role or self._infer_role(values, span)
        encoding = encoding or self._infer_encoding(values)
        self.fields[name] = Field(
            name, values, role=role, span=span, state=state, encoding=encoding,
            coverage=coverage, directed=directed, weighted=weighted, subtype=subtype,
            uncertainty=uncertainty, provenance=provenance or {})
        return self.fields[name]

    def field(self, name):
        return self.fields[name]

    def fields_over(self, axis_name):
        return [f for f in self.fields.values() if f.span and axis_name in f.span]

    # ---- inference: only `values` is required; the rest is resolved ----
    def _infer_span(self, values, span):
        if span is not None:
            return list(span)
        shape = _shape(values)
        chosen = []
        for dim in shape:
            cands = [n for n, a in self.axes.items() if len(a) == dim]
            if len(cands) != 1:
                raise ValueError(
                    "cannot infer axis for a dimension of length %d among axes %s; "
                    "pass span=[...]" % (dim, list(self.axes)))
            chosen.append(cands[0])
        return chosen

    @staticmethod
    def _infer_role(values, span):
        if len(span) == 1:
            arr = np.asarray(values)
            return "label" if arr.dtype.kind in ("U", "S", "O") else "measure"
        return "measure"

    @staticmethod
    def _infer_encoding(values):
        import scipy.sparse as sp
        if sp.issparse(values):
            return values.getformat()
        return "dense"

    def __repr__(self):
        return "Dataset(kind=%r, axes=%s, fields=%s)" % (
            self.kind, list(self.axes), list(self.fields))


def _shape(values):
    import scipy.sparse as sp
    if sp.issparse(values):
        return tuple(values.shape)
    return tuple(np.asarray(values).shape)
