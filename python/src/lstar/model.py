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


class Categorical:
    """A categorical / factor value: integer codes into an ordered set of category labels.

    `codes[i] == k` means element i is `categories[k]`; `-1` is **missing**. This is the in-memory form
    of the `categorical` encoding -- the one mechanism behind dtype-faithful labels (the category set,
    order, and missingness all survive a round trip) and, via induction, factor axes. `np.asarray(cat)`
    decodes to label strings (missing -> "").
    """

    def __init__(self, codes, categories, ordered=False):
        self.codes = np.asarray(codes).astype(np.int64, copy=False)
        self.categories = np.asarray(categories, dtype=str)
        self.ordered = bool(ordered)

    @property
    def shape(self):
        return (int(self.codes.shape[0]),)

    def __len__(self):
        return int(self.codes.shape[0])

    def __array__(self, dtype=None):
        safe = np.clip(self.codes, 0, max(len(self.categories) - 1, 0))
        out = np.where(self.codes >= 0, self.categories[safe], "")
        return out.astype(dtype) if dtype is not None else out

    def __repr__(self):
        return "Categorical(n=%d, k=%d, ordered=%s)" % (len(self), len(self.categories), self.ordered)


def _is_categorical(v):
    """True for an L* Categorical or a duck-typed pandas.Categorical (.codes/.categories, not ndarray)."""
    return isinstance(v, Categorical) or (
        not isinstance(v, np.ndarray) and hasattr(v, "codes") and hasattr(v, "categories"))


def as_categorical(v, ordered=None):
    """Coerce a Categorical / pandas.Categorical into an L* `Categorical` (codes, categories, ordered)."""
    if isinstance(v, Categorical):
        return v if ordered is None else Categorical(v.codes, v.categories, ordered)
    codes = np.asarray(v.codes)                                # pandas.Categorical: -1 is missing
    cats = np.asarray(v.categories, dtype=str)
    return Categorical(codes, cats, bool(v.ordered if ordered is None else ordered))


@dataclass
class Axis:
    name: str
    labels: Any
    origin: str = OBSERVED
    role: Optional[str] = None          # observation | feature | coordinate | factor | None
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
    mask: Any = None                     # optional uint8 validity mask, 1 == missing (nullable Int/bool/string)
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
                  uncertainty=None, mask=None, provenance=None):
        if _is_categorical(values):
            values = as_categorical(values)               # normalize pandas.Categorical -> L* Categorical
        span = self._infer_span(values, span)
        role = role or self._infer_role(values, span)
        encoding = encoding or self._infer_encoding(values)
        if mask is not None:
            mask = np.asarray(mask).astype(np.uint8, copy=False)   # 1 == missing
        self.fields[name] = Field(
            name, values, role=role, span=span, state=state, encoding=encoding,
            coverage=coverage, directed=directed, weighted=weighted, subtype=subtype,
            uncertainty=uncertainty, mask=mask, provenance=provenance or {})
        if _is_categorical(values):            # a categorical label induces its factor axis (data-driven)
            self._auto_induce(name)
        return self.fields[name]

    def field(self, name):
        return self.fields[name]

    # ---- induction (model.md "Three induction rules") ----
    def induce(self, field_name):
        """Derive the axis a field induces, idempotently, recording the link (`induced_by`).

        A categorical `label` induces a **factor axis**: its categories, in category order, become an
        ordered set of axis labels (`role="factor"`, `origin="derived"`). Per-group results -- DE,
        pseudobulk, PAGA -- are then ordinary fields *over this axis* rather than special cases. The
        categories live once (on the field); the axis carries the same labels and an `induced_by` back
        link, so `validate()` can confirm they never drift.

        **Canonical identity = the field's name (bare) + its ordered label set.** Inducing a label whose
        canonical axis already exists *with the same labels* returns it (so DE from scanpy and markers
        from pagoda2 land on one axis and align); a name clash with *different* labels is an error the
        writer must resolve -- never a silent merge of two distinct label sets. Returns the axis name.
        """
        fl = self.fields[field_name]
        if not _is_categorical(fl.values):
            raise ValueError(
                "induce(%r): only a categorical label induces an axis (encoding=%r)"
                % (field_name, fl.encoding))
        labels = np.asarray(as_categorical(fl.values).categories, dtype=str)
        name = field_name                          # canonical name = the inducing field's name, bare
        existing = self.axes.get(name)
        if existing is not None:
            ex = np.asarray(existing.labels, dtype=str)
            if ex.shape == labels.shape and bool(np.all(ex == labels)):
                if existing.induced_by is None:    # adopt the link if reusing a hand-declared axis
                    existing.induced_by = field_name
                return name                        # same identity -> reuse (independent results align)
            raise ValueError(
                "induce(%r): axis '%s' already exists with different labels (%d vs %d); rename the "
                "field or the axis -- collisions are never silently merged"
                % (field_name, name, len(ex), len(labels)))
        self.add_axis(name, labels, origin=DERIVED, role="factor", induced_by=field_name)
        return name

    def _auto_induce(self, field_name):
        """Eager auto-induce on `add_field` of a categorical: create the factor axis, but never abort a
        field add over a name clash -- a genuine clash (a categorical named like an existing axis with
        different labels) is left for `validate()` to surface rather than raised here."""
        try:
            self.induce(field_name)
        except ValueError:
            pass

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
        if _is_categorical(values):
            return "label"
        if len(span) == 1:
            arr = np.asarray(values)
            return "label" if arr.dtype.kind in ("U", "S", "O") else "measure"
        return "measure"

    @staticmethod
    def _infer_encoding(values):
        import scipy.sparse as sp
        if _is_categorical(values):
            return "categorical"
        if _is_stream_source(values):       # a chunked sparse source (LazyCSX / backed h5ad)
            return values.fmt
        if sp.issparse(values):
            return values.getformat()
        return "dense"

    def __repr__(self):
        return "Dataset(kind=%r, axes=%s, fields=%s)" % (
            self.kind, list(self.axes), list(self.fields))


def _is_stream_source(v):
    """A streaming sparse source: yields scipy blocks without materializing the whole matrix.
    Used so `add_field` can hold a backed/lazy source and `write` can stream it (low-memory writes)."""
    return hasattr(v, "blocks") and hasattr(v, "fmt") and hasattr(v, "shape")


def _shape(values):
    import scipy.sparse as sp
    if _is_categorical(values):
        return (len(values),)
    if sp.issparse(values) or _is_stream_source(values):
        return tuple(values.shape)
    return tuple(np.asarray(values).shape)
