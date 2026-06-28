"""Structural validation of an L* Dataset / store.

Checks the invariants any conforming writer (Python, C++, R) must satisfy: spans reference
existing axes, arities/shapes match axis lengths, relations span two axes, etc. The role and
state vocabularies are *open*, so unrecognized values are warnings, not errors.

    issues = lstar.validate(ds)            # list of "ERROR: ..."/"WARN: ..." strings
    lstar.validate(ds, strict=True)        # raises ValueError on any ERROR
"""
import re

import numpy as np
import scipy.sparse as sp

CORE_ROLES = {"measure", "embedding", "loading", "relation", "label", "sequence",
              "design", "transform", "vector", "probability", "coordinate", "factor"}
CORE_STATES = {None, "", "raw", "lognorm", "scaled", "clr", "permutation"}


def validate(ds, strict=False):
    issues = []
    err = lambda m: issues.append("ERROR: " + m)
    warn = lambda m: issues.append("WARN: " + m)

    axlen = {n: len(a) for n, a in ds.axes.items()}

    for name, f in ds.fields.items():
        span = list(f.span or [])
        missing = [ax for ax in span if ax not in ds.axes]
        if missing:
            err("field '%s' span references unknown axes %s" % (name, missing))
            continue

        # Partial coverage: the field covers only a subset of `index_axis`, keyed by `index` (int
        # positions into that axis). The effective length along that axis is len(index), so shape checks
        # compare against `flen` (axis lengths with the index_axis overridden), and the index itself must
        # be in range. A modality measured on a cell subset (multiome / CITE-seq dropout) uses this.
        flen = axlen
        if getattr(f, "index", None) is not None:
            idx = np.asarray(f.index); iax = f.index_axis
            if iax not in span:
                err("field '%s' partial index_axis '%s' not in span %s" % (name, iax, span))
                continue
            if idx.size and (idx.min() < 0 or idx.max() >= axlen[iax]):
                err("field '%s' partial index out of range [0,%d) for axis '%s'"
                    % (name, axlen[iax], iax))
            if idx.size and np.unique(idx).size != idx.size:
                warn("field '%s' partial index has duplicate positions" % name)
            flen = dict(axlen); flen[iax] = int(idx.size)

        enc = f.encoding
        if enc in ("csr", "csc", "coo"):
            if len(span) != 2:
                err("field '%s' (%s) must span 2 axes, got %s" % (name, enc, span))
            elif sp.issparse(f.values):
                shp, exp = tuple(f.values.shape), (flen[span[0]], flen[span[1]])
                if shp != exp:
                    err("field '%s' sparse shape %s != axis lengths %s" % (name, shp, exp))
        elif enc == "categorical":
            cat = f.values
            n = len(cat)
            if len(span) == 1 and n != flen[span[0]]:
                err("field '%s' categorical length %d != axis '%s' length %d"
                    % (name, n, span[0], flen[span[0]]))
            k = len(getattr(cat, "categories", []))
            codes = np.asarray(getattr(cat, "codes", []))
            if codes.size and (codes.min() < -1 or codes.max() >= k):
                err("field '%s' categorical codes out of range [-1, %d)" % (name, k))
            ax = ds.axes.get(name)            # name clash that suppressed auto-induce (never silent)
            if ax is not None and ax.induced_by != name:
                albls = np.asarray(ax.labels, dtype=str)
                acats = np.asarray(getattr(cat, "categories", []), dtype=str)
                if albls.shape != acats.shape or not bool(np.all(albls == acats)):
                    warn("categorical field '%s' shares a name with axis '%s' (different labels) but "
                         "did not induce it -- a name clash; rename one to give it a factor axis" % (name, name))
        elif enc == "utf8" or (f.role == "label" and len(span) == 1):
            arr = np.asarray(f.values)
            if len(span) == 1 and (arr.ndim != 1 or arr.shape[0] != flen[span[0]]):
                err("field '%s' label shape %s != axis '%s' length %d"
                    % (name, arr.shape, span[0], flen[span[0]]))
        else:  # dense
            arr = np.asarray(f.values)
            if arr.ndim != len(span):
                err("field '%s' dense ndim %d != span arity %d" % (name, arr.ndim, len(span)))
            else:
                for i, ax in enumerate(span):
                    if arr.shape[i] != flen[ax]:
                        err("field '%s' dim %d (%d) != axis '%s' length %d"
                            % (name, i, arr.shape[i], ax, flen[ax]))

        if f.role == "relation" and len(span) != 2:
            err("field '%s' relation must span 2 axes, got %s" % (name, span))

        if getattr(f, "mask", None) is not None:          # nullable validity mask: 1 == missing
            mk = np.asarray(f.mask)
            if mk.dtype.kind not in ("u", "i", "b"):
                err("field '%s' mask dtype %s is not integer/bool" % (name, mk.dtype))
            if len(span) == 1 and mk.shape[0] != flen[span[0]]:
                err("field '%s' mask length %d != axis '%s' length %d"
                    % (name, mk.shape[0], span[0], flen[span[0]]))

        # open vocabularies -> warnings only
        if f.role and f.role not in CORE_ROLES and not f.role.startswith("x-"):
            warn("field '%s' uses non-core role '%s'" % (name, f.role))
        if f.role == "measure" and f.state not in CORE_STATES and not str(f.state).startswith("x-"):
            warn("field '%s' uses non-core state '%s'" % (name, f.state))

    # induced factor axes must stay consistent with their inducing field -- induction is *checkable*,
    # not merely conventional (model.md / induction_design.md §4): the axis labels ARE the field's
    # categories, so any drift between them is a writer bug.
    for name, ax in ds.axes.items():
        if ax.role != "factor" or not ax.induced_by:
            continue
        f = ds.fields.get(ax.induced_by)
        if f is None:
            err("factor axis '%s' induced_by '%s' but no such field" % (name, ax.induced_by))
            continue
        cats = np.asarray(getattr(f.values, "categories", []), dtype=str)
        albls = np.asarray(ax.labels, dtype=str)
        if cats.shape != albls.shape or not bool(np.all(cats == albls)):
            err("factor axis '%s' labels disagree with inducing field '%s' categories "
                "(induced-axis drift)" % (name, ax.induced_by))

    # Profile contract: `viewer@0.1` guarantees a set of precomputed fields with a fixed orientation
    # (docs/format.md "The viewer profile"). A reader relies on the tag, so a stamped-but-incomplete
    # or wrongly-oriented store is an ERROR, not a convention.
    if "viewer@0.1" in (getattr(ds, "profiles", None) or []):
        _check_viewer_profile(ds, err, warn)

    if strict and any(i.startswith("ERROR") for i in issues):
        raise ValueError("L* validation failed:\n  " + "\n  ".join(
            i for i in issues if i.startswith("ERROR")))
    return issues


def _check_viewer_profile(ds, err, warn):
    """Enforce the `viewer@0.1` contract (docs/format.md). Required fields must exist for >=1 grouping,
    with the canonical orientation: stats group-major [<g>, genes], markers gene-major [genes, <g>]."""
    fields = ds.fields
    counts = fields.get("counts")
    gene_axis = counts.span[1] if counts is not None and len(counts.span or []) > 1 else "genes"
    cell_axis = counts.span[0] if counts is not None and (counts.span or []) else "cells"

    cm = fields.get("counts_cellmajor")
    if cm is None:
        err("viewer@0.1: required field 'counts_cellmajor' is missing")
    elif cm.encoding != "csr" or list(cm.span or []) != [cell_axis, gene_axis]:
        err("viewer@0.1: 'counts_cellmajor' must be csr over [%s, %s], got %s over %s"
            % (cell_axis, gene_axis, cm.encoding, list(cm.span or [])))

    od = fields.get("od_score")
    if od is None:
        err("viewer@0.1: required field 'od_score' is missing")
    elif list(od.span or []) != [gene_axis]:
        err("viewer@0.1: 'od_score' must span [%s], got %s" % (gene_axis, list(od.span or [])))

    groupings = sorted(m.group(1) for m in (re.match(r"^stats_(.+)_sum$", n) for n in fields) if m)
    if not groupings:
        err("viewer@0.1: no grouping found (expected stats_<g>_sum / markers_<g>_* for some <g>)")

    for g in groupings:
        for stat in ("sum", "sumsq", "nexpr"):
            nm = "stats_%s_%s" % (g, stat)
            f = fields.get(nm)
            if f is None:
                err("viewer@0.1: required field '%s' is missing" % nm)
            elif len(f.span or []) != 2 or f.span[1] != gene_axis or f.span[0] == gene_axis:
                err("viewer@0.1: '%s' must be group-major [<%s>, %s]; got %s"
                    % (nm, g, gene_axis, list(f.span or [])))
        for m in ("lfc", "padj"):
            nm = "markers_%s_%s" % (g, m)
            f = fields.get(nm)
            if f is None:
                err("viewer@0.1: required field '%s' is missing" % nm)
            elif len(f.span or []) != 2 or f.span[0] != gene_axis or f.span[1] == gene_axis:
                err("viewer@0.1: '%s' must be gene-major [%s, <%s>]; got %s  "
                    "(markers are gene-major; do not transpose)" % (nm, gene_axis, g, list(f.span or [])))

    order = fields.get("counts_cellmajor_order")
    if order is not None and list(order.span or []) != [cell_axis]:
        warn("viewer@0.1: 'counts_cellmajor_order' should span [%s], got %s"
             % (cell_axis, list(order.span or [])))
