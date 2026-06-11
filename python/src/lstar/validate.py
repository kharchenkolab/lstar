"""Structural validation of an L* Dataset / store.

Checks the invariants any conforming writer (Python, C++, R) must satisfy: spans reference
existing axes, arities/shapes match axis lengths, relations span two axes, etc. The role and
state vocabularies are *open*, so unrecognized values are warnings, not errors.

    issues = lstar.validate(ds)            # list of "ERROR: ..."/"WARN: ..." strings
    lstar.validate(ds, strict=True)        # raises ValueError on any ERROR
"""
import numpy as np
import scipy.sparse as sp

CORE_ROLES = {"measure", "embedding", "loading", "relation", "label", "sequence",
              "design", "transform", "vector", "probability", "coordinate"}
CORE_STATES = {None, "", "raw", "lognorm", "scaled", "clr"}


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

        enc = f.encoding
        if enc in ("csr", "csc", "coo"):
            if len(span) != 2:
                err("field '%s' (%s) must span 2 axes, got %s" % (name, enc, span))
            elif sp.issparse(f.values):
                shp, exp = tuple(f.values.shape), (axlen[span[0]], axlen[span[1]])
                if shp != exp:
                    err("field '%s' sparse shape %s != axis lengths %s" % (name, shp, exp))
        elif enc == "utf8" or (f.role == "label" and len(span) == 1):
            arr = np.asarray(f.values)
            if len(span) == 1 and (arr.ndim != 1 or arr.shape[0] != axlen[span[0]]):
                err("field '%s' label shape %s != axis '%s' length %d"
                    % (name, arr.shape, span[0], axlen[span[0]]))
        else:  # dense
            arr = np.asarray(f.values)
            if arr.ndim != len(span):
                err("field '%s' dense ndim %d != span arity %d" % (name, arr.ndim, len(span)))
            else:
                for i, ax in enumerate(span):
                    if arr.shape[i] != axlen[ax]:
                        err("field '%s' dim %d (%d) != axis '%s' length %d"
                            % (name, i, arr.shape[i], ax, axlen[ax]))

        if f.role == "relation" and len(span) != 2:
            err("field '%s' relation must span 2 axes, got %s" % (name, span))

        # open vocabularies -> warnings only
        if f.role and f.role not in CORE_ROLES and not f.role.startswith("x-"):
            warn("field '%s' uses non-core role '%s'" % (name, f.role))
        if f.role == "measure" and f.state not in CORE_STATES and not str(f.state).startswith("x-"):
            warn("field '%s' uses non-core state '%s'" % (name, f.state))

    if strict and any(i.startswith("ERROR") for i in issues):
        raise ValueError("L* validation failed:\n  " + "\n  ".join(
            i for i in issues if i.startswith("ERROR")))
    return issues
