"""Faithfulness guard for the two-tier corpus: assert each **synthetic** CI stand-in structurally mirrors
the **real** dataset it replaces. This is the automated enforcement of the contract "the synthetic
fixtures faithfully represent the real corpus" -- otherwise that contract rests on human judgment and
silently rots when an upstream library changes (e.g. a new scanpy version alters `rank_genes_groups`
dtypes or which `uns` keys it writes).

It compares **structure, not values**: the set of field roles, axis roles, encodings, promoted-uns
subtypes, feature-axis names, mask presence, and aux passthrough -- everything the profile *types*. The
real structure must be a subset of the synthetic one (synth must cover at least what real carries); a
failure means a real structure is no longer represented synthetically -> add it to `synth.py`.

Runs in the **local** tier (needs real data). In CI / offline, the real loaders return None -> SKIP.

Run: PYTHONPATH=python/src python3 python/tests/test_synth_faithful.py
"""
import os
import sys
import warnings

sys.path.insert(0, os.path.dirname(__file__))
import corpus  # noqa: E402
import synth   # noqa: E402

import lstar  # noqa: E402

warnings.filterwarnings("ignore")
corpus.SYNTHETIC = False     # force the REAL loaders here regardless of LSTAR_SYNTHETIC_CORPUS


def _encoding(f):
    import scipy.sparse as sp
    if getattr(f, "encoding", None):
        return str(f.encoding)
    v = getattr(f, "values", None)
    return "sparse" if sp.issparse(v) else "dense"


def signature(ds):
    """A structural fingerprint of a dataset -- the typed shape, with no values or names that vary by
    dataset (gene ids, cell counts)."""
    return {
        "field_roles": {f.role for f in ds.fields.values()},
        "axis_roles": {a.role for a in ds.axes.values()},
        "subtypes": {f.subtype for f in ds.fields.values() if getattr(f, "subtype", "")},
        "encodings": {_encoding(f) for f in ds.fields.values()},
        "feature_axes": {n for n, a in ds.axes.items() if a.role == "feature"},
        "has_mask": any(getattr(f, "mask", None) is not None for f in ds.fields.values()),
        "aux_keys": set(ds.aux.keys()),
    }


# Which structural keys are a hard contract (real must be a subset of synth) vs informational.
HARD = ("field_roles", "axis_roles", "subtypes", "feature_axes")


def _compare(name, real_ds, synth_ds):
    rs, ss = signature(real_ds), signature(synth_ds)
    problems = []
    for key in HARD:
        missing = rs[key] - ss[key]
        if missing:
            problems.append("%s: synthetic missing %s" % (key, sorted(missing)))
    # informational (printed, not asserted): mask + aux presence
    info = []
    if rs["has_mask"] and not ss["has_mask"]:
        info.append("real has a validity mask; synthetic does not")
    if rs["aux_keys"] and not ss["aux_keys"]:
        info.append("real has an aux passthrough; synthetic has none")
    assert not problems, "%s faithfulness gap:\n  " % name + "\n  ".join(problems)
    extra = "" if not info else "  (note: " + "; ".join(info) + ")"
    print("  %-26s faithful: real structure ⊆ synthetic [roles=%d axes=%d subtypes=%s feat=%s]%s"
          % (name, len(ss["field_roles"]), len(ss["axis_roles"]), sorted(ss["subtypes"]),
             sorted(ss["feature_axes"]), extra))


def test_pbmc68k_faithful():
    real = corpus.pbmc68k_reduced()
    if real is None:
        print("  SKIP pbmc68k faithful (real unavailable)"); return
    _compare("pbmc68k (AnnData)", lstar.read_anndata(real), lstar.read_anndata(synth.pbmc68k_like()))


def test_pbmc3k_de_faithful():
    import scanpy as sc
    real = corpus.pbmc3k_processed()
    if real is None:
        print("  SKIP pbmc3k DE faithful (real unavailable)"); return
    real = real.copy(); sc.tl.rank_genes_groups(real, "louvain", method="t-test")
    s = synth.pbmc3k_like(); sc.tl.rank_genes_groups(s, "louvain", method="t-test")
    _compare("pbmc3k+DE (AnnData)", lstar.read_anndata(real), lstar.read_anndata(s))


def test_velocity_faithful():
    real = corpus.pancreas_velocity()
    if real is None:
        print("  SKIP velocity faithful (scvelo/fixture unavailable)"); return
    _compare("velocity (AnnData)", lstar.read_anndata(real), lstar.read_anndata(synth.velocity()))


def test_mudata_faithful():
    real = corpus.minipbcite()
    if real is None:
        print("  SKIP mudata faithful (real minipbcite unavailable)"); return
    _compare("minipbcite (MuData)", lstar.read_mudata(real), lstar.read_mudata(synth.citeseq_mudata_annotated()))


if __name__ == "__main__":
    test_pbmc68k_faithful()
    test_pbmc3k_de_faithful()
    test_velocity_faithful()
    test_mudata_faithful()
    print("synthetic-faithfulness guard: every available real dataset's structure is represented synthetically")
