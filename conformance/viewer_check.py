"""Conformance driver for the `viewer@0.1` profile (docs/format.md "The viewer profile").

Subcommands:
  canonical                 build a synthetic sample, extend_for_viewer, assert validate() clean,
                            print the produced viewer-field shapes/spans.
  validate  <store>         read a store and assert it satisfies the viewer@0.1 contract (exit 1 on
                            any viewer ERROR) -- used to check what *other* surfaces wrote.
  equiv     <a> <b>         assert two viewer-extended stores agree on ALL viewer fields, for every
                            grouping. Tolerances (generous for f4-vs-f8 across surfaces, tight enough to
                            catch a method/orientation drift): stats sum/sumsq rtol 1e-4, nexpr 1e-5,
                            markers lfc 2e-3, od_score 5e-3; counts_cellmajor[_order] exact; spans/shapes
                            identical. The looser od_score reflects the naive-vs-stable variance across
                            surfaces (see audit); tighten once od variance is single-sourced.

Exit code 0 == pass. Used by conformance/viewer.sh; importable pieces kept tiny on purpose.
"""
import sys

import numpy as np
import scipy.sparse as sp

import lstar


def _synthetic(nc=200, ng=50, K=5, seed=0, fmt="csc", multi=False):
    rng = np.random.default_rng(seed)
    X = sp.random(nc, ng, density=0.2, format="csc", random_state=seed)
    X.data = (rng.poisson(3, size=X.data.shape) + 1).astype(np.int32)
    X = X.tocsr() if fmt == "csr" else X.tocsc()           # counts may arrive in either encoding (A1)
    ds = lstar.Dataset(kind="sample")
    ds.add_axis("cells", [f"c{i}" for i in range(nc)])
    ds.add_axis("genes", [f"g{j}" for j in range(ng)])
    ds.add_axis("umap", ["umap1", "umap2"])
    ds.add_field("counts", X, role="measure", span=["cells", "genes"], state="raw")
    ds.add_field("umap", rng.normal(size=(nc, 2)), role="embedding", span=["cells", "umap"])
    ds.add_field("leiden", rng.integers(0, K, size=nc).astype(str), role="label", span=["cells"])
    if multi:                                              # competing groupings -> exercises detection (#4)
        ds.add_field("louvain", rng.integers(0, 4, size=nc).astype(str), role="label", span=["cells"])
        ds.add_field("annotation", rng.integers(0, 6, size=nc).astype(str), role="label", span=["cells"])
        ds.add_field("phase", np.array(["G1", "S", "G2M"])[rng.integers(0, 3, size=nc)], role="label", span=["cells"])
        # a GENE-axis label (like AnnData's highly_variable) -- must NOT be picked as a CELL grouping (T1.3):
        # every surface must exclude it (JS previously accepted any 1-D label -> ngenes codes -> crash).
        ds.add_field("highly_variable", np.array(["True", "False"])[rng.integers(0, 2, size=ng)], role="label", span=["genes"])
    return ds


def _viewer_errors(ds):
    return [i for i in lstar.validate(ds) if i.startswith("ERROR") and "viewer@0.1" in i]


def _fail(msg):
    print("  FAIL: " + msg)
    sys.exit(1)


def cmd_make_base(out, fmt="csc", multi=""):
    """Write a bare (un-prepped) synthetic store — the shared input both preps extend. `fmt` sets the
    counts encoding (csc|csr) to exercise A1 normalization; `multi`="multi" adds competing groupings
    (louvain/annotation/phase) to exercise grouping detection (#4) across surfaces."""
    lstar.write(_synthetic(fmt=fmt, multi=(multi == "multi")), out)
    print(f"  OK: wrote base store {out} (counts={fmt}{', multi-grouping' if multi == 'multi' else ''})")


def cmd_prep_lstar(base, out, basis=None):
    """lstar's own prep: read the base store, extend_for_viewer, write it out. `basis` ("lognorm" or
    None=raw) forwards to extend_for_viewer so corpus data with no raw counts preps from a log measure."""
    ds = lstar.read(base)
    lstar.extend_for_viewer(ds, basis=basis)
    lstar.write(ds, out)
    print(f"  OK: lstar extend_for_viewer -> {out}" + (f" (basis={basis})" if basis else ""))


def cmd_canonical():
    ds = _synthetic()
    lstar.extend_for_viewer(ds)
    errs = _viewer_errors(ds)
    if errs:
        _fail("canonical extend_for_viewer is not viewer@0.1-clean:\n    " + "\n    ".join(errs))
    for nm in sorted(ds.fields):
        if nm.startswith(("stats_", "markers_", "counts_cellmajor", "od_score")):
            f = ds.field(nm)
            shp = tuple(np.asarray(f.values).shape) if not sp.issparse(f.values) else tuple(f.values.shape)
            print(f"  [py ] {nm:32s} {f.encoding:5s} {list(f.span or [])} {shp}")
    print("  OK: canonical Python prep satisfies viewer@0.1")


def cmd_validate(store):
    ds = lstar.read(store)
    if "viewer@0.1" not in (ds.profiles or []):
        _fail(f"{store} does not declare the viewer@0.1 profile (profiles={list(ds.profiles)})")
    errs = _viewer_errors(ds)
    if errs:
        _fail(f"{store} violates the viewer@0.1 contract:\n    " + "\n    ".join(errs))
    print(f"  OK: {store} satisfies viewer@0.1")


def _groupings_of(ds):
    return sorted(nm[len("stats_"):-len("_sum")] for nm in ds.fields
                  if nm.startswith("stats_") and nm.endswith("_sum"))


def _cmp_dense(da, db, nm, rtol, atol):
    va, vb = np.asarray(da.field(nm).values), np.asarray(db.field(nm).values)
    sa, sb = list(da.field(nm).span or []), list(db.field(nm).span or [])
    if sa != sb:
        _fail(f"{nm}: span differs {sa} != {sb} (orientation mismatch)")
    if va.shape != vb.shape:
        _fail(f"{nm}: shape differs {va.shape} != {vb.shape}")
    if not np.allclose(va, vb, rtol=rtol, atol=atol):
        _fail(f"{nm}: values differ (max abs {float(np.nanmax(np.abs(va - vb))):.3g} > rtol {rtol})")
    print(f"  OK: {nm} agrees ({sa} {va.shape})")


def cmd_equiv(a, b):
    da, db = lstar.read(a), lstar.read(b)
    # ALL groupings must match (not just the first) -- corpus data carries several (louvain/phase/...);
    # a divergent grouping-detection would surface here as a differing set or a permuted primary.
    grpa, grpb = _groupings_of(da), _groupings_of(db)
    if not grpa:
        _fail("no stats_<g>_sum field found")
    if grpa != grpb:
        _fail(f"grouping SET differs: {grpa} != {grpb} (detection diverged)")
    # tolerances are f4-vs-f8 generous (pagoda3 stores float32, lstar float64) but tight enough to
    # catch any method/orientation drift (a wrong od method or transpose differs by orders of magnitude).
    for g in grpa:
        for nm, rtol, atol in [(f"stats_{g}_sum", 1e-4, 1e-2), (f"stats_{g}_sumsq", 1e-4, 1e-2),
                               (f"stats_{g}_nexpr", 1e-5, 1e-3), (f"markers_{g}_lfc", 2e-3, 1e-2)]:
            if nm in da.fields and nm in db.fields:
                _cmp_dense(da, db, nm, rtol, atol)
    _cmp_dense(da, db, "od_score", 5e-3, 5e-2)

    # physical-layout equivalence -- the fields the old equiv OMITTED, which let Python's cluster+Hilbert
    # reorder and R's cluster-only reorder (and JS's identity stub) all pass as "equivalent". The target
    # contract is a byte-identical store across surfaces, so these must match exactly.
    oa = np.rint(np.asarray(da.field("counts_cellmajor_order").values)).astype("i8")
    ob = np.rint(np.asarray(db.field("counts_cellmajor_order").values)).astype("i8")
    if oa.shape != ob.shape:
        _fail(f"counts_cellmajor_order: shape differs {oa.shape} != {ob.shape}")
    if not np.array_equal(oa, ob):
        ndiff = int(np.sum(oa != ob))
        _fail(f"counts_cellmajor_order: physical row order differs ({ndiff}/{oa.size} cells) -- the "
              f"cell reorder is not identical across surfaces")
    print(f"  OK: counts_cellmajor_order matches ({oa.size} cells)")

    ca, cb = da.field("counts_cellmajor").values, db.field("counts_cellmajor").values
    A = ca.toarray() if sp.issparse(ca) else np.asarray(ca)
    B = cb.toarray() if sp.issparse(cb) else np.asarray(cb)
    if A.shape != B.shape:
        _fail(f"counts_cellmajor: shape differs {A.shape} != {B.shape}")
    if not np.array_equal(A, B):
        _fail(f"counts_cellmajor: physical cell-major payload differs (max abs "
              f"{float(np.max(np.abs(A.astype('f8') - B.astype('f8')))):.3g})")
    print(f"  OK: counts_cellmajor matches ({A.shape})")

    print(f"  OK: {a} and {b} agree on the viewer fields ({len(grpa)} grouping(s): {', '.join(grpa)})")


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else ""
    if cmd == "canonical":
        cmd_canonical()
    elif cmd == "validate":
        cmd_validate(sys.argv[2])
    elif cmd == "make-base":
        cmd_make_base(sys.argv[2], sys.argv[3] if len(sys.argv) > 3 else "csc",
                      sys.argv[4] if len(sys.argv) > 4 else "")
    elif cmd == "prep-lstar":
        cmd_prep_lstar(sys.argv[2], sys.argv[3], sys.argv[4] if len(sys.argv) > 4 else None)
    elif cmd == "equiv":
        cmd_equiv(sys.argv[2], sys.argv[3])
    else:
        print(__doc__)
        sys.exit(2)
