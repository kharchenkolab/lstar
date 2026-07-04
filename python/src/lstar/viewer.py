"""Extend an L* dataset with the precomputed fields the **lstar-viewer** web app needs for fast,
latency-cheap browsing of a `.lstar.zarr` store.

A plain converted store (counts + an embedding + categorical labels) is enough to *render* an
embedding, color by a gene or a metadatum, and crosstab labels. But the viewer's heavy interactions
-- differential expression / variable-gene ranking / dotplots, and low-latency reads over a remote
store -- want extra precomputed fields plus a locality-friendly cell row order. :func:`extend_for_viewer`
adds exactly those, natively, so the app loads them instead of recomputing in the browser:

  * ``counts_cellmajor``          -- a CSR (cell-major) copy of ``counts`` (the substrate for on-the-fly
                                     scope compute and per-cell range reads).
  * ``stats_<g>_{sum,sumsq,nexpr}`` -- per-(group, gene) sufficient statistics over ``log1p(counts)``
                                     for each categorical grouping ``g`` (over an induced ``groups_<g>``
                                     axis); DE / dotplots are built from these.
  * ``markers_<g>_{lfc,padj}``    -- a 1-vs-rest marker table per group, derived from the stats.
  * ``od_score``                  -- a global per-gene overdispersion residual (variable-gene ranking).
  * ``counts_cellmajor_order``    -- (when ``order="hybrid"``) each cell's PHYSICAL row in the reordered
                                     ``counts_cellmajor``, after a cluster-contiguous + Hilbert-local
                                     reorder; the reader keys on the ``_order`` sibling field so a
                                     cluster/lasso selection coalesces into a few byte-range reads.

The computations mirror the viewer's JS store-prep (``lstar-viewer/prep/prep.ts`` and
``prep/reorder.mjs``) so a natively-extended store is byte-equivalent to the JS-prepped one on the
fields that matter (the stats match exactly; the marker ``lfc``/``padj`` match to ~1e-3).
"""
import numpy as np
import scipy.sparse as sp

from .kernels import col_sum_by_group, markers_one_vs_rest, overdispersion, cell_order, _N_GRID
from .model import as_categorical, _is_categorical

VIEWER_PROFILE = "viewer@0.1"
N_GRID = _N_GRID                                       # single grid source lives in kernels (was a 2nd copy)


# ---------------------------------------------------------------------------------------------------
# auto-detection of the grouping labels and the embedding to extend over
# ---------------------------------------------------------------------------------------------------

# Canonical grouping-detection policy -- the single source R (.VIEWER_PREFERRED_GROUPINGS) and JS
# (policy.ts) must match; enforced against conformance/viewer_policy.json by conformance/policy_linter.py.
_PREFERRED_GROUPINGS = ("leiden", "cluster", "clusters", "cell_type", "celltype", "cell_types",
                        "louvain", "seurat_clusters", "annotation", "cluster_label")
_MIN_GROUPS, _MAX_GROUPS = 2, 60                       # single-sourced with viewer_policy.json (policy_linter)
_LOGNORM_NAMES = ("X", "data", "logcounts")            # lognorm measure-name fallback (no state=="lognorm" match)
_PREFERRED_EMBEDDINGS = ("umap",)


def _label_codes(ds, name):
    """The (sorted-unique categories, per-cell int codes) for a categorical-or-utf8 cell label.

    Categories are the sorted unique label set -- the same ordering the viewer's JS prep uses
    (``[...new Set(labels)].sort()``) -- so the induced ``groups_<g>`` axis aligns field-for-field
    with a JS-prepped store. Returns ``(groups, codes)`` with ``codes`` an int32 array over cells.
    """
    fl = ds.field(name)
    if _is_categorical(fl.values):
        labels = np.asarray(as_categorical(fl.values), dtype=str)   # decode to strings (missing -> "")
    else:
        labels = np.asarray(fl.values, dtype=str)
    groups, codes = np.unique(labels, return_inverse=True)
    return np.asarray(groups, dtype=str), np.asarray(codes, dtype=np.int32)


def _detect_groupings(ds, min_groups=_MIN_GROUPS, max_groups=_MAX_GROUPS):
    """Categorical cell labels usable as groupings: a ``label``-role field over the cell axis with
    2..~max_groups distinct values. Names that look like a clustering / cell-type annotation are
    preferred (and sorted to the front); near-unique id-like labels are skipped.
    """
    cell_axis = _cell_axis(ds)
    out = []
    for name, fl in ds.fields.items():
        if fl.role != "label" or not fl.span or list(fl.span) != [cell_axis]:
            continue
        if not _is_categorical(fl.values):                 # a string label is fine; numeric/sparse is not
            arr = np.asarray(fl.values)
            if arr.ndim != 1 or arr.dtype.kind not in ("U", "S", "O"):
                continue
        groups, _ = _label_codes(ds, name)
        if min_groups <= len(groups) <= max_groups:
            out.append(name)

    def _rank(nm):
        low = nm.lower()
        for i, p in enumerate(_PREFERRED_GROUPINGS):
            if p in low:
                return (0, i, nm)
        return (1, 0, nm)

    return sorted(out, key=_rank)


def _detect_embedding(ds):
    """The primary embedding field name (``role="embedding"`` over the cell axis), preferring ``umap``.

    A two-dimensional embedding is required for the Hilbert reorder; the first such is used.
    """
    cell_axis = _cell_axis(ds)
    cands = []
    for name, fl in ds.fields.items():
        if fl.role != "embedding" or not fl.span or fl.span[0] != cell_axis:
            continue
        arr = np.asarray(fl.values)
        if arr.ndim == 2 and arr.shape[1] >= 2:
            cands.append(name)
    if not cands:
        return None
    cands.sort(key=lambda nm: (next((i for i, p in enumerate(_PREFERRED_EMBEDDINGS) if p in nm.lower()),
                                    len(_PREFERRED_EMBEDDINGS)), nm))
    return cands[0]


def _cell_axis(ds):
    """The observation (cell) axis name -- the first span axis of ``counts`` (defaults to ``cells``)."""
    if "counts" in ds.fields and ds.field("counts").span:
        return ds.field("counts").span[0]
    return "cells"


def _gene_axis(ds):
    """The feature (gene) axis name -- the second span axis of ``counts`` (defaults to ``genes``)."""
    if "counts" in ds.fields and ds.field("counts").span and len(ds.field("counts").span) > 1:
        return ds.field("counts").span[1]
    return "genes"


def _select_counts_basis(ds, counts=None, basis=None):
    """Choose ``(field_name, apply_log1p)`` the viewer navigators are built from.

    Selected by *content/state*, not a magic field name (a converter that named its raw matrix ``X``
    or a modality is still viewer-preppable). A **raw** measure is preferred and the kernels apply
    ``log1p`` to it. Pass ``counts=<field>`` to force a measure, or ``basis="lognorm"`` to prep
    (approximately) from an already log-normalized measure -- then ``apply_log1p`` is False, so the
    stats are var-of-lognorm rather than var-of-log1p(counts)."""
    twod = [n for n, f in ds.fields.items()
            if f.role == "measure" and f.span and len(f.span) == 2
            and (f.span[0] == "cells" or str(f.span[0]).startswith("cells"))]
    present = ", ".join("%s[%s]" % (n, ds.field(n).state) for n in twod) or "(none)"
    if counts is not None:
        if counts not in ds.fields:
            raise ValueError("extend_for_viewer: counts=%r is not a measure (present cells x genes "
                             "measures: %s)" % (counts, present))
        return counts, (basis != "lognorm" and ds.field(counts).state != "lognorm")
    if basis == "lognorm":
        pick = next((n for n in twod if ds.field(n).state == "lognorm"), None) \
            or next((n for n in twod if n in _LOGNORM_NAMES), None)
        if pick is None:
            raise ValueError("extend_for_viewer: basis='lognorm' but no log-normalized measure "
                             "found (present: %s)" % present)
        return pick, False
    pick = ("counts" if "counts" in twod else None) \
        or next((n for n in twod if ds.field(n).state == "raw"), None)
    if pick is not None:
        return pick, True
    raise ValueError(
        "extend_for_viewer: no raw counts measure found (present cells x genes measures: %s). "
        "Viewer prep needs raw counts; pass counts=<field>, provide a raw-counts measure, or pass "
        "basis='lognorm' to prep (approximately) from a log-normalized measure." % present)


# ---------------------------------------------------------------------------------------------------
# the public entry point
# ---------------------------------------------------------------------------------------------------

def extend_for_viewer(ds, groupings=None, order="hybrid", embedding=None, markers=True,
                      counts=None, basis=None, primary=None):
    """Add the viewer's precomputed fields to ``ds`` in place (and return it).

    Parameters
    ----------
    ds : Dataset
        A sample dataset with a raw ``counts`` measure (CSC, cells x genes), at least one categorical
        cell label, and (for ``order="hybrid"``) a 2-D embedding.
    groupings : list[str] | None
        Categorical label field names to build stats/markers for. ``None`` auto-detects (labels with
        2..~60 distinct values, clustering/cell-type names preferred).
    primary : str | None
        The grouping the *viewer opens on* (its default grouping). Hoisted to the front of the prepared
        groupings, so it is the ``counts_cellmajor`` locality-reorder key AND its stats/markers are
        computed first — the eager-prepare a fast launch waits on. Unlike ordering ``groupings`` by
        hand, ``primary`` COMPOSES with auto-detect: ``primary="cell_type"`` with ``groupings=None``
        preps *every* detected grouping but keys the reorder on ``cell_type``. This is the hook the
        viewer's ``view()`` launcher uses to align the prep with what the first workspace shows (the
        auto-detect policy prefers clusterings, but the viewer may open on a cell-type annotation).
        ``None`` (default) keeps the current behavior: the first detected grouping is primary.
    order : "hybrid" | "none"
        ``"hybrid"`` (default) physically reorders ``counts_cellmajor`` rows by (first grouping's
        cluster code, then a Hilbert index over the embedding) and records each cell's physical row in
        ``counts_cellmajor_order``. ``"none"`` skips the reorder (rows stay in cell order).
    embedding : str | None
        The embedding field used for the Hilbert key; ``None`` auto-detects (prefers ``umap``).
    markers : bool
        Also compute the 1-vs-rest ``markers_<g>_{lfc,padj}`` tables (default ``True``).
    counts : str | None
        Name of the count measure to build from. ``None`` (default) auto-detects by state: a measure
        named ``counts``, else any measure with ``state == "raw"``. Raises a clear error listing the
        present measures if none is found (e.g. a scaled ``X`` + lognorm ``raw`` with no counts).
    basis : {None, "lognorm"}
        ``None`` = raw basis (``log1p``-transformed). ``"lognorm"`` preps -- approximately -- from an
        already log-normalized measure (values used as-is; stats are var-of-lognorm, not
        var-of-log1p(counts)).

    Returns
    -------
    Dataset
        The same dataset, with ``counts_cellmajor``, ``stats_<g>_*``, ``markers_<g>_*``, ``od_score``
        and (for ``order="hybrid"``) ``counts_cellmajor_order`` added.
    """
    # Select the count basis by content/state (not the literal name "counts"): a raw measure is
    # preferred and log1p'd; `counts=`/`basis=` override. Clear error if nothing usable.
    counts_field, use_lognorm = _select_counts_basis(ds, counts=counts, basis=basis)
    basis_state = ds.field(counts_field).state
    cf_span = ds.field(counts_field).span
    cell_axis, gene_axis = cf_span[0], cf_span[1]

    X = ds.field(counts_field).values
    X = X.tocsc() if sp.issparse(X) else sp.csc_matrix(X)
    ncells, ngenes = X.shape

    if groupings is None:
        groupings = _detect_groupings(ds)
    groupings = [g for g in groupings if g in ds.fields]
    # Hoist the viewer's primary grouping to the front (guaranteed present): it becomes the reorder key
    # + is summarized first. Composes with auto-detect above — the rest of the groupings are still prepped.
    if primary is not None:
        if primary not in ds.fields:
            raise ValueError("extend_for_viewer: primary=%r is not a field in the dataset" % primary)
        # must be a 1-D grouping over the CELL axis, else the reorder crashes cryptically (a gene-axis label
        # gives ngenes codes; an embedding is 2-D). span==[cell_axis] is the check identical across Py/R/JS.
        if list(ds.field(primary).span or []) != [cell_axis]:
            raise ValueError("extend_for_viewer: primary=%r must be a grouping over the cell axis %r "
                             "(a 1-D label spanning [%s]), not span=%r" % (primary, cell_axis, cell_axis, ds.field(primary).span))
        groupings = [primary] + [g for g in groupings if g != primary]
    if not groupings:
        raise ValueError("extend_for_viewer: no categorical grouping found (pass groupings=[...])")

    if embedding is None:
        embedding = _detect_embedding(ds)

    # 1) cell-major (CSR) copy of counts -- the substrate for scope compute and per-cell range reads.
    #    Compact raw integer counts to int32; a lognorm-basis measure is float and MUST stay float (an
    #    int32 cast there truncates the values to garbage -- and diverges from R, which keeps float).
    Xr = X.tocsr()
    if use_lognorm:                                    # use_lognorm==True is the RAW basis (log1p applied downstream)
        Xr.data = Xr.data.astype(np.int32, copy=False)
    Xr.indices = Xr.indices.astype(np.int32, copy=False)
    Xr.indptr = Xr.indptr.astype(np.int32, copy=False)

    # 2) global overdispersion residual (a single "group" over all cells -> per-gene mean/var(log1p)).
    od = _od_score(X, ncells, ngenes, lognorm=use_lognorm)
    ds.add_field("od_score", od, role="measure", span=[gene_axis], state=None, encoding="dense",
                 provenance={"cache": VIEWER_PROFILE, "method": "viewer.od",
                             "basis": "log1p" if use_lognorm else "lognorm-input", "trend": "lowess"})

    # 3) per-grouping sufficient stats + (optional) marker tables, each over an induced groups_<g> axis.
    for g in groupings:
        groups, codes = _label_codes(ds, g)
        K = len(groups)
        gaxis = "groups_%s" % g
        if gaxis not in ds.axes:
            ds.add_axis(gaxis, groups, origin="derived", role="feature")
        S, SS, NE = col_sum_by_group(X, codes, K, lognorm=use_lognorm)  # (K, ngenes); log1p iff raw basis
        _cache = {"cache": VIEWER_PROFILE}
        ds.add_field("stats_%s_sum" % g, S, role="measure", span=[gaxis, gene_axis], encoding="dense", provenance=_cache)
        ds.add_field("stats_%s_sumsq" % g, SS, role="measure", span=[gaxis, gene_axis], encoding="dense", provenance=_cache)
        ds.add_field("stats_%s_nexpr" % g, NE, role="measure", span=[gaxis, gene_axis], encoding="dense", provenance=_cache)
        if markers:
            lfc, padj = _markers(S, NE, codes, ncells, K, ngenes)      # dense (ngenes, K) -- genes x groups
            _mk = {"cache": VIEWER_PROFILE, "method": "viewer.markers", "test": "1-vs-rest"}
            ds.add_field("markers_%s_lfc" % g, lfc, role="measure", span=[gene_axis, gaxis], encoding="dense", provenance=_mk)
            ds.add_field("markers_%s_padj" % g, padj, role="measure", span=[gene_axis, gaxis], encoding="dense", provenance=_mk)

    # 4) hybrid cell order: reorder counts_cellmajor rows, record each cell's physical row.
    if order == "hybrid":
        primary = groupings[0]
        _, primary_codes = _label_codes(ds, primary)
        emb = ds.field(embedding).values if (embedding is not None and embedding in ds.fields) else None
        # pos_of[cell] = physical row, from the SHARED core reorder (cluster code, then Hilbert(embedding)
        # when present) -- identical across Python/R/JS. perm is its inverse (physical row -> cell) for
        # the CSR row gather.
        pos_of = cell_order(primary_codes, emb, grid=N_GRID)
        perm = np.empty(ncells, dtype=np.int64)
        perm[pos_of] = np.arange(ncells)
        Xr = _reorder_csr_rows(Xr, perm)
        ds.add_field("counts_cellmajor_order", pos_of.astype(np.float64), role="measure", span=[cell_axis],
                     state="permutation", encoding="dense",
                     provenance={"cache": VIEWER_PROFILE, "method": "viewer.reorder",
                                 "curve": "hilbert" if emb is not None else "cluster",
                                 "grid": N_GRID, "group": primary})

    ds.add_field("counts_cellmajor", Xr, role="measure", span=[cell_axis, gene_axis],
                 state=(basis_state or "raw"), encoding="csr", provenance={"cache": VIEWER_PROFILE})

    if VIEWER_PROFILE not in ds.profiles:
        ds.profiles.append(VIEWER_PROFILE)
    return ds


# ---------------------------------------------------------------------------------------------------
# computations (faithful ports of the viewer's JS store-prep)
# ---------------------------------------------------------------------------------------------------

def _od_score(X, ncells, ngenes, lognorm=True):
    """Per-gene overdispersion score (pagoda2 lowess + F-test) over ``log1p(counts)`` across all
    cells. Delegates to the shared :func:`lstar.kernels.overdispersion` (core C++ / numpy fallback).
    ``lognorm=False`` when the basis is already log-normalized (the values are used as-is)."""
    S, SS, NE = col_sum_by_group(X, np.zeros(ncells, dtype=np.int32), 1, lognorm=lognorm)
    mean = S[0] / ncells
    var = np.maximum(SS[0] / ncells - mean * mean, 0.0)
    nobs = NE[0].astype("i8")                          # expressing cells per gene (the F-test dof)
    return overdispersion(mean, var, nobs)


def _markers(S, NE, codes, ncells, K, ngenes):
    """1-vs-rest marker table from the sufficient stats; returns ``(lfc, padj)`` each dense
    ``(ngenes, K)`` (gene-major). Delegates to the shared :func:`lstar.kernels.markers_one_vs_rest`."""
    nper = np.bincount(codes, minlength=K)
    return markers_one_vs_rest(S, NE, nper, ncells)


def _reorder_csr_rows(Xr, perm):
    """Return a CSR matrix whose physical row p is row ``perm[p]`` of ``Xr`` (an int32 CSR copy)."""
    Xr = Xr.tocsr()
    indptr = Xr.indptr.astype(np.int64, copy=False)
    counts = np.diff(indptr)[perm].astype(np.int64)
    new_indptr = np.empty(len(perm) + 1, dtype=np.int32)
    new_indptr[0] = 0
    new_indptr[1:] = np.cumsum(counts)
    # build the gather index: physical row p contributes source positions [start, start+count) of
    # source row perm[p]. `repeat(start - dest_start, count) + arange(nnz)` lays those slices end-to-end.
    src_starts = indptr[perm]
    nnz = int(counts.sum())
    dest_starts = new_indptr[:-1].astype(np.int64)
    offset = np.repeat(src_starts - dest_starts, counts)
    gather = offset + np.arange(nnz, dtype=np.int64)
    new_data = Xr.data[gather]                          # preserve dtype (int32 raw counts / float lognorm)
    new_indices = Xr.indices[gather].astype(np.int32, copy=False)
    out = sp.csr_matrix((new_data, new_indices, new_indptr), shape=Xr.shape)
    return out
