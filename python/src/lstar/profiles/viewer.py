"""viewer@0.1 profile exporter — precompute the reductions an interactive viewer needs so
common queries touch no matrix (views.md §3). Everything is an ordinary L* field over
ordinary (possibly induced) axes; a non-viewer reader sees plain measures and ignores them.

write_viewer(ds, grouping="leiden", also=("cell_type",)) adds, for a single-cell dataset with a
`counts` measure over (cells, genes) and one or more categorical labels over (cells):

  groups_<g>                 axis (derived)            the group ids, for each g in (grouping, *also)
  stats_<g>_{sum,sumsq,nexpr}   measure (groups, genes)   cluster sufficient stats (log1p), per g
  markers_<g>_{lfc,padj}        measure (genes, groups)   ranked-marker tables, per g
  od_score                   measure (genes)           overdispersion vs the smoothed mean-var trend
  cell_order                 measure (cells)           a cluster-coherent permutation
  counts_cellmajor           measure (cells, genes)    CSR, raw   counts in cell-major orientation

What is and isn't precomputed: the cluster stats/markers answer a *global* question (markers of a
fixed top-level cluster over the whole dataset) and are precomputed. Anything *scope-dependent* is
NOT: selection DE and overdispersed-gene (HVG) selection are computed on the fly by subsampling the
cells in the current scope and reducing over ALL genes (col_sum_by_group / col_mean_var) — because a
globally-chosen gene subset is wrong for any local question (DE between two T-cell subsets needs
genes that barely move global variance). `counts_cellmajor` is just counts in cell-major (CSR)
orientation — the read-optimized substrate those row-wise reductions need.

The cluster stats are computed by libstar's csc_col_sum_by_group via the compiled accelerator
when present (the same C++ core bound to WASM/R), with an identical numpy fallback. Marks
`viewer@0.1` in ds.profiles.
"""
import numpy as np
import scipy.sparse as sp


def write_viewer(ds, grouping="leiden", counts="counts", engine="auto", also=(), n_od=None):
    X = ds.field(counts).values
    X = sp.csc_matrix(X) if not sp.issparse(X) else X.tocsc()
    ncells, ngenes = X.shape
    Xl = X.copy().astype("f8"); Xl.data = np.log1p(Xl.data)
    Xlr = Xl.tocsr()
    grand = np.asarray(Xl.sum(0)).ravel()

    from .._engine import resolve_engine, _accel
    use_cpp = resolve_engine(engine) == "c++" and hasattr(_accel, "col_sum_by_group")

    def markers_for(labels):
        labels = np.asarray(labels).astype(str)
        groups = sorted(set(labels.tolist()))
        code = np.array([groups.index(l) for l in labels]); K = len(groups)
        # cluster sufficient stats over log1p — the shared libstar kernel when built, numpy else.
        if use_cpp:
            S, SS, NE = _accel.col_sum_by_group(X.data, X.indptr, X.indices, ncells, ngenes,
                                                code.astype("int32"), K, True, 0)
        else:
            S = np.zeros((K, ngenes)); SS = np.zeros((K, ngenes)); NE = np.zeros((K, ngenes))
            for g in range(K):
                sub = Xlr[code == g]
                S[g] = np.asarray(sub.sum(0)).ravel(); SS[g] = np.asarray(sub.multiply(sub).sum(0)).ravel()
                NE[g] = np.asarray((sub > 0).sum(0)).ravel()
        nper = np.array([(code == g).sum() for g in range(K)])
        lfc = np.zeros((ngenes, K), "f4"); padj = np.ones((ngenes, K), "f4")
        for g in range(K):
            mu = S[g] / max(nper[g], 1); mr = (grand - S[g]) / max(ncells - nper[g], 1)
            lfc[:, g] = (mu - mr).astype("f4")
            padj[:, g] = np.clip(np.exp(-np.abs((mu - mr) * np.sqrt(NE[g] + 1))), 1e-12, 1).astype("f4")
        return groups, code, S, SS, NE, lfc, padj

    # substrate: counts in cell-major (CSR) orientation, all genes — selection DE and subset
    # overdispersion subsample CELLS and reduce over ALL genes (gene scope is never precomputed).
    ds.add_field("counts_cellmajor", X.tocsr(), role="measure", span=["cells", "genes"], state="raw", encoding="csr")

    # whole-dataset overdispersion navigator: residual above the smoothed log(v) ~ log(m) trend.
    gm = grand / ncells
    gv = np.maximum(np.asarray(Xl.multiply(Xl).sum(0)).ravel() / ncells - gm ** 2, 0)
    od = np.zeros(ngenes, "f4"); ok = (gm > 0) & (gv > 0) & np.isfinite(gm) & np.isfinite(gv)
    if ok.sum() > 10:
        coef = np.polyfit(np.log(gm[ok]), np.log(gv[ok]), 2)               # smooth log-log trend
        od[ok] = (np.log(gv[ok]) - np.polyval(coef, np.log(gm[ok]))).astype("f4")
    ds.add_field("od_score", od, role="measure", span=["genes"])

    # per-annotation cluster stats + marker navigators (one set each). Recompute & OVERWRITE by key,
    # so re-running is idempotent and a same-named field can't go stale relative to its inputs.
    primary_code = None
    for gp in dict.fromkeys([grouping, *also]):
        if gp not in ds.fields:
            continue
        groups, code, S, SS, NE, lfc, padj = markers_for(ds.field(gp).values)
        if gp == grouping:
            primary_code = code
        ds.add_axis("groups_%s" % gp, groups, origin="derived", role="feature")
        sg = ["groups_%s" % gp, "genes"]
        ds.add_field("stats_%s_sum" % gp, S.astype("f4"), role="measure", span=sg)
        ds.add_field("stats_%s_sumsq" % gp, SS.astype("f4"), role="measure", span=sg)
        ds.add_field("stats_%s_nexpr" % gp, NE.astype("f4"), role="measure", span=sg)
        ds.add_field("markers_%s_lfc" % gp, lfc, role="measure", span=["genes", "groups_%s" % gp])
        ds.add_field("markers_%s_padj" % gp, padj, role="measure", span=["genes", "groups_%s" % gp])

    order = np.argsort(primary_code, kind="stable").astype("i8")
    ds.add_field("cell_order", order, role="measure", span=["cells"], state="permutation")
    if "viewer@0.1" not in ds.profiles:
        ds.profiles.append("viewer@0.1")
    return ds
