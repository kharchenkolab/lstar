"""viewer@0.1 profile exporter — precompute the reductions an interactive viewer needs so
common queries touch no matrix (views.md §3). Everything is an ordinary L* field over
ordinary (possibly induced) axes; a non-viewer reader sees plain measures and ignores them.

write_viewer(ds, grouping="leiden") adds, for a single-cell dataset with a `counts` measure
over (cells, genes) and a categorical grouping label over (cells):

  groups_<grouping>          axis (derived)            the group ids
  stats_<grouping>_{sum,sumsq,nexpr}   measure (groups, genes)   cluster sufficient stats (log1p)
  markers_<grouping>_{lfc,padj}        measure (genes, groups)   ranked-marker tables
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


def write_viewer(ds, grouping="leiden", counts="counts", n_od=300, engine="auto"):
    X = ds.field(counts).values
    X = sp.csc_matrix(X) if not sp.issparse(X) else X.tocsc()
    ncells, ngenes = X.shape
    gene_labels = list(ds.axis("genes").labels)

    labels = np.asarray(ds.field(grouping).values).astype(str)
    groups = sorted(set(labels.tolist()))
    gidx = {g: i for i, g in enumerate(groups)}
    code = np.array([gidx[l] for l in labels])
    K = len(groups)

    Xl = X.copy().astype("f8"); Xl.data = np.log1p(Xl.data)
    Xlr = Xl.tocsr()

    # cluster sufficient stats over log1p — the shared libstar kernel (csc_col_sum_by_group, the
    # same C++ bound to WASM/R) when the accelerator is built, else an identical numpy reference.
    from .._engine import resolve_engine, _accel
    use_cpp = resolve_engine(engine) == "c++" and hasattr(_accel, "col_sum_by_group")
    if use_cpp:
        S, SS, NE = _accel.col_sum_by_group(X.data, X.indptr, X.indices, ncells, ngenes,
                                            code.astype("int32"), K, True, 0)
    else:
        S = np.zeros((K, ngenes)); SS = np.zeros((K, ngenes)); NE = np.zeros((K, ngenes))
        for g in range(K):
            sub = Xlr[code == g]
            S[g] = np.asarray(sub.sum(0)).ravel()
            SS[g] = np.asarray(sub.multiply(sub).sum(0)).ravel()
            NE[g] = np.asarray((sub > 0).sum(0)).ravel()
    nper = np.array([(code == g).sum() for g in range(K)])
    grand = np.asarray(Xl.sum(0)).ravel()

    # marker tables: lfc = group mean(log1p) - rest mean(log1p); padj a monotone proxy
    lfc = np.zeros((ngenes, K), "f4"); padj = np.ones((ngenes, K), "f4")
    for g in range(K):
        mu = S[g] / max(nper[g], 1); mr = (grand - S[g]) / max(ncells - nper[g], 1)
        lfc[:, g] = (mu - mr).astype("f4")
        padj[:, g] = np.clip(np.exp(-np.abs((mu - mr) * np.sqrt(NE[g] + 1))), 1e-12, 1).astype("f4")

    order = np.argsort(code, kind="stable").astype("i8")          # cluster-coherent cell permutation

    # counts in cell-major (CSR) orientation, all genes — the substrate for on-the-fly, scope-correct
    # compute. Selection DE and overdispersion subsample CELLS and read their rows over ALL genes;
    # the gene scope is never precomputed, because a global od subset is wrong for any local question
    # (a T-cell-only DE needs genes that don't move global variance). `n_od` is retained but unused.
    panel = X.tocsr()

    # Recompute and OVERWRITE the whole profile from the current (counts, grouping); add_axis/
    # add_field overwrite by key, so re-running is idempotent and a same-named field can't go stale.
    ds.add_axis("groups_%s" % grouping, groups, origin="derived", role="feature")
    span_gg = ["groups_%s" % grouping, "genes"]
    ds.add_field("stats_%s_sum" % grouping, S.astype("f4"), role="measure", span=span_gg)
    ds.add_field("stats_%s_sumsq" % grouping, SS.astype("f4"), role="measure", span=span_gg)
    ds.add_field("stats_%s_nexpr" % grouping, NE.astype("f4"), role="measure", span=span_gg)
    ds.add_field("markers_%s_lfc" % grouping, lfc, role="measure", span=["genes", "groups_%s" % grouping])
    ds.add_field("markers_%s_padj" % grouping, padj, role="measure", span=["genes", "groups_%s" % grouping])
    ds.add_field("cell_order", order, role="measure", span=["cells"], state="permutation")
    ds.add_field("counts_cellmajor", panel, role="measure", span=["cells", "genes"], state="raw", encoding="csr")
    if "viewer@0.1" not in ds.profiles:
        ds.profiles.append("viewer@0.1")
    return ds
