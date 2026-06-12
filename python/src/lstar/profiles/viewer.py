"""viewer@0.1 profile exporter — precompute the reductions an interactive viewer needs so
common queries touch no matrix (views.md §3). Everything is an ordinary L* field over
ordinary (possibly induced) axes; a non-viewer reader sees plain measures and ignores them.

write_viewer(ds, grouping="leiden") adds, for a single-cell dataset with a `counts` measure
over (cells, genes) and a categorical grouping label over (cells):

  groups_<grouping>          axis (derived)            the group ids
  od_genes                   axis (derived)            the overdispersed-gene subset
  stats_<grouping>_{sum,sumsq,nexpr}   measure (groups, genes)   cluster sufficient stats (log1p)
  markers_<grouping>_{lfc,padj}        measure (genes, groups)   ranked-marker tables
  cell_order                 measure (cells)           a cluster-coherent permutation
  de_panel                   measure (cells, od_genes) CSR, log1p  the cell-major DE panel
                                                                    (subsample DE at O(rows))

The cluster stats mirror libstar's csc_col_sum_by_group; this numpy reference is the exporter
(offline) path. Marks `viewer@0.1` in ds.profiles.
"""
import numpy as np
import scipy.sparse as sp


def write_viewer(ds, grouping="leiden", counts="counts", n_od=300):
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

    # cluster sufficient stats over log1p
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

    # overdispersed-gene subset (top variance) and a cluster-coherent cell order
    gm = grand / ncells
    gvar = np.asarray(Xl.multiply(Xl).sum(0)).ravel() / ncells - gm ** 2
    od_idx = np.sort(np.argsort(-gvar)[: min(n_od, ngenes)])
    order = np.argsort(code, kind="stable").astype("i8")

    # cell-major DE panel: (cells, od_genes), CSR, log1p — read a few hundred rows for subsample DE
    panel = Xlr[:, od_idx].tocsr().astype("f4")

    if ("groups_%s" % grouping) not in ds.axes:
        ds.add_axis("groups_%s" % grouping, groups, origin="derived", role="feature")
    if "od_genes" not in ds.axes:
        ds.add_axis("od_genes", [gene_labels[i] for i in od_idx], origin="derived", role="feature")
    span_gg = ["groups_%s" % grouping, "genes"]
    ds.add_field("stats_%s_sum" % grouping, S.astype("f4"), role="measure", span=span_gg)
    ds.add_field("stats_%s_sumsq" % grouping, SS.astype("f4"), role="measure", span=span_gg)
    ds.add_field("stats_%s_nexpr" % grouping, NE.astype("f4"), role="measure", span=span_gg)
    ds.add_field("markers_%s_lfc" % grouping, lfc, role="measure", span=["genes", "groups_%s" % grouping])
    ds.add_field("markers_%s_padj" % grouping, padj, role="measure", span=["genes", "groups_%s" % grouping])
    ds.add_field("cell_order", order, role="measure", span=["cells"], state="permutation")
    ds.add_field("de_panel", panel, role="measure", span=["cells", "od_genes"], state="lognorm", encoding="csr")
    if "viewer@0.1" not in ds.profiles:
        ds.profiles.append("viewer@0.1")
    return ds
