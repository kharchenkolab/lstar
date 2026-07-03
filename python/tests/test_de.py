"""Differential expression as a factor-axis bundle -- grounded in **real** scanpy `rank_genes_groups`
output (pbmc3k/pbmc68k, real pipeline, real data), not fabricated structures. A one-vs-rest result is
typed into measures over `(factor, genes)`, round-trips through L*, and regenerates `rank_genes_groups`
on write-back; pairwise / reference-group DE stays verbatim in the passthrough; `markers()` is the tidy
view. Covers the real-world variants: t-test (full names/scores/lfc/pvals/pvals_adj) and logreg
(names/scores only).

Run: PYTHONPATH=python/src python3 python/tests/test_de.py
"""
import os
import sys
import tempfile
import warnings

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
import corpus  # noqa: E402

import lstar  # noqa: E402

warnings.filterwarnings("ignore")


def _store():
    return os.path.join(tempfile.mkdtemp(), "de.lstar.zarr")


def test_de_bundle_roundtrip_ttest():
    a = corpus.pbmc3k_with_de("t-test")                      # REAL one-vs-rest DE, full statistics
    if a is None:
        print("  SKIP test_de_bundle_roundtrip_ttest (corpus unavailable)"); return
    rgg = a.uns["rank_genes_groups"]
    groups = list(rgg["names"].dtype.names)

    ds = lstar.read_anndata(a)
    sf = ds.field("de.louvain.score")
    gene_axis = sf.span[1]                                   # genes or genes_raw, per use_raw
    assert sf.span[0] == "louvain" and sf.subtype == "de"
    for st in ("lfc", "pval", "padj"):                       # t-test yields the full bundle
        assert "de.louvain.%s" % st in ds.fields, st
    assert "rank_genes_groups" not in ds.aux.get("anndata.uns", {})
    assert not lstar.validate(ds)

    genes = list(np.asarray(ds.axis(gene_axis).labels))
    forder = list(np.asarray(ds.axis("louvain").labels))
    for grp in groups:                                       # scattered scores match scanpy exactly
        top_gene = str(rgg["names"][grp][0]); top_score = float(rgg["scores"][grp][0])
        v = np.asarray(sf.values)[forder.index(grp), genes.index(top_gene)]
        assert np.isclose(v, top_score, rtol=1e-5), (grp, v, top_score)

    p = _store(); lstar.write(ds, p)
    a2 = lstar.write_anndata(lstar.read(p))                  # regenerate rank_genes_groups
    rgg2 = a2.uns["rank_genes_groups"]
    assert set(rgg2["names"].dtype.names) == set(groups)
    for grp in groups:
        m1 = dict(zip(rgg["names"][grp].astype(str), rgg["scores"][grp].astype(float)))
        m2 = dict(zip(rgg2["names"][grp].astype(str), rgg2["scores"][grp].astype(float)))
        assert set(m1) == set(m2)
        for gene in list(m1)[:50]:
            assert np.isclose(m1[gene], m2[gene], rtol=1e-4), (grp, gene)
    print("DE (real t-test): full bundle typed over (louvain,%s), round-trips, regenerates exactly" % gene_axis)


def test_de_logreg_scoreonly_variant():
    a = corpus.pbmc68k_reduced()                             # REAL DE: method=logreg -> names+scores only
    if a is None:
        print("  SKIP test_de_logreg_scoreonly_variant (corpus unavailable)"); return
    assert str(a.uns["rank_genes_groups"]["params"]["method"]) == "logreg"
    ds = lstar.read_anndata(a)
    assert "de.bulk_labels.score" in ds.fields                # the only stat logreg provides
    assert "de.bulk_labels.lfc" not in ds.fields              # no lfc/pval for logreg -> not invented
    assert not lstar.validate(ds)
    a2 = lstar.write_anndata(lstar.read(_w(ds)))
    assert set(a2.uns["rank_genes_groups"]["names"].dtype.names) == \
        set(a.uns["rank_genes_groups"]["names"].dtype.names)
    print("DE (real logreg): score-only variant typed faithfully (no invented lfc/pvals); regenerates")


def test_markers_tidy_view():
    a = corpus.pbmc3k_with_de("t-test")
    if a is None:
        print("  SKIP test_markers_tidy_view (corpus unavailable)"); return
    ds = lstar.read_anndata(a)
    m = lstar.markers(ds, "louvain", top=5, sort_by="score")
    assert list(m.columns[:2]) == ["group", "gene"]
    assert set(m["group"]) == set(np.asarray(ds.axis("louvain").labels))
    assert (m.groupby("group").size() == 5).all()
    rgg = a.uns["rank_genes_groups"]
    for grp, sub in m.groupby("group"):
        top = sub.sort_values("score", ascending=False).iloc[0]["gene"]
        assert top == str(rgg["names"][str(grp)][0])
    print("markers(): tidy (group, gene, score, lfc, pval, padj) top-N matches scanpy's ranking")


def test_pseudobulk_bundle():
    a = corpus.pbmc68k_reduced()                             # real expression + real louvain grouping
    if a is None:
        print("  SKIP test_pseudobulk_bundle (corpus unavailable)"); return
    import scipy.sparse as sp
    src = a.raw if a.raw is not None else a            # real COUNTS live in .raw (X is scaled)
    ds = lstar.Dataset(kind="sample")
    genes = np.asarray(src.var_names, dtype=str)
    ds.add_axis("cells", np.asarray(a.obs_names, dtype=str)); ds.add_axis("genes", genes)
    X = src.X; X = X.tocsr() if sp.issparse(X) else sp.csr_matrix(np.asarray(X))
    ds.add_field("counts", X, role="measure", span=["cells", "genes"], state="raw")
    codes = np.asarray(a.obs["louvain"].cat.codes)
    ds.add_field("louvain", lstar.Categorical(codes, np.asarray(a.obs["louvain"].cat.categories, dtype=str)),
                 span=["cells"])

    pb = lstar.pseudobulk(ds, "louvain", field="counts")
    assert ds.field("pb.louvain.mean").span == ["louvain", "genes"]
    rows0 = np.nonzero(codes == 0)[0]
    assert np.allclose(pb["mean"][0], np.asarray(X[rows0].mean(axis=0)).ravel(), rtol=1e-4)  # f32 data
    assert np.all((pb["frac"] >= 0) & (pb["frac"] <= 1))
    assert not lstar.validate(ds)
    assert np.allclose(np.asarray(lstar.read(_w(ds)).field("pb.louvain.frac").values), pb["frac"], rtol=1e-4)
    print("pseudobulk (real pbmc68k counts): (factor,genes) mean+frac match manual reduction; round-trips")


def test_pseudobulk_kernel_reduction():
    """pseudobulk routes its per-(group,gene) reduction through the SHARED col_sum_by_group kernel (was a
    numpy per-group loop that duplicated it). Its mean/frac must equal a direct reduction, both bases."""
    import scipy.sparse as sp
    rng = np.random.default_rng(0); nc, ng, K = 150, 30, 4
    X = sp.random(nc, ng, density=0.3, format="csr", random_state=0)
    X.data = (rng.poisson(3, X.data.shape) + 1).astype(np.float64)
    codes = rng.integers(0, K, size=nc)
    ds = lstar.Dataset(kind="sample")
    ds.add_axis("cells", [f"c{i}" for i in range(nc)]); ds.add_axis("genes", [f"g{j}" for j in range(ng)])
    ds.add_field("counts", X, role="measure", span=["cells", "genes"], state="raw")
    ds.add_field("grp", lstar.Categorical(codes.astype(np.int32), np.array([str(k) for k in range(K)])), span=["cells"])
    for lognorm in (False, True):
        pb = lstar.pseudobulk(ds, "grp", field="counts", lognorm=lognorm, add=False)
        for k in range(K):
            rows = np.nonzero(codes == k)[0]
            sub = X[rows]; vals = sub.copy()
            if lognorm:
                vals.data = np.log1p(vals.data)
            assert np.allclose(pb["mean"][k], np.asarray(vals.sum(0)).ravel() / len(rows))
            assert np.allclose(pb["frac"][k], np.asarray((sub > 0).sum(0)).ravel() / len(rows))
    print("pseudobulk kernel-reduction: mean/frac == direct reduction (raw + lognorm)")


def test_pairwise_de_stays_passthrough():
    a = corpus.pbmc3k_processed()
    if a is None:
        print("  SKIP test_pairwise_de_stays_passthrough (corpus unavailable)"); return
    import scanpy as sc
    g0 = str(a.obs["louvain"].cat.categories[0])
    a = a.copy(); sc.tl.rank_genes_groups(a, "louvain", method="t-test", reference=g0)  # REAL pairwise
    ds = lstar.read_anndata(a)
    assert not any(f.subtype == "de" for f in ds.fields.values())     # reference-group -> not typed
    assert "rank_genes_groups" in ds.aux["anndata.uns"]               # kept verbatim in passthrough
    a2 = lstar.write_anndata(lstar.read(_w(ds)))
    assert a2.uns["rank_genes_groups"]["params"]["reference"] == g0
    grp = a.uns["rank_genes_groups"]["names"].dtype.names[1]
    assert list(a2.uns["rank_genes_groups"]["names"][grp].astype(str)) == \
        list(a.uns["rank_genes_groups"]["names"][grp].astype(str))
    print("pairwise/reference-group DE (real): survives verbatim in the passthrough (not typed)")


def _w(ds):
    p = _store(); lstar.write(ds, p); return p


if __name__ == "__main__":
    test_de_bundle_roundtrip_ttest()
    test_de_logreg_scoreonly_variant()
    test_markers_tidy_view()
    test_pseudobulk_bundle()
    test_pseudobulk_kernel_reduction()
    test_pairwise_de_stays_passthrough()
