"""Differential expression as a factor-axis bundle. A one-vs-rest `uns['rank_genes_groups']` is typed
into measures over `(factor, genes)` (the induced clustering axis), round-trips through L*, and
regenerates `rank_genes_groups` on write-back. Pairwise / reference-group DE stays in the lossless
passthrough (not typed). `markers()` gives the tidy long-form view.

Run: PYTHONPATH=python/src python3 python/tests/test_de.py
"""
import os
import tempfile

import numpy as np

import lstar


def _store():
    return os.path.join(tempfile.mkdtemp(), "de.lstar.zarr")


def _adata_with_de(reference="rest"):
    import anndata as ad
    import pandas as pd
    import scanpy as sc
    import scipy.sparse as sp
    rng = np.random.default_rng(0)
    n, g = 120, 30
    X = sp.csr_matrix(rng.poisson(0.4, size=(n, g)).astype("float32"))
    obs = pd.DataFrame({"leiden": pd.Categorical([str(i % 3) for i in range(n)])},
                       index=[f"cell{i}" for i in range(n)])
    a = ad.AnnData(X=X, obs=obs, var=pd.DataFrame(index=[f"g{j}" for j in range(g)]))
    sc.pp.normalize_total(a, target_sum=1e4); sc.pp.log1p(a)
    sc.tl.rank_genes_groups(a, "leiden", method="t-test", reference=reference)
    return a


def test_de_bundle_roundtrip():
    import warnings; warnings.filterwarnings("ignore")
    a = _adata_with_de("rest")
    rgg = a.uns["rank_genes_groups"]
    groups = list(rgg["names"].dtype.names)

    ds = lstar.read_anndata(a)
    # typed into a (factor, genes) bundle over the induced 'leiden' factor axis
    sc_field = ds.field("de.leiden.score")
    assert sc_field.span == ["leiden", "genes"] and sc_field.subtype == "de"
    assert "de.leiden.lfc" in ds.fields and "de.leiden.padj" in ds.fields
    assert not lstar.validate(ds)
    # rank_genes_groups removed from the passthrough (typed, not double-stored)
    assert "rank_genes_groups" not in ds.aux.get("anndata.uns", {})

    # the scattered scores match scanpy: group's top-ranked gene has its score at the gene-keyed slot
    genes = list(np.asarray(ds.axis("genes").labels))
    forder = list(np.asarray(ds.axis("leiden").labels))
    for grp in groups:
        top_gene = str(rgg["names"][grp][0]); top_score = float(rgg["scores"][grp][0])
        v = np.asarray(sc_field.values)[forder.index(grp), genes.index(top_gene)]
        assert np.isclose(v, top_score, rtol=1e-5), (grp, v, top_score)

    # round-trip through the store, then regenerate rank_genes_groups on write-back
    p = _store(); lstar.write(ds, p)
    a2 = lstar.write_anndata(lstar.read(p))
    rgg2 = a2.uns["rank_genes_groups"]
    assert set(rgg2["names"].dtype.names) == set(groups)
    for grp in groups:
        m1 = dict(zip(rgg["names"][grp].astype(str), rgg["scores"][grp].astype(float)))
        m2 = dict(zip(rgg2["names"][grp].astype(str), rgg2["scores"][grp].astype(float)))
        assert set(m1) == set(m2)
        for gene in m1:
            assert np.isclose(m1[gene], m2[gene], rtol=1e-4), (grp, gene, m1[gene], m2[gene])
    print("DE one-vs-rest: typed to (factor,genes), round-trips, regenerates rank_genes_groups exactly")


def test_markers_tidy_view():
    import warnings; warnings.filterwarnings("ignore")
    a = _adata_with_de("rest")
    ds = lstar.read_anndata(a)
    m = lstar.markers(ds, "leiden", top=5, sort_by="score")
    assert list(m.columns[:2]) == ["group", "gene"]
    assert set(m["group"]) == set(np.asarray(ds.axis("leiden").labels))
    assert (m.groupby("group").size() == 5).all()                      # top-5 per group
    # the top gene per group (by score) matches scanpy's first-ranked gene
    rgg = a.uns["rank_genes_groups"]
    for grp, sub in m.groupby("group"):
        top = sub.sort_values("score", ascending=False).iloc[0]["gene"]
        assert top == str(rgg["names"][str(grp)][0])
    print("markers(): tidy (group, gene, score, lfc, pval, padj) with top-N per group")


def test_pairwise_de_stays_passthrough():
    import warnings; warnings.filterwarnings("ignore")
    a = _adata_with_de(reference="0")                                  # reference group -> NOT one-vs-rest
    ds = lstar.read_anndata(a)
    assert not any(f.subtype == "de" for f in ds.fields.values())      # not typed
    assert "rank_genes_groups" in ds.aux["anndata.uns"]                # kept verbatim in passthrough
    p = _store(); lstar.write(ds, p)
    a2 = lstar.write_anndata(lstar.read(p))
    rgg, rgg2 = a.uns["rank_genes_groups"], a2.uns["rank_genes_groups"]
    assert rgg2["params"]["reference"] == "0"
    g0 = rgg["names"].dtype.names[1]
    assert list(rgg2["names"][g0].astype(str)) == list(rgg["names"][g0].astype(str))  # verbatim
    print("pairwise/reference-group DE survives verbatim in the passthrough (not typed)")


if __name__ == "__main__":
    test_de_bundle_roundtrip()
    test_markers_tidy_view()
    test_pairwise_de_stays_passthrough()
