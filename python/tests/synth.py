"""Synthetic-but-faithful fixtures for CI -- the synthetic half of the two-tier corpus.

Philosophy (per project decision): **CI runs synthetic-only** (no real datasets committed to or
downloaded by github -- they are large and slow); the **extensive real corpus runs locally** (see
`corpus.py` + `conformance/sweep/`) and is what proves the profiles against real-world structure. The
job of this module is to make the CI fixtures *structurally faithful* to that real corpus, so a passing
CI run is meaningful.

The faithfulness trick: we don't hand-fabricate the tricky structures (categoricals, `*_colors`,
`uns['pca']`, the `neighbors` OverloadedDict, `rank_genes_groups` structured arrays, RNA+ADT modalities).
We generate **synthetic counts with real cluster/marker structure** and then run the **real scanpy /
mudata pipeline** over them -- so the structures are produced by the same library code that produced the
real objects, just on synthetic input. That is exactly "synthetically generated examples that properly
represent the corpus": same classes, same uns/obsm/layers shapes, same dtypes -- deterministic, offline,
and dependency-light (clusters come from the synthetic ground truth, so no louvain/leiden/umap packages
are needed in CI).

Used by `corpus.py` when `LSTAR_SYNTHETIC_CORPUS=1` (set in CI). Locally, unset, the real loaders run.
"""
import warnings

import numpy as np


def _palette(n):
    """The exact hex palette scanpy assigns to `uns['<key>_colors']` (so colors are byte-real)."""
    try:
        from scanpy.plotting.palettes import default_20, default_28
        pal = list(default_20) if n <= 20 else list(default_28)
    except Exception:                                   # palettes always ship with scanpy, but be safe
        pal = ["#%02x%02x%02x" % (37 * i % 256, 89 * i % 256, 151 * i % 256) for i in range(max(n, 1))]
    if n > len(pal):
        pal = (pal * (n // len(pal) + 1))
    return np.array(pal[:n], dtype=object).astype(str)


def _counts(n_obs, n_genes, n_clusters, seed, marker_fold=6.0):
    """A synthetic count matrix with real structure: each cluster up-regulates its own marker-gene block
    (so PCA separates the clusters, the kNN graph is meaningful, and DE finds distinct markers) over a
    gamma-distributed baseline with per-cell depth variation. Returns (counts f32 [obs x genes], labels)."""
    rng = np.random.default_rng(seed)
    labels = rng.integers(0, n_clusters, n_obs)
    for k in range(n_clusters):                         # guarantee every cluster is non-empty
        labels[k] = k
    gene_base = rng.gamma(0.4, 1.5, n_genes) + 0.05
    mu = np.tile(gene_base, (n_obs, 1))
    block = max(1, n_genes // (n_clusters + 1))
    for k in range(n_clusters):
        rows = np.nonzero(labels == k)[0]
        mu[np.ix_(rows, np.arange(k * block, (k + 1) * block))] *= marker_fold
    depth = rng.uniform(0.5, 1.6, n_obs)[:, None]
    return rng.poisson(mu * depth).astype(np.float32), labels


def _pbmc(n_obs, n_genes, n_clusters, seed, n_comps=20, subset_hvg=False):
    """Build a fully-processed AnnData the way scanpy would: normalize_total -> log1p -> `.raw` -> (HVG
    subset) -> scale -> pca -> neighbors. Returns (adata, labels). Categoricals / colors / DE are layered
    on by the callers. `subset_hvg`: when True, X is restricted to the HVG subset so `.raw` keeps a
    *larger* gene set (the real `pbmc3k_processed` shape: X=1838 HVG, raw=13714); when False, X and `.raw`
    share the gene set (the real `pbmc68k_reduced` shape -- so an X<->raw chain through Seurat/SCE, which
    carry one feature space per assay, doesn't see two divergent gene axes)."""
    warnings.filterwarnings("ignore")
    import anndata as ad
    import pandas as pd
    import scanpy as sc
    import scipy.sparse as sp

    counts, labels = _counts(n_obs, n_genes, n_clusters, seed)
    obs = pd.DataFrame(index=[f"Cell{i}" for i in range(n_obs)])
    var = pd.DataFrame(index=[f"Gene{j}" for j in range(n_genes)])
    a = ad.AnnData(sp.csr_matrix(counts), obs=obs, var=var)
    # real-looking QC numerics (a plain float `percent_mito` exercises the "NaN, no mask" path)
    a.obs["n_genes"] = (counts > 0).sum(1).astype(np.int64)
    a.obs["n_counts"] = counts.sum(1).astype(np.float32)
    a.obs["percent_mito"] = (counts[:, :max(1, n_genes // 50)].sum(1) / np.maximum(counts.sum(1), 1)).astype(np.float32)

    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)
    a.raw = a                                           # log-normalized full gene set kept in `.raw`
    sc.pp.highly_variable_genes(a, n_top_genes=max(n_comps + 5, int(n_genes * 0.6)))
    if subset_hvg:
        a = a[:, a.var["highly_variable"]].copy()       # X = HVG subset (fewer genes than .raw: divergent)
    sc.pp.scale(a, max_value=10)
    sc.tl.pca(a, n_comps=min(n_comps, a.n_obs - 1, a.n_vars - 1))
    sc.pp.neighbors(a, n_neighbors=10)                  # -> obsp distances/connectivities + uns neighbors
    a.obsm["X_umap"] = (a.obsm["X_pca"][:, :2] +        # set umap directly (no umap-learn needed in CI)
                        np.random.default_rng(seed + 1).normal(0, 0.05, (a.n_obs, 2))).astype(np.float32)
    return a, labels


def pbmc68k_like(seed=0):
    """Stand-in for `pbmc68k_reduced`: real-pipeline AnnData with `bulk_labels`/`phase`/`louvain`
    categoricals + their `*_colors`, real `uns['pca']`/`neighbors`, divergent `.raw`, and a **logreg**
    one-vs-rest `rank_genes_groups` grouped by `bulk_labels` (the score-only DE variant)."""
    import pandas as pd
    import scanpy as sc

    a, labels = _pbmc(n_obs=300, n_genes=400, n_clusters=5, seed=seed, subset_hvg=False)
    names = np.array(["CD4+", "CD8+", "B", "NK", "Mono", "DC", "Mega", "pDC"])[: labels.max() + 1]
    a.obs["bulk_labels"] = pd.Categorical(names[labels])
    a.obs["louvain"] = pd.Categorical([str(l) for l in labels])
    phase = np.array(["G1", "S", "G2M"])[np.random.default_rng(seed + 2).integers(0, 3, a.n_obs)]
    a.obs["phase"] = pd.Categorical(phase)
    for col in ("bulk_labels", "louvain", "phase"):
        a.uns[col + "_colors"] = _palette(len(a.obs[col].cat.categories))
    sc.tl.rank_genes_groups(a, "bulk_labels", method="logreg")   # names+scores only (no lfc/pvals)
    return a


def pbmc3k_like(seed=1):
    """Stand-in for `pbmc3k_processed`: real-pipeline AnnData with a `louvain` categorical (+colors),
    `n_genes`/`percent_mito`, divergent `.raw` -- but **no** baked DE (the DE tests add it themselves)."""
    import pandas as pd

    a, labels = _pbmc(n_obs=360, n_genes=500, n_clusters=6, seed=seed, subset_hvg=True)
    a.obs["louvain"] = pd.Categorical([str(l) for l in labels])
    a.uns["louvain_colors"] = _palette(len(a.obs["louvain"].cat.categories))
    return a


def velocity(seed=7):
    """Stand-in for the scVelo pancreas fixture: `spliced`/`unspliced` count layers, a `clusters`
    categorical (+colors), and -- the point -- `uns['velocity_graph']` / `uns['velocity_graph_neg']` as
    sparse **cell x cell** transition graphs (scVelo's canonical uns location, NOT obsp)."""
    warnings.filterwarnings("ignore")
    import anndata as ad
    import pandas as pd
    import scanpy as sc
    import scipy.sparse as sp

    n_obs, n_genes, n_clusters = 150, 250, 4
    spliced, labels = _counts(n_obs, n_genes, n_clusters, seed)
    rng = np.random.default_rng(seed + 3)
    unspliced = rng.poisson(spliced * 0.4 + 0.2).astype(np.float32)      # introns: a fraction of mature
    a = ad.AnnData(sp.csr_matrix(spliced.copy()),
                   obs=pd.DataFrame(index=[f"Cell{i}" for i in range(n_obs)]),
                   var=pd.DataFrame(index=[f"Gene{j}" for j in range(n_genes)]))
    a.layers["spliced"] = sp.csr_matrix(spliced)
    a.layers["unspliced"] = sp.csr_matrix(unspliced)
    sc.pp.normalize_total(a); sc.pp.log1p(a)
    # real scVelo output carries the standard pipeline products too (moments runs pca+neighbors) -- run
    # them so the synthetic stand-in mirrors the real object's full structure (pca/neighbors/umap), not
    # just the velocity-specific layers/graph (enforced by test_synth_faithful).
    sc.pp.pca(a, n_comps=10); sc.pp.neighbors(a, n_neighbors=10)
    a.obsm["X_umap"] = (a.obsm["X_pca"][:, :2] + rng.normal(0, 0.05, (n_obs, 2))).astype(np.float32)
    a.obs["clusters"] = pd.Categorical(np.array(["Ductal", "Ngn3", "Pre", "Beta"])[labels])
    a.uns["clusters_colors"] = _palette(n_clusters)
    # velocity graph: each cell transitions to ~12 neighbours; cosine-like signed weights, pos/neg split
    deg = 12
    rows = np.repeat(np.arange(n_obs), deg)
    cols = rng.integers(0, n_obs, n_obs * deg)
    w = rng.uniform(-1, 1, n_obs * deg).astype(np.float32)
    pos = w > 0
    a.uns["velocity_graph"] = sp.csr_matrix((w[pos], (rows[pos], cols[pos])), shape=(n_obs, n_obs))
    a.uns["velocity_graph_neg"] = sp.csr_matrix((w[~pos], (rows[~pos], cols[~pos])), shape=(n_obs, n_obs))
    return a


# ---- multimodal (MuData): RNA + ADT, the same L* shape as a Seurat multi-assay object ---------------

def _citeseq_arrays(n_obs=80, n_genes=27, n_prot=29, seed=11):
    """Synthetic CITE-seq counts: an RNA block with cluster/marker structure + an ADT block whose protein
    levels track the cell clusters (so the two modalities are correlated, like real CITE-seq). Dims match
    the historical real fixture (80 cells x 27 genes + 29 proteins) so the dimension asserts are stable."""
    rna, labels = _counts(n_obs, n_genes, min(5, n_genes), seed)
    rng = np.random.default_rng(seed + 4)
    prot_base = rng.gamma(1.0, 4.0, n_prot) + 1.0
    mu = np.tile(prot_base, (n_obs, 1))
    block = max(1, n_prot // (labels.max() + 2))
    for k in range(labels.max() + 1):
        rows = np.nonzero(labels == k)[0]
        mu[np.ix_(rows, np.arange(k * block, (k + 1) * block))] *= 5.0
    adt = rng.poisson(mu).astype(np.float32)
    cells = [f"Cell{i}" for i in range(n_obs)]
    genes = [f"Gene{j}" for j in range(n_genes)]
    prots = [f"ADT{j}" for j in range(n_prot)]
    return rna, adt, labels, cells, genes, prots


def citeseq_matrices(n_obs=80, n_genes=27, n_prot=29, seed=11):
    """(rna, adt, cells, genes, proteins) -- the language-agnostic CITE-seq arrays (CSR float32)."""
    import scipy.sparse as sp
    rna, adt, _labels, cells, genes, prots = _citeseq_arrays(n_obs, n_genes, n_prot, seed)
    return sp.csr_matrix(rna), sp.csr_matrix(adt), cells, genes, prots


def citeseq_mudata(seed=11):
    """A minimal RNA+ADT MuData (27 genes + 29 proteins) -- the multimodal round-trip fixture."""
    try:
        import anndata as ad
        import mudata
        import pandas as pd
    except Exception:
        return None
    rna, adt, cells, genes, prots = citeseq_matrices(seed=seed)
    arna = ad.AnnData(rna, obs=pd.DataFrame(index=cells), var=pd.DataFrame(index=genes))
    aadt = ad.AnnData(adt, obs=pd.DataFrame(index=cells), var=pd.DataFrame(index=prots))
    return mudata.MuData({"rna": arna, "prot": aadt})


def write_citeseq_mtx(outdir, seed=11):
    """Materialize the CITE-seq arrays as language-agnostic Matrix-Market + label files (cells x features,
    matching the historical fixture layout) so the R conformance scripts (Seurat/SCE multimodal) read the
    *same* synthetic RNA+ADT the Python MuData test uses -- one generator, both languages, nothing
    committed. Returns the output dir."""
    import os

    import scipy.io as sio
    rna, adt, cells, genes, prots = citeseq_matrices(seed=seed)
    os.makedirs(outdir, exist_ok=True)
    sio.mmwrite(os.path.join(outdir, "rna.mtx"), rna)            # cells x genes
    sio.mmwrite(os.path.join(outdir, "adt.mtx"), adt)            # cells x proteins
    for fn, items in (("cells.txt", cells), ("genes.txt", genes), ("proteins.txt", prots)):
        with open(os.path.join(outdir, fn), "w") as fh:
            fh.write("\n".join(items) + "\n")
    return outdir


def citeseq_mudata_annotated(seed=12):
    """A richer RNA+ADT MuData standing in for the downloaded minipbcite: per-modality processing
    (own PCA in `obsm` + `uns`), a **global** `celltype` categorical (-> factor axis), a global `obsm`
    (WNN-style joint embedding), and per-modality `uns` -- exercising the global-obs / global-obsm /
    per-mod-uns paths the real minipbcite test asserts on."""
    try:
        import anndata as ad
        import mudata
        import pandas as pd
        import scanpy as sc
    except Exception:
        return None
    warnings.filterwarnings("ignore")
    n_obs = 411
    rna, _l = _counts(n_obs, 120, 5, seed)
    adt, labels = _counts(n_obs, 15, 5, seed + 1)
    cells = [f"Cell{i}" for i in range(n_obs)]
    arna = ad.AnnData(__import__("scipy.sparse", fromlist=["csr_matrix"]).csr_matrix(rna),
                      obs=pd.DataFrame(index=cells), var=pd.DataFrame(index=[f"Gene{j}" for j in range(120)]))
    aadt = ad.AnnData(__import__("scipy.sparse", fromlist=["csr_matrix"]).csr_matrix(adt),
                      obs=pd.DataFrame(index=cells), var=pd.DataFrame(index=[f"ADT{j}" for j in range(15)]))
    for a in (arna, aadt):
        sc.pp.normalize_total(a); sc.pp.log1p(a); sc.pp.pca(a, n_comps=10)   # per-mod obsm['X_pca'] + uns['pca']
    md = mudata.MuData({"rna": arna, "prot": aadt})
    names = np.array(["CD4+", "CD8+", "B", "NK", "Mono"])
    md.obs["celltype"] = pd.Categorical(names[labels])                       # global obs -> factor axis
    md.obs["leiden"] = pd.Categorical([str(l) for l in labels])
    md.obsm["X_wnn_umap"] = (arna.obsm["X_pca"][:, :2]).astype(np.float32)   # global joint embedding
    return md


if __name__ == "__main__":                                       # `python synth.py citeseq <dir>` for R
    import sys
    if len(sys.argv) >= 3 and sys.argv[1] == "citeseq":
        print(write_citeseq_mtx(sys.argv[2]))
    else:
        raise SystemExit("usage: python synth.py citeseq <outdir>")
