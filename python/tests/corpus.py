"""Example datasets the test suite is grounded in -- the **two-tier corpus**.

*Locally* (the default), these loaders fetch/derive **real** datasets (cached in the gitignored
`testdata/` dir; large/local-only atlases used when present) so the profiles are exercised against
structures a real tool actually produced -- categoricals, color palettes, `rank_genes_groups`
(t-test/wilcoxon/logreg/pairwise), PCA variance, graphs, RNA+ADT modalities. The extensive real breadth
lives in `conformance/sweep/`.

*In CI* (`LSTAR_SYNTHETIC_CORPUS=1`), each loader instead returns a **synthetic-but-faithful** stand-in
from `synth.py` -- synthetic counts run through the *real* scanpy/mudata pipeline, so the same library
code produces the same structure (no real datasets committed to or downloaded by github, which are large
and slow). Keeping the synthetic fixtures structurally representative of the real corpus is the explicit
contract; the local real runs + the sweep are what verify it.

Loaders return `None` only when a dataset is genuinely unavailable (offline / missing optional dep) --
callers then print a SKIP note rather than silently passing.
"""
import os
import warnings

TESTDATA = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "testdata"))
FIXTURES = os.path.abspath(os.path.join(os.path.dirname(__file__), "fixtures"))

# CI sets this -> serve the synthetic, download-free corpus (faithful to the real one, see synth.py).
SYNTHETIC = os.environ.get("LSTAR_SYNTHETIC_CORPUS") == "1"
if SYNTHETIC:
    import synth


def _sc():
    warnings.filterwarnings("ignore")
    import scanpy as sc
    os.makedirs(TESTDATA, exist_ok=True)
    sc.settings.datasetdir = TESTDATA
    sc.settings.verbosity = 0
    return sc


def pbmc68k_reduced():
    """700x765 real PBMC, fully processed: obs categoricals (bulk_labels/phase/louvain), `*_colors`,
    `uns['pca']`, a real one-vs-rest `rank_genes_groups` (method=logreg -> names+scores only), and the
    `neighbors` OverloadedDict. Downloads ~ a few MB; cached."""
    if SYNTHETIC:
        return synth.pbmc68k_like()
    try:
        return _sc().datasets.pbmc68k_reduced()
    except Exception as e:                     # pragma: no cover - only on a genuine fetch failure
        print("  [corpus] pbmc68k_reduced unavailable:", e)
        return None


def pbmc3k_processed():
    """2638x1838 real PBMC, fully processed (louvain, `*_colors`, pca, neighbors, tsne/umap). ~23 MB."""
    if SYNTHETIC:
        return synth.pbmc3k_like()
    try:
        return _sc().datasets.pbmc3k_processed()
    except Exception as e:                     # pragma: no cover
        print("  [corpus] pbmc3k_processed unavailable:", e)
        return None


def pbmc3k_with_de(method="t-test", reference="rest", groupby="louvain"):
    """pbmc3k_processed with a **real** `rank_genes_groups` computed by the real scanpy pipeline --
    method in {t-test, wilcoxon, logreg, ...}; `reference="rest"` (one-vs-rest) or a group name
    (pairwise). t-test/wilcoxon yield the full names/scores/logfoldchanges/pvals/pvals_adj bundle."""
    sc = _sc()
    a = pbmc3k_processed()
    if a is None:
        return None
    a = a.copy()
    try:
        sc.tl.rank_genes_groups(a, groupby, method=method, reference=reference)
    except Exception as e:                     # pragma: no cover
        print("  [corpus] rank_genes_groups(%s,%s) failed:" % (method, reference), e)
        return None
    return a


# A real, realistically-sized atlas (Tabula Muris Senis droplet, Marrow; official annotations).
# Local-only (1.2 GB); used for a realistic-size smoke when present, else noted.
MARROW = "/home/pkharchenko/cacoa/age/tab.muris/" \
         "tabula-muris-senis-droplet-processed-official-annotations-Marrow.h5ad"


def marrow_backed():
    """The real Marrow atlas (40220x20138) opened backed, or None if not on this machine."""
    import anndata as ad
    if not os.path.exists(MARROW):
        return None
    try:
        return ad.read_h5ad(MARROW, backed="r")
    except Exception as e:                     # pragma: no cover
        print("  [corpus] Marrow unavailable:", e)
        return None


CITESEQ = os.path.join(TESTDATA, "citeseq")     # local cache (gitignored), derived from real minipbcite


def citeseq_matrices():
    """CITE-seq RNA+ADT arrays (80 cells × 27 genes + 29 proteins) as language-agnostic Matrix-Market, so
    the MuData (Python) and Seurat/SCE (R) multimodal tests share one source. CI: synthetic (`synth`).
    Local: subsampled from the **real** downloaded minipbcite, cached under `testdata/citeseq/`. Returns
    (rna, adt, cells, genes, proteins) or None."""
    if SYNTHETIC:
        return synth.citeseq_matrices()
    import scipy.io as sio
    if not os.path.exists(os.path.join(CITESEQ, "rna.mtx")):
        if _derive_real_citeseq(CITESEQ) is None:
            return None
    rd = lambda f: open(os.path.join(CITESEQ, f)).read().split()
    rna = sio.mmread(os.path.join(CITESEQ, "rna.mtx")).tocsr().astype("float32")
    adt = sio.mmread(os.path.join(CITESEQ, "adt.mtx")).tocsr().astype("float32")
    return rna, adt, rd("cells.txt"), rd("genes.txt"), rd("proteins.txt")


def _derive_real_citeseq(outdir):
    """Subsample the real minipbcite (80 cells × 27 genes + 29 proteins) to Matrix-Market + label files
    -- the local real source for the multimodal tests, replacing a committed fixture. None if unavailable."""
    md = minipbcite()
    if md is None:
        return None
    try:
        import numpy as np
        import scipy.io as sio
        import scipy.sparse as sp
        rna_a, adt_a = md.mod["rna"], md.mod["prot"]
        ci = np.arange(min(80, md.n_obs)); gi = np.arange(min(27, rna_a.n_vars)); pi = np.arange(min(29, adt_a.n_vars))
        dense = lambda X: X.toarray() if sp.issparse(X) else np.asarray(X)
        rna = sp.csr_matrix(dense(rna_a.X)[np.ix_(ci, gi)].astype("float32"))
        adt = sp.csr_matrix(dense(adt_a.X)[np.ix_(ci, pi)].astype("float32"))
        os.makedirs(outdir, exist_ok=True)
        sio.mmwrite(os.path.join(outdir, "rna.mtx"), rna)
        sio.mmwrite(os.path.join(outdir, "adt.mtx"), adt)
        for fn, items in (("cells.txt", md.obs_names[ci]), ("genes.txt", rna_a.var_names[gi]),
                          ("proteins.txt", adt_a.var_names[pi])):
            open(os.path.join(outdir, fn), "w").write("\n".join(map(str, items)) + "\n")
        return outdir
    except Exception as e:                         # pragma: no cover
        print("  [corpus] could not derive real citeseq:", e); return None


def citeseq_mudata():
    """A CITE-seq MuData (RNA + ADT, 27 genes + 29 proteins). None if mudata absent."""
    if SYNTHETIC:
        return synth.citeseq_mudata()
    try:
        import anndata as ad
        import mudata
        import pandas as pd
    except Exception:
        return None
    got = citeseq_matrices()
    if got is None:
        return None
    rna, adt, cells, genes, prots = got
    arna = ad.AnnData(rna, obs=pd.DataFrame(index=cells), var=pd.DataFrame(index=genes))
    aadt = ad.AnnData(adt, obs=pd.DataFrame(index=cells), var=pd.DataFrame(index=prots))
    return mudata.MuData({"rna": arna, "prot": aadt})


def minipbcite():
    """Real CITE-seq MuData (411 cells, RNA+ADT) -- downloaded + cached locally. CI: a richer synthetic
    annotated MuData (`synth`). Returns None if mudata absent or the fetch fails."""
    if SYNTHETIC:
        return synth.citeseq_mudata_annotated()
    try:
        import mudata
    except Exception:
        return None
    os.makedirs(TESTDATA, exist_ok=True)
    p = os.path.join(TESTDATA, "minipbcite.h5mu")
    if not os.path.exists(p):
        try:
            import urllib.request
            urllib.request.urlretrieve(
                "https://github.com/gtca/h5xx-datasets/raw/main/datasets/minipbcite.h5mu", p)
        except Exception as e:                     # pragma: no cover
            print("  [corpus] minipbcite download failed:", e); return None
    try:
        return mudata.read_h5mu(p)
    except Exception as e:                         # pragma: no cover
        print("  [corpus] minipbcite read failed:", e); return None


def pancreas_velocity():
    """RNA-velocity object: `spliced`/`unspliced` layers, a `clusters` categorical + `*_colors`, and a
    `uns['velocity_graph']` (+ `_neg`) -- the canonical scVelo location (cell x cell, NOT obsp). CI:
    synthetic (`synth`). Local: regenerated by running the **real** scVelo pipeline on scVelo's pancreas
    dataset (subsampled), cached under `testdata/`. None if scvelo is unavailable."""
    if SYNTHETIC:
        return synth.velocity()
    import anndata as ad
    p = os.path.join(TESTDATA, "pancreas_velocity_small.h5ad")
    if os.path.exists(p):
        try:
            return ad.read_h5ad(p)
        except Exception:                          # pragma: no cover
            pass
    try:                                           # regenerate from real scVelo (heavy; local only)
        import numpy as np
        import scvelo as scv
        a = scv.datasets.pancreas()
        a = a[np.random.default_rng(0).choice(a.n_obs, 150, replace=False)].copy()
        scv.pp.filter_and_normalize(a, n_top_genes=250, min_shared_counts=0)
        scv.pp.moments(a); scv.tl.velocity(a); scv.tl.velocity_graph(a)
        os.makedirs(TESTDATA, exist_ok=True)
        a.write_h5ad(p)
        return a
    except Exception as e:                         # pragma: no cover
        print("  [corpus] pancreas_velocity unavailable (needs scvelo):", e); return None
