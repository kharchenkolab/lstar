"""Collection-level grouped reducer (`collection_pseudobulk`): aggregate a collection's per-sample
`counts.<s>` measures into `(joint-group x genes)`, grouped by a factor over the **union** cells axis,
**streamed one sample at a time** (the joint matrix is never materialized) with float64 accumulation.
This is the storage-backed Conos / pagoda2 pseudobulk path. The test proves the streamed result equals a
materialized reference, exercises gene alignment across samples (a gene absent in one sample contributes
0), and validates the output fields.

Run: PYTHONPATH=python/src python3 python/tests/test_collection_reduce.py
"""
import os
import sys
import tempfile

import numpy as np
import scipy.sparse as sp

sys.path.insert(0, os.path.dirname(__file__))
import lstar  # noqa: E402


def _collection(per_sample_genes=False):
    """A synthetic 3-sample collection: per-sample cells/measures, a derived union cells axis, and a
    joint clustering over the union. `per_sample_genes` gives each sample its own gene axis with a
    slightly different gene set (to exercise label-alignment + 0-fill)."""
    rng = np.random.default_rng(0)
    ng = 12
    ds = lstar.Dataset(kind="collection")
    samples = {"S1": 30, "S2": 40, "S3": 20}
    ds.add_axis("samples", list(samples), role="sample")
    if not per_sample_genes:
        ds.add_axis("genes", [f"g{i}" for i in range(ng)], role="feature")
    union, mats, gene_lists = [], [], []
    for si, (s, nc) in enumerate(samples.items()):
        cells_s = [f"{s}_c{i}" for i in range(nc)]
        ds.add_axis("cells.%s" % s, cells_s, role="observation")
        if per_sample_genes:                                   # sample s drops gene (ng-1-si): different sets
            gl = [f"g{i}" for i in range(ng) if i != ng - 1 - si]
            ds.add_axis("genes.%s" % s, gl, role="feature"); gax = "genes.%s" % s
        else:
            gl = [f"g{i}" for i in range(ng)]; gax = "genes"
        M = sp.random(nc, len(gl), density=0.4, format="csr", random_state=si + 1).astype("float32")
        ds.add_field("counts.%s" % s, M, role="measure", span=["cells.%s" % s, gax], state="raw")
        union += cells_s; mats.append((gl, M.toarray()))
    ds.add_axis("cells", union, origin="derived", role="observation")          # the union axis
    codes = rng.integers(0, 4, len(union))
    ds.add_field("cluster", lstar.Categorical(codes, np.array(["A", "B", "C", "D"])), span=["cells"])
    return ds, codes, mats


def _materialized_ref(mats, codes, K=4):
    """Reference: assemble the full union x (canonical genes) dense matrix and reduce per group."""
    gidx, genes = {}, []
    for gl, _ in mats:
        for g in gl:
            if g not in gidx:
                gidx[g] = len(genes); genes.append(g)
    U = np.zeros((len(codes), len(genes)), dtype=np.float64); r = 0
    for gl, A in mats:
        cols = [gidx[g] for g in gl]
        U[r:r + A.shape[0]][:, cols] = A; r += A.shape[0]
    mean = np.zeros((K, len(genes))); frac = np.zeros((K, len(genes)))
    for k in range(K):
        rows = np.nonzero(codes == k)[0]
        mean[k] = U[rows].mean(0); frac[k] = (U[rows] > 0).mean(0)
    return mean, frac, genes


def _check(per_sample_genes):
    ds, codes, mats = _collection(per_sample_genes)
    out = lstar.collection_pseudobulk(ds, "cluster", field="counts", lognorm=False)
    rmean, rfrac, genes = _materialized_ref(mats, codes)
    assert out["mean"].shape == (4, len(genes))
    assert np.allclose(out["mean"], rmean, rtol=1e-5, atol=1e-7), "streamed mean != materialized"
    assert np.allclose(out["frac"], rfrac, rtol=1e-5, atol=1e-7), "streamed frac != materialized"
    assert ds.field("pb.cluster.mean").span == ["cluster", "genes"]
    assert ds.field("pb.cluster.mean").provenance.get("collection") is True
    assert not lstar.validate(ds)
    # round-trips through the store
    p = os.path.join(tempfile.mkdtemp(), "c.lstar.zarr"); lstar.write(ds, p)
    assert not lstar.validate(lstar.read(p))
    return len(genes)


def test_collection_pseudobulk_shared_genes():
    _check(per_sample_genes=False)
    print("collection pseudobulk (shared genes): streamed (cluster x genes) == materialized; round-trips")


def test_collection_pseudobulk_aligns_differing_gene_sets():
    ng = _check(per_sample_genes=True)
    print("collection pseudobulk (per-sample gene sets, %d-gene union): label-aligned, 0-filled where "
          "absent; streamed == materialized" % ng)


if __name__ == "__main__":
    test_collection_pseudobulk_shared_genes()
    test_collection_pseudobulk_aligns_differing_gene_sets()
