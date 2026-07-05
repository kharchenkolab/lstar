"""Corpus-driven viewer@0.1 base builder: load a curated corpus dataset (REAL locally; synthetic-but-
faithful in CI via LSTAR_SYNTHETIC_CORPUS=1) and convert it to an L* store -- no prep. The shell driver
(viewer_corpus.sh) then preps it on each surface and cross-checks ALL viewer fields.

Corpus data is where the divergences actually bite: real categoricals in non-alphabetical stored order
(#3), several competing groupings like bulk_labels/phase/louvain (#4), real embeddings, native CSR.
Skips a dataset cleanly (exit 7) when its loader is unavailable (offline / missing dep). No real data is
committed or downloaded by CI -- the synthetic stand-ins (python/tests/synth.py) run the same pipeline.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "python", "src"))
sys.path.insert(0, os.path.join(HERE, "..", "python", "tests"))

import lstar                                            # noqa: E402
from lstar.profiles.anndata import read_anndata         # noqa: E402

DATASETS = ["pbmc68k", "pbmc3k", "dense_primary"]        # curated corpora + a synthetic DENSE-primary case


def _load(name):
    import corpus
    return {"pbmc68k": corpus.pbmc68k_reduced, "pbmc3k": corpus.pbmc3k_processed}[name]()


def _build_dense_primary(out):
    """A store whose PRIMARY measure is DENSE (encoding="dense", on-disk /fields/X/values) — the shape an
    SCE logcounts assay or a scaled/dense AnnData X takes. The corpora above all carry a SPARSE (csc) counts
    basis, so without this case the DENSE measure-read path is exercised on NO surface through the viewer
    parity — which is exactly how a JS extendForViewer that hardcoded the sparse /data path shipped: dense
    primary -> NotFoundError -> viewer-opt silently skipped. Always available (no corpus/scanpy dependency)."""
    import numpy as np
    rng = np.random.default_rng(0)
    n, g, npc = 120, 30, 5
    ds = lstar.Dataset(kind="sample")
    ds.add_axis("cells", ["c%d" % i for i in range(n)])
    ds.add_axis("genes", ["g%d" % i for i in range(g)])
    ds.add_axis("pca", ["PC%d" % i for i in range(npc)], origin="derived")
    X = (rng.poisson(0.4, size=(n, g)) * rng.random((n, g))).astype("float32")   # dense, many exact zeros
    ds.add_field("logcounts", X, role="measure", span=["cells", "genes"], state="lognorm")   # DENSE primary
    ds.add_field("cluster", np.array(["k%d" % (i % 4) for i in range(n)]), role="label", span=["cells"])
    ds.add_field("phase", np.array(["G1", "S", "G2M"])[np.arange(n) % 3], role="label", span=["cells"])
    ds.add_field("pca", rng.standard_normal((n, npc)).astype("float32"), role="embedding", span=["cells", "pca"])
    lstar.write(ds, out)
    print("  [synthetic] dense_primary -> L* base: DENSE primary measure `logcounts` (encoding=dense), "
          "groupings=[cluster, phase], embedding=[pca]")
    print("BASIS=lognorm")                                # dense measure is log-normalized -> prep from it as lognorm


def cmd_base(name, out):
    if name == "dense_primary":
        _build_dense_primary(out)
        return
    try:
        a = _load(name)
    except Exception as e:                               # a broken/missing optional dep -> clean skip
        print(f"  SKIP {name} (corpus load error: {type(e).__name__}: {str(e)[:120]})")
        sys.exit(7)
    if a is None:
        print(f"  SKIP {name} (unavailable)")
        sys.exit(7)
    ds = read_anndata(a)
    lstar.write(ds, out)
    meas = [(n, f.state) for n, f in ds.fields.items() if f.role == "measure" and f.span and len(f.span) == 2]
    labels = [n for n, f in ds.fields.items() if f.role == "label"]
    embs = [n for n, f in ds.fields.items() if f.role == "embedding"]
    has_raw = any(f.state == "raw" for _, f in ds.fields.items() if f.role == "measure")
    print(f"  [corpus] {name} -> L* base: measures={meas}; labels={labels}; embeddings={embs}")
    print("BASIS=" + ("" if has_raw else "lognorm"))    # machine-readable line for the shell driver


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else ""
    if cmd == "base":
        cmd_base(sys.argv[2], sys.argv[3])
    elif cmd == "datasets":
        print(" ".join(DATASETS))
    else:
        print(__doc__)
        sys.exit(2)
