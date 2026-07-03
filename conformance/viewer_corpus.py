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

DATASETS = ["pbmc68k", "pbmc3k"]                         # curated: multiple categoricals + an embedding


def _load(name):
    import corpus
    return {"pbmc68k": corpus.pbmc68k_reduced, "pbmc3k": corpus.pbmc3k_processed}[name]()


def cmd_base(name, out):
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
