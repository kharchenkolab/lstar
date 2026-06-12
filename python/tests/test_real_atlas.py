"""Realistic-size grounding on a REAL atlas (Tabula Muris Senis droplet, Marrow; 40220 x 20138, official
annotations). Local-only (1.2 GB) -- skips with a note where the file isn't present (the downloadable
corpus, pbmc3k/pbmc68k, grounds CI). Reads backed (matrices stay on disk), so the profile, factor-axis
induction, color/PCA promotion, and validation are exercised against a real, realistically-sized object
with many real annotations -- not a toy.

Run: PYTHONPATH=python/src python3 python/tests/test_real_atlas.py
"""
import os
import sys
import warnings

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
import corpus  # noqa: E402

import lstar  # noqa: E402

warnings.filterwarnings("ignore")


def test_marrow_realistic_size():
    a = corpus.marrow_backed()
    if a is None:
        print("  SKIP test_marrow_realistic_size (TMS Marrow atlas not on this machine; "
              "local-only realistic-size smoke -- CI grounds on the downloadable corpus)"); return
    ds = lstar.read_anndata(a)                                  # backed: X stays on disk

    assert len(ds.axis("cells")) == 40220 and len(ds.axis("genes")) == 20138
    # real categorical annotations (age, sex, cell_ontology_class, leiden, louvain, ...) -> factor axes
    facs = [n for n, x in ds.axes.items() if x.role == "factor"]
    assert len(facs) >= 10, facs
    assert ds.axis("cell_ontology_class").role == "factor"
    # real color palette promoted, bound to its factor axis
    cc = ds.field("cell_ontology_class_colors")
    assert cc.subtype == "color" and cc.span == ["cell_ontology_class"]
    assert len(np.asarray(cc.values)) == len(ds.axis("cell_ontology_class"))
    # real PCA variance promoted over the pca axis
    assert ds.field("pca_variance_ratio").span == ["pca"]
    # the untyped tail (leiden/louvain/neighbors params) preserved in the passthrough, not dropped
    assert {"leiden", "louvain", "neighbors"} & set(ds.aux.get("anndata.uns", {}))

    errs = [i for i in lstar.validate(ds) if i.startswith("ERROR")]
    assert not errs, errs
    print("real atlas (TMS Marrow 40220x20138): read backed + validate clean; %d factor axes from real "
          "annotations; colors + pca variance promoted" % len(facs))


if __name__ == "__main__":
    test_marrow_realistic_size()
