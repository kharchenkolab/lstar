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


def test_synthetic_backed_promotion():
    """The real-atlas test above is local-only (skips in CI). This guards the **backed read + promotion**
    path in CI on a small synthetic atlas-like object: write a synth pbmc-like `.h5ad`, open it *backed*
    (X stays on disk), and confirm factor axes + `*_colors` + pca-variance still promote and validate."""
    import tempfile

    import anndata as ad
    sys.path.insert(0, os.path.dirname(__file__))
    import synth  # noqa: E402

    a = synth.pbmc68k_like()
    p = os.path.join(tempfile.mkdtemp(), "synth_atlas.h5ad"); a.write_h5ad(p)
    ab = ad.read_h5ad(p, backed="r")                            # X stays on disk
    ds = lstar.read_anndata(ab)
    facs = [n for n, x in ds.axes.items() if x.role == "factor"]
    assert "bulk_labels" in facs and ds.axis("bulk_labels").role == "factor"
    assert ds.field("bulk_labels_colors").subtype == "color"
    assert ds.field("pca_variance_ratio").span == ["pca"]
    assert not [i for i in lstar.validate(ds) if i.startswith("ERROR")]
    print("backed promotion (synthetic atlas): factor axes + colors + pca variance promote on a backed read")


if __name__ == "__main__":
    test_marrow_realistic_size()
    test_synthetic_backed_promotion()
