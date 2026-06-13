"""Fetch REAL imaging-based spatial datasets from `squidpy.datasets` into testdata/spatial/ as .h5ad.

Run with the squidpy side-venv (py>=3.10), NOT lstar's system python:
    python3 -m venv /tmp/sq && /tmp/sq/bin/pip install squidpy
    /tmp/sq/bin/python3 conformance/sweep/sweep_spatial_fetch_sq.py

Each `sq.datasets.*()` returns a preprocessed AnnData with `obsm['spatial']` (cell/bead coords) and,
for some, image keys under `uns['spatial']`. These are imaging-based platforms with TARGETED gene panels
(small feature axis) and -- for MERFISH/seqFISH -- molecule/coord provenance distinct from spot grids:

  - merfish      : Moffitt hypothalamus MERFISH (~73k cells, ~150-gene panel, has a `Bregma` z-section
                   factor -> a multi-section imaging case)
  - seqfish      : Lohoff mouse embryo seqFISH (~57k cells, 3D-ish coords)
  - slideseqv2   : mouse hippocampus Slide-seqV2 (~41k beads -- high-res bead array, not a spot grid)
  - imc          : imaging mass cytometry (PROTEIN spatial -- a proteins-style feature axis, no RNA)

These are written by squidpy's anndata (0.11) -> a newer .h5ad schema than lstar's system anndata 0.8,
so the sweep also exercises cross-version .h5ad recognition (format-variety bonus).

LOCAL-ONLY: testdata/ is gitignored; not committed, not fetched by CI.
"""
import warnings, os, sys
warnings.filterwarnings("ignore")
import squidpy as sq

OUT = "testdata/spatial"
os.makedirs(OUT, exist_ok=True)

JOBS = [
    ("sq_merfish", sq.datasets.merfish),
    ("sq_seqfish", sq.datasets.seqfish),
    ("sq_slideseqv2", sq.datasets.slideseqv2),
    ("sq_imc", sq.datasets.imc),
]

for name, fn in JOBS:
    dest = os.path.join(OUT, name + ".h5ad")
    if os.path.exists(dest):
        print("skip (cached):", name)
        continue
    try:
        a = fn(path=os.path.join(OUT, "_sqcache", name + ".h5ad"))
        try:
            a.var_names_make_unique()
        except Exception:
            pass
        a.write_h5ad(dest)
        print("OK  %-16s shape=%s obsm=%s uns=%s"
              % (name, a.shape, list(a.obsm.keys()), list(a.uns.keys())))
    except Exception as e:
        print("FAIL %-16s %s" % (name, str(e)[:140]))
        sys.stdout.flush()
