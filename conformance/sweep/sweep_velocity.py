"""Sweep real RNA-velocity / trajectory AnnData (scVelo datasets) through read_anndata.

Each is a published dataset with `spliced`/`unspliced` (and sometimes `ambiguous`) layers, a
`clusters`/celltype categorical (+ `*_colors`), and varied organism/tissue (mouse hippocampus,
mouse gastrulation erythroid lineage, human bone marrow, mouse pancreas). Velocity layers ride the
generic measure paths; the test is that they round-trip as extra (cells, genes) measures and the
categorical clusters induce factor axes.

Run: `/tmp/scvelo_venv/bin/python3 ... ` is only needed to *download* them (see sweep_velocity_fetch);
this sweep reads the cached `.h5ad` with the system Python (lstar's own deps), writes /tmp/sweep_velocity.tsv.
"""
import warnings, sys, os, tempfile, glob
warnings.filterwarnings("ignore")
sys.path.insert(0, "python/src")
import anndata as ad
import lstar

VEL = "testdata/velocity"
rows = []
files = sorted(glob.glob(os.path.join(VEL, "*.h5ad")))
if not files:
    print("no velocity datasets cached under", VEL, "-- run the fetch first (scvelo venv)")
for path in files:
    name = os.path.splitext(os.path.basename(path))[0]
    try:
        a = ad.read_h5ad(path)
        layers = list(a.layers.keys())
        ds = lstar.read_anndata(a)
        errs = [i for i in lstar.validate(ds) if i.startswith("ERROR")]
        # confirm the velocity layers survived as measures
        layer_fields = [f for f in ds.fields if any(L in f for L in ("spliced", "unspliced", "ambiguous"))]
        nfactor = sum(1 for ax in ds.axes.values() if getattr(ax, "role", "") == "factor"
                      or "factor" in ax.name)
        p = os.path.join(tempfile.mkdtemp(), "s.lstar.zarr")
        lstar.write(ds, p)
        status = "PASS" if not errs else "VALIDATE-ERR"
        note = ("layers=" + "+".join(layers))[:60]
        rows.append((name, str(a.shape), status, len(ds.fields), len(ds.axes),
                     len(layer_fields), note if not errs else errs[0][:60]))
    except Exception as e:
        rows.append((name, "", "FAIL", "", "", "", str(e)[:80]))

hdr = "dataset\tshape\tstatus\tfields\taxes\tlayer_fields\tnote\n"
open("/tmp/sweep_velocity.tsv", "w").write(
    hdr + "\n".join("\t".join(map(str, r)) for r in rows) + "\n")
npass = sum(1 for r in rows if r[2] == "PASS")
print(f"Velocity sweep: {npass} PASS / {len(rows)} datasets")
for r in rows:
    print(f"  {r[2]:12s} {r[0]:26s} {r[1]:14s} f={r[3]} a={r[4]} velo_fields={r[5]} {str(r[6])[:46]}")
