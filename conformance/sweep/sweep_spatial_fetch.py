"""Fetch a spread of REAL 10x Visium samples (scanpy `visium_sge`) into testdata/spatial/ as .h5ad.

`sc.datasets.visium_sge(sample_id)` downloads a Space Ranger output (filtered feature-barcode .h5 +
the `spatial/` dir with tissue image + scalefactors) and returns an AnnData carrying
`obsm['spatial']` (spot pixel coords) and `uns['spatial'][library_id]` (images/scalefactors/metadata).
We cache the assembled AnnData as a single `.h5ad` per sample so the sweep (sweep_spatial.py) reads it
with lstar's own deps -- no scanpy/download needed at sweep time.

LOCAL-ONLY: testdata/ is gitignored; nothing here is committed or fetched by CI.

Sample selection (a breadth spread across the 27 available, not all 27):
  - tissues: breast cancer, heart, lymph node, kidney, brain (human + mouse)
  - a MULTI-SECTION pair: mouse brain sagittal posterior section 1 + 2 (a spatial *collection*)
  - a TARGETED-vs-PARENT pair: human cerebellum (panel-restriction case, spaceranger 1.2.0)

Run:  python3 conformance/sweep/sweep_spatial_fetch.py
"""
import warnings, os, sys
warnings.filterwarnings("ignore")
import scanpy as sc

OUT = "testdata/spatial"
os.makedirs(OUT, exist_ok=True)
sc.settings.datasetdir = OUT          # download cache lives alongside the .h5ad
sc.settings.verbosity = 1

SAMPLES = [
    "V1_Breast_Cancer_Block_A_Section_1",       # human breast cancer
    "V1_Human_Heart",                            # human heart
    "V1_Human_Lymph_Node",                       # human lymph node
    "V1_Mouse_Kidney",                           # mouse kidney
    "V1_Adult_Mouse_Brain",                      # mouse brain (coronal)
    "V1_Mouse_Brain_Sagittal_Posterior",         # multi-section pair (1/2)
    "V1_Mouse_Brain_Sagittal_Posterior_Section_2",
    "Targeted_Visium_Human_Cerebellum_Neuroscience",   # targeted/parent pair (panel restriction)
    "Parent_Visium_Human_Cerebellum",
]

for sid in SAMPLES:
    dest = os.path.join(OUT, sid + ".h5ad")
    if os.path.exists(dest):
        print("skip (cached):", sid)
        continue
    try:
        a = sc.datasets.visium_sge(sample_id=sid)
        a.var_names_make_unique()
        a.write_h5ad(dest)
        print("OK  %-50s shape=%s obsm=%s uns_spatial=%s"
              % (sid, a.shape, list(a.obsm.keys()), list(a.uns.get("spatial", {}).keys())))
    except Exception as e:
        print("FAIL %-50s %s" % (sid, str(e)[:120]))
        sys.stdout.flush()
