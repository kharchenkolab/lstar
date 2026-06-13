import warnings, sys, os, tempfile, glob; warnings.filterwarnings("ignore")
sys.path.insert(0, "python/src")
import scanpy as sc, anndata as ad, lstar
sc.settings.datasetdir = "testdata"; sc.settings.verbosity = 0
rows = []
# scanpy's REAL bundled/downloadable datasets + local atlases
jobs = [("pbmc3k_processed", lambda: sc.datasets.pbmc3k_processed()),
        ("pbmc68k_reduced", lambda: sc.datasets.pbmc68k_reduced()),
        ("paul15", lambda: sc.datasets.paul15()),
        ("moignard15", lambda: sc.datasets.moignard15()),
        ("burczynski06", lambda: sc.datasets.burczynski06()),
        ("TMS_Marrow[local]", lambda: ad.read_h5ad("/home/pkharchenko/cacoa/age/tab.muris/tabula-muris-senis-droplet-processed-official-annotations-Marrow.h5ad", backed="r")),
        ("micropatterns[local]", lambda: ad.read_h5ad("/home/pkharchenko/igor/micropatterns/luis/micropatterns_reanalysis.h5ad", backed="r"))]
# real scVelo RNA-velocity / trajectory datasets cached under testdata/velocity (spliced/unspliced layers)
for vp in sorted(glob.glob("testdata/velocity/*.h5ad")):
    nm = "velocity:" + os.path.splitext(os.path.basename(vp))[0]
    jobs.append((nm, (lambda p=vp: ad.read_h5ad(p))))
# real 10x CITE-seq read as a single AnnData (genes+proteins mixed in var via feature_types) -- the raw
# 10x form (the MuData-split form is in sweep_mudata.py); checks the mixed-feature AnnData path.
for cp in sorted(glob.glob("testdata/citeseq_10x/*.h5")):
    nm = "10xcite:" + os.path.splitext(os.path.basename(cp))[0]
    jobs.append((nm, (lambda p=cp: (lambda a: (a.var_names_make_unique() or a))(sc.read_10x_h5(p, gex_only=False)))))
for name, load in jobs:
    try:
        a = load(); ds = lstar.read_anndata(a)
        errs = [i for i in lstar.validate(ds) if i.startswith("ERROR")]
        p = os.path.join(tempfile.mkdtemp(), "s.lstar.zarr")
        if not a.isbacked: lstar.write(ds, p)
        rows.append((name, str(a.shape), "PASS" if not errs else "VALIDATE-ERR", len(ds.fields), len(ds.axes), ";".join(ds.dropped[:2]) or (errs[0][:60] if errs else "")))
    except Exception as e:
        rows.append((name, "", "FAIL", "", "", str(e)[:80]))
open("/tmp/sweep_anndata.tsv","w").write("dataset\tshape\tstatus\tfields\taxes\tnote\n" + "\n".join("\t".join(map(str,r)) for r in rows))
npass = sum(1 for r in rows if r[2]=="PASS")
print(f"AnnData sweep: {npass} PASS / {len(rows)} datasets")
for r in rows: print(f"  {r[2]:12s} {r[0]:30s} {r[1]:14s} f={r[3]} a={r[4]} {r[5][:50]}")
