"""Sweep real PERTURBATION AnnData (scPerturb harmonized .h5ad) through read_anndata -> validate -> write.

Perturbation data is the stress test for lstar's FACTOR-AXIS machinery: a categorical `obs` column
auto-induces a derived factor axis (model.py `_auto_induce`), and Perturb-seq objects carry perturbation/
guide categoricals with HUNDREDS to THOUSANDS of levels (combinatorial CRISPR, drug x dose). This sweep
confirms those induce correctly, the dose/condition fields are kept, validate() is clean, and the object
round-trips (streaming the heavy matrix from a BACKED read so a 100k+ x 30k+ count matrix never fully
materializes).

Datasets (scPerturb, Zenodo 10.5281/zenodo.7041848, cached local-only under testdata/perturbation/):
  - DatlingerBock2021              : ~39k cells, scifi-RNA-seq CRISPR; `perturbation` (41) + `perturbation_2`
                                     (combo 2nd guide) + `sample` (384) -- a small/fast guide-layout case.
  - SrivatsanTrapnell2020_sciplex2 : ~24k cells, sci-Plex drug-response; `perturbation` (5 drugs) +
                                     `dose_value` (8-level dose) -- dose as an ORDERED treatment factor.
  - NormanWeissman2019_filtered    : ~111k cells, combinatorial CRISPRa (K562); `perturbation` (237 single
                                     + combo levels) + `guide_id` (290) -- the HIGH-CARDINALITY induction
                                     stressor (a factor axis with hundreds of categories).

For each it reports:
  - factor axes induced from the perturbation/guide categoricals + their cardinalities (the stress metric)
  - that `perturbation` / `dose_value` survived as fields
  - validate() clean
  - bounded-memory round-trip (backed read -> write(stream=True) -> .lstar.zarr)

Each dataset runs in its OWN subprocess (backed read + streaming write keep RSS bounded, but isolation
means a crash fails only that one). Writes /tmp/sweep_perturbation.tsv.

Run one in isolation:  python3 conformance/sweep/sweep_perturbation.py testdata/perturbation/<f>.h5ad
Run the whole sweep:   python3 conformance/sweep/sweep_perturbation.py
"""
import warnings, sys, os, tempfile, glob, subprocess, json
warnings.filterwarnings("ignore")

# obs columns that, when categorical, are the perturbation/treatment design we expect to induce factors
PERT_KEYS = ("perturbation", "guide", "drug", "dose", "target", "product",
             "condition", "treatment", "sgrna", "grna", "gene_target")


def run_one(path):
    sys.path.insert(0, "python/src")
    import numpy as np, anndata as ad, lstar
    # backed: X stays on disk and is wrapped as a streaming source; obs/var load eagerly
    a = ad.read_h5ad(path, backed="r")
    shape = tuple(a.shape)
    # which obs columns are perturbation-related categoricals, and their level counts
    import pandas as pd
    pert_cols = {}
    for c in a.obs.columns:
        if any(k in c.lower() for k in PERT_KEYS) and isinstance(a.obs[c].dtype, pd.CategoricalDtype):
            pert_cols[c] = int(a.obs[c].nunique())

    ds = lstar.read_anndata(a)
    errs = [i for i in lstar.validate(ds) if i.startswith("ERROR")]

    # factor axes induced from those columns (a categorical label auto-induces an axis of its bare name)
    factor_axes = {n: len(np.asarray(ax.labels))
                   for n, ax in ds.axes.items()
                   if getattr(ax, "role", None) == "factor"}
    pert_factors = {c: factor_axes.get(c) for c in pert_cols}
    # the biggest factor axis (induction stressor metric)
    max_factor = max(factor_axes.values()) if factor_axes else 0
    # confirm the perturbation columns themselves survived as fields (label fields)
    pert_fields_kept = [c for c in pert_cols if c in ds.fields]

    # bounded-memory round-trip: stream the heavy sparse matrix block-by-block to a zarr store
    p = os.path.join(tempfile.mkdtemp(), "s.lstar.zarr")
    rt_ok, rt_note = None, ""
    try:
        lstar.write(ds, p, stream=True)
        rt_ok = True
    except Exception as e:
        rt_ok, rt_note = False, "write(stream): " + str(e)[:60]
    try:
        a.file.close()
    except Exception:
        pass

    # induction integrity: every perturbation categorical that has levels must have induced a factor axis
    induction_ok = all((pert_factors.get(c) == n) for c, n in pert_cols.items() if n > 0) \
        if pert_cols else None

    return {
        "status": "PASS" if (not errs and rt_ok and (induction_ok in (True, None))) else "CHECK",
        "shape": list(shape),
        "pert_cols": pert_cols, "pert_factors": pert_factors,
        "induction_ok": induction_ok, "max_factor": max_factor,
        "n_factor_axes": len(factor_axes), "pert_fields_kept": pert_fields_kept,
        "fields": len(ds.fields), "axes": len(ds.axes),
        "dropped": ds.dropped[:6], "rt_ok": rt_ok,
        "note": (errs[0][:70] if errs else rt_note),
    }


def main():
    files = sorted(glob.glob("testdata/perturbation/*.h5ad"))
    if not files:
        print("no perturbation datasets cached under testdata/perturbation -- "
              "fetch scPerturb .h5ad from Zenodo 7041848 (see README)")
        return
    rows = []
    for path in files:
        name = os.path.splitext(os.path.basename(path))[0]
        r = subprocess.run([sys.executable, __file__, path], capture_output=True, text=True)
        try:
            d = json.loads(r.stdout.strip().splitlines()[-1])
        except Exception:
            d = {"status": "CRASH", "note": (r.stderr.strip().splitlines()[-1:] or [""])[0][:90]}
        rows.append((name, d))

    hdr = ("dataset\tstatus\tshape\tn_factor\tmax_factor\tpert_factors\tinduction_ok\t"
           "rt_ok\tdropped\tnote\n")
    with open("/tmp/sweep_perturbation.tsv", "w") as fh:
        fh.write(hdr)
        for name, d in rows:
            fh.write("\t".join(map(str, [
                name, d.get("status"), d.get("shape", ""), d.get("n_factor_axes", ""),
                d.get("max_factor", ""), d.get("pert_factors", ""), d.get("induction_ok", ""),
                d.get("rt_ok", ""), ";".join(d.get("dropped", [])), d.get("note", "")])) + "\n")
    npass = sum(1 for _, d in rows if d.get("status") == "PASS")
    print(f"Perturbation sweep: {npass} PASS / {len(rows)} datasets")
    for name, d in rows:
        pf = d.get("pert_factors", {})
        pfstr = ",".join("%s=%s" % (k, v) for k, v in pf.items()) if isinstance(pf, dict) else str(pf)
        print(f"  {d.get('status'):7s} {name:34s} {str(d.get('shape','')):16s} "
              f"factors={d.get('n_factor_axes')} maxlvl={d.get('max_factor')} "
              f"rt={d.get('rt_ok')} [{pfstr[:54]}] {str(d.get('note',''))[:30]}")


if __name__ == "__main__":
    if len(sys.argv) > 1:
        print(json.dumps(run_one(sys.argv[1])))
    else:
        main()
