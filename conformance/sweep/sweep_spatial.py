"""Sweep real SPATIAL AnnData (10x Visium + imaging-based) through read_anndata -> validate -> write.

lstar's spatial support is *conceptual* (deliberate, see misc/format_coverage.md Tier 3): `obsm['spatial']`
becomes a NAMED OBSERVED coordinate axis (subtype `spatial`) that round-trips back to `obsm['spatial']`;
`uns['spatial']` (tissue images, scalefactors, vendor metadata) stays in the LOSSLESS PASSTHROUGH
(aux/anndata.uns), NOT typed. This sweep validates that handling on REAL objects and -- importantly --
RECORDS what spatial structure is (correctly) deferred vs what silently disappears.

For each dataset it checks:
  1. `obsm['spatial']` -> an `embedding` field named `spatial`, subtype `spatial`, over an OBSERVED
     coordinate axis `spatial` (not a derived embedding axis).
  2. round-trip: write_anndata(read_anndata(a)) puts the coords back at `obsm['spatial']` (shape + values).
  3. `uns['spatial']` survives in the passthrough (ds.aux['anndata.uns']['spatial']) with its library
     keys/scalefactors -- and the image arrays under it (the deferred tier) are reported, not asserted.
  4. validate() is clean.

Datasets (cached local-only under testdata/spatial/, gitignored):
  - Visium .h5ad from sweep_spatial_fetch.py (breast/heart/lymph/kidney/brain, human+mouse, a
    multi-section pair, a targeted-vs-parent pair).
  - imaging-based squidpy datasets (merfish/seqfish/slideseqv2/imc) cached by sweep_spatial_fetch_sq.py
    (if the squidpy side-venv was available); these carry a targeted panel + molecule/coord tables.

Each dataset runs in its OWN subprocess so a big-read OOM/crash fails only that one. Writes
/tmp/sweep_spatial.tsv.

Run one in isolation:  python3 conformance/sweep/sweep_spatial.py testdata/spatial/<f>.h5ad
Run the whole sweep:   python3 conformance/sweep/sweep_spatial.py
"""
import warnings, sys, os, tempfile, glob, subprocess, json
warnings.filterwarnings("ignore")


def run_one(path):
    sys.path.insert(0, "python/src")
    import numpy as np, anndata as ad, lstar
    a = ad.read_h5ad(path)
    has_spatial_obsm = "spatial" in a.obsm
    uns_spatial = a.uns.get("spatial", None)
    # enumerate what's IN uns['spatial'] so we can report deferred (image) structure
    uns_spatial_keys, img_keys, scalefactor_keys = [], [], []
    if isinstance(uns_spatial, dict):
        uns_spatial_keys = list(uns_spatial.keys())
        for lib, sub in uns_spatial.items():
            if isinstance(sub, dict):
                if "images" in sub and isinstance(sub["images"], dict):
                    img_keys += ["%s/%s" % (lib, ik) for ik in sub["images"].keys()]
                if "scalefactors" in sub and isinstance(sub["scalefactors"], dict):
                    scalefactor_keys = list(sub["scalefactors"].keys())

    ds = lstar.read_anndata(a)
    errs = [i for i in lstar.validate(ds) if i.startswith("ERROR")]

    # (1) spatial coord axis: observed coordinate axis + embedding field with subtype 'spatial'
    sax = ds.axes.get("spatial")
    sfield = ds.fields.get("spatial")
    axis_observed = bool(sax is not None and getattr(sax, "origin", None) == lstar.OBSERVED
                         and getattr(sax, "role", None) == "coordinate")
    field_ok = bool(sfield is not None and getattr(sfield, "subtype", None) == "spatial"
                    and "spatial" in (sfield.span or []))

    # (2) round-trip: spatial coords land back at obsm['spatial']
    p = os.path.join(tempfile.mkdtemp(), "s.lstar.zarr")
    lstar.write(ds, p)
    rt_ok, rt_note = None, ""
    try:
        a2 = lstar.write_anndata(ds)
        if has_spatial_obsm:
            if "spatial" not in a2.obsm:
                rt_ok, rt_note = False, "spatial absent from obsm after round-trip"
            else:
                s0, s1 = np.asarray(a.obsm["spatial"]), np.asarray(a2.obsm["spatial"])
                rt_ok = (s0.shape == s1.shape) and bool(np.allclose(np.nan_to_num(s0),
                                                                    np.nan_to_num(s1)))
                rt_note = "shape %s->%s" % (s0.shape, s1.shape)
        else:
            rt_ok = True
    except Exception as e:
        rt_ok, rt_note = False, "write_anndata: " + str(e)[:60]

    # (3) uns['spatial'] survives in the passthrough
    aux_uns = ds.aux.get("anndata.uns", {})
    uns_passthrough = "spatial" in aux_uns
    # do the library/scalefactor keys survive verbatim?
    sf_survived = None
    if uns_passthrough and scalefactor_keys:
        kept = aux_uns.get("spatial", {})
        try:
            lib0 = next(iter(kept.values())) if isinstance(kept, dict) else {}
            sf_survived = set(scalefactor_keys).issubset(set((lib0 or {}).get("scalefactors", {}).keys()))
        except Exception:
            sf_survived = None

    return {
        "status": "PASS" if (not errs and axis_observed and field_ok and rt_ok and
                             (uns_passthrough or not isinstance(uns_spatial, dict))) else "CHECK",
        "shape": list(a.shape),
        "axis_observed": axis_observed, "field_ok": field_ok, "rt_ok": rt_ok,
        "uns_passthrough": uns_passthrough, "sf_survived": sf_survived,
        "uns_spatial_libs": uns_spatial_keys,
        "deferred_images": img_keys,            # NOT represented (deferred tier) -- recorded, not asserted
        "n_obsm": list(a.obsm.keys()),
        "fields": len(ds.fields), "axes": len(ds.axes),
        "dropped": ds.dropped[:6],
        "note": (errs[0][:70] if errs else rt_note),
    }


def main():
    files = sorted(glob.glob("testdata/spatial/*.h5ad"))
    if not files:
        print("no spatial datasets cached under testdata/spatial -- run sweep_spatial_fetch.py first")
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

    hdr = ("dataset\tstatus\tshape\taxis_obs\tfield_ok\trt_ok\tuns_pass\tsf_ok\t"
           "uns_libs\tdeferred_images\tdropped\tnote\n")
    with open("/tmp/sweep_spatial.tsv", "w") as fh:
        fh.write(hdr)
        for name, d in rows:
            fh.write("\t".join(map(str, [
                name, d.get("status"), d.get("shape", ""), d.get("axis_observed", ""),
                d.get("field_ok", ""), d.get("rt_ok", ""), d.get("uns_passthrough", ""),
                d.get("sf_survived", ""), "+".join(d.get("uns_spatial_libs", [])),
                "+".join(d.get("deferred_images", [])), ";".join(d.get("dropped", [])),
                d.get("note", "")])) + "\n")
    npass = sum(1 for _, d in rows if d.get("status") == "PASS")
    print(f"Spatial sweep: {npass} PASS / {len(rows)} datasets")
    for name, d in rows:
        print(f"  {d.get('status'):7s} {name:46s} {str(d.get('shape','')):14s} "
              f"axisObs={d.get('axis_observed')} field={d.get('field_ok')} "
              f"rt={d.get('rt_ok')} unsPass={d.get('uns_passthrough')} "
              f"imgs={len(d.get('deferred_images', []))} {str(d.get('note',''))[:38]}")


if __name__ == "__main__":
    if len(sys.argv) > 1:
        print(json.dumps(run_one(sys.argv[1])))
    else:
        main()
