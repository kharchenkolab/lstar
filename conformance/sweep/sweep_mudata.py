"""Sweep real multi-modal MuData (`.h5mu`) through read_mudata -> write -> write_mudata.

Datasets (all real, cached local-only under testdata/):
  - minipbcite.h5mu            : 411-cell CITE-seq (RNA+ADT), annotated (celltype, per-mod PCA)
  - 5k_pbmc_protein.h5mu       : 10x 5k PBMC CITE-seq (RNA 33538 + ADT 32), built from the public .h5
  - pbmc_multiome_3k.h5mu      : 10x PBMC multiome (RNA 36601 + ATAC 98319 peaks), built from .h5
  - (any other *.h5mu under testdata/)

The MuData profile turns each modality into a feature axis (canonical genes/proteins/peaks), so this
exercises >1 feature axis over one shared cells axis -- the multimodal L* shape -- and the ATAC modality
exercises a *peaks* feature axis arriving from MuData (not Signac). Writes /tmp/sweep_mudata.tsv.

Run a single dataset in isolation:  python3 conformance/sweep/sweep_mudata.py <path.h5mu>
Run the whole sweep (subprocess-isolated, so a big read OOM/crash only fails that one): no args.
"""
import warnings, sys, os, tempfile, glob, subprocess, json
warnings.filterwarnings("ignore")


def run_one(path):
    sys.path.insert(0, "python/src")
    import mudata
    import lstar
    md = mudata.read_h5mu(path)
    mods = {k: tuple(v.shape) for k, v in md.mod.items()}
    ds = lstar.read_mudata(md)
    errs = [i for i in lstar.validate(ds) if i.startswith("ERROR")]
    feat_axes = [ax.name for ax in ds.axes.values()
                 if ax.name in ("genes", "proteins", "peaks")]
    p = os.path.join(tempfile.mkdtemp(), "s.lstar.zarr")
    lstar.write(ds, p)
    md2 = lstar.write_mudata(ds)                 # round-trip back to MuData
    ok_rt = set(md2.mod) == set(md.mod) and all(
        md2.mod[k].shape == md.mod[k].shape for k in md.mod)
    return {
        "status": "PASS" if (not errs and ok_rt) else ("VALIDATE-ERR" if errs else "ROUNDTRIP-MISMATCH"),
        "shape": list(md.shape), "mods": mods, "fields": len(ds.fields),
        "axes": len(ds.axes), "feat_axes": feat_axes,
        "dropped": ds.dropped[:4], "note": (errs[0][:70] if errs else ""),
    }


def main():
    files = sorted(glob.glob("testdata/minipbcite.h5mu") +
                   glob.glob("testdata/mudata_examples/*.h5mu"))
    if not files:
        print("no .h5mu cached; nothing to sweep")
        return
    rows = []
    for path in files:
        name = os.path.splitext(os.path.basename(path))[0]
        # isolate each dataset in its own subprocess (big reads can OOM the parent)
        r = subprocess.run([sys.executable, __file__, path], capture_output=True, text=True)
        try:
            d = json.loads(r.stdout.strip().splitlines()[-1])
        except Exception:
            d = {"status": "CRASH", "note": (r.stderr.strip().splitlines()[-1:] or [""])[0][:80]}
        rows.append((name, d))

    hdr = "dataset\tstatus\tshape\tmods\tfields\taxes\tfeat_axes\tnote\n"
    with open("/tmp/sweep_mudata.tsv", "w") as fh:
        fh.write(hdr)
        for name, d in rows:
            fh.write("\t".join(map(str, [
                name, d.get("status"), d.get("shape", ""),
                d.get("mods", ""), d.get("fields", ""), d.get("axes", ""),
                "+".join(d.get("feat_axes", [])), d.get("note", "")])) + "\n")
    npass = sum(1 for _, d in rows if d.get("status") == "PASS")
    print(f"MuData sweep: {npass} PASS / {len(rows)} datasets")
    for name, d in rows:
        mods = d.get("mods", {})
        modstr = "+".join(f"{k}{v}" for k, v in mods.items()) if isinstance(mods, dict) else str(mods)
        print(f"  {d.get('status'):16s} {name:24s} {modstr[:46]:46s} "
              f"feat_axes={'+'.join(d.get('feat_axes', []))} {str(d.get('note',''))[:40]}")


if __name__ == "__main__":
    if len(sys.argv) > 1:
        print(json.dumps(run_one(sys.argv[1])))     # single-dataset worker: emit one JSON line
    else:
        main()
