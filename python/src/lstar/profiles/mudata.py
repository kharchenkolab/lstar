"""MuData profile (scverse multimodal): read_mudata / write_mudata.

A MuData (`.h5mu`) holds several **modalities** (RNA, ADT/protein, ATAC/peaks) over a **shared cell
set**. That is exactly L*'s "multiple feature spaces, one observation axis" shape -- the *same*
representation the Seurat multi-assay profile produces -- so a CITE-seq object round-trips MuData <-> L*
<-> Seurat. Each modality becomes a **feature axis** (canonical `genes`/`proteins`/`peaks`) with its
measures over `(cells, <feature-axis>)`; per-modality `var` -> feature fields; global `obs` -> cell
fields (categoricals induce factor axes); global + per-modality `obsm` -> embeddings; `uns` -> the
lossless passthrough. A modality covering only a subset of cells (partial overlap, `obsmap` with 0s)
gets its own `cells.<mod>` axis (faithful, not zero-filled).

mudata is imported lazily, so `import lstar` never requires it.
"""
import numpy as np

from ..model import Dataset, OBSERVED, DERIVED
from .anndata import _by_dtype_series, _coord_axis, _uniq, _to_pandas_col

PROFILE = "mudata@0.1"

# canonical feature-axis names (the shared multimodal vocabulary) -- so a modality lands on the same
# axis regardless of source format / capitalisation.
_MODALITY = {"rna": "genes", "gex": "genes", "gene": "genes", "genes": "genes",
             "prot": "proteins", "adt": "proteins", "protein": "proteins", "proteins": "proteins",
             "antibody": "proteins", "atac": "peaks", "peak": "peaks", "peaks": "peaks",
             "chromatin": "peaks"}
_AXIS_MODALITY = {"genes": "rna", "proteins": "prot", "peaks": "atac"}  # canonical axis -> mudata mod name


def _modality_axis(m, taken):
    fax = _MODALITY.get(str(m).lower(), str(m).lower())
    return fax if fax not in taken else str(m).lower()


def _mudata_version():
    try:
        import mudata
        return "mudata@%s" % mudata.__version__
    except Exception:
        return "mudata@?"


def read_mudata(md, kind="sample"):
    """Read a MuData object into an L* Dataset (one shared `cells` axis + a feature axis per modality)."""
    import scipy.sparse as sp

    ds = Dataset(kind=kind)
    ds.profiles = [PROFILE, _mudata_version()]
    cells = np.asarray(md.obs_names, dtype=str)
    ds.add_axis("cells", cells, origin=OBSERVED, role="observation")
    cell_pos = {c: i for i, c in enumerate(cells)}

    for m in md.mod:
        am = md.mod[m]
        fax = _modality_axis(m, ds.axes)
        ds.add_axis(fax, np.asarray(am.var_names, dtype=str), origin=OBSERVED, role="feature")

        mcells = np.asarray(am.obs_names, dtype=str)
        aligned = len(mcells) == len(cells) and np.array_equal(mcells, cells)
        cax = "cells"; midx = None; iax = None
        if not aligned:
            if all(c in cell_pos for c in mcells):   # a cell subset of the shared union -> partial coverage
                midx = np.array([cell_pos[c] for c in mcells], dtype=np.int64); iax = "cells"
            else:                                    # genuinely disjoint cells -> its own observation axis
                cax = _uniq(ds, "cells.%s" % m, "obs")
                ds.add_axis(cax, mcells, origin=OBSERVED, role="observation")

        def _add_measure(name, X, state):
            ds.add_field(name, X, role="measure", span=[cax, fax], state=state, index=midx, index_axis=iax,
                         provenance={"mudata": "%s/%s" % (m, name.split(".")[-1])})
        _add_measure(_uniq(ds, m, "mod"), am.X, _guess_state_mod(am))
        for lk in list(am.layers.keys()):
            from .anndata import _infer_state
            _add_measure(_uniq(ds, "%s.%s" % (m, lk), "mod"), am.layers[lk],
                         _infer_state(am.layers[lk], name=str(lk)))

        for col in am.var.columns:                   # per-modality feature metadata over the feature axis
            vals, role, mask = _by_dtype_series(am.var[col])
            ds.add_field(_uniq(ds, "%s.%s" % (m, str(col)), fax), vals, role=role, span=[fax], mask=mask,
                         provenance={"mudata": "%s/var/%s" % (m, col)})
        for k in list(am.obsm.keys()):               # per-modality embeddings (own PCA/UMAP)
            v = np.asarray(am.obsm[k])
            if v.ndim != 2:
                continue
            cname = _uniq(ds, "%s_%s" % (m, k[2:] if str(k).startswith("X_") else k), "coord")
            _coord_axis(ds, cname, v.shape[1])
            ds.add_field(_uniq(ds, "%s_%s" % (m, k), fax), v, role="embedding", span=[cax, cname],
                         index=midx, index_axis=iax, provenance={"mudata": "%s/obsm/%s" % (m, k)})
        for k in list(am.varm.keys()):               # per-modality loadings (MOFA LFs, totalVI, ...): a
            v = np.asarray(am.varm[k])               # loading over (feature axis, factor coord). When the
            if v.ndim != 2:                          # factor coord is shared (same stripped name as the
                continue                             # global factor scores' obsm), they land on ONE axis.
            base = k[2:] if str(k).startswith("X_") else str(k)
            # Reuse a same-named coordinate axis only when its length matches (the genuine shared-factor
            # case: MOFA/totalVI factors load from several modalities onto ONE axis). If a same-named axis
            # exists with a DIFFERENT length, this is a per-modality reduction (e.g. each modality's own
            # `PCs`: rna 50 comps, prot 31) -- namespace it so the two don't collide on a length-mismatched
            # axis (else the loadings span a wrong-length axis -> validate error). Caught on real minipbcite.
            if base in ds.axes and len(ds.axes[base]) != v.shape[1]:
                base = "%s_%s" % (m, base)
            cname = _coord_axis(ds, base, v.shape[1])    # reuses the factor axis if scores already made it
            ds.add_field(_uniq(ds, "%s_%s_loadings" % (m, k), fax), v, role="loading", span=[fax, cname],
                         provenance={"mudata": "%s/varm/%s" % (m, k)})
        uns = {k: v for k, v in am.uns.items() if not str(k).startswith("lstar/")}
        if uns:
            ds.aux["mudata.%s.uns" % m] = uns

    # global (MuData-level) obs / obsm / uns over the shared cells axis
    for col in md.obs.columns:
        vals, role, mask = _by_dtype_series(md.obs[col])
        ds.add_field(_uniq(ds, str(col), "cells"), vals, role=role, span=["cells"], mask=mask,
                     provenance={"mudata": "obs/%s" % col})
    for k in list(md.obsm.keys()):
        v = np.asarray(md.obsm[k])
        if v.ndim != 2:
            continue
        cname = _uniq(ds, k[2:] if str(k).startswith("X_") else str(k), "coord")
        _coord_axis(ds, cname, v.shape[1])
        ds.add_field(_uniq(ds, str(k), "cells"), v, role="embedding", span=["cells", cname],
                     provenance={"mudata": "obsm/%s" % k})
    for k in list(getattr(md, "obsp", {}) or {}):    # global joint graph (WNN connectivities) -> relation
        mm = md.obsp[k]
        if sp.issparse(mm) and mm.shape == (len(cells), len(cells)):
            ds.add_field(_uniq(ds, str(k), "cells"), mm, role="relation", span=["cells", "cells"],
                         directed=False, weighted=True, provenance={"mudata": "obsp/%s" % k})
    # Facet-set provenance on joint products (WNN/MOFA/totalVI): a joint embedding records the input
    # **feature axes** that fed it (`provenance["input_axes"]`). For factor models we can infer it -- the
    # contributing feature axes are exactly those carrying a `loading` over the embedding's coordinate
    # axis (MOFA: genes + proteins both load on `factors`). A producer that knows its facets (pagoda2
    # records `facets=c("RNA","ADT")` on a reduction) sets the same key; it round-trips (provenance is
    # preserved across Py/C++/R).
    for nm, f in list(ds.fields.items()):
        if f.role != "embedding" or len(f.span or []) != 2:
            continue
        coord = f.span[1]
        feats = sorted({g.span[0] for g in ds.fields.values()
                        if g.role == "loading" and len(g.span or []) == 2 and g.span[1] == coord
                        and ds.axes.get(g.span[0]) is not None and ds.axes[g.span[0]].role == "feature"})
        if feats:
            f.provenance = dict(f.provenance or {}, input_axes=feats)
    guns = {k: v for k, v in md.uns.items() if not str(k).startswith("lstar/")}
    if guns:
        ds.aux["mudata.uns"] = guns
    return ds


def _guess_state_mod(am):
    # content-based (shared with the anndata reader): neg -> scaled, non-neg int -> raw, else lognorm.
    # (Was a `max>30 & integer` heuristic that mislabeled low-count raw as lognorm and missed scaled.)
    from .anndata import _infer_state
    try:
        return _infer_state(am.X)
    except Exception:
        return None


def write_mudata(ds):
    """Write an L* Dataset back to a MuData object: one AnnData per feature axis (modality)."""
    import anndata as ad
    import mudata
    import numpy as np
    import pandas as pd

    cells = np.asarray(ds.axis("cells").labels, dtype=str)
    feat_axes = [n for n, a in ds.axes.items() if a.role == "feature"]
    mods = {}
    for fax in feat_axes:
        feats = np.asarray(ds.axis(fax).labels, dtype=str)
        mname = _AXIS_MODALITY.get(fax, fax)
        # measures over (cells, fax) -> this modality's X / layers (field name "<mod>" or "<mod>.<layer>")
        meas = [(n, f) for n, f in ds.fields.items()
                if f.role == "measure" and f.span and len(f.span) == 2 and f.span[1] == fax
                and f.span[0] == "cells"]
        if not meas:
            continue
        prov_mod = (meas[0][1].provenance or {}).get("mudata", "/").split("/")[0]
        primary = next((f for n, f in meas if f.state == "raw"), meas[0][1])
        # partial coverage: this modality was measured on only a subset of the shared cells (cells[index])
        midx = getattr(primary, "index", None)
        mod_cells = cells[np.asarray(midx)] if midx is not None else cells
        layers = {}
        for n, f in meas:
            if f is primary:
                continue
            lk = n.split(".", 1)[1] if "." in n else n
            layers[lk] = f.values
        var = {}
        for n, f in ds.fields.items():
            if f.span == [fax] and f.role in ("label", "measure"):
                key = n.split(".", 1)[1] if n.startswith(fax + ".") or "." in n else n
                var[key] = _to_pandas_col(f)
        am = ad.AnnData(X=primary.values, var=pd.DataFrame(var, index=feats) if var else pd.DataFrame(index=feats),
                        obs=pd.DataFrame(index=mod_cells), layers=layers or None)
        mod = prov_mod or mname
        for n, f in ds.fields.items():               # per-mod embeddings -> obsm, per-mod loadings -> varm
            prov = str((f.provenance or {}).get("mudata", ""))
            if f.role == "embedding" and prov.startswith(mod + "/obsm/"):
                am.obsm[prov.split("/", 2)[2]] = np.asarray(f.values)
            elif f.role == "loading" and prov.startswith(mod + "/varm/"):
                am.varm[prov.split("/", 2)[2]] = np.asarray(f.values)
        mods[mod] = am

    md = mudata.MuData(mods)
    # global products: obs (cell fields, no modality provenance), obsm (joint embeddings), obsp (joint graphs)
    gobs = {}
    for n, f in ds.fields.items():
        prov = str((f.provenance or {}).get("mudata", ""))
        if f.span == ["cells"] and not prov.count("/") and f.role in ("label", "measure"):
            gobs[n] = _to_pandas_col(f)
        elif f.role == "embedding" and prov.startswith("obsm/"):
            md.obsm[prov.split("/", 1)[1]] = np.asarray(f.values)
        elif f.role == "relation" and f.span == ["cells", "cells"] and prov.startswith("obsp/"):
            md.obsp[prov.split("/", 1)[1]] = f.values
    if gobs:
        gdf = pd.DataFrame(gobs, index=cells)
        for c in gdf.columns:
            md.obs[c] = gdf[c].values
    return md


def convert_h5mu(h5mu_path, store_path, **write_kwargs):
    """Convert a `.h5mu` to an L* store."""
    import mudata

    from ..zarr_io import write as _write
    return _write(read_mudata(mudata.read_h5mu(h5mu_path)), store_path, **write_kwargs)
