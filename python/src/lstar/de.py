"""Differential-expression / marker results as ordinary fields over a **factor axis**.

A DE result is a *bundle* of optional measures over `(factor, gene-axis)` in one canonical orientation
(rows = factor groups, cols = genes), one field per statistic (`score`/`lfc`/`pval`/`padj`), sharing the
gene axis. Because the factor axis is induced from the clustering label (see `model.md` induction), a
per-group result needs no special slot -- it is just a measure, and `markers()` gives the tidy view.

Field naming: ``de.<factor>.<stat>`` (role=measure, span=[<factor>, <gene-axis>], subtype="de"). The
inverse-orientation question is settled: **factor-first**, so DE/pseudobulk/PAGA all share it.
"""
import numpy as np

DE_STATS = ("score", "lfc", "pval", "padj")          # the recognized statistics (all optional)


def de_field_name(factor, stat):
    return "de.%s.%s" % (factor, stat)


def de_factors(ds):
    """The set of factors that carry a DE bundle in `ds`."""
    out = []
    for f in ds.fields.values():
        if f.subtype == "de" and (f.provenance or {}).get("de_factor"):
            fac = f.provenance["de_factor"]
            if fac not in out:
                out.append(fac)
    return out


def de_bundle(ds, factor):
    """The `{stat: Field}` DE bundle for `factor` (empty dict if none)."""
    out = {}
    for f in ds.fields.values():
        if f.subtype == "de" and (f.provenance or {}).get("de_factor") == factor:
            st = f.provenance.get("de_stat")
            if st:
                out[st] = f
    return out


def markers(ds, factor, top=None, sort_by="score", descending=True):
    """A tidy long-form marker table for `factor`: one row per (group, gene) with whichever of
    `score`/`lfc`/`pval`/`padj` the bundle carries. `top` keeps the top-N genes per group (by
    `sort_by`). Returns a pandas DataFrame. The ergonomic counterpart to the (factor, genes) bundle.
    """
    import pandas as pd

    bundle = de_bundle(ds, factor)
    if not bundle:
        raise ValueError("no DE bundle for factor %r (have: %s)" % (factor, de_factors(ds)))
    any_field = next(iter(bundle.values()))
    gene_axis = any_field.span[1]
    groups = np.asarray(ds.axis(factor).labels, dtype=str)
    genes = np.asarray(ds.axis(gene_axis).labels, dtype=str)
    cols = {st: np.asarray(f.values, dtype=float) for st, f in bundle.items()}

    frames = []
    for gi, grp in enumerate(groups):
        row = {st: v[gi] for st, v in cols.items()}
        keep = ~np.isnan(next(iter(row.values())))             # genes present for this group
        idx = np.nonzero(keep)[0]
        if sort_by in row:
            order = np.argsort(row[sort_by][idx])
            if descending:
                order = order[::-1]
            idx = idx[order]
        if top is not None:
            idx = idx[:top]
        d = {"group": np.repeat(grp, len(idx)), "gene": genes[idx]}
        for st, v in row.items():
            d[st] = v[idx]
        frames.append(pd.DataFrame(d))
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(
        columns=["group", "gene", *bundle.keys()])
