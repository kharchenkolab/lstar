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


def pseudobulk(ds, factor, field="counts", lognorm=False, add=True):
    """Per-group **pseudobulk** over `(factor, genes)`: each group's mean expression and fraction of
    cells expressing, from a measure (`field`, cells x genes) grouped by the per-cell factor label. The
    symmetric companion to a DE bundle -- ordinary `pb.<factor>.<stat>` measures over the same factor
    axis. With `add=True` the fields are added to `ds`; returns the `{stat: array}` dict either way.
    """
    import scipy.sparse as sp

    from .model import as_categorical

    cat = as_categorical(ds.field(factor).values)
    codes = np.asarray(cat.codes)
    groups = np.asarray(cat.categories, dtype=str)
    M = ds.field(field).values
    M = M.tocsr() if sp.issparse(M) else sp.csr_matrix(M)
    gene_axis = ds.field(field).span[1]
    K, ng = len(groups), M.shape[1]
    mean = np.zeros((K, ng)); frac = np.zeros((K, ng))
    for k in range(K):
        rows = np.nonzero(codes == k)[0]
        if not len(rows):
            continue
        sub = M[rows]
        vals = sub.copy()
        if lognorm:
            vals.data = np.log1p(vals.data)
        mean[k] = np.asarray(vals.sum(axis=0)).ravel() / len(rows)       # mean (log1p) expression
        frac[k] = np.asarray((sub > 0).sum(axis=0)).ravel() / len(rows)  # fraction expressing
    out = {"mean": mean, "frac": frac}
    if add:
        for stat, arr in out.items():
            ds.add_field("pb.%s.%s" % (factor, stat), arr, role="measure", span=[factor, gene_axis],
                         subtype="pseudobulk", state=("lognorm" if lognorm else None),
                         provenance={"pb_factor": factor, "pb_stat": stat, "from": field})
    return out


def collection_pseudobulk(ds, factor, field="counts", lognorm=False, add=True, union_axis="cells"):
    """Pseudobulk over a **collection**: aggregate the per-sample `<field>.<s>` measures into
    `(factor-group x genes)`, grouped by `factor` — a categorical over the **union** `cells` axis (a joint
    clustering over a Conos-style integration) — **streamed one sample at a time** so the joint matrix is
    never materialized, with **float64 accumulation** over float32 storage. The collection generalization
    of `pseudobulk()`: it walks the per-sample measures and accumulates into the induced factor axis.

    Genes are aligned by label across samples (a gene absent from a sample contributes 0 there); the
    result is `pb.<factor>.{mean,frac}` measures over `(factor, genes)`. Returns the `{stat: array}` dict.
    """
    import scipy.sparse as sp

    from .model import as_categorical

    cat = as_categorical(ds.field(factor).values)
    groups = np.asarray(cat.categories, dtype=str); K = len(groups)
    union_codes = np.asarray(cat.codes)                              # group index per union cell (-1 = none)
    cell_pos = {c: i for i, c in enumerate(np.asarray(ds.axis(union_axis).labels, dtype=str))}

    persamp = [(nm, f) for nm, f in ds.fields.items()               # `<field>.<sample>` over (cells.<s>, genes)
               if f.role == "measure" and nm.startswith(field + ".")
               and f.span and len(f.span) == 2 and f.span[0].startswith("cells.")]
    if not persamp:
        raise ValueError("no per-sample '%s.<sample>' measures over (cells.<sample>, genes)" % field)

    gene_index, genes = {}, []                                       # canonical gene set (union, first-seen)
    for _nm, f in persamp:
        for g in np.asarray(ds.axis(f.span[1]).labels, dtype=str):
            if g not in gene_index:
                gene_index[g] = len(genes); genes.append(g)
    ng = len(genes)

    sum64 = np.zeros((K, ng), dtype=np.float64)                      # float64 accumulation
    nz = np.zeros((K, ng), dtype=np.float64)
    cnt = np.zeros(K, dtype=np.int64)
    for _nm, f in persamp:                                           # bounded memory: one sample at a time
        M = f.values; M = M.tocsr() if sp.issparse(M) else sp.csr_matrix(M)
        if lognorm:
            M = M.copy(); M.data = np.log1p(M.data)
        gcols = np.array([gene_index[g] for g in np.asarray(ds.axis(f.span[1]).labels, dtype=str)])
        rows_union = np.array([cell_pos[c] for c in np.asarray(ds.axis(f.span[0]).labels, dtype=str)])
        grp = union_codes[rows_union]                               # group of each sample row
        for k in range(K):
            rows = np.nonzero(grp == k)[0]
            if not len(rows):
                continue
            sub = M[rows]
            sum64[k, gcols] += np.asarray(sub.sum(axis=0)).ravel()
            nz[k, gcols] += np.asarray((sub > 0).sum(axis=0)).ravel()
            cnt[k] += len(rows)
    denom = np.where(cnt > 0, cnt, 1)[:, None]
    out = {"mean": sum64 / denom, "frac": nz / denom}
    if "genes" not in ds.axes:
        ds.add_axis("genes", genes, origin="derived", role="feature")
    if add:
        for stat, arr in out.items():
            ds.add_field("pb.%s.%s" % (factor, stat), arr, role="measure", span=[factor, "genes"],
                         subtype="pseudobulk", state=("lognorm" if lognorm else None),
                         provenance={"pb_factor": factor, "pb_stat": stat, "from": field, "collection": True})
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
