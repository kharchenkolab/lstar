---
name: lstar-format-conversion
description: Convert / interchange single-cell & spatial omics data between formats with the lstar package — AnnData (.h5ad), Seurat (.rds), SingleCellExperiment, Conos, pagoda2, and .lstar.zarr stores — via the one-command `lstar convert` CLI or the in-language read_X/write_X profiles, bridged by lstar's uniform L* data model. Also covers heterogeneous-sample COLLECTIONS (Conos), what each conversion preserves vs. drops, version-robust format detection, and lazy/streamed per-gene stats.
when_to_use: Use to move one single-cell object or file from one toolkit/format to another (AnnData<->Seurat<->SCE<->Conos<->pagoda2<->zarr), to round-trip through .lstar.zarr, to flatten/export a multi-sample collection (e.g. a Conos integration) to Seurat v5 or AnnData, or to compute per-gene/per-cell stats over a large store without loading it. lstar is the "glue", not an analysis package.
avoid_when: Do not use to ANALYZE data (cluster, find markers, integrate) — that is pagoda2 / conos / Seurat / scanpy. lstar only moves data between them. For a Conos integration itself, use the conos recipe; for single-dataset analysis, the pagoda2 recipe.
invocation: interactive+batch
requires_tools: [run_r, run_python, Bash]
capabilities_needed: [lstar]
keywords: [lstar, L*, convert, conversion, interchange, glue, h5ad, AnnData, Seurat, SingleCellExperiment, SCE, Conos, pagoda2, zarr, collection, export, import, format detection, read_anndata, write_anndata, read_seurat, write_seurat, write_conos, lstar convert, backend, package-free, lazy, streaming]
produces: [converted.h5ad, converted.rds, store.lstar.zarr, fidelity_report.json]
domain: genomics
source: "lstar (private repo ~/p21/lstar) — source-verified R profiles (R/R/profile_*.R), Python CLI (python/src/lstar/cli.py), and conformance/."
---

# Converting single-cell data between formats with lstar

lstar is a lightweight "glue" library: a uniform data model (**L\***) plus a **Zarr**
interchange format, with bidirectional converters for **AnnData, Seurat,
SingleCellExperiment, Conos, and pagoda2** in **Python, R, and C++**. Use it to move a
dataset between toolkits, round-trip through `.lstar.zarr`, or flatten a multi-sample
collection — without writing bespoke per-pair conversion code, and without fabricating
data a target format can't hold (what doesn't fit goes to `dropped`, not silent loss).

The mental model: every conversion is `convert(X -> Y) = write_Y(read_X(obj))`, with an
L* dataset (in memory, or an on-disk `.lstar.zarr` store) as the universal intermediate.
lstar is **not** an analysis package — it does not cluster, integrate, or call markers.
For that, hand the converted object to pagoda2 / conos / Seurat / scanpy.

## Bundled references — load on demand

This SKILL.md covers the common conversions. Load a reference for the full matrix, the
data model, a specific language API, or collection details:

- `references/conversions.md` — the `lstar convert` CLI in full, the read_X/write_X
  matrix, what each route preserves vs. drops, version recognition, the deterministic
  role->slot mapping, and `--check` native-acceptance.
- `references/model.md` — the L* model: axes, fields (role/span/encoding/state),
  collections (heterogeneous samples are a collection, NOT one aligned tensor), store
  layout. Read this when a conversion's structure is confusing.
- `references/r.md` — full R API: `lstar_read`/`lstar_write`, the `read_*`/`write_*`
  profiles, `write_conos`/`read_conos`, backed readers.
- `references/python.md` — full Python API: `read_anndata`/`write_anndata`, `read`/`write`,
  lazy reads, `stream_col_stats`, the `lstar convert` CLI internals.
- `references/cpp.md` — the libstar C++ core (structs, Zarr IO, the OpenMP accelerator).
- `references/recipes.md` — task-oriented recipes (conversions, collections, round-trips,
  per-gene stats at scale).

## Install

lstar lives in a (private) repo with an **R package under `R/`** and a **Python package
under `python/`**; it is not on CRAN/PyPI. Install whichever language you need from the
repo. The `lstar convert` CLI needs the Python package; the `read_*`/`write_*` profiles
need the R package. The domain packages (anndata / SeuratObject / SingleCellExperiment)
are only needed for the formats you actually touch — lstar has a package-free `--backend
direct` fallback for `.h5ad` (h5py) and Seurat/SCE `.rds` (base R).

```bash
# Python (for the `lstar convert` CLI and AnnData routes):
pip install -e /path/to/lstar               # the repo root (editable); released as `pip install lstar-sc`
python3 -c "import lstar; print(lstar.__version__)"
```

```r
# R (for the read_*/write_* profiles and Conos/SCE/Seurat routes):
if (!requireNamespace("lstar", quietly = TRUE)) {
  # local source checkout (the R package is the repo's R/ subdir):
  install.packages("/path/to/lstar/R", repos = NULL, type = "source")
}
library(lstar)
```

## Decisions to surface up front

1. **Direction and formats.** Name the source and target explicitly (e.g. AnnData ->
   Seurat). Detection is by file extension; override with `--from`/`--to` if a path is
   ambiguous.
2. **Single object vs. collection.** A heterogeneous multi-sample dataset (e.g. a Conos
   panel: per-sample counts + a joint graph, divergent gene sets) is a **collection**, NOT
   one `cells x genes` matrix. lstar preserves that structure (per-sample axes + a joint
   layer); flattening it to a single AnnData unions the genes and is lossy in a documented
   way. See `references/model.md`.
3. **What gets dropped.** Conversions are lossless on the shared core (counts, data/X,
   pca + loadings, umap, labels, metadata). Anything a target can't hold lands in
   `dropped` (or `uns['lstar/dropped']`), never silently. Use `--report` to see it.
4. **Backend.** `--backend auto` (default) uses the format's native package when present,
   else lstar's package-free codec (`direct`). Force `--backend direct` to convert without
   installing anndata/SeuratObject/SingleCellExperiment.
5. **No analysis here.** lstar moves data; it does not analyze it. After converting, hand
   off to the appropriate analysis tool.

---

## Step 1 — Convert with the one-command CLI (the quickest path)

`lstar convert SRC DST` detects both formats from their paths, reads SRC into the L*
model, and writes DST. Extensions: `.h5ad` -> anndata, `.rds` -> Seurat/SCE (sniffed),
`.lstar.zarr` / `.zarr` -> store.

```bash
lstar convert sample.h5ad sample.rds                 # AnnData -> Seurat (.rds)
lstar convert sample.h5ad sample.rds --to sce        # ... -> SingleCellExperiment instead
lstar convert sample.h5ad sample.lstar.zarr --report # -> store + a fidelity report (fields kept + dropped)
lstar inspect sample.h5ad                            # read + print its L* structure, no write
```

`--check` (default on) opens the result in its native library and runs a canonical-ops
smoke test (verifies the destination's own tools accept it, not just that bytes
round-tripped); `--strict` makes a failed check non-zero. `--backend direct` converts
without the domain package (`.h5ad` via h5py; Seurat/SCE `.rds` via base R).

**Report:** print the `--report` summary (source/target format + version, fields
preserved, anything in `dropped`). Confirm the output file exists and `lstar inspect DST`
shows the expected axes/fields.

For the full CLI, the conversion matrix, and the role->slot mapping, read
`references/conversions.md`.

---

## Step 2 — Convert in R with the read_X / write_X profiles

In an R session the same conversions are `write_<target>(read_<source>(obj))`, with the L*
dataset as the intermediate. `read_*` = object -> L* dataset; `write_*` = L* dataset ->
object.

```r
library(lstar)
sce <- write_sce(read_seurat(seurat_obj))    # Seurat -> SingleCellExperiment (in memory)
so  <- write_seurat(read_sce(sce))            # ... and back
ds  <- read_seurat(seurat_obj)                # -> an L* dataset
lstar_write(ds, "sample.lstar.zarr")          # persist to a zarr store
ds2 <- lstar_read("sample.lstar.zarr")        # ... and read it back
```

`read_seurat_backed("file.h5ad")` reads an AnnData file into a (disk-backed) Seurat
object. Seurat legacy v2 -> v5 and the pagoda2 accessor-vs-slot differences are handled.

**Report:** `class()` + `dim()` of the converted object; `ds$dropped` for anything not
representable. For the full R API, read `references/r.md`.

---

## Step 3 — Convert in Python (AnnData routes)

```python
from lstar import read_anndata, write_anndata, write, read
write(read_anndata(adata), "a.lstar.zarr")     # AnnData -> store
adata2 = write_anndata(read("a.lstar.zarr"))    # store -> AnnData (X/obs/obsm/obsp restored; uns -> dropped)
```

For lazy/streamed per-gene stats over a big store without densifying:

```python
ds = lstar.read("big.lstar.zarr", lazy=True)
mean, var, nnz = lstar.stream_col_stats(ds.field("counts").values, lognorm=True, n_threads=0)
```

**Report:** `adata2.shape`, and `adata2.uns.get("lstar/dropped")`. For the full Python API
+ the C++ accelerator, read `references/python.md` and `references/cpp.md`.

---

## Step 4 — Collections (multi-sample, e.g. a Conos panel)

A Conos integration is a **collection**: per-sample raw counts + a joint graph / embedding
/ clustering over the union of cells (conos integrates in graph space — there is **no**
corrected expression matrix). lstar round-trips it and flattens it to native collection
formats without fabricating a corrected matrix.

```r
ds <- write_conos(con)                          # Conos -> L* collection
con2 <- read_conos(ds)                          # ... and back to a live Conos (graph restored)
so  <- write_seurat(ds)                         # -> Seurat v5 split assay (per-sample layers + joint Graph + DimReduc)
lstar_write(ds, "panel.lstar.zarr")             # persist the whole collection
```

```bash
lstar convert panel.lstar.zarr panel.h5ad       # -> single AnnData: X = raw joint counts (no corrected
                                                #    matrix), joint graph in obsp/connectivities, embedding
                                                #    in obsm, sample/clusters in obs (a documented flattening)
```

**Report:** for a collection, `ds$kind == "collection"`, the per-sample `cells.<s>` axes,
and the joint `graph`/`embedding`/cluster fields. For the collection model and what
flattening unions, read `references/model.md` and `references/conversions.md`.

---

## Batch variant

For `args == "batch"` (an orchestrator converting many files): loop the `lstar convert`
call (or the R/Python profile) over the inputs, suppress per-file `--report` prints, and
emit ONE summary line (`converted N files: <src-fmt> -> <dst-fmt>, M with dropped fields`).
Still write each output file named per the caller's pattern. Reserve `--check` for a final
spot-check rather than every file unless the caller asked for strict validation.

---

## Final response checklist

When the conversion is complete, summarize:

- Source format/version and target format; single object vs. collection.
- The route used (CLI / R profile / Python) and backend (native vs direct).
- What was preserved and what landed in `dropped` (from `--report` / `ds$dropped`).
- Output file(s) produced, and a native-acceptance confirmation (`--check` or
  `lstar inspect`).
- Caveats: collection flattening is lossy by construction; no analysis was performed.
