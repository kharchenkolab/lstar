"""Native-acceptance checks — after a conversion, open the produced object in its OWN ecosystem's
library and run a short canonical-ops smoke, so we know the native analysis tools won't choke on it.

This is the "third reader" — the destination toolchain (scanpy / Seurat / scran-scater) — that round-trip
fidelity never exercises: ``X → store → X`` proves we preserved *our* representation, not that the object
we handed back is *canonical* and that native tools accept it. Each check degrades gracefully: the heavy
analysis libraries are optional; absent → the check falls back to *open + structural invariants* and
reports ``ops skipped``, never a hard failure. The library needed merely to *open* the target being
absent yields ``skip`` (we couldn't check at all).

``check(dst, fmt) -> {"format", "status": "pass"|"fail"|"skip", "detail"}``.
"""
from __future__ import annotations

import os
import subprocess
import tempfile


def _res(fmt, status, detail):
    return {"format": fmt, "status": status, "detail": detail}


def check(dst: str, fmt: str) -> dict:
    if fmt == "anndata":
        return _check_anndata(dst)
    if fmt == "mudata":
        return _check_mudata(dst)
    if fmt in ("seurat", "sce", "rds"):
        return _check_r(dst, "sce" if fmt == "sce" else "seurat")
    if fmt == "store":
        return _res("store", "skip", "store is L*-native — nothing to native-check")
    return _res(fmt, "skip", "no native check defined for this format")


def _check_anndata(dst: str) -> dict:
    try:
        import anndata as ad
    except ImportError:
        return _res("anndata", "skip", "anndata not installed — cannot open the target")
    a = ad.read_h5ad(dst)
    if not (a.n_obs and a.n_vars):
        return _res("anndata", "fail", f"empty AnnData {a.shape}")
    if not a.obs_names.is_unique:
        return _res("anndata", "fail", "obs_names are not unique (scanpy/anndata require unique indices)")
    detail = f"opened {a.shape}, obs_names unique"
    try:
        import scanpy as sc
    except ImportError:
        return _res("anndata", "pass", detail + "; scanpy absent → ops skipped")
    try:
        c = a.copy()
        if c.X is None:                                   # scanpy ops read .X; populate from a layer
            if c.layers:
                c.X = next(iter(c.layers.values()))
            else:
                return _res("anndata", "pass", detail + "; no X/layers → ops skipped")
        sc.pp.normalize_total(c)
        sc.pp.log1p(c)
        sc.pp.pca(c, n_comps=max(2, min(5, min(c.shape) - 1)))
        ran = "normalize_total,log1p,pca"
        cats = [col for col in c.obs.columns
                if str(c.obs[col].dtype) == "category" and c.obs[col].nunique() > 1]
        if cats:
            sc.tl.rank_genes_groups(c, cats[0], method="t-test")
            ran += f",rank_genes_groups[{cats[0]}]"
        return _res("anndata", "pass", detail + f"; scanpy {ran} OK")
    except Exception as e:
        return _res("anndata", "fail", f"scanpy op failed: {type(e).__name__}: {e}")


def _check_mudata(dst: str) -> dict:
    try:
        import mudata as md
    except ImportError:
        return _res("mudata", "skip", "mudata not installed — cannot open the target")
    m = md.read_h5mu(dst)
    if not (m.n_obs and len(m.mod)):
        return _res("mudata", "fail", "empty MuData (no cells or no modalities)")
    return _res("mudata", "pass", f"opened {m.shape}, {len(m.mod)} modalities ({','.join(m.mod)})")


# R driver: open the .rds in its native class, assert structural invariants the native tools rely on,
# then (only if the full analysis package is present) run a canonical-ops smoke. Always exits 0 and
# reports via a single `NATIVE_CHECK<TAB>status<TAB>detail` line, so a Python caller reads one line.
_NATIVE_CHECK_R = r'''
args <- commandArgs(trailingOnly = TRUE); fmt <- args[1]; path <- args[2]
rlib <- Sys.getenv("LSTAR_RLIB", ""); if (nzchar(rlib)) .libPaths(c(rlib, .libPaths()))
emit <- function(status, detail) {
  cat(sprintf("NATIVE_CHECK\t%s\t%s\n", status, detail)); quit(save = "no", status = 0) }
obj <- tryCatch(readRDS(path), error = function(e) emit("fail", paste("readRDS:", conditionMessage(e))))

if (identical(fmt, "sce")) {
  if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) emit("skip", "SingleCellExperiment not installed")
  suppressMessages(library(SingleCellExperiment))
  if (!methods::is(obj, "SingleCellExperiment")) emit("fail", paste("not a SingleCellExperiment (got", class(obj)[1], ")"))
  d <- dim(obj); if (any(d == 0)) emit("fail", "empty SCE")
  an <- SummarizedExperiment::assayNames(obj)
  for (r in reducedDimNames(obj)) {                                   # reducedDim rowname alignment
    m <- reducedDim(obj, r)
    if (!is.null(rownames(m)) && !identical(rownames(m), colnames(obj)))
      emit("fail", sprintf("reducedDim %s rownames != cells", r)) }
  detail <- sprintf("opened %dx%d; assays={%s}; reducedDims={%s}; invariants OK",
                    d[1], d[2], paste(an, collapse = ","), paste(reducedDimNames(obj), collapse = ","))
  if (requireNamespace("scran", quietly = TRUE) && requireNamespace("scater", quietly = TRUE)) {
    res <- tryCatch({ suppressMessages({ library(scran); library(scater) })
      o <- obj; if (!("logcounts" %in% an)) o <- scater::logNormCounts(o)
      scran::modelGeneVar(o); o <- scater::runPCA(o, ncomponents = min(5L, ncol(o) - 1L)); "scran/scater logNorm/modelGeneVar/PCA OK" },
      error = function(e) paste("OP_FAIL:", conditionMessage(e)))
    if (grepl("^OP_FAIL", res)) emit("fail", res) else detail <- paste(detail, ";", res)
  } else detail <- paste(detail, "; scran/scater absent -> ops skipped")
  emit("pass", detail)
} else {
  if (!requireNamespace("SeuratObject", quietly = TRUE)) emit("skip", "SeuratObject not installed")
  suppressMessages(library(SeuratObject))
  if (!methods::is(obj, "Seurat")) emit("fail", paste("not a Seurat object (got", class(obj)[1], ")"))
  a <- DefaultAssay(obj)
  m <- tryCatch(GetAssayData(obj, assay = a, layer = "counts"),
                error = function(e) tryCatch(GetAssayData(obj, assay = a), error = function(e2) NULL))
  if (is.null(m) || any(dim(m) == 0)) emit("fail", "empty / unreadable default assay")
  if (!identical(rownames(obj@meta.data), colnames(obj))) emit("fail", "meta.data rownames != cell names")
  for (r in Reductions(obj)) {                                        # DimReduc keys must end in "_"
    k <- Key(obj[[r]]); if (!grepl("_$", k)) emit("fail", sprintf("reduction %s key '%s' not _-terminated", r, k)) }
  detail <- sprintf("opened %dx%d (assay %s); reductions={%s}; invariants OK",
                    nrow(m), ncol(m), a, paste(Reductions(obj), collapse = ","))
  if (requireNamespace("Seurat", quietly = TRUE)) {
    res <- tryCatch({ suppressMessages(library(Seurat))
      o <- NormalizeData(obj, verbose = FALSE)
      o <- FindVariableFeatures(o, verbose = FALSE, nfeatures = min(50L, nrow(o)))
      o <- ScaleData(o, verbose = FALSE)
      o <- RunPCA(o, verbose = FALSE, npcs = min(5L, ncol(o) - 1L)); "Seurat normalize/HVG/scale/PCA OK" },
      error = function(e) paste("OP_FAIL:", conditionMessage(e)))
    if (grepl("^OP_FAIL", res)) emit("fail", res) else detail <- paste(detail, ";", res)
  } else detail <- paste(detail, "; Seurat absent -> ops skipped")
  emit("pass", detail)
}
'''


def _check_r(dst: str, fmt: str) -> dict:
    rscript = os.environ.get("LSTAR_RSCRIPT", "Rscript")
    fh = tempfile.NamedTemporaryFile("w", suffix=".R", delete=False)
    fh.write(_NATIVE_CHECK_R)
    fh.close()
    try:
        proc = subprocess.run([rscript, fh.name, fmt, dst],
                              stdin=subprocess.DEVNULL, capture_output=True, text=True)
    except FileNotFoundError:
        return _res(fmt, "skip", "Rscript not found — cannot run the native check")
    finally:
        os.unlink(fh.name)
    line = next((ln for ln in proc.stdout.splitlines() if ln.startswith("NATIVE_CHECK")), None)
    if not line:
        tail = "\n".join((proc.stderr or proc.stdout).strip().splitlines()[-4:])
        return _res(fmt, "fail", f"native check did not run: {tail}")
    parts = (line.split("\t") + ["", ""])[:3]
    return _res(fmt, parts[1], parts[2])
