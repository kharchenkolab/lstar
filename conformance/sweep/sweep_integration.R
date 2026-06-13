# Sweep real multi-sample INTEGRATION datasets read AS COLLECTIONS (not flattened): ifnb (ctrl/stim),
# panc8 (8 datasets / 5 technologies), pbmcsca (multi-method). Each is one Seurat object spanning several
# samples; instead of reading it as one aligned matrix we SPLIT it by its sample column into per-sample
# `cells.<s>` / `genes.<s>` axes + `counts.<s>` measures + a `samples` axis + a union `cells` axis with a
# `sample` label -- the L* collection shape ([[feedback-collection-not-tensor]]) -- then write + validate.
# This is the integration scenario as a heterogeneous collection, exercising the same model the Conos
# profile builds, on real published integration data.
#
# Run: Rscript conformance/sweep/sweep_integration.R   (writes /tmp/sweep_integration.tsv)
.libPaths(c("/home/pkharchenko/p21/lstar/.Rlib", .libPaths()))
suppressMessages({library(Seurat); library(SeuratObject); library(SeuratData); library(lstar); library(Matrix)})

# dataset -> object name -> the metadata column that identifies the sample/batch
jobs <- list(
  list(d = "ifnb",    obj = "ifnb",    by = "stim"),       # ctrl vs stim
  list(d = "panc8",   obj = "panc8",   by = "dataset"),    # 8 datasets across 5 technologies
  list(d = "pbmcsca", obj = "pbmcsca", by = "Method")      # multi-method comparison
)

# Build an L* collection dataset by splitting a Seurat object by a sample column.
build_collection <- function(so, by) {
  meta <- so[[]]
  if (!(by %in% colnames(meta))) stop("no sample col: ", by)
  samp <- as.character(meta[[by]])
  samp <- gsub("[^A-Za-z0-9]+", "_", samp)                 # sanitize sample names for axis names
  usamp <- unique(samp)
  counts <- GetAssayData(so, assay = DefaultAssay(so), layer = "counts")
  if (is.null(counts) || nrow(counts) == 0)
    counts <- GetAssayData(so, assay = DefaultAssay(so), layer = "data")
  genes <- rownames(counts)
  ds <- list(kind = "collection", spec_version = "0.1",
             profiles = paste0("integration-split@", as.character(packageVersion("Seurat"))),
             dropped = character(0), axes = list(), fields = list())
  ucells <- character(0); ulabel <- character(0)
  for (sn in usamp) {
    idx <- which(samp == sn)
    cn <- colnames(counts)[idx]
    sub <- counts[, idx, drop = FALSE]
    # drop genes with zero total in this sample so per-sample gene sets genuinely differ (heterogeneous)
    keep <- Matrix::rowSums(sub) > 0
    sub <- sub[keep, , drop = FALSE]
    ds$axes[[paste0("cells.", sn)]] <- list(labels = cn, origin = "observed", role = "observation")
    ds$axes[[paste0("genes.", sn)]] <- list(labels = rownames(sub), origin = "observed", role = "feature")
    # L* measure is (cells, genes): transpose the Seurat (genes x cells) matrix
    ds$fields[[paste0("counts.", sn)]] <- list(role = "measure",
        span = c(paste0("cells.", sn), paste0("genes.", sn)), state = "raw", values = Matrix::t(sub))
    ucells <- c(ucells, cn); ulabel <- c(ulabel, rep(sn, length(cn)))
  }
  ds$axes[["samples"]] <- list(labels = usamp, origin = "observed", role = "sample")
  ds$axes[["cells"]]   <- list(labels = ucells, origin = "derived", role = "observation")
  ds$fields[["n_cells"]] <- list(role = "measure", span = "samples",
      values = as.numeric(table(factor(samp, levels = usamp))))
  ds$fields[["sample"]]  <- list(role = "label", span = "cells", values = ulabel)
  class(ds) <- "lstar_dataset"
  list(ds = ds, nsamp = length(usamp))
}

rep <- file("/tmp/sweep_integration.tsv", "w")
cat("dataset\tby\tstatus\tsamples\taxes\tfields\tncells\tnote\n", file = rep)
ok <- 0; fail <- 0
for (j in jobs) {
  cat(sprintf("[integration] %-10s by=%-8s ... ", j$d, j$by)); flush.console()
  r <- tryCatch({
    pkg <- paste0(j$d, ".SeuratData")
    suppressWarnings(suppressMessages(data(list = j$obj, package = pkg)))
    so <- UpdateSeuratObject(get(j$obj))
    bc <- build_collection(so, j$by)
    ds <- bc$ds
    errs <- character(0)
    p <- tempfile(fileext = ".zarr"); lstar_write(ds, p)        # write IS the validation in R (cores check spans)
    # round-trip: read back and confirm the collection structure survived
    ds2 <- lstar_read(p)
    unlink(p, recursive = TRUE)
    n_csamp <- sum(grepl("^counts\\.", names(ds2$fields)))
    if (n_csamp != bc$nsamp) errs <- c(errs, sprintf("ERROR: read-back has %d per-sample counts, expected %d", n_csamp, bc$nsamp))
    list(s = if (length(errs)) "VALIDATE-ERR" else "PASS", nsamp = bc$nsamp,
         a = length(ds$axes), fld = length(ds$fields), nc = ncol(so),
         n = if (length(errs)) substr(errs[1], 1, 90) else "collection round-trips")
  }, error = function(e) list(s = "FAIL", nsamp = "", a = "", fld = "", nc = "",
                              n = substr(conditionMessage(e), 1, 100)))
  if (r$s == "PASS") ok <- ok + 1 else fail <- fail + 1
  cat(r$s, "\n"); flush.console()
  cat(sprintf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", j$d, j$by, r$s, r$nsamp, r$a, r$fld, r$nc, r$n),
      file = rep); flush(rep)
}
close(rep)
cat(sprintf("integration-collection sweep: %d PASS / %d FAIL\n", ok, fail))
