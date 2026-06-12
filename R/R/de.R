# Differential-expression / marker results as fields over a factor axis (see docs/model.md induction).
# A DE bundle is the set of `de.<factor>.<stat>` measures over (factor, genes); `lstar_markers()` is the
# tidy long-form view, the R counterpart of `lstar.markers()` in Python.

#' Tidy marker table for a factor's DE bundle.
#'
#' Gathers the `de.<factor>.<stat>` fields (score/lfc/pval/padj over the factor x genes axes) into a
#' long-form data frame: one row per (group, gene) with whichever statistics are present.
#'
#' @param ds an `lstar_dataset` (e.g. from [lstar_read()]).
#' @param factor the factor-axis name the DE was computed over (e.g. `"leiden"`).
#' @param top optional: keep only the top-N genes per group (by `sort_by`).
#' @param sort_by statistic to rank within a group (default `"score"`).
#' @param descending sort descending (default `TRUE`).
#' @return a data frame with columns `group`, `gene`, and the available statistics.
#' @export
lstar_markers <- function(ds, factor, top = NULL, sort_by = "score", descending = TRUE) {
  prefix <- paste0("de.", factor, ".")
  nms <- names(ds$fields)
  de_fields <- nms[startsWith(nms, prefix)]
  if (!length(de_fields)) stop(sprintf("no DE bundle for factor '%s' in dataset", factor))
  stats <- substring(de_fields, nchar(prefix) + 1L)
  gene_axis <- as.character(ds$fields[[de_fields[1]]]$span)[2]
  groups <- as.character(ds$axes[[factor]]$labels)
  genes <- as.character(ds$axes[[gene_axis]]$labels)
  mats <- stats::setNames(lapply(de_fields, function(nm) as.matrix(ds$fields[[nm]]$values)), stats)

  parts <- lapply(seq_along(groups), function(gi) {
    keep <- which(!is.na(mats[[1]][gi, ]))                      # genes ranked for this group
    if (sort_by %in% stats) keep <- keep[order(mats[[sort_by]][gi, keep], decreasing = descending)]
    if (!is.null(top)) keep <- utils::head(keep, top)
    df <- data.frame(group = groups[gi], gene = genes[keep], stringsAsFactors = FALSE)
    for (st in stats) df[[st]] <- mats[[st]][gi, keep]
    df
  })
  do.call(rbind, parts)
}
