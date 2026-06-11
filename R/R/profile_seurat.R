# Seurat profile: build a Seurat object from an L* dataset and back, across Seurat versions.
# L* measures over (cells, genes) are transposed to Seurat's (genes x cells) orientation.
#
# Version recognition. Seurat's object layout changed between major versions and the profile
# adapts rather than assuming one shape:
#   - v5 (Assay5): per-layer access via Layers()/LayerData(); a split assay is a collection.
#   - v3/v4 (Assay): three fixed slots (counts/data/scale.data); no layers, no collection.
#   - SeuratObject < 5 lacks the Layers()/LayerData() API entirely -> fall back to GetAssayData().
# The detected versions are recorded in ds$profiles ("seurat@0.1", "SeuratObject@<v>", "assay@<v>").

.seurat_versions <- function(so, assay) {
  sov <- tryCatch(as.character(utils::packageVersion("SeuratObject")), error = function(e) "?")
  acl <- tryCatch(class(so[[assay]])[1], error = function(e) "Assay")
  av <- if (identical(acl, "Assay5")) "v5" else if (identical(acl, "Assay")) "v3" else acl
  c("seurat@0.1", paste0("SeuratObject@", sov), paste0("assay@", av))
}

# A version-agnostic view of an assay's layers: $names, $get(layer)->matrix, $cells(layer).
.seurat_layer_access <- function(so, assay) {
  have_layers <- "Layers" %in% getNamespaceExports("SeuratObject")
  if (have_layers) {
    list(names = SeuratObject::Layers(so, assay = assay),
         get = function(L) SeuratObject::LayerData(so, assay = assay, layer = L),
         cells = function(L, m) {
           lc <- colnames(m)
           if (is.null(lc)) lc <- tryCatch(SeuratObject::Cells(so[[assay]], layer = L),
                                           error = function(e) colnames(so))
           lc
         })
  } else {                                  # SeuratObject < 5: fixed slots, no layers/collections
    nz <- function(s) { m <- tryCatch(SeuratObject::GetAssayData(so, slot = s, assay = assay),
                                      error = function(e) NULL); !is.null(m) && length(m) > 0 && nrow(m) > 0 }
    list(names = Filter(nz, c("counts", "data", "scale.data")),
         get = function(L) SeuratObject::GetAssayData(so, slot = L, assay = assay),
         cells = function(L, m) colnames(so))
  }
}

.dimreduc_key <- function(nm) {
  switch(tolower(nm), pca = "PC_", umap = "UMAP_", tsne = "tSNE_", paste0(toupper(nm), "_"))
}

.fields_over <- function(ds, span, role = NULL) {
  keep <- character(0)
  for (nm in names(ds$fields)) {
    f <- ds$fields[[nm]]
    if (identical(as.character(f$span), as.character(span)) &&
        (is.null(role) || identical(f$role, role))) keep <- c(keep, nm)
  }
  keep
}

.pick <- function(ds, cands, prefer_state = NULL, prefer_name = NULL, exclude = NULL) {
  cands <- setdiff(cands, exclude)
  if (length(cands) == 0) return(NULL)
  if (!is.null(prefer_state)) {
    hit <- cands[vapply(cands, function(n) identical(ds$fields[[n]]$state, prefer_state), logical(1))]
    if (length(hit)) return(hit[[1]])
  }
  if (!is.null(prefer_name) && prefer_name %in% cands) return(prefer_name)
  cands[[1]]
}

#' Build a Seurat object from an L* dataset.
#'
#' Measures over `(cells, genes)` become assay layers (transposed to Seurat's genes x cells
#' orientation); embeddings become `DimReduc`s; arity-1 cell fields become `meta.data`.
#'
#' @param ds an `lstar_dataset`
#' @return a `Seurat` object.
#' @seealso [read_seurat()]
#' @export
write_seurat <- function(ds) {
  if (!requireNamespace("SeuratObject", quietly = TRUE)) stop("SeuratObject is required")
  cells <- as.character(ds$axes$cells$labels)
  genes <- as.character(ds$axes$genes$labels)

  meas <- .fields_over(ds, c("cells", "genes"), role = "measure")
  counts_nm <- .pick(ds, meas, prefer_state = "raw", prefer_name = "counts")
  data_nm <- .pick(ds, meas, prefer_state = "lognorm", prefer_name = "X", exclude = counts_nm)
  if (is.null(counts_nm) && is.null(data_nm)) stop("no (cells x genes) measure to build an assay")

  gxc <- function(nm) {
    m <- Matrix::t(ds$fields[[nm]]$values)        # (cells x genes) -> (genes x cells)
    dimnames(m) <- list(genes, cells)
    m
  }

  primary <- if (!is.null(counts_nm)) counts_nm else data_nm
  so <- SeuratObject::CreateSeuratObject(counts = gxc(primary), assay = "RNA")
  if (!is.null(counts_nm) && !is.null(data_nm)) {
    SeuratObject::LayerData(so, assay = "RNA", layer = "data") <- gxc(data_nm)
  }

  # meta.data: arity-1 fields over cells
  for (nm in names(ds$fields)) {
    f <- ds$fields[[nm]]
    if (identical(as.character(f$span), "cells") && length(f$values) == length(cells)) {
      so[[nm]] <- stats::setNames(as.vector(f$values), cells)
    }
  }

  # reductions: embedding fields over (cells, <coord>), with loadings if present
  for (nm in names(ds$fields)) {
    f <- ds$fields[[nm]]
    if (identical(f$role, "embedding") && length(f$span) == 2 && f$span[[1]] == "cells") {
      coord <- f$span[[2]]
      emb <- as.matrix(f$values)
      rownames(emb) <- cells
      key <- .dimreduc_key(nm)
      colnames(emb) <- paste0(key, seq_len(ncol(emb)))
      load_nm <- paste0(coord, "_loadings")
      if (!is.null(ds$fields[[load_nm]])) {
        L <- as.matrix(ds$fields[[load_nm]]$values)
        rownames(L) <- genes
        colnames(L) <- colnames(emb)
        dr <- SeuratObject::CreateDimReducObject(embeddings = emb, loadings = L,
                                                 key = key, assay = "RNA")
      } else {
        dr <- SeuratObject::CreateDimReducObject(embeddings = emb, key = key, assay = "RNA")
      }
      so[[nm]] <- dr
    }
  }
  so
}

#' Read a Seurat object into an L* dataset.
#'
#' Handles Seurat v3/v4 (`Assay`) and v5 (`Assay5`); a v5 assay split by sample
#' (`split(assay, f = ...)`) is read as an L* collection. The detected versions are recorded in
#' `ds$profiles`.
#'
#' @param so a `Seurat` object
#' @param assay assay to read (default: the default assay)
#' @return an `lstar_dataset` (of kind `"sample"`, or `"collection"` for a split assay).
#' @seealso [write_seurat()]
#' @export
read_seurat <- function(so, assay = SeuratObject::DefaultAssay(so)) {
  cells <- colnames(so)
  genes <- rownames(so[[assay]])
  ds <- list(kind = "sample", spec_version = "0.1",
             profiles = .seurat_versions(so, assay),
             dropped = character(0), axes = list(), fields = list())
  ds$axes$cells <- list(labels = cells, origin = "observed", role = "observation")
  ds$axes$genes <- list(labels = genes, origin = "observed", role = "feature")
  lac <- .seurat_layer_access(so, assay)

  add <- function(nm, values, role, span, state = "") {
    ds$fields[[nm]] <<- list(role = role, span = span, state = state, subtype = "", values = values)
  }

  add2 <- function(nm, values, role, span, state = "", subtype = "") {
    ds$fields[[nm]] <<- list(role = role, span = span, state = state, subtype = subtype,
                             values = values)
  }
  state_of <- function(L) switch(L, counts = "raw", data = "lognorm", scale.data = "scaled", "")
  name_of <- function(L) switch(L, counts = "counts", data = "X", L)
  # A Seurat v5 assay split for integration (split(assay, f = sample)) holds a *collection*:
  # layers named "<root>.<sample>" each cover only their sample's cells. Joined layers cover all
  # cells (one aligned matrix). Distinguish by whether a layer's cells are a strict subset.
  samples_seen <- character(0)
  sample_of_cell <- stats::setNames(rep(NA_character_, length(cells)), cells)
  for (L in lac$names) {
    m <- lac$get(L)                                              # genes x (layer cells)
    lc <- lac$cells(L, m)
    parts <- regmatches(L, regexec("^(counts|data|scale\\.data)\\.(.+)$", L))[[1]]
    if (length(parts) == 3 && length(lc) < length(cells)) {     # a per-sample collection layer
      root <- parts[2]; sn <- parts[3]; ca <- paste0("cells.", sn)
      if (!ca %in% names(ds$axes))
        ds$axes[[ca]] <- list(labels = lc, origin = "observed", role = "observation")
      add2(paste0(name_of(root), ".", sn), Matrix::t(m), "measure", c(ca, "genes"),
           state = state_of(root))
      samples_seen <- union(samples_seen, sn); sample_of_cell[lc] <- sn
    } else {                                                    # a joined (aligned) layer
      add(name_of(L), Matrix::t(m), "measure", c("cells", "genes"), state = state_of(L))
    }
  }

  md <- so[[]]
  for (col in colnames(md)) {
    v <- md[[col]]
    if (is.numeric(v)) add(col, as.numeric(v), "measure", "cells")
    else add(col, as.character(v), "label", "cells")
  }

  # If the assay was split, record the collection structure: a samples axis + a per-cell
  # sample label (the design grouping the layers encode).
  if (length(samples_seen)) {
    ds$kind <- "collection"
    ds$axes$samples <- list(labels = samples_seen, origin = "observed", role = "sample")
    add2("sample", unname(sample_of_cell), "label", "cells", subtype = "design")
  }

  for (rn in SeuratObject::Reductions(so)) {
    dr <- so[[rn]]
    emb <- SeuratObject::Embeddings(dr)
    ds$axes[[rn]] <- list(labels = colnames(emb), origin = "derived", role = "coordinate")
    add(rn, unname(as.matrix(emb)), "embedding", c("cells", rn))
    ld <- SeuratObject::Loadings(dr)
    if (length(ld) > 0 && nrow(ld) > 0)
      add(paste0(rn, "_loadings"), unname(as.matrix(ld)), "loading", c("genes", rn))
  }

  class(ds) <- "lstar_dataset"
  ds
}
