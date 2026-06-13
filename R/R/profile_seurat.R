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
  pov <- tryCatch(as.character(utils::packageVersion("SeuratObject")), error = function(e) "?")
  # the OBJECT's serialized version (a v3/v4 object loaded under Seurat 5 still reports its own version),
  # plus EVERY assay's class -- a multimodal object mixes classes (RNA Assay5 + ADT Assay5, or a v3
  # Assay + an SCTAssay). This is the per-object/per-assay version tracking the corpus relies on.
  ov  <- tryCatch(as.character(SeuratObject::Version(so)), error = function(e) pov)
  assays <- tryCatch(SeuratObject::Assays(so), error = function(e) assay)
  cls <- vapply(assays, function(a) tryCatch(class(so[[a]])[1], error = function(e) "Assay"), character(1))
  c("seurat@0.1", paste0("SeuratObject@", pov), paste0("object@", ov),
    paste0("assay@", assays, ":", cls))
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

  # Collection (a Seurat v5 *split*/integration object reads as per-sample measures over
  # (cells.<s>, genes)). Join them into the union (cells, genes) measure(s) so a single assay can be
  # built, recording each cell's sample so the assay can be re-split into v5 layers below.
  split_by <- NULL
  persamp <- Filter(function(nm) {
    f <- ds$fields[[nm]]; sp <- as.character(f$span)
    identical(f$role, "measure") && length(sp) == 2 && sp[2] == "genes" && startsWith(sp[1], "cells.")
  }, names(ds$fields))
  if (length(persamp)) {
    soc <- stats::setNames(rep(NA_character_, length(cells)), cells)
    roots <- sub("\\.[^.]+$", "", persamp)
    for (root in unique(roots)) {
      M <- Matrix::Matrix(0, nrow = length(cells), ncol = length(genes), sparse = TRUE,
                          dimnames = list(cells, genes))
      st <- "raw"
      for (nm in persamp[roots == root]) {
        ca <- as.character(ds$fields[[nm]]$span)[1]
        sc <- as.character(ds$axes[[ca]]$labels)
        M[sc, ] <- as(ds$fields[[nm]]$values, "CsparseMatrix")
        soc[sc] <- sub("^cells\\.", "", ca)
        st <- ds$fields[[nm]]$state %||% st
        ds$fields[[nm]] <- NULL
      }
      ds$fields[[root]] <- list(values = M, role = "measure", span = c("cells", "genes"), state = st)
    }
    split_by <- soc
  }

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
      val <- f$values                                     # keep a factor a factor (don't strip to character)
      so[[nm]] <- stats::setNames(if (is.factor(val)) val else as.vector(val), cells)
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
      sd_nm <- paste0(coord, "_stdev")
      sdv <- if (!is.null(ds$fields[[sd_nm]])) as.numeric(ds$fields[[sd_nm]]$values) else numeric(0)
      if (!is.null(ds$fields[[load_nm]])) {
        L <- as.matrix(ds$fields[[load_nm]]$values)
        rownames(L) <- genes
        colnames(L) <- colnames(emb)
        dr <- SeuratObject::CreateDimReducObject(embeddings = emb, loadings = L,
                                                 key = key, assay = "RNA", stdev = sdv)
      } else {
        dr <- SeuratObject::CreateDimReducObject(embeddings = emb, key = key, assay = "RNA", stdev = sdv)
      }
      so[[nm]] <- dr
    }
  }

  if (!is.null(ds$fields[["ident"]])) {            # restore the active identity (Idents)
    iv <- ds$fields[["ident"]]$values
    if (!is.factor(iv)) iv <- as.factor(iv)
    SeuratObject::Idents(so) <- stats::setNames(iv, cells)
  }

  # Multimodal: rebuild every non-`genes` feature axis as its own assay (ADT, ATAC, ...) from its
  # measures `<axis>.<layer>` over (cells, <axis>).
  for (fax in names(ds$axes)) {
    ax <- ds$axes[[fax]]
    if (!identical(ax$role, "feature") || fax == "genes") next
    fm <- Filter(function(nm) { f <- ds$fields[[nm]]; sp <- as.character(f$span)
      identical(f$role, "measure") && length(sp) == 2 && sp[1] == "cells" && sp[2] == fax }, names(ds$fields))
    if (!length(fm)) next
    feats <- as.character(ax$labels)
    fxc <- function(nm) { mm <- Matrix::t(as(ds$fields[[nm]]$values, "CsparseMatrix")); dimnames(mm) <- list(feats, cells); mm }
    states <- vapply(fm, function(nm) ds$fields[[nm]]$state %||% "", character(1))
    primary <- fm[match("raw", states)]; if (is.na(primary)) primary <- fm[1]
    aobj <- SeuratObject::CreateAssay5Object(counts = fxc(primary))
    for (nm in setdiff(fm, primary)) {
      lyr <- sub(paste0("^", fax, "\\."), "", nm)               # "ADT.data" -> "data"
      SeuratObject::LayerData(aobj, layer = lyr) <- fxc(nm)
    }
    so[[fax]] <- aobj
  }

  if (!is.null(split_by) && inherits(so[["RNA"]], "Assay5")) {   # re-split a v5 integration object
    so[["RNA"]] <- split(so[["RNA"]], f = unname(split_by[cells]))
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
  others <- setdiff(tryCatch(SeuratObject::Assays(so), error = function(e) assay), assay)
  lac <- .seurat_layer_access(so, assay)

  add <- function(nm, values, role, span, state = "") {
    ds$fields[[nm]] <<- list(role = role, span = span, state = state, subtype = "", values = values)
  }

  add2 <- function(nm, values, role, span, state = "", subtype = "") {
    ds$fields[[nm]] <<- list(role = role, span = span, state = state, subtype = subtype,
                             values = values)
  }
  add_factor <- function(nm, v, span) {                   # a factor column -> categorical field that
    v <- as.factor(v)                                     # induces a bare-named `factor` axis (its levels)
    ds$fields[[nm]] <<- list(role = "label", span = span, state = "", subtype = "",
                             encoding = "categorical", values = v)
    if (is.null(ds$axes[[nm]]))                           # don't clobber an existing axis of that name
      ds$axes[[nm]] <<- list(labels = levels(v), origin = "derived", role = "factor", induced_by = nm)
  }
  # A Signac ChromatinAssay (scATAC) carries the defining peak **genomic ranges** (a GRanges over its
  # features) plus external **fragment** files. Type the ranges as arity-1 feature fields over the peak
  # axis (chromosome -> factor; start/end -> measures) so the coordinates survive the round-trip, and
  # record the fragment files (external, can't inline). Without this they're silently lost -- a real
  # corpus bug (a pure-ATAC object, where the ChromatinAssay is the *default* assay, never hit the
  # other-assay recording path). Slot-accessed + tryCatch so the profile never hard-depends on Signac.
  capture_chromatin <- function(a_name, axis, n) {
    a <- tryCatch(so[[a_name]], error = function(e) NULL)
    if (is.null(a) || !methods::is(a, "ChromatinAssay")) return(invisible())
    gr <- tryCatch(methods::slot(a, "ranges"), error = function(e) NULL)
    got <- FALSE
    if (!is.null(gr) && length(gr) == n && requireNamespace("GenomicRanges", quietly = TRUE)) {
      sn <- tryCatch(as.character(GenomicRanges::seqnames(gr)), error = function(e) NULL)
      st <- tryCatch(as.numeric(GenomicRanges::start(gr)), error = function(e) NULL)
      en <- tryCatch(as.numeric(GenomicRanges::end(gr)), error = function(e) NULL)
      if (length(sn) == n) { add_factor(paste0(axis, "_seqnames"), sn, axis); got <- TRUE }
      if (length(st) == n) add2(paste0(axis, "_start"), st, "measure", axis, subtype = "genomic_pos")
      if (length(en) == n) add2(paste0(axis, "_end"), en, "measure", axis, subtype = "genomic_pos")
    }
    if (!got) ds$dropped <<- c(ds$dropped, sprintf("assay/%s/ranges (Signac ChromatinAssay; uncaptured)", a_name))
    frag <- tryCatch(methods::slot(a, "fragments"), error = function(e) list())
    if (length(frag)) ds$dropped <<- c(ds$dropped,
      sprintf("assay/%s/fragments (%d external file(s), not inlined)", a_name, length(frag)))
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
    } else if (nrow(m) != length(genes)) {
      # a layer over a *feature subset* (real objects keep scale.data over variable features only, e.g.
      # pbmc3k.final: 2000 HVGs of 13714). Partial coverage isn't typed yet -> record the loss rather
      # than mis-span it over the full gene axis (which crashes write-back) or silently drop it.
      ds$dropped <- c(ds$dropped, sprintf("layer/%s (%d of %d features)", L, nrow(m), length(genes)))
    } else {                                                    # a joined (aligned) layer
      add(name_of(L), Matrix::t(m), "measure", c("cells", "genes"), state = state_of(L))
    }
  }
  capture_chromatin(assay, "genes", length(genes))              # default assay may itself be ATAC (peaks)

  # Multimodal: every *other* assay is its own feature space over the shared `cells` axis (CITE-seq
  # RNA+ADT, 10x multiome RNA+ATAC). Capture it as a feature axis named after the assay + measures
  # `<assay>.<layer>` over (cells, <assay>), rather than dropping it. (A Signac ChromatinAssay's genomic
  # ranges/fragments aren't typed yet -> recorded.)
  for (a in others) {
    feats <- rownames(so[[a]])
    if (a %in% names(ds$axes)) { ds$dropped <- c(ds$dropped, paste0("assay/", a, " (axis-name clash)")); next }
    ds$axes[[a]] <- list(labels = feats, origin = "observed", role = "feature")
    laca <- .seurat_layer_access(so, a)
    for (L in laca$names) {
      m <- laca$get(L)
      if (nrow(m) != length(feats)) {                           # subset layer -> partial coverage
        ds$dropped <- c(ds$dropped, sprintf("assay/%s layer/%s (%d of %d features)", a, L, nrow(m), length(feats)))
        next
      }
      add(paste0(a, ".", L), Matrix::t(m), "measure", c("cells", a), state = state_of(L))
    }
    capture_chromatin(a, a, length(feats))                      # ATAC as a non-default assay (peaks axis)
  }

  md <- so[[]]
  for (col in colnames(md)) {
    v <- md[[col]]
    if (is.numeric(v)) add(col, as.numeric(v), "measure", "cells")
    else if (is.factor(v)) add_factor(col, v, "cells")    # preserve levels/order + induce a factor axis
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
    sd <- tryCatch(SeuratObject::Stdev(dr), error = function(e) numeric(0))   # per-dim stdev -> measure
    if (length(sd) == ncol(emb)) add(paste0(rn, "_stdev"), as.numeric(sd), "measure", rn)
    ld <- SeuratObject::Loadings(dr)
    if (length(ld) > 0 && nrow(ld) > 0) {
      if (nrow(ld) == length(genes))                     # loadings over all genes -> a loading field
        add(paste0(rn, "_loadings"), unname(as.matrix(ld)), "loading", c("genes", rn))
      else                                               # real PCA loadings are over variable features
        ds$dropped <- c(ds$dropped,                      # only (e.g. pbmc3k.final: 2000/13714) -> record
                        sprintf("loadings/%s (%d of %d features)", rn, nrow(ld), length(genes)))
    }
  }

  # active identity (Idents): the active-vs-stored distinction is otherwise lost. Capture it as a
  # categorical 'ident' field (inducing its factor axis) flagged `active_ident`; restored on write.
  id <- tryCatch(SeuratObject::Idents(so), error = function(e) NULL)
  if (!is.null(id) && length(id) == length(cells)) {
    idf <- droplevels(as.factor(id))
    ds$fields[["ident"]] <- list(role = "label", span = "cells", state = "", subtype = "active_ident",
                                 encoding = "categorical", values = idf)
    if (is.null(ds$axes[["ident"]]))
      ds$axes[["ident"]] <- list(labels = levels(idf), origin = "derived", role = "factor",
                                 induced_by = "ident")
  }

  class(ds) <- "lstar_dataset"
  ds
}
