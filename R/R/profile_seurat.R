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
# A profile's `cache` fields (provenance$cache, e.g. the viewer@0.1 navigators) are regenerable from
# primary data, so a converter to a foreign object drops them (recorded in `dropped`) rather than
# carrying a redundant/mis-aligned copy. Read-by-name is unaffected; only this format mapping skips them.
.lstar_drop_cache <- function(ds) {
  cache <- names(Filter(function(f) !is.null(f$provenance$cache), ds$fields))
  if (length(cache)) { ds$dropped <- unique(c(ds$dropped, cache)); ds$fields[cache] <- NULL }
  ds
}

# Canonical multimodal feature-axis vocabulary -- the SAME mapping mudata uses (profiles/mudata.py),
# so a modality lands on the same axis (`genes`/`proteins`/`peaks`) regardless of source format. A
# Seurat assay name maps by its lowercased name; an unknown assay keeps a sanitized name. The original
# assay name is preserved in the measure's `provenance$assay` so write_seurat restores it exactly.
.SEURAT_MODALITY_AXIS <- c(rna = "genes", gex = "genes", "gene expression" = "genes",
  adt = "proteins", prot = "proteins", protein = "proteins", proteins = "proteins",
  antibody = "proteins", "antibody capture" = "proteins", cite = "proteins",
  atac = "peaks", peak = "peaks", peaks = "peaks", "chromatin accessibility" = "peaks")
.modality_axis <- function(name, taken = character(0)) {
  ax <- unname(.SEURAT_MODALITY_AXIS[tolower(name)])
  if (is.na(ax)) ax <- make.names(tolower(name))
  base <- ax; i <- 1L; while (ax %in% taken) { ax <- paste0(base, i); i <- i + 1L }
  ax
}

write_seurat <- function(ds) {
  if (!requireNamespace("SeuratObject", quietly = TRUE)) stop("SeuratObject is required")
  ds <- .lstar_drop_cache(ds)
  cells <- as.character(ds$axes$cells$labels)
  genes <- as.character(ds$axes$genes$labels)

  # Collection -> a Seurat v5 *split* assay. Per-sample measures live over (cells.<s>, genes[.<s>]): a
  # Seurat-v5-split origin shares the `genes` axis across samples, whereas a Conos collection gives each
  # sample its OWN genes.<s> axis (gene sets that overlap, differ, or are entirely disjoint). Either way,
  # union the gene sets, place each sample's counts into the union columns BY NAME (so divergent genes are
  # tolerated, absent genes left 0), and record each cell's sample so the assay re-splits into v5 layers.
  # No corrected/integrated expression is fabricated -- the joint layer is the graph + embedding + clusters.
  split_by <- NULL
  persamp <- Filter(function(nm) {
    f <- ds$fields[[nm]]; sp <- as.character(f$span)
    identical(f$role, "measure") && length(sp) == 2 && startsWith(sp[1], "cells.") &&
      (sp[2] == "genes" || startsWith(sp[2], "genes."))
  }, names(ds$fields))
  if (length(persamp)) {
    if (!length(genes))                                  # no shared genes axis (Conos): union per-sample sets
      genes <- Reduce(union, lapply(persamp, function(nm)
        as.character(ds$axes[[as.character(ds$fields[[nm]]$span)[2]]]$labels)))
    soc <- stats::setNames(rep(NA_character_, length(cells)), cells)
    cell_pos <- stats::setNames(seq_along(cells), cells)            # name -> union row/col, vectorized
    gene_pos <- stats::setNames(seq_along(genes), genes)
    roots <- sub("\\.[^.]+$", "", persamp)
    for (root in unique(roots)) {
      # Assemble the union (cells x genes) matrix from triplets in ONE sparseMatrix() call. Each sample's
      # block is remapped into the union rows/cols BY NAME -- divergent genes tolerated, absent left 0.
      # (Avoid `M[sc, sg] <- block` subset-assignment into a zeroed sparse matrix: it is O(panel size) per
      # write and hangs on real panels with ~30k genes.)
      ii <- integer(0); jj <- integer(0); xx <- numeric(0); st <- "raw"
      for (nm in persamp[roots == root]) {
        sp <- as.character(ds$fields[[nm]]$span)
        sc <- as.character(ds$axes[[sp[1]]]$labels)
        sg <- as.character(ds$axes[[sp[2]]]$labels)
        tr <- as(as(ds$fields[[nm]]$values, "CsparseMatrix"), "TsparseMatrix")   # 0-based (i, j, x) triplets
        ii <- c(ii, cell_pos[sc[tr@i + 1L]]); jj <- c(jj, gene_pos[sg[tr@j + 1L]]); xx <- c(xx, tr@x)
        soc[sc] <- sub("^cells\\.", "", sp[1])
        st <- ds$fields[[nm]]$state %||% st
        ds$fields[[nm]] <- NULL
      }
      M <- Matrix::sparseMatrix(i = ii, j = jj, x = xx, dims = c(length(cells), length(genes)),
                                dimnames = list(cells, genes))
      ds$fields[[root]] <- list(values = M, role = "measure", span = c("cells", "genes"), state = st)
    }
    split_by <- soc
    # per-sample latent spaces (pca.<s>) are sample-local rotations, not a joint reduction: not
    # representable as one Seurat DimReduc over all cells, so drop them rather than fabricate a joint space.
    for (nm in names(ds$fields)) {
      sp <- as.character(ds$fields[[nm]]$span)
      if (identical(ds$fields[[nm]]$role, "embedding") && length(sp) == 2 && startsWith(sp[1], "cells.")) {
        ds$dropped <- c(ds$dropped, nm); ds$fields[[nm]] <- NULL
      }
    }
  }

  meas <- .fields_over(ds, c("cells", "genes"), role = "measure")
  full_meas <- Filter(function(nm) is.null(ds$fields[[nm]]$index), meas)   # partial measures (subset
  counts_nm <- .pick(ds, full_meas, prefer_state = "raw", prefer_name = "counts")  # scale.data) can't be
  data_nm <- .pick(ds, full_meas, prefer_state = "lognorm", prefer_name = "X", exclude = counts_nm)  # the
  # full primary/data layer -> restored separately below (scale.data) over their covered genes.
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
  # scale.data (state "scaled"), possibly over a *gene subset* (partial coverage) -> restore the layer on
  # its covered genes (real Seurat keeps scale.data over the variable features only).
  scaled_nm <- NULL
  for (nm in setdiff(meas, c(counts_nm, data_nm)))
    if (identical(ds$fields[[nm]]$state, "scaled")) { scaled_nm <- nm; break }
  if (!is.null(scaled_nm)) {
    sf <- ds$fields[[scaled_nm]]
    sgenes <- if (!is.null(sf$index)) genes[as.integer(sf$index) + 1L] else genes
    sm <- Matrix::t(as.matrix(sf$values)); dimnames(sm) <- list(sgenes, cells)
    SeuratObject::LayerData(so, assay = "RNA", layer = "scale.data") <- sm
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
    if (identical(f$subtype, "spatial")) next        # spatial coords aren't a DimReduc (and may be partial)
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
        lspan <- as.character(ds$fields[[load_nm]]$span)      # loadings may be over a *subset* feature
        lfeat <- if (lspan[1] == "genes" || is.null(ds$axes[[lspan[1]]])) genes   # axis (variable feats)
                 else as.character(ds$axes[[lspan[1]]]$labels)
        rownames(L) <- lfeat
        colnames(L) <- colnames(emb)
        dr <- SeuratObject::CreateDimReducObject(embeddings = emb, loadings = L,
                                                 key = key, assay = "RNA", stdev = sdv)
      } else {
        dr <- SeuratObject::CreateDimReducObject(embeddings = emb, key = key, assay = "RNA", stdev = sdv)
      }
      so[[nm]] <- dr
    }
  }

  # joint cell-cell graphs (relation over (cells, cells)) -> Seurat Graph objects. For a Conos collection
  # this is the integration graph -- the substance of the joint analysis, carried natively (no corrected
  # expression needed; Seurat stores graphs independently of any assay matrix).
  for (nm in names(ds$fields)) {
    f <- ds$fields[[nm]]
    if (identical(f$role, "relation") && length(f$span) == 2 && all(as.character(f$span) == "cells")) {
      A <- as(f$values, "CsparseMatrix"); dimnames(A) <- list(cells, cells)
      g <- SeuratObject::as.Graph(A)
      SeuratObject::DefaultAssay(g) <- "RNA"
      so[[nm]] <- g
    }
  }

  if (!is.null(ds$fields[["ident"]])) {            # restore the active identity (Idents)
    iv <- ds$fields[["ident"]]$values
    if (!is.factor(iv)) iv <- as.factor(iv)
    SeuratObject::Idents(so) <- stats::setNames(iv, cells)
  }

  # Multimodal: rebuild every non-default assay from its measures over a canonical feature axis
  # (`proteins`/`peaks`/...). Group by the source assay name in `provenance$assay` (so an ADT assay on
  # the `proteins` axis comes back as "ADT"); fall back to the axis name for older stores that lack it.
  other_meas <- Filter(function(nm) { f <- ds$fields[[nm]]; sp <- as.character(f$span)
    identical(f$role, "measure") && length(sp) == 2 && sp[1] == "cells" && sp[2] != "genes" &&
      !is.null(ds$axes[[sp[2]]]) && identical(ds$axes[[sp[2]]]$role, "feature") }, names(ds$fields))
  if (length(other_meas)) {
    assay_of <- vapply(other_meas, function(nm)
      ds$fields[[nm]]$provenance$assay %||% as.character(ds$fields[[nm]]$span)[2], character(1))
    for (aname in unique(assay_of)) {
      fm <- other_meas[assay_of == aname]
      fax <- as.character(ds$fields[[fm[1]]]$span)[2]
      feats <- as.character(ds$axes[[fax]]$labels)
      fxc <- function(nm) { mm <- Matrix::t(as(ds$fields[[nm]]$values, "CsparseMatrix")); dimnames(mm) <- list(feats, cells); mm }
      states <- vapply(fm, function(nm) ds$fields[[nm]]$state %||% "", character(1))
      primary <- fm[match("raw", states)]; if (is.na(primary)) primary <- fm[1]
      aobj <- SeuratObject::CreateAssay5Object(counts = fxc(primary))
      for (nm in setdiff(fm, primary)) {
        lyr <- ds$fields[[nm]]$provenance$layer %||% sub(paste0("^", aname, "\\."), "", nm)  # "ADT.data" -> "data"
        SeuratObject::LayerData(aobj, layer = lyr) <- fxc(nm)
      }
      so[[aname]] <- aobj
    }
  }

  if (!is.null(split_by) && inherits(so[["RNA"]], "Assay5")) {   # re-split a v5 integration object
    so[["RNA"]] <- split(so[["RNA"]], f = unname(split_by[cells]))
  }
  if (length(ds$dropped)) so@misc$lstar_dropped <- unique(as.character(ds$dropped))  # never silently lose
  so
}

# ---- Seurat v2 (pre-Assay) ---------------------------------------------------------------------------
# A "very old" serialized Seurat object is the lowercase S4 class `seurat` (NOT `Seurat`): it predates the
# Assay/Assay5 classes and the SeuratObject package, so none of the SeuratObject accessors work on it. Its
# data live in fixed slots -- genes x cells matrices `raw.data`/`data`/`scale.data`, a `dr` list of
# `dim.reduction` objects (each `cell.embeddings`/`gene.loadings`/`sdev`/`key`), `meta.data`, an `ident`
# factor, an `snn` graph, an `assay` list for multimodal -- read here via attr() because S4 slots ARE
# stored as attributes, so they survive `readRDS()` even when the ancient `seurat` class isn't defined in a
# modern R. Read-only: write_seurat always emits a modern object (converting old -> new is the whole point).
.is_seurat_v2 <- function(so) inherits(so, "seurat") && !inherits(so, "Seurat")

.read_seurat_v2 <- function(so) {
  S <- function(nm) tryCatch(attr(so, nm), error = function(e) NULL)
  raw <- S("raw.data"); dat <- S("data"); scl <- S("scale.data")
  ref <- if (!is.null(dat) && length(dat)) dat else if (!is.null(raw) && length(raw)) raw else scl
  if (is.null(ref)) stop("read_seurat: Seurat v2 object carries no raw.data/data/scale.data matrix")
  cells <- S("cell.names"); if (is.null(cells) || !length(cells)) cells <- colnames(ref)
  cells <- as.character(cells); genes <- as.character(rownames(ref))
  ver <- tryCatch(as.character(S("version")), error = function(e) NA_character_)
  if (length(ver) != 1 || is.na(ver)) ver <- "2"
  ds <- list(kind = "sample", spec_version = "0.1",
             profiles = c("seurat@0.1", "object@seurat", paste0("Seurat_pkg@", ver)),
             dropped = character(0), axes = list(), fields = list())
  ds$axes$cells <- list(labels = cells, origin = "observed", role = "observation")
  ds$axes$genes <- list(labels = genes, origin = "observed", role = "feature")

  add <- function(nm, values, role, span, state = "", subtype = "") {
    ds$fields[[nm]] <<- list(role = role, span = span, state = state, subtype = subtype, values = values)
  }
  add_factor <- function(nm, v, span) {
    v <- as.factor(v)
    ds$fields[[nm]] <<- list(role = "label", span = span, state = "", subtype = "",
                             encoding = "categorical", values = v)
    if (is.null(ds$axes[[nm]]))
      ds$axes[[nm]] <<- list(labels = levels(v), origin = "derived", role = "factor", induced_by = nm)
  }
  # a genes x cells measure -> (cells, genes); a gene SUBSET (scale.data over var.genes only) is kept as a
  # typed partial-coverage measure (index into the gene axis), mirroring the v3/v5 scale.data path.
  add_measure <- function(nm, m, state) {
    if (is.null(m) || !length(m) || !nrow(m)) return(invisible())
    if (nrow(m) == length(genes)) {
      add(nm, Matrix::t(m), "measure", c("cells", "genes"), state = state)
    } else {
      gi <- match(as.character(rownames(m)), genes)
      if (!anyNA(gi)) {
        ds$fields[[nm]] <<- list(role = "measure", span = c("cells", "genes"), state = state, subtype = "",
                                 values = Matrix::t(m), coverage = "partial",
                                 index = as.integer(gi - 1L), index_axis = "genes")
      } else ds$dropped <<- c(ds$dropped, sprintf("%s (%d of %d genes, unnamed)", nm, nrow(m), length(genes)))
    }
  }
  add_measure("counts", raw, "raw")
  add_measure("X", dat, "lognorm")
  add_measure("scale.data", scl, "scaled")

  # dr: a named list of `dim.reduction` (pca/tsne/umap/ica/...) -> an embedding axis + scores, the per-dim
  # stdev as a measure, and gene loadings (over all genes, or over a subset feature axis like real v2 PCA).
  dr <- S("dr")
  if (is.list(dr)) for (rn in names(dr)) {
    d <- dr[[rn]]
    emb <- tryCatch(attr(d, "cell.embeddings"), error = function(e) NULL)
    if (is.null(emb) || !nrow(emb)) next
    emb <- as.matrix(emb); dimn <- colnames(emb)
    if (is.null(dimn)) dimn <- paste0(toupper(rn), seq_len(ncol(emb)))
    er <- rownames(emb); ord <- if (!is.null(er) && setequal(er, cells)) match(cells, er) else seq_len(nrow(emb))
    ds$axes[[rn]] <- list(labels = dimn, origin = "derived", role = "coordinate")
    add(rn, unname(emb[ord, , drop = FALSE]), "embedding", c("cells", rn))
    sdv <- tryCatch(as.numeric(attr(d, "sdev")), error = function(e) numeric(0))
    if (length(sdv) == ncol(emb)) add(paste0(rn, "_stdev"), sdv, "measure", rn)
    ld <- tryCatch(attr(d, "gene.loadings"), error = function(e) NULL)
    if (!is.null(ld) && nrow(ld) > 0 && ncol(ld) == ncol(emb)) {
      ld <- as.matrix(ld)
      if (nrow(ld) == length(genes)) {
        add(paste0(rn, "_loadings"), unname(ld), "loading", c("genes", rn))
      } else {
        hvg <- rownames(ld)
        if (!is.null(hvg) && length(hvg) == nrow(ld)) {
          fax <- paste0(rn, "_features")
          ds$axes[[fax]] <- list(labels = as.character(hvg), origin = "derived", role = "feature")
          add(paste0(rn, "_loadings"), unname(ld), "loading", c(fax, rn))
        } else ds$dropped <- c(ds$dropped, sprintf("loadings/%s (%d features, unnamed)", rn, nrow(ld)))
      }
    }
  }

  # meta.data (cells x columns): numeric -> measure, factor -> categorical (induces a factor axis), else label
  md <- S("meta.data")
  if (is.data.frame(md) && nrow(md)) {
    mdr <- rownames(md); ord <- if (!is.null(mdr) && setequal(mdr, cells)) match(cells, mdr) else seq_len(nrow(md))
    for (col in colnames(md)) {
      v <- md[[col]][ord]
      if (is.numeric(v)) add(col, as.numeric(v), "measure", "cells")
      else if (is.factor(v)) add_factor(col, v, "cells")
      else add(col, as.character(v), "label", "cells")
    }
  }

  # ident: the active identity (a factor over cells) -> a categorical `ident` label, flagged active_ident
  id <- S("ident")
  if (!is.null(id) && length(id)) {
    idr <- names(id); ord <- if (!is.null(idr) && setequal(idr, cells)) match(cells, idr) else seq_along(id)
    idf <- droplevels(as.factor(id[ord]))
    if (length(idf) == length(cells)) {
      ds$fields[["ident"]] <- list(role = "label", span = "cells", state = "", subtype = "active_ident",
                                   encoding = "categorical", values = idf)
      if (is.null(ds$axes[["ident"]]))
        ds$axes[["ident"]] <- list(labels = levels(idf), origin = "derived", role = "factor", induced_by = "ident")
    }
  }

  # snn (dgCMatrix, cells x cells) -> a cell-cell relation; var.genes -> a per-gene 0/1 measure
  snn <- S("snn")
  if (!is.null(snn) && length(dim(snn)) == 2 && all(dim(snn) == length(cells)))
    add("snn", snn, "relation", c("cells", "cells"))
  vg <- S("var.genes")
  if (!is.null(vg) && length(vg)) add("variable_features", as.numeric(genes %in% as.character(vg)), "measure", "genes")

  # multimodal: an `assay` list (each a v2 `assay` object: raw.data/data over its own features) -> a second
  # feature space over the shared cells axis (the v2 analogue of a CITE-seq ADT assay).
  asy <- S("assay")
  if (is.list(asy)) for (an in names(asy)) {
    a <- asy[[an]]
    araw <- tryCatch(attr(a, "raw.data"), error = function(e) NULL)
    adat <- tryCatch(attr(a, "data"), error = function(e) NULL)
    aref <- if (!is.null(adat) && length(adat)) adat else araw
    if (is.null(aref) || !nrow(aref)) next
    if (an %in% names(ds$axes)) { ds$dropped <- c(ds$dropped, paste0("assay/", an, " (axis-name clash)")); next }
    ds$axes[[an]] <- list(labels = as.character(rownames(aref)), origin = "observed", role = "feature")
    if (!is.null(araw) && nrow(araw) == nrow(aref)) add(paste0(an, ".counts"), Matrix::t(araw), "measure", c("cells", an), state = "raw")
    if (!is.null(adat) && nrow(adat) == nrow(aref)) add(paste0(an, ".data"), Matrix::t(adat), "measure", c("cells", an), state = "lognorm")
  }

  # analysis cruft with no cross-format home -> recorded in `dropped`, never silently lost
  for (slt in c("cluster.tree", "calc.params", "kmeans", "hvg.info", "imputed", "spatial", "misc"))
    if (length(S(slt))) ds$dropped <- c(ds$dropped, paste0("seurat-v2/", slt))

  class(ds) <- "lstar_dataset"
  ds
}

# ---- Seurat v3/v4/v5 PACKAGE-FREE read (base R, no SeuratObject) --------------------------------------
# The `--backend direct` fallback used by `lstar convert` when SeuratObject isn't installed: read a modern
# Seurat .rds by walking its S4 slots via attr() (the same mechanism the v2 path uses), so only base R +
# Matrix are needed. Produces the same core L* dataset `read_seurat` builds (counts/data/scale.data
# measures, reductions, meta.data, active ident). v5 layer names come from the `cells`/`features` LogMaps
# (`rownames(attr(assay,"cells"))`). A package-backed matrix slot (BPCells/DelayedArray) can't be decoded
# packagelessly -> a clear error names the package. Verified value-equal to native by convert_cli.sh.
.read_seurat_direct <- function(so) {
  S <- function(obj, nm) tryCatch(attr(obj, nm), error = function(e) NULL)
  assays <- S(so, "assays")
  if (is.null(assays) || !length(assays))
    stop("read_seurat (direct): object has no 'assays' slot -- not a Seurat object ",
         "(a SingleCellExperiment? install SingleCellExperiment and re-run for the native path)")
  aname <- S(so, "active.assay")
  if (is.null(aname) || !length(aname) || !nzchar(aname[1])) aname <- names(assays)[1]
  a <- assays[[aname]]
  asparse <- function(m) {
    if (is.null(m) || !length(m)) return(NULL)
    if (methods::is(m, "externalptr") || !(is.matrix(m) || methods::is(m, "Matrix")))
      stop(sprintf("read_seurat (direct): assay '%s' holds a package-backed matrix (e.g. BPCells / ",
                   aname), "DelayedArray) that can't be decoded without its package -- install it and ",
           "re-run (lstar will then use the native path)")
    m
  }
  layers <- S(a, "layers")
  if (!is.null(layers) && length(layers)) {            # v5 Assay5: names from the cells/features LogMaps
    cl <- rownames(S(a, "cells")); ft <- rownames(S(a, "features"))
    getL <- function(L) {
      m <- asparse(layers[[L]]); if (is.null(m)) return(NULL)
      if (is.null(rownames(m)) && nrow(m) == length(ft) && ncol(m) == length(cl)) dimnames(m) <- list(ft, cl)
      m
    }
    raw <- getL("counts"); dat <- getL("data"); scl <- getL("scale.data"); cells <- cl; genes <- ft
  } else {                                             # v3/v4 Assay
    raw <- asparse(S(a, "counts")); dat <- asparse(S(a, "data")); scl <- asparse(S(a, "scale.data"))
    ref <- if (!is.null(dat)) dat else if (!is.null(raw)) raw else scl
    if (is.null(ref)) stop("read_seurat (direct): assay carries no counts/data/scale.data matrix")
    cells <- colnames(ref); genes <- rownames(ref)
  }
  cells <- as.character(cells); genes <- as.character(genes)
  ver <- tryCatch(as.character(S(so, "version")), error = function(e) NA_character_)
  if (!length(ver) || is.na(ver[1])) ver <- "5"

  ds <- list(kind = "sample", spec_version = "0.1",
             profiles = c("seurat@0.1", "object@seurat", paste0("Seurat_pkg@", ver[1])),
             dropped = character(0), axes = list(), fields = list())
  ds$axes$cells <- list(labels = cells, origin = "observed", role = "observation")
  ds$axes$genes <- list(labels = genes, origin = "observed", role = "feature")
  add <- function(nm, values, role, span, state = "", subtype = "") {
    ds$fields[[nm]] <<- list(role = role, span = span, state = state, subtype = subtype, values = values)
  }
  add_factor <- function(nm, v, span) {
    v <- as.factor(v)
    ds$fields[[nm]] <<- list(role = "label", span = span, state = "", subtype = "",
                             encoding = "categorical", values = v)
    if (is.null(ds$axes[[nm]]))
      ds$axes[[nm]] <<- list(labels = levels(v), origin = "derived", role = "factor", induced_by = nm)
  }
  add_measure <- function(nm, m, state) {
    if (is.null(m) || !length(m) || !nrow(m)) return(invisible())
    m <- as(m, "CsparseMatrix")
    if (nrow(m) == length(genes)) {
      add(nm, Matrix::t(m), "measure", c("cells", "genes"), state = state)
    } else {
      gi <- match(as.character(rownames(m)), genes)
      if (!anyNA(gi)) ds$fields[[nm]] <<- list(role = "measure", span = c("cells", "genes"), state = state,
        subtype = "", values = Matrix::t(m), coverage = "partial", index = as.integer(gi - 1L), index_axis = "genes")
      else ds$dropped <<- c(ds$dropped, sprintf("%s (%d of %d genes, unnamed)", nm, nrow(m), length(genes)))
    }
  }
  add_measure("counts", raw, "raw"); add_measure("X", dat, "lognorm"); add_measure("scale.data", scl, "scaled")

  dr <- S(so, "reductions")                            # v3-v5 slot names: cell.embeddings / stdev / feature.loadings
  if (is.list(dr)) for (rn in names(dr)) {
    d <- dr[[rn]]; emb <- S(d, "cell.embeddings")
    if (is.null(emb) || !nrow(emb)) next
    emb <- as.matrix(emb); dimn <- colnames(emb)
    if (is.null(dimn)) dimn <- paste0(toupper(rn), seq_len(ncol(emb)))
    er <- rownames(emb); ord <- if (!is.null(er) && setequal(er, cells)) match(cells, er) else seq_len(nrow(emb))
    ds$axes[[rn]] <- list(labels = dimn, origin = "derived", role = "coordinate")
    add(rn, unname(emb[ord, , drop = FALSE]), "embedding", c("cells", rn))
    sdv <- tryCatch(as.numeric(S(d, "stdev")), error = function(e) numeric(0))
    if (length(sdv) == ncol(emb)) add(paste0(rn, "_stdev"), sdv, "measure", rn)
    ld <- S(d, "feature.loadings")
    if (!is.null(ld) && nrow(ld) > 0 && ncol(ld) == ncol(emb)) {
      ld <- as.matrix(ld)
      if (nrow(ld) == length(genes)) add(paste0(rn, "_loadings"), unname(ld), "loading", c("genes", rn))
      else { hvg <- rownames(ld); if (!is.null(hvg) && length(hvg) == nrow(ld)) {
        fax <- paste0(rn, "_features"); ds$axes[[fax]] <- list(labels = as.character(hvg), origin = "derived", role = "feature")
        add(paste0(rn, "_loadings"), unname(ld), "loading", c(fax, rn)) } }
    }
  }

  md <- S(so, "meta.data")                             # meta.data + active ident -> cell fields (same as v2)
  if (is.data.frame(md) && nrow(md)) {
    mdr <- rownames(md); ord <- if (!is.null(mdr) && setequal(mdr, cells)) match(cells, mdr) else seq_len(nrow(md))
    for (col in colnames(md)) {
      v <- md[[col]][ord]
      if (is.numeric(v)) add(col, as.numeric(v), "measure", "cells")
      else if (is.factor(v)) add_factor(col, v, "cells")
      else add(col, as.character(v), "label", "cells")
    }
  }
  id <- S(so, "active.ident")
  if (!is.null(id) && length(id)) {
    idr <- names(id); ord <- if (!is.null(idr) && setequal(idr, cells)) match(cells, idr) else seq_along(id)
    idf <- droplevels(as.factor(id[ord]))
    if (length(idf) == length(cells)) {
      ds$fields[["ident"]] <- list(role = "label", span = "cells", state = "", subtype = "active_ident",
                                   encoding = "categorical", values = idf)
      if (is.null(ds$axes[["ident"]]))
        ds$axes[["ident"]] <- list(labels = levels(idf), origin = "derived", role = "factor", induced_by = "ident")
    }
  }
  class(ds) <- "lstar_dataset"
  ds
}

# ---- Seurat PACKAGE-FREE write (base R, no SeuratObject) ----------------------------------------------
# The `--backend direct` write fallback: build a native-valid Seurat object from an L* dataset using a
# PINNED v4 SeuratObject schema (slot names/types lifted verbatim from SeuratObject::getSlots), `new()`d in
# base R, with the S4 class identity forged to "SeuratObject" so a real SeuratObject session reconstructs
# it as a genuine Seurat object and its tools accept it (verified by the native-acceptance check). Only the
# WRITE direction is hard packagelessly (reading just needs attr()); this is the pinned-version answer.
.forge_s4 <- function(obj) {                            # claim SeuratObject ownership of the S4 class
  cl <- attr(obj, "class"); attr(cl, "package") <- "SeuratObject"; attr(obj, "class") <- cl; obj
}

# A private environment to hold the pinned-schema S4 class defs. We do NOT use globalenv() (CRAN forbids
# modifying it) and cannot use the package namespace (locked at runtime); this env is ours and unlocked.
# S4 classes register in the global class table regardless of `where`, so new()/isClass still find them.
.lstar_seurat_classes <- new.env(parent = emptyenv())

.seurat_pinned_classes <- function() {
  # Guarded: only define when the class isn't already present (so we never clash with a loaded
  # SeuratObject). Direct-write only runs when SeuratObject is absent, so these are normally fresh.
  if (!methods::isClass("JackStrawData"))
    methods::setClass("JackStrawData", where = .lstar_seurat_classes, methods::representation(empirical.p.values = "matrix",
      fake.reduction.scores = "matrix", empirical.p.values.full = "matrix", overall.p.values = "matrix"))
  if (!methods::isClass("DimReduc"))
    methods::setClass("DimReduc", where = .lstar_seurat_classes, methods::representation(cell.embeddings = "matrix", feature.loadings = "matrix",
      feature.loadings.projected = "matrix", assay.used = "character", global = "logical",
      stdev = "numeric", jackstraw = "JackStrawData", misc = "list", key = "character"))
  if (!methods::isClass("Assay"))
    methods::setClass("Assay", where = .lstar_seurat_classes, methods::representation(counts = "ANY", data = "ANY", scale.data = "matrix",
      assay.orig = "ANY", var.features = "vector", meta.features = "data.frame", misc = "ANY",
      key = "character"))
  if (!methods::isClass("LogMap"))                       # v5: a logical (entities x layers) membership map
    methods::setClass("LogMap", where = .lstar_seurat_classes, contains = "matrix")
  if (!methods::isClass("Assay5"))
    methods::setClass("Assay5", where = .lstar_seurat_classes, methods::representation(layers = "list", cells = "LogMap",
      features = "LogMap", default = "integer", assay.orig = "character", meta.data = "data.frame",
      misc = "list", key = "character"))
  if (!methods::isClass("Seurat"))
    methods::setClass("Seurat", where = .lstar_seurat_classes, methods::representation(assays = "list", meta.data = "data.frame",
      active.assay = "character", active.ident = "factor", graphs = "list", neighbors = "list",
      reductions = "list", images = "list", project.name = "character", misc = "list",
      version = "ANY", commands = "list", tools = "list"))
}

.build_seurat_direct <- function(ds) {
  .seurat_pinned_classes()
  em <- function() matrix(numeric(0), 0, 0)
  cells <- as.character(ds$axes$cells$labels); genes <- as.character(ds$axes$genes$labels); nc <- length(cells)
  gxc <- function(nm) { m <- Matrix::t(as(ds$fields[[nm]]$values, "CsparseMatrix")); dimnames(m) <- list(genes, cells); m }
  is_full_meas <- function(nm) { f <- ds$fields[[nm]]; identical(f$role, "measure") &&
      length(f$span) == 2 && all(c("cells", "genes") == as.character(f$span)) && is.null(f$index) }
  meas <- Filter(is_full_meas, names(ds$fields))
  pick <- function(state, name) { for (nm in meas) { f <- ds$fields[[nm]]
      if (identical(f$state, state) || identical(nm, name)) return(nm) }; NULL }
  counts_nm <- pick("raw", "counts"); data_nm <- pick("lognorm", "X")
  # Build a v5 Assay5 with explicit layers (counts / data / scale.data), so a counts-only assay has NO
  # data layer -- matching the native writer and avoiding a spurious lognorm `X` (a v4 Assay always
  # surfaces a `data` layer). Layer membership is carried by the cells/features LogMaps.
  layers <- list(); lnames <- character(0); used <- character(0)
  # a real Assay5 layer is a bare matrix (NULL dimnames); cell/feature names live in the cells/features LogMaps
  addL <- function(lname, nm) { m <- gxc(nm); dimnames(m) <- NULL; layers[[lname]] <<- m
    lnames <<- c(lnames, lname); used <<- c(used, nm) }
  primary <- if (!is.null(counts_nm)) counts_nm else data_nm   # primary -> the counts layer (as CreateSeuratObject does)
  if (is.null(primary)) stop("write_seurat (direct): no full (cells x genes) measure to build an assay")
  addL("counts", primary)
  if (!is.null(counts_nm) && !is.null(data_nm)) addL("data", data_nm)   # a distinct lognorm layer only if both exist
  scaled_nm <- NULL
  for (nm in setdiff(meas, used)) if (identical(ds$fields[[nm]]$state, "scaled") && is.null(ds$fields[[nm]]$index)) { scaled_nm <- nm; break }
  if (!is.null(scaled_nm)) addL("scale.data", scaled_nm)
  cnt <- gxc(primary)   # for the QC bookkeeping below (genes x cells, dimnamed)
  mklm <- function(nms) .forge_s4(methods::new("LogMap", matrix(TRUE, length(nms), length(lnames), dimnames = list(nms, lnames))))
  assay <- .forge_s4(methods::new("Assay5", layers = layers, cells = mklm(cells), features = mklm(genes),
    default = 1L, assay.orig = character(0), meta.data = data.frame(row.names = genes),
    misc = list(), key = "rna_"))

  reductions <- list()
  embs <- Filter(function(nm) { f <- ds$fields[[nm]]
    identical(f$role, "embedding") && length(f$span) == 2 && as.character(f$span)[1] == "cells" }, names(ds$fields))
  for (nm in embs) {
    f <- ds$fields[[nm]]; rn <- as.character(f$span)[2]
    key <- .dimreduc_key(rn)                            # SeuratObject requires embedding colnames to start with the key
    emb <- as.matrix(f$values); rownames(emb) <- cells; colnames(emb) <- paste0(key, seq_len(ncol(emb)))
    ld_nm <- paste0(nm, "_loadings"); sd_nm <- paste0(nm, "_stdev")
    ld <- em()
    if (!is.null(ds$fields[[ld_nm]])) { ld <- as.matrix(ds$fields[[ld_nm]]$values)
      fax <- as.character(ds$fields[[ld_nm]]$span)[1]; rownames(ld) <- as.character(ds$axes[[fax]]$labels); colnames(ld) <- colnames(emb) }
    sdv <- if (!is.null(ds$fields[[sd_nm]])) as.numeric(ds$fields[[sd_nm]]$values) else numeric(0)
    js <- .forge_s4(methods::new("JackStrawData", empirical.p.values = em(), fake.reduction.scores = em(),
      empirical.p.values.full = em(), overall.p.values = em()))
    reductions[[nm]] <- .forge_s4(methods::new("DimReduc", cell.embeddings = emb, feature.loadings = ld,
      feature.loadings.projected = em(), assay.used = "RNA", global = FALSE, stdev = sdv,
      jackstraw = js, misc = list(), key = key))
    used <- c(used, nm, ld_nm, sd_nm)
  }

  md <- data.frame(row.names = cells); ident <- NULL
  for (nm in setdiff(names(ds$fields), used)) {
    f <- ds$fields[[nm]]
    if (length(f$span) == 1 && as.character(f$span)[1] == "cells" && f$role %in% c("label", "measure")) {
      v <- f$values
      if (identical(f$subtype, "active_ident")) { ident <- as.factor(v); next }
      md[[nm]] <- if (is.factor(v)) v else if (identical(f$role, "measure")) as.numeric(v) else as.character(v)
    }
  }
  # CreateSeuratObject bookkeeping (match the native writer): default identity + per-cell QC totals,
  # added only when the store didn't already carry them (a Seurat -> store -> Seurat round-trip would).
  if (is.null(md[["orig.ident"]])) md[["orig.ident"]] <- factor(rep("SeuratProject", nc))
  if (is.null(md[["nCount_RNA"]])) md[["nCount_RNA"]] <- as.numeric(Matrix::colSums(cnt))
  if (is.null(md[["nFeature_RNA"]])) md[["nFeature_RNA"]] <- as.numeric(Matrix::colSums(cnt > 0))
  if (is.null(ident)) ident <- as.factor(md[["orig.ident"]])
  names(ident) <- cells

  .forge_s4(methods::new("Seurat", assays = list(RNA = assay), meta.data = md, active.assay = "RNA",
    active.ident = ident, graphs = list(), neighbors = list(), reductions = reductions, images = list(),
    project.name = "lstar", misc = list(), version = as.package_version("5.0.0"),
    commands = list(), tools = list()))
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
  if (.is_seurat_v2(so)) return(.read_seurat_v2(so))   # pre-Assay v2 `seurat` class -> dedicated slot reader
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
      # a layer over a *feature subset* (real objects keep scale.data over the variable features only,
      # e.g. pbmc3k.final: 2000 HVGs of 13714) -> a typed **partial-coverage** measure over (cells, genes)
      # keyed by an index into the gene axis (0-based), not dropped and not mis-spanned over all genes.
      gi <- match(rownames(m), genes)
      if (!anyNA(gi)) {
        ds$fields[[name_of(L)]] <- list(role = "measure", span = c("cells", "genes"), state = state_of(L),
                                        subtype = "", values = Matrix::t(m), coverage = "partial",
                                        index = as.integer(gi - 1L), index_axis = "genes")
      } else {                                                  # unnamed subset -> can't place it; record
        ds$dropped <- c(ds$dropped, sprintf("layer/%s (%d of %d features)", L, nrow(m), length(genes)))
      }
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
    cax <- .modality_axis(a, names(ds$axes))                    # ADT->proteins, ATAC->peaks (mudata vocab)
    if (cax %in% names(ds$axes)) { ds$dropped <- c(ds$dropped, paste0("assay/", a, " (axis-name clash)")); next }
    ds$axes[[cax]] <- list(labels = feats, origin = "observed", role = "feature")
    laca <- .seurat_layer_access(so, a)
    for (L in laca$names) {
      m <- laca$get(L)
      if (nrow(m) != length(feats)) {                           # subset layer -> partial coverage
        ds$dropped <- c(ds$dropped, sprintf("assay/%s layer/%s (%d of %d features)", a, L, nrow(m), length(feats)))
        next
      }
      ds$fields[[paste0(a, ".", L)]] <- list(role = "measure", span = c("cells", cax), state = state_of(L),
        subtype = "", values = Matrix::t(m), provenance = list(assay = a, layer = L, feature_axis = cax))
    }
    capture_chromatin(a, cax, length(feats))                    # ATAC as a non-default assay (peaks axis)
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
      if (nrow(ld) == length(genes)) {                   # loadings over all genes -> a loading field
        add(paste0(rn, "_loadings"), unname(as.matrix(ld)), "loading", c("genes", rn))
      } else {
        # Real PCA loadings are over the variable features only (e.g. pbmc3k.final: 2000 of 13714).
        # Capture them faithfully over a *subset feature axis* (the HVG names) instead of dropping --
        # the values are real and worth keeping. write_seurat reads the loading's own feature axis back.
        hvg <- rownames(ld)
        if (!is.null(hvg) && length(hvg) == nrow(ld)) {
          fax <- paste0(rn, "_features")
          ds$axes[[fax]] <- list(labels = as.character(hvg), origin = "derived", role = "feature")
          add(paste0(rn, "_loadings"), unname(as.matrix(ld)), "loading", c(fax, rn))
        } else {                                         # no feature names -> can't place it; record loss
          ds$dropped <- c(ds$dropped,
                          sprintf("loadings/%s (%d of %d features, unnamed)", rn, nrow(ld), length(genes)))
        }
      }
    }
  }

  # Neighbor objects (FindNeighbors(return.neighbor=TRUE)): a kNN result as nn.idx (cells x k) + nn.dist.
  # Type each as a weighted cell-cell relation (distance-weighted), the analogue of a stored Graph, so
  # it isn't silently lost.
  for (nn_nm in tryCatch(SeuratObject::Neighbors(so), error = function(e) character(0))) {
    nn <- tryCatch(so[[nn_nm]], error = function(e) NULL)
    idx <- tryCatch(SeuratObject::Indices(nn), error = function(e) NULL)
    if (is.null(idx) || nrow(idx) != length(cells)) next
    dst <- tryCatch(SeuratObject::Distances(nn), error = function(e) NULL)
    x <- if (!is.null(dst) && all(dim(dst) == dim(idx))) as.numeric(dst) else rep(1, length(idx))
    m <- Matrix::sparseMatrix(i = rep(seq_len(nrow(idx)), times = ncol(idx)), j = as.integer(idx),
                              x = x, dims = c(length(cells), length(cells)))
    add(paste0("nn_", nn_nm), m, "relation", c("cells", "cells"))
  }

  # Stored cell-cell graphs (SeuratObject::Graphs(): FindNeighbors() nn/snn graphs, or an externally
  # attached integration graph such as Conos's joint kNN written by write_seurat()). Each is a cells x
  # cells sparse adjacency; carry it as a `relation` field (the inverse of write_seurat's Graph emission)
  # so a joint/integration graph survives a Seurat round-trip instead of being silently dropped.
  for (gn in tryCatch(SeuratObject::Graphs(so), error = function(e) character(0))) {
    g <- tryCatch(so[[gn]], error = function(e) NULL)
    A <- tryCatch(methods::as(g, "CsparseMatrix"), error = function(e) NULL)
    if (is.null(A) || nrow(A) != length(cells) || ncol(A) != length(cells)) next
    rn <- rownames(A)
    if (!is.null(rn) && !identical(rn, cells) && all(cells %in% rn)) A <- A[cells, cells, drop = FALSE]
    add(gn, methods::as(A, "CsparseMatrix"), "relation", c("cells", "cells"))
  }

  # Spatial images (Visium / Slide-seq / FOV): tissue coordinates live in `so@images`, NOT in Reductions,
  # so without this they'd be silently lost. Capture each image's coordinates as a `spatial` *observed*
  # coordinate axis (mirroring the AnnData obsm['spatial'] path); coordinates over a cell *subset* (a
  # multi-section object) use partial coverage. Conceptual spatial support only -- the pixel images
  # themselves are a deferred tier and recorded in `dropped`, not typed.
  imgs <- tryCatch(SeuratObject::Images(so), error = function(e) character(0))
  for (img in imgs) {
    co <- tryCatch(as.data.frame(SeuratObject::GetTissueCoordinates(so[[img]])), error = function(e) NULL)
    if (is.null(co) || !nrow(co)) { ds$dropped <- c(ds$dropped, paste0("image/", img, " (coords unreadable)")); next }
    cc <- intersect(c("cell", "cells", "barcode"), colnames(co))
    icells <- if (length(cc)) as.character(co[[cc[1]]]) else rownames(co)
    xy <- as.matrix(co[, vapply(co, is.numeric, logical(1)), drop = FALSE])     # x / y (/ z) numeric cols
    pos <- match(icells, cells)
    if (is.null(icells) || ncol(xy) == 0 || anyNA(pos)) {
      ds$dropped <- c(ds$dropped, paste0("image/", img, " (coords unmapped)")); next
    }
    sax <- if (length(imgs) > 1) paste0("spatial.", img) else "spatial"
    ds$axes[[sax]] <- list(labels = paste0(sax, seq_len(ncol(xy))), origin = "observed", role = "coordinate")
    fld <- if (length(imgs) > 1) paste0("spatial_", img) else "spatial"
    if (length(icells) == length(cells) && all(pos == seq_along(cells))) {
      add2(fld, unname(xy), "embedding", c("cells", sax), subtype = "spatial")
    } else {                                                                     # subset -> partial coverage
      ds$fields[[fld]] <- list(role = "embedding", span = c("cells", sax), state = "", subtype = "spatial",
                               values = unname(xy), coverage = "partial",
                               index = as.integer(pos - 1L), index_axis = "cells")
    }
    ds$dropped <- c(ds$dropped, paste0("image/", img, "/pixels (deferred spatial tier)"))   # images not typed
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
