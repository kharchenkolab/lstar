# Sweep REAL 10x Visium Seurat spatial objects (SeuratData::stxBrain, 4 sections) through
# read_seurat -> lstar_write -> write_seurat, and RECORD what spatial structure survives on the
# Seurat side.
#
# stxBrain ships 4 sections (anterior1/2, posterior1/2) = a multi-section spatial COLLECTION. Each is a
# Seurat object whose Visium coordinates + tissue images live in `so@images[[<slice>]]` (a VisiumV1/V2
# SpatialImage object: GetTissueCoordinates() x/y + ScaleFactors() spot/fiducial/hires/lowres), NOT in
# `so@reductions`.
#
# As of fix #6 (14b0225), the lstar Seurat profile reads `so@images` and -- mirroring the AnnData path
# (obsm['spatial'] -> typed observed coordinate axis) -- captures the coordinates as a `spatial` observed
# coordinate axis (subtype `spatial`; multi-section subsets use partial coverage), recording the pixel
# images in `ds$dropped`. This sweep CONFIRMS that capture on real objects: it loads each section and
# records that a `spatial` coordinate axis / field appears in the dataset and that the image loss is
# recorded (not silent). (Before #6 the Seurat side silently dropped these coords; this sweep was the
# gap detector that surfaced it.)
#
# Datasets are SeuratData packages (installed into .Rlib, local-only -- not committed, not in CI).
# Run:  R_LIBS=.Rlib Rscript conformance/sweep/sweep_spatial.R   (writes /tmp/sweep_spatial_seurat.tsv)

.libPaths(c("/home/pkharchenko/p21/lstar/.Rlib", .libPaths()))
suppressMessages({library(Seurat); library(SeuratObject); library(SeuratData); library(lstar)})

# install stxBrain / ssHippo on demand (local-only)
for (d in c("stxBrain", "ssHippo")) {
  if (!(d %in% tryCatch(InstalledData()$Dataset, error = function(e) character(0)))) {
    cat("[install]", d, "...\n"); flush.console()
    tryCatch(InstallData(d), error = function(e) cat("install err:", conditionMessage(e), "\n"))
  }
}

# (dataset-package, item, platform) tuples to sweep
jobs <- list(
  c("stxBrain.SeuratData", "anterior1",  "Visium"),
  c("stxBrain.SeuratData", "anterior2",  "Visium"),
  c("stxBrain.SeuratData", "posterior1", "Visium"),
  c("stxBrain.SeuratData", "posterior2", "Visium"),
  c("ssHippo.SeuratData",  "ssHippo",    "SlideSeqV2")
)

rep <- file("/tmp/sweep_spatial_seurat.tsv", "w")
cat("dataset\titem\tplatform\tstatus\tdim\tn_images\timg_class\tcoords_in_ds\tspatial_axis\trt_img\taxes\tdropped\tnote\n",
    file = rep)
ok <- 0; total <- 0

for (j in jobs) {
  pkg <- j[1]; on <- j[2]; platform <- j[3]; total <- total + 1
  cat(sprintf("[spatial] %-22s %-12s ... ", pkg, on)); flush.console()
  r <- tryCatch({
    suppressWarnings(suppressMessages(data(list = on, package = pkg)))
    so <- UpdateSeuratObject(get(on))
    imgs <- tryCatch(Images(so), error = function(e) character(0))
    img_cls <- if (length(imgs)) class(so@images[[imgs[1]]])[1] else ""
    dimstr <- paste(dim(so), collapse = "x")

    ds <- read_seurat(so)
    # does the dataset carry the spatial coordinates in any form?
    has_spatial_axis <- "spatial" %in% names(ds$axes)
    coord_fields <- Filter(function(nm) {
      f <- ds$fields[[nm]]
      identical(f$role, "embedding") &&
        (isTRUE(identical(f$subtype, "spatial")) || nm == "spatial")
    }, names(ds$fields))
    coords_in_ds <- length(coord_fields) > 0 || has_spatial_axis

    p <- tempfile(fileext = ".zarr"); lstar_write(ds, p)
    so2 <- write_seurat(ds)
    rt_img <- length(tryCatch(Images(so2), error = function(e) character(0)))
    unlink(p, recursive = TRUE)

    # whether the profile recorded the image loss (vs dropping it silently)
    img_dropped <- any(grepl("image|images|spatial|tissue|VisiumV|SlideSeq", ds$dropped, ignore.case = TRUE))
    note <- if (length(imgs) && !coords_in_ds) {
      if (img_dropped) "coords in so@images NOT captured (recorded in dropped)"
      else "coords in so@images NOT captured (SILENTLY -- not in dropped)"
    } else if (coords_in_ds) "spatial coords captured" else ""

    list(s = "PASS", dim = dimstr, ni = length(imgs), ic = img_cls,
         cds = coords_in_ds, sa = has_spatial_axis, rti = rt_img,
         a = length(ds$axes), d = paste(head(ds$dropped, 3), collapse = ";"), n = note)
  }, error = function(e) list(s = "FAIL", dim = "", ni = "", ic = "", cds = "", sa = "",
                              rti = "", a = "", d = "", n = substr(conditionMessage(e), 1, 110)))
  if (identical(r$s, "PASS")) ok <- ok + 1
  cat(r$s, "-", r$n, "\n"); flush.console()
  cat(sprintf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
              pkg, on, platform, r$s, r$dim, r$ni, r$ic, r$cds, r$sa, r$rti, r$a, r$d, r$n),
      file = rep); flush(rep)
}
close(rep)
cat(sprintf("SPATIAL (Seurat) SWEEP: %d/%d loaded; coordinate handling recorded above\n", ok, total))
