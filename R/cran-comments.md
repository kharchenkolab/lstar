## Submission

This is a new submission. `lstar` provides a uniform data model (L\*) and a 'Zarr' interchange format
for single-cell / spatial omics, with a header-only C++ core and bidirectional converters for 'Seurat',
'SingleCellExperiment', 'Conos' and 'pagoda2' objects (the same on-disk store is also readable from
Python and C++).

## R CMD check results

`R CMD check --as-cran` gives 0 ERRORs and 0 WARNINGs. NOTEs:

* "New submission" (expected).
* "checking for future file timestamps ... unable to verify current time" — environment-specific (a
  clock/network check on the build host); it does not reflect the package and does not appear on CRAN's
  builders.

The compiled C++ core's shared object is stripped of debug symbols in `src/Makevars` (`strip -S`,
guarded), so the installed package is small (~2 MB) and there is no "installed package size" NOTE.

The vignette (`converting-formats.Rmd`, knitr/rmarkdown) requires pandoc to build; it builds on CRAN's
builders. (Local checks without pandoc were run with `--no-build-vignettes`, which produces a benign
"no files in inst/doc" WARNING that does not occur on CRAN.)

## Dependencies

All hard dependencies (Matrix, methods, stats, utils) and the LinkingTo (cpp11) are on CRAN. The format
converters use heavy single-cell packages only conditionally (declared in Suggests, all on CRAN /
Bioconductor: SeuratObject, Seurat, SingleCellExperiment, SummarizedExperiment, S4Vectors,
GenomicRanges, conos, pagoda2, ...) and the package functions, examples and tests degrade gracefully
when they are absent.

## Test environments

* local: Ubuntu Linux, R 4.4.1
* (planned before submission: win-builder devel + release, R-hub, macOS)
