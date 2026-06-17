## Submission

This is a new submission. `lstar` provides a uniform data model (L\*) and a 'Zarr' interchange format
for single-cell / spatial omics, with a header-only C++ core and bidirectional converters for 'Seurat',
'SingleCellExperiment', 'Conos' and 'pagoda2' objects (the same on-disk store is also readable from
Python and C++).

## R CMD check results

On CRAN's builders `R CMD check --as-cran` is expected to give 0 ERRORs / 0 WARNINGs / 1 NOTE
("New submission"). A local check (Ubuntu, R 4.4.1) is clean apart from items that reflect missing
optional *tools* on the local host, not the package, and which are present on CRAN's builders:

* WARNING "'qpdf' is needed for checks on size reduction of PDFs" — qpdf is not installed locally.
* NOTE "checking HTML version of manual ... no command 'tidy' found" — the HTML validator is not installed
  locally.
* NOTE "checking for future file timestamps ... unable to verify current time" — a clock/network check on
  the build host.

The compiled C++ core's shared object is stripped of debug symbols in `src/Makevars` (`strip -S`,
guarded), so the installed package is small (~2 MB) and there is no "installed package size" NOTE. The
vignette (`converting-formats.Rmd`, knitr/rmarkdown) builds cleanly with pandoc.

## Dependencies

All hard dependencies (Matrix, methods, stats, utils) and the LinkingTo (cpp11) are on CRAN. The format
converters use heavy single-cell packages only conditionally (declared in Suggests, all on CRAN /
Bioconductor: SeuratObject, Seurat, SingleCellExperiment, SummarizedExperiment, S4Vectors,
GenomicRanges, conos, pagoda2, ...) and the package functions, examples and tests degrade gracefully
when they are absent.

## Test environments

* local: Ubuntu Linux, R 4.4.1
* (planned before submission: win-builder devel + release, R-hub, macOS)
