## Submission

This is an update to the CRAN package `lstar` (0.1.0 -> 0.1.6). `lstar` provides a uniform data model
(L\*) and a 'Zarr' interchange format for single-cell / spatial omics, with a header-only C++ core
(vendored under `inst/include/lstar`) and bidirectional converters for 'Seurat', 'SingleCellExperiment',
'Conos' and 'pagoda2' objects (the same on-disk store is also readable from Python and C++).

The version is 0.1.6 (not 0.1.0). The R package version is aligned with the companion Python package
(`lstar-sc` on PyPI, 0.1.6) and the shared on-disk format; the R DESCRIPTION had lagged at 0.1.0
while the shared release line moved to 0.1.6. See `NEWS.md`.

## R CMD check results

`R CMD check --no-manual --as-cran` on the local host (Ubuntu 20.04, R 4.4.1) gives:

    Status: 1 WARNING, 3 NOTEs

Every item is either a local-host artifact (absent on CRAN's builders) or an expected, justified
NOTE. There are 0 ERRORs and no code/build WARNINGs.

### NOTE: checking CRAN incoming feasibility

* The local check reports "New maintainer / Old maintainer" only because an older 0.1.0 copy is
  installed in the local check library; this is an update of the existing CRAN package (0.1.0 ->
  0.1.6) and the maintainer is unchanged: Peter Kharchenko <pk.restricted@gmail.com>. This artifact
  does not arise on CRAN's builders, which compare against the published version.

All Suggests are on CRAN or Bioconductor. Two optional integrations target packages that are not on
a mainstream repository — the interactive viewer (`pagoda3`) and a disk-backed Seurat reader
(`BPCells`, used only by `read_seurat_backed()`). Neither is declared in DESCRIPTION; both are resolved
at call time (the package name held in a variable, dispatched via `getExportedValue()`), so the package
installs, checks and runs without them and each gives a clear install hint when invoked and absent.

### NOTE: checking pragmas in C/C++ headers and code

The vendored core header `inst/include/lstar/lstar.hpp` contains a single, tightly scoped
diagnostic-suppression block:

    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wdeprecated-declarations"
    #include "nlohmann/json.hpp"
    #pragma GCC diagnostic pop

It suppresses one warning **only** around the bundled third-party header `nlohmann/json.hpp`, which
instantiates `std::char_traits<unsigned char>` via its `std::basic_string<unsigned char>` output
adapters. libc++ on recent toolchains (Xcode 26.5+) deprecates `char_traits<T>` for non-standard
`T`, which turns into a build failure under `-Werror`-style configurations. The suppression is not
used to hide warnings in our own code (which names no such instantiation) and cannot be removed
without risking the build of the vendored JSON dependency on those toolchains. It is restored with a
matching `pop` immediately after the include. There are no other diagnostic-suppression pragmas in
the package (the remaining `#pragma` directives are `#pragma once` and OpenMP `#pragma omp`).

### NOTE / WARNING: local-tool artifacts (absent on CRAN's builders)

* NOTE "checking for future file timestamps ... unable to verify current time" — a clock/network
  check on the local build host; no package content is involved.
* WARNING "'qpdf' is needed for checks on size reduction of PDFs" — `qpdf` is not installed on the
  local host. CRAN's builders provide it, so this does not arise there.

The compiled core's shared object is stripped of debug symbols in `src/Makevars` (`strip -S`,
guarded), so the installed package is small and there is no "installed package size" NOTE. The
vignette (`converting-formats.Rmd`, knitr/rmarkdown) builds cleanly.

## Dependencies

All hard dependencies (Matrix, methods, stats, utils) and the LinkingTo (cpp11) are on CRAN. The
format converters use heavy single-cell packages only conditionally (declared in Suggests, all on
CRAN / Bioconductor: SeuratObject, Seurat, SingleCellExperiment, SummarizedExperiment, S4Vectors,
GenomicRanges, igraph, conos, pagoda2, HDF5Array). Package functions, examples and tests degrade
gracefully when they are absent. Two optional integrations target non-mainstream packages that are
therefore intentionally NOT declared in DESCRIPTION and referenced fully dynamically (via
`getExportedValue()`), each degrading to a clear install hint when absent: `view()` forwards to the
separate `pagoda3` viewer package, and `read_seurat_backed()` uses `BPCells` for a disk-backed matrix.

## Test environments

* local: Ubuntu 20.04 Linux, R 4.4.1
* win-builder: Windows, R-release (R 4.6.1) and R-devel — Status: 1 NOTE (the vendored `nlohmann/json`
  diagnostic pragma documented above); "checking CRAN incoming feasibility ... OK".
