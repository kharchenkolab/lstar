## Resubmission

This is a resubmission. Two incoming-pretest issues have been fixed:

* Removed an orphaned `Additional_repositories: https://bnprks.r-universe.dev` entry from `DESCRIPTION` —
  no declared dependency used it (`BPCells` is referenced fully dynamically, not a `Suggests`), so the
  pretest could not verify it (reported as `? ?`). All `Suggests` are on CRAN or Bioconductor.
* Removed a `#pragma GCC diagnostic ignored "-Wdeprecated-declarations"` that had wrapped the vendored
  `nlohmann/json.hpp` include. It silenced a libc++ (Xcode 26.5+) deprecation of
  `char_traits<unsigned char>` — a warning that arises only under a `-Werror`-style local toolchain and on
  no CRAN builder. The package compiles cleanly without it (as 0.1.0, which shipped no such pragma, did).
  No compiler diagnostics are suppressed in the package sources.

There are no other changes.

## Submission

This is an update to the CRAN package `lstar` (0.1.0 -> 0.2.1). `lstar` provides a uniform data model
(L\*) and a 'Zarr' interchange format for single-cell / spatial omics, with a header-only C++ core
(vendored under `inst/include/lstar`) and bidirectional converters for 'Seurat', 'SingleCellExperiment',
'Conos' and 'pagoda2' objects (the same on-disk store is also readable from Python and C++). Since the
0.1.0 CRAN release the on-disk store defaults to the 'Zarr' v3 layout (v2 remains available; both are read
transparently), with optional 'Zstd' compression and sharded writes; `lstar_read()`/`lstar_write()` and the
converters keep their signatures. See `NEWS.md` for the full history.

The version is 0.2.1 (not 0.1.0). The R package version tracks the companion Python package (`lstar-sc` on
PyPI, 0.2.1) and the shared on-disk format; the R DESCRIPTION had lagged at 0.1.0 while the shared release
line moved ahead. See `NEWS.md`.

## R CMD check results

`R CMD check --no-manual --as-cran` on the local host (Ubuntu 20.04, R 4.4.1) gives:

    Status: 1 WARNING, 2 NOTEs

Every item is either a local-host artifact (absent on CRAN's builders) or an expected, justified
NOTE. There are 0 ERRORs and no code/build WARNINGs.

### NOTE: checking CRAN incoming feasibility

* The local check reports "New maintainer / Old maintainer" only because an older 0.1.0 copy is
  installed in the local check library; this is an update of the existing CRAN package (0.1.0 ->
  0.2.1) and the maintainer is unchanged: Peter Kharchenko <pk.restricted@gmail.com>. This artifact
  does not arise on CRAN's builders, which compare against the published version.

All Suggests are on CRAN or Bioconductor. Two optional integrations target packages that are not on
a mainstream repository — the interactive viewer (`pagoda3`) and a disk-backed Seurat reader
(`BPCells`, used only by `read_seurat_backed()`). Neither is declared in DESCRIPTION; both are resolved
at call time (the package name held in a variable, dispatched via `getExportedValue()`), so the package
installs, checks and runs without them and each gives a clear install hint when invoked and absent.

### NOTE: GNU make is a SystemRequirements

`src/Makevars` uses a GNU make conditional (`ifeq` + `$(shell pkg-config --exists libzstd)`) to detect
`libzstd` at build time: when present, the vendored core is compiled with 'Zstd' decompression enabled (so
it reads 'Zarr' v3 stores written with the Zstd codec); when absent it degrades to gzip-only with a clear
run-time message. `GNU make` is therefore declared in `SystemRequirements`, raising this expected NOTE.
GNU make is the default on all CRAN build platforms, and nothing depends on a specific make beyond that
single conditional.

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

* local: Ubuntu 20.04 Linux, R 4.4.1 — `R CMD check --no-manual --as-cran`, Status as above.
* GitHub Actions (r-lib/actions, `--as-cran` with error-on = warning): ubuntu-latest, windows-latest and
  macOS-latest (Xcode 26.5), each on R-release AND R-devel — all six pass with 0 ERRORs and 0 WARNINGs.
  This is the primary cross-platform evidence for this submission.
* win-builder (R-release + R-devel) was used during preparation; the only items it surfaced were the
  incoming-pretest issues since removed (the orphaned `Additional_repositories` and the diagnostic pragma).
