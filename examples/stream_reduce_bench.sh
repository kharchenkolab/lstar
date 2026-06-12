#!/usr/bin/env bash
# Benchmark the blocked (bounded-memory) per-gene reduction against a full read, on a real-sized
# gene-major (CSC) measure: peak resident memory (via /usr/bin/time -v, which captures the C++
# allocations that live outside R's heap) and OpenMP thread scaling. The point of the C++/R blocked
# reader is that a per-gene mean/variance pass (HVG selection, variance modeling) runs over an atlas
# too large to hold -- so we show it computes the same statistics in a fraction of the memory.
#
# Usage: bash examples/stream_reduce_bench.sh [path/to/file.h5ad]
#   Builds a chunked CSC L* store from the given .h5ad (or a local default, or synthetic), then
#   benchmarks. Needs the lstar R package installed in ./.Rlib (conformance/run.sh does that).
set -e
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
RLIB="$ROOT/.Rlib"
STORE=/tmp/stream_reduce_bench.lstar.zarr
DEFAULT="/home/pkharchenko/cacoa/age/tab.muris/tabula-muris-senis-droplet-processed-official-annotations-Marrow.h5ad"
H5="${1:-$DEFAULT}"

# (1) Build a chunked, gene-major (CSC) store. Stream the .h5ad in (bounded), then recompress the
#     measure as CSC so a per-gene reduction reads contiguous gene columns.
PYTHONPATH="$ROOT/python/src" python3 - "$H5" "$STORE" <<'PY'
import os, sys, warnings; warnings.filterwarnings("ignore")
import numpy as np, scipy.sparse as sp
import lstar
h5, store = sys.argv[1], sys.argv[2]
if os.path.exists(h5):
    import anndata as ad
    X = sp.csc_matrix(ad.read_h5ad(h5).X)               # cells x genes, gene-major
    src = os.path.basename(h5)
else:
    X = sp.random(200_000, 20_000, density=0.05, format="csc", random_state=0)
    X.data = np.round(X.data * 9 + 1).astype("f4"); src = "synthetic"
ds = lstar.Dataset(kind="sample")
ds.add_axis("cells", [f"c{i}" for i in range(X.shape[0])])
ds.add_axis("genes", [f"g{i}" for i in range(X.shape[1])])
ds.add_field("counts", X.astype("f4"), role="measure", span=["cells", "genes"], state="raw")
lstar.write(ds, store, chunk_elems=2_000_000)           # chunked -> blocked reader streams
print(f"source: {src}  ({X.shape[0]}x{X.shape[1]}, {X.nnz:,} nonzeros, CSC, chunked)")
PY

run_rss() {   # $1=label  $2=R expression -> prints wall time + peak RSS
  local out; out=$(/usr/bin/time -v Rscript -e ".libPaths(c('$RLIB', .libPaths())); suppressMessages(library(lstar)); t<-proc.time(); $2; cat(sprintf('  wall %.2fs\n',(proc.time()-t)[3]))" 2>&1)
  local rss; rss=$(echo "$out" | awk -F': ' '/Maximum resident/{print $2}')
  echo "$out" | grep -E '^  '
  printf "  peak RSS: %.0f MB\n" "$(echo "$rss/1024" | bc -l)"
}

echo
echo "== memory: full read+reduce  vs  blocked reduce =="
echo "-- full read (whole measure into R, then reduce) --"
run_rss full "m<-as(field_value(lstar_read('$STORE'),'counts'),'CsparseMatrix'); nr<-nrow(m); X<-m; X@x<-log1p(X@x); cs<-Matrix::colSums(X); v<-(Matrix::colSums(X*X)-cs*cs/nr)/(nr-1); cat(sprintf('  genes=%d mean(var)=%.4f\n', length(v), mean(v)))"
echo "-- blocked reduce (bounded: one block at a time) --"
run_rss blocked "s<-stream_col_stats('$STORE','counts', block=2048L, n_threads=1L, lognorm=TRUE); cat(sprintf('  genes=%d mean(var)=%.4f\n', length(s\$var), mean(s\$var)))"

echo
echo "== thread scaling of the blocked reduce (log1p, block=2048) =="
Rscript -e ".libPaths(c('$RLIB', .libPaths())); suppressMessages(library(lstar))
for (nt in c(1L,2L,4L,8L)) {
  t <- proc.time()
  s <- stream_col_stats('$STORE','counts', block=2048L, n_threads=nt, lognorm=TRUE)
  el <- (proc.time()-t)[3]
  cat(sprintf('  %2d threads: %5.2fs\n', nt, el))
}" 2>&1 | grep -E '^  '
echo "  (note: this reduction is dominated by reading+decoding chunks, not arithmetic, so"
echo "   intra-block OpenMP threads barely move the wall clock -- the win here is memory, not"
echo "   cores. Threading pays off for the compute-bound kernels, e.g. csc_col_sum_by_group.)"
