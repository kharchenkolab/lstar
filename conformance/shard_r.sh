#!/usr/bin/env bash
# Sharded v3 WRITES on the R surface. R's lstar_write(shard_elems=) threads through the cpp11 glue to the
# shared C++ core (libzarr ArraySpec.shards), so R can produce a sharded v3 store -- the hosting
# optimization (fewer objects) the pagoda3 viewer reads. On the maximal store:
#   1. R writes it UNSHARDED v3 and SHARDED v3 (chunked); the two are value-identical (sharding changes
#      the object layout, never the data) -- compared R-to-R so R's own nullable representation is
#      consistent (a Python-written seed differs at masked positions, which is expected + benign).
#   2. zarr-python confirms the sharded store genuinely uses the sharding_indexed codec.
#   3. lstar_write rejects shard_elems without format="v3" / without chunk_elems (clear errors).
# Needs a zarr>=3.1 python, the installed R package, and the cmake C++ core (test_v3 comparator).
set -e
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
RLIB="$ROOT/.Rlib"; BIN="$ROOT/core/build/test_v3"
export PYTHONPATH="$ROOT/python/src"
SEED=/tmp/shr_seed.lstar.zarr; RUN=/tmp/shr_r_unsharded.lstar.zarr; RSH=/tmp/shr_r_sharded.lstar.zarr

if [ ! -x "$BIN" ]; then
  cmake -S "$ROOT/core" -B "$ROOT/core/build" -DCMAKE_BUILD_TYPE=Release >/tmp/shr_cmake.log 2>&1
  cmake --build "$ROOT/core/build" --target test_v3 -j4 >>/tmp/shr_cmake.log 2>&1
fi

python3 "$ROOT/conformance/v3_gen.py" "$SEED" | sed 's/^/  [py] /'
rm -rf "$RUN" "$RSH"
Rscript -e '.libPaths(c("'"$RLIB"'", .libPaths())); suppressMessages(library(lstar))
ds <- lstar_read("'"$SEED"'")
lstar_write(ds, "'"$RUN"'", chunk_elems = 1000, format = "v3")                     # 1. R unsharded v3
lstar_write(ds, "'"$RSH"'", chunk_elems = 1000, format = "v3", shard_elems = 4000) # 1. R sharded v3
# 3. validation: shard_elems is v3-only and needs chunk_elems
r1 <- tryCatch({lstar_write(ds, tempfile(), chunk_elems = 1000, format = "v2", shard_elems = 4000); ""},
               error = function(e) conditionMessage(e))
r2 <- tryCatch({lstar_write(ds, tempfile(), format = "v3", shard_elems = 4000); ""},
               error = function(e) conditionMessage(e))
stopifnot(grepl("v3", r1), grepl("chunk_elems", r2))
cat("  [R ] wrote unsharded + sharded v3; rejects shard+v2 and shard-without-chunk\n")' \
  2>&1 | grep -vE "^Warning|deprecat|masked|following object|Attaching|^$"

"$BIN" compare "$RUN" "$RSH"                                       # 1. R sharded == R unsharded (value-identical)
python3 "$ROOT/conformance/v3_verify.py" shardcheck "$RSH" "$RUN"  # 2. sharding_indexed + zarr-python == R unsharded
echo "  sharded v3 writes on the R surface: R produces a conformant sharded store == its unsharded output"
