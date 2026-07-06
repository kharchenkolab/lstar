#!/usr/bin/env bash
# Zstd read parity on the R surface (R + zarr-python). zarr-python 3's DEFAULT v3 compressor is zstd, so a
# "wild" v3 store is zstd-encoded. The C++ core reads it when built with libzstd (CMake find_library); this
# leg proves the R binding does too, once R is built with libzstd (via R/configure's probe → -lzstd).
#   1. Python writes the maximal store as a ZSTD-compressed v3 store (zarr-python's own default codec) and
#      confirms the store genuinely uses the zstd codec.
#   2. R reads the zstd store and it is value-identical to the (gzip/uncompressed) seed.
# If the installed R package was built WITHOUT libzstd (libzstd-dev absent on the runner), R errors clearly
# on the zstd codec and this leg SKIPS rather than fails — the build is gzip-only-by-design there.
# Self-contained; needs a zarr>=3.1 python and the installed R package.
set -e
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
RLIB="$ROOT/.Rlib"
export PYTHONPATH="$ROOT/python/src"

SEED=/tmp/zstdr_seed.lstar.zarr
ZST=/tmp/zstdr_zstd.lstar.zarr

python3 "$ROOT/conformance/v3_gen.py" "$SEED" | sed 's/^/  [py] /'

# 1. Python re-writes the seed as a zstd-compressed v3 store, and asserts the codec is really zstd.
python3 - "$SEED" "$ZST" <<'PY'
import sys, json, numcodecs, lstar
seed, out = sys.argv[1], sys.argv[2]
lstar.write(lstar.read(seed), out, format="v3", compressor=numcodecs.Zstd(level=5))
codecs = [c.get("name") for c in json.load(open(out + "/fields/counts/data/zarr.json"))["codecs"]]
assert "zstd" in codecs, f"expected zstd codec, got {codecs}"
print(f"  [py] wrote zstd v3 store (counts/data codecs: {codecs})")
PY

# 2. R reads it == the seed (or SKIPS if this R build has no zstd).
Rscript -e '.libPaths(c("'"$RLIB"'", .libPaths()))
suppressMessages(library(lstar))
z <- tryCatch(lstar_read("'"$ZST"'"), error = function(e) {
  if (grepl("zstd", conditionMessage(e), ignore.case = TRUE)) {
    cat("  [R ] SKIP: this R build lacks libzstd (gzip-only); zstd read errors clearly:\n        ",
        conditionMessage(e), "\n"); quit(status = 0)
  }
  stop(e)
})
s <- lstar_read("'"$SEED"'")
stopifnot(identical(sort(names(z$fields)), sort(names(s$fields))),
          identical(sort(names(z$axes)),   sort(names(s$axes))))
zc <- z$fields[["counts"]]; sc <- s$fields[["counts"]]
stopifnot(isTRUE(all.equal(as.numeric(zc$data), as.numeric(sc$data))),
          identical(as.integer(zc$indices), as.integer(sc$indices)))
cat(sprintf("  [R ] read zstd v3 store == seed (%d fields, %d axes; counts data+indices identical)\n",
            length(z$fields), length(z$axes)))' \
  2>&1 | grep -vE "^Warning|deprecat|masked|following object|Attaching|^$"

echo "  zstd read parity on the R surface: R reads zarr-python's default (zstd) v3 store == reference"
