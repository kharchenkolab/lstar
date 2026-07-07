# Vendored zstd decoder (for the WASM reader)

Emscripten has no zstd port, so the browser/Node reader (`js/dist/lstar_io.*`) compiles a vendored zstd
**decompress-only** amalgamation to decode Zstd-compressed Zarr v3 stores (zarr-python 3's default v3
compressor). The reader never encodes, and libzarr is built with `-DLIBZARR_ZSTD_DECODE_ONLY`, so only the
decompressor is linked (no `ZSTD_compress*` symbols).

## Files

- `zstddeclib.c` — decompress-only single-file amalgamation.
- `zstd.h`, `zstd_errors.h` — public headers (libzarr's `codecs_zstd.hpp` includes `<zstd.h>`).
- `LICENSE` — zstd's BSD/GPLv2 dual license.

## Provenance / regeneration

Generated from **zstd v1.5.6** with its own tool:

```sh
# from a zstd v1.5.6 source tree:
cd build/single_file_libs && sh create_single_file_decoder.sh   # -> zstddeclib.c
cp zstddeclib.c lib/zstd.h lib/zstd_errors.h LICENSE  <this dir>
```

Wired into `js/build.sh` (the `lstar_io` build): compiled as C to an object, then linked into the C++
reader with `-DLIBZARR_HAS_ZSTD -DLIBZARR_ZSTD_DECODE_ONLY`.
