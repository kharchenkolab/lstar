# Vendored zstd (for the WASM reader + writer)

Emscripten has no zstd port, so the WASM modules vendor zstd from source. Two amalgamations, because the
reader and writer have opposite needs:

- The **reader** (`js/dist/lstar_io.*`) only ever *decodes* Zstd-compressed Zarr v3 stores (zarr-python 3's
  default v3 compressor). libzarr is built with `-DLIBZARR_ZSTD_DECODE_ONLY`, so only the small
  decompress-only amalgamation is linked (no `ZSTD_compress*` symbols).
- The **writer** (`js/dist/lstar_writer.*`) *encodes* chunks with zstd (full JS-writer parity), so it links
  the FULL amalgamation with `-DLIBZARR_HAS_ZSTD` (no `DECODE_ONLY`). Loaded only when writing a store.

## Files

- `zstddeclib.c` — decompress-only single-file amalgamation (reader).
- `zstd.c` — full single-file amalgamation, encode + decode (writer).
- `zstd.h`, `zstd_errors.h` — public headers (libzarr's `codecs_zstd.hpp` includes `<zstd.h>`).
- `LICENSE` — zstd's BSD/GPLv2 dual license.

## Provenance / regeneration

Generated from **zstd v1.5.6** with its own tools:

```sh
# from a zstd v1.5.6 source tree:
cd build/single_file_libs
sh create_single_file_decoder.sh                                 # -> zstddeclib.c (reader)
sh create_single_file_library.sh                                 # -> zstd.c       (writer)
cp zstddeclib.c zstd.c lib/zstd.h lib/zstd_errors.h ../../LICENSE  <this dir>
```

Wired into `js/build.sh`: `zstddeclib.c` -> the `lstar_io` reader build (`-DLIBZARR_ZSTD_DECODE_ONLY`);
`zstd.c` -> the `lstar_writer` build (`-DLIBZARR_HAS_ZSTD`, full encode). Each is compiled as C to an object
first, then linked into its C++ module (a `.c` source can't share the `.cpp`'s `-std=c++17`).
