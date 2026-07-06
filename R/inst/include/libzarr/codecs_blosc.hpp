// SPDX-License-Identifier: MIT

#ifndef LIBZARR_CODECS_BLOSC_HPP
#define LIBZARR_CODECS_BLOSC_HPP

#ifndef LIBZARR_HAS_BLOSC
#error "libzarr/codecs_blosc.hpp requires c-blosc: define LIBZARR_HAS_BLOSC and link blosc"
#endif

#include <cstddef>
#include <cstdint>
#include <optional>
#include <string>

#include <blosc.h>

#include "libzarr/detail/common.hpp"
#include "libzarr/types.hpp"

/// \file codecs_blosc.hpp
/// c-blosc-backed codec. Blosc frames are self-describing, so decode needs no
/// parameters; the configuration only drives encoding.

namespace zarr::detail {

struct BloscParams {
  std::string cname = "lz4";
  int clevel = 5;
  int shuffle = BLOSC_SHUFFLE;
  std::size_t typesize = 1;
  std::size_t blocksize = 0;  // 0 = automatic
};

inline Bytes blosc_decompress_bytes(const Bytes& src, std::optional<std::uint64_t> expected,
                                    const char* what) {
  std::size_t nbytes = 0;
  if (blosc_cbuffer_validate(src.data(), src.size(), &nbytes) < 0) {
    throw error(std::string(what) + ": corrupt blosc frame");
  }
  if (expected && nbytes != *expected) {
    throw error(std::string(what) + ": blosc frame decodes to " + std::to_string(nbytes) +
                " bytes, expected " + std::to_string(*expected));
  }
  Bytes out(nbytes);
  const int n = blosc_decompress_ctx(src.data(), out.data(), out.size(), /*numinternalthreads=*/1);
  if (n < 0 || static_cast<std::size_t>(n) != nbytes) {
    throw error(std::string(what) + ": blosc decompression failed (" + std::to_string(n) + ")");
  }
  return out;
}

inline Bytes blosc_compress_bytes(const Bytes& src, const BloscParams& params, const char* what) {
  Bytes out(src.size() + BLOSC_MAX_OVERHEAD);
  const int n = blosc_compress_ctx(params.clevel, params.shuffle, params.typesize, src.size(),
                                   src.data(), out.data(), out.size(), params.cname.c_str(),
                                   params.blocksize, /*numinternalthreads=*/1);
  if (n <= 0) {
    throw error(std::string(what) + ": blosc compression failed (" + std::to_string(n) + ")");
  }
  out.resize(static_cast<std::size_t>(n));
  return out;
}

}  // namespace zarr::detail

#endif  // LIBZARR_CODECS_BLOSC_HPP
