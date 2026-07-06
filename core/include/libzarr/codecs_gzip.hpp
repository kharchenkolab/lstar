// SPDX-License-Identifier: MIT

#ifndef LIBZARR_CODECS_GZIP_HPP
#define LIBZARR_CODECS_GZIP_HPP

#ifndef LIBZARR_HAS_ZLIB
#error "libzarr/codecs_gzip.hpp requires zlib: define LIBZARR_HAS_ZLIB and link against zlib"
#endif

#include <algorithm>
#include <climits>
#include <cstdint>
#include <optional>
#include <string>

#include <zlib.h>

#include "libzarr/detail/common.hpp"
#include "libzarr/types.hpp"

/// \file codecs_gzip.hpp
/// zlib-backed byte codecs. Two distinct framings share one algorithm:
/// v2 "zlib" is RFC 1950, v2/v3 "gzip" is RFC 1952 — the difference is only
/// the header/trailer, selected via zlib's windowBits.

namespace zarr::detail {

class ZlibDeflateGuard {
 public:
  explicit ZlibDeflateGuard(z_stream* zs) : zs_(zs) {}
  ZlibDeflateGuard(const ZlibDeflateGuard&) = delete;
  ZlibDeflateGuard& operator=(const ZlibDeflateGuard&) = delete;
  ZlibDeflateGuard(ZlibDeflateGuard&&) = delete;
  ZlibDeflateGuard& operator=(ZlibDeflateGuard&&) = delete;
  ~ZlibDeflateGuard() { deflateEnd(zs_); }

 private:
  z_stream* zs_;
};

class ZlibInflateGuard {
 public:
  explicit ZlibInflateGuard(z_stream* zs) : zs_(zs) {}
  ZlibInflateGuard(const ZlibInflateGuard&) = delete;
  ZlibInflateGuard& operator=(const ZlibInflateGuard&) = delete;
  ZlibInflateGuard(ZlibInflateGuard&&) = delete;
  ZlibInflateGuard& operator=(ZlibInflateGuard&&) = delete;
  ~ZlibInflateGuard() { inflateEnd(zs_); }

 private:
  z_stream* zs_;
};

/// Compresses `src` at `level` (0-9). gzip_framing selects RFC 1952 vs 1950.
inline Bytes deflate_bytes(const Bytes& src, int level, bool gzip_framing, const char* what) {
  z_stream zs{};
  const int window_bits = gzip_framing ? 15 + 16 : 15;
  if (deflateInit2(&zs, level, Z_DEFLATED, window_bits, 8, Z_DEFAULT_STRATEGY) != Z_OK) {
    throw error(std::string(what) + ": deflateInit2 failed");
  }
  const ZlibDeflateGuard guard(&zs);
  Bytes out(deflateBound(&zs, static_cast<uLong>(src.size())));
  std::size_t in_pos = 0;
  std::size_t out_pos = 0;
  int ret = Z_OK;
  do {
    const std::size_t in_step = std::min<std::size_t>(src.size() - in_pos, UINT_MAX);
    // zlib's C API takes a non-const next_in but never writes through it.
    // NOLINTNEXTLINE(cppcoreguidelines-pro-type-const-cast)
    zs.next_in = const_cast<Bytef*>(src.data() + in_pos);
    zs.avail_in = static_cast<uInt>(in_step);
    zs.next_out = out.data() + out_pos;
    zs.avail_out = static_cast<uInt>(std::min<std::size_t>(out.size() - out_pos, UINT_MAX));
    const bool last_input = in_pos + in_step == src.size();
    ret = deflate(&zs, last_input ? Z_FINISH : Z_NO_FLUSH);
    if (ret != Z_OK && ret != Z_STREAM_END && ret != Z_BUF_ERROR) {
      throw error(std::string(what) + ": deflate failed (" + std::to_string(ret) + ")");
    }
    in_pos += in_step - zs.avail_in;
    out_pos = static_cast<std::size_t>(zs.next_out - out.data());
    if (out_pos == out.size() && ret != Z_STREAM_END) {
      out.resize(out.size() + out.size() / 2 + 64);
    }
  } while (ret != Z_STREAM_END);
  out.resize(out_pos);
  return out;
}

/// Decompresses `src`. When `expected_size` is known (a decoded chunk's byte
/// count) the output is allocated exactly and any mismatch is an error —
/// this doubles as the decompression-bomb guard. Reads auto-detect zlib vs
/// gzip framing (windowBits 15+32): the two ids get confused in the wild,
/// and accepting both on read is harmless.
inline Bytes inflate_bytes(const Bytes& src, std::optional<std::uint64_t> expected_size,
                           const char* what) {
  z_stream zs{};
  if (inflateInit2(&zs, 15 + 32) != Z_OK) {
    throw error(std::string(what) + ": inflateInit2 failed");
  }
  const ZlibInflateGuard guard(&zs);
  Bytes out;
  if (expected_size) {
    out.resize(checked_size(*expected_size, what));
  } else {
    out.resize(std::max<std::size_t>(src.size() * 3, 64));
  }
  std::size_t in_pos = 0;
  std::size_t out_pos = 0;
  int ret = Z_OK;
  while (true) {
    const std::size_t in_step = std::min<std::size_t>(src.size() - in_pos, UINT_MAX);
    // zlib's C API takes a non-const next_in but never writes through it.
    // NOLINTNEXTLINE(cppcoreguidelines-pro-type-const-cast)
    zs.next_in = const_cast<Bytef*>(src.data() + in_pos);
    zs.avail_in = static_cast<uInt>(in_step);
    zs.next_out = out.data() + out_pos;
    zs.avail_out = static_cast<uInt>(std::min<std::size_t>(out.size() - out_pos, UINT_MAX));
    ret = inflate(&zs, Z_NO_FLUSH);
    if (ret != Z_OK && ret != Z_STREAM_END && ret != Z_BUF_ERROR) {
      throw error(std::string(what) + ": corrupt compressed data (zlib error " +
                  std::to_string(ret) + ")");
    }
    in_pos += in_step - zs.avail_in;
    out_pos = static_cast<std::size_t>(zs.next_out - out.data());
    if (ret == Z_STREAM_END) {
      break;
    }
    if (out_pos == out.size()) {
      if (expected_size) {
        throw error(std::string(what) + ": decompressed data exceeds expected " +
                    std::to_string(*expected_size) + " bytes");
      }
      out.resize(out.size() * 2);
    } else if (in_pos == src.size()) {
      throw error(std::string(what) + ": compressed data is truncated");
    }
  }
  if (in_pos != src.size()) {
    throw error(std::string(what) + ": trailing garbage after compressed data");
  }
  if (expected_size && out_pos != *expected_size) {
    throw error(std::string(what) + ": decompressed to " + std::to_string(out_pos) +
                " bytes, expected " + std::to_string(*expected_size));
  }
  out.resize(out_pos);
  return out;
}

}  // namespace zarr::detail

#endif  // LIBZARR_CODECS_GZIP_HPP
