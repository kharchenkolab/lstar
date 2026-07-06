// SPDX-License-Identifier: MIT

#ifndef LIBZARR_CODECS_ZSTD_HPP
#define LIBZARR_CODECS_ZSTD_HPP

#ifndef LIBZARR_HAS_ZSTD
#error "libzarr/codecs_zstd.hpp requires libzstd: define LIBZARR_HAS_ZSTD and link zstd"
#endif

#include <cstdint>
#include <optional>
#include <string>

#include <zstd.h>

#include "libzarr/detail/common.hpp"
#include "libzarr/types.hpp"

/// \file codecs_zstd.hpp
/// zstd codec (zarr-python 3's default compressor, and numcodecs Zstd in v2
/// stores). Frames written by libzarr carry the content size; reads handle
/// frames without one via streaming decompression.

namespace zarr::detail {

class ZstdCctxGuard {
 public:
  ZstdCctxGuard() : cctx_(ZSTD_createCCtx()) {}
  ZstdCctxGuard(const ZstdCctxGuard&) = delete;
  ZstdCctxGuard& operator=(const ZstdCctxGuard&) = delete;
  ZstdCctxGuard(ZstdCctxGuard&&) = delete;
  ZstdCctxGuard& operator=(ZstdCctxGuard&&) = delete;
  ~ZstdCctxGuard() { ZSTD_freeCCtx(cctx_); }
  [[nodiscard]] ZSTD_CCtx* get() const { return cctx_; }

 private:
  ZSTD_CCtx* cctx_;
};

class ZstdDctxGuard {
 public:
  ZstdDctxGuard() : dctx_(ZSTD_createDCtx()) {}
  ZstdDctxGuard(const ZstdDctxGuard&) = delete;
  ZstdDctxGuard& operator=(const ZstdDctxGuard&) = delete;
  ZstdDctxGuard(ZstdDctxGuard&&) = delete;
  ZstdDctxGuard& operator=(ZstdDctxGuard&&) = delete;
  ~ZstdDctxGuard() { ZSTD_freeDCtx(dctx_); }
  [[nodiscard]] ZSTD_DCtx* get() const { return dctx_; }

 private:
  ZSTD_DCtx* dctx_;
};

/// Compresses `src` at `level` (0 = zstd default). The frame records the
/// content size; `checksum` appends the xxhash frame checksum.
inline Bytes zstd_compress_bytes(const Bytes& src, int level, bool checksum, const char* what) {
  const ZstdCctxGuard cctx;
  if (cctx.get() == nullptr) {
    throw error(std::string(what) + ": ZSTD_createCCtx failed");
  }
  ZSTD_CCtx_setParameter(cctx.get(), ZSTD_c_compressionLevel, level);
  ZSTD_CCtx_setParameter(cctx.get(), ZSTD_c_checksumFlag, checksum ? 1 : 0);
  Bytes out(ZSTD_compressBound(src.size()));
  const std::size_t n = ZSTD_compress2(cctx.get(), out.data(), out.size(), src.data(), src.size());
  if (ZSTD_isError(n) != 0) {
    throw error(std::string(what) + ": zstd compression failed (" + ZSTD_getErrorName(n) + ")");
  }
  out.resize(n);
  return out;
}

/// Decompresses `src`. Frames carrying a content size are validated against
/// `expected` (the decompression-bomb guard); frames without one fall back
/// to streaming decompression.
inline Bytes zstd_decompress_bytes(const Bytes& src, std::optional<std::uint64_t> expected,
                                   const char* what) {
  const unsigned long long content = ZSTD_getFrameContentSize(src.data(), src.size());
  if (content == ZSTD_CONTENTSIZE_ERROR) {
    throw error(std::string(what) + ": corrupt zstd frame");
  }
  if (content != ZSTD_CONTENTSIZE_UNKNOWN) {
    if (expected && content != *expected) {
      throw error(std::string(what) + ": zstd frame decodes to " + std::to_string(content) +
                  " bytes, expected " + std::to_string(*expected));
    }
    Bytes out(checked_size(content, what));
    const std::size_t n = ZSTD_decompress(out.data(), out.size(), src.data(), src.size());
    if (ZSTD_isError(n) != 0 || n != out.size()) {
      throw error(std::string(what) + ": zstd decompression failed");
    }
    return out;
  }

  // No recorded content size (seen from streaming writers): decompress
  // incrementally, growing the output — capped by `expected` when known.
  const ZstdDctxGuard dctx;
  if (dctx.get() == nullptr) {
    throw error(std::string(what) + ": ZSTD_createDCtx failed");
  }
  Bytes out;
  if (expected) {
    out.resize(checked_size(*expected, what));
  } else {
    out.resize(std::max<std::size_t>(src.size() * 3, 64));
  }
  ZSTD_inBuffer in{src.data(), src.size(), 0};
  ZSTD_outBuffer ob{out.data(), out.size(), 0};
  while (true) {
    const std::size_t ret = ZSTD_decompressStream(dctx.get(), &ob, &in);
    if (ZSTD_isError(ret) != 0) {
      throw error(std::string(what) + ": corrupt zstd data (" + ZSTD_getErrorName(ret) + ")");
    }
    if (ret == 0) {
      break;  // frame complete
    }
    if (ob.pos == ob.size) {
      if (expected) {
        throw error(std::string(what) + ": decompressed data exceeds expected " +
                    std::to_string(*expected) + " bytes");
      }
      out.resize(out.size() * 2);
      ob.dst = out.data();
      ob.size = out.size();
    } else if (in.pos == in.size) {
      throw error(std::string(what) + ": zstd data is truncated");
    }
  }
  if (in.pos != in.size) {
    throw error(std::string(what) + ": trailing garbage after zstd data");
  }
  if (expected && ob.pos != *expected) {
    throw error(std::string(what) + ": decompressed to " + std::to_string(ob.pos) +
                " bytes, expected " + std::to_string(*expected));
  }
  out.resize(ob.pos);
  return out;
}

}  // namespace zarr::detail

#endif  // LIBZARR_CODECS_ZSTD_HPP
