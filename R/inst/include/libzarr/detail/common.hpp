// SPDX-License-Identifier: MIT

#ifndef LIBZARR_DETAIL_COMMON_HPP
#define LIBZARR_DETAIL_COMMON_HPP

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <string>
#include <string_view>
#include <vector>

#include "libzarr/types.hpp"

/// \file common.hpp
/// Internal helpers: endianness, checked arithmetic, base64, and N-dimensional
/// box copies. Everything here is implementation detail, not public API.

namespace zarr::detail {

inline bool host_is_little_endian() {
  const std::uint16_t probe = 1;
  std::uint8_t first = 0;
  std::memcpy(&first, &probe, 1);
  return first == 1;
}

/// Reverses the bytes of each `width`-byte element in place. For complex
/// types callers pass the component width (a complex is two floats, each
/// independently byte-ordered).
inline void byteswap_inplace(std::uint8_t* data, std::uint64_t count, std::uint32_t width) {
  for (std::uint64_t i = 0; i < count; ++i) {
    std::reverse(data + i * width, data + (i + 1) * width);
  }
}

inline std::uint64_t ceil_div(std::uint64_t a, std::uint64_t b) {
  return a / b + static_cast<std::uint64_t>(a % b != 0);
}

/// Product of `dims`, throwing zarr::error on uint64 overflow. `what` names
/// the quantity in the error message.
inline std::uint64_t checked_product(const std::vector<std::uint64_t>& dims, const char* what) {
  std::uint64_t product = 1;
  for (const std::uint64_t d : dims) {
    if (d != 0 && product > std::numeric_limits<std::uint64_t>::max() / d) {
      throw error(std::string(what) + ": size overflows uint64");
    }
    product *= d;
  }
  return product;
}

/// Narrow a uint64 byte count to std::size_t (which is 32-bit on wasm32),
/// throwing if the platform cannot address it.
inline std::size_t checked_size(std::uint64_t bytes, const char* what) {
  if (bytes > std::numeric_limits<std::size_t>::max()) {
    throw error(std::string(what) + ": " + std::to_string(bytes) +
                " bytes exceeds this platform's addressable size");
  }
  return static_cast<std::size_t>(bytes);
}

inline bool starts_with(std::string_view text, std::string_view prefix) {
  return text.size() >= prefix.size() && text.compare(0, prefix.size(), prefix) == 0;
}

inline bool ends_with(std::string_view text, std::string_view suffix) {
  return text.size() >= suffix.size() &&
         text.compare(text.size() - suffix.size(), suffix.size(), suffix) == 0;
}

/// CRC-32C (Castagnoli, reflected poly 0x82F63B78, per RFC 3720 §B.4) — the
/// checksum required by the v3 crc32c codec and the sharding index. This is
/// NOT the zip/zlib CRC-32 (IEEE), which uses a different polynomial.
///
/// Portable byte-at-a-time table implementation; the runtime-dispatched
/// public crc32c() below prefers the SSE4.2 instruction where available.
inline std::uint32_t crc32c_table(const std::uint8_t* data, std::size_t size) {
  static const std::array<std::uint32_t, 256> table = [] {
    std::array<std::uint32_t, 256> t{};
    for (std::uint32_t i = 0; i < 256; ++i) {
      std::uint32_t c = i;
      for (int k = 0; k < 8; ++k) {
        c = ((c & 1U) != 0) ? 0x82F63B78U ^ (c >> 1U) : c >> 1U;
      }
      t[i] = c;
    }
    return t;
  }();
  std::uint32_t c = 0xFFFFFFFFU;
  for (std::size_t i = 0; i < size; ++i) {
    c = table[(c ^ data[i]) & 0xFFU] ^ (c >> 8U);
  }
  return c ^ 0xFFFFFFFFU;
}

// The hardware path uses the SSE4.2 CRC32 instruction (`crc32` — Castagnoli,
// the same polynomial), ~8x the table method. It exists only on x86 under
// GCC/Clang, where target-attributed code can emit the instruction without
// building the whole TU for SSE4.2 and __builtin_cpu_supports picks it at
// run time. Everywhere else (ARM, wasm, MSVC) the table path is used, so the
// core stays portable and WASM-clean.
#if (defined(__x86_64__) || defined(__i386__)) && (defined(__GNUC__) || defined(__clang__))
#define LIBZARR_CRC32C_X86 1
#include <nmmintrin.h>

__attribute__((target("sse4.2"))) inline std::uint32_t crc32c_sse42(const std::uint8_t* data,
                                                                    std::size_t size) {
  std::uint64_t crc = 0xFFFFFFFFU;
  std::size_t i = 0;
  for (; i + 8 <= size; i += 8) {
    std::uint64_t word = 0;
    std::memcpy(&word, data + i, 8);  // x86 is little-endian; matches byte order
    crc = _mm_crc32_u64(crc, word);
  }
  auto c = static_cast<std::uint32_t>(crc);
  for (; i < size; ++i) {
    c = _mm_crc32_u8(c, data[i]);
  }
  return c ^ 0xFFFFFFFFU;
}
#endif

/// True when crc32c() uses the SSE4.2 instruction on this CPU.
inline bool crc32c_uses_hardware() {
#ifdef LIBZARR_CRC32C_X86
  static const bool hw = __builtin_cpu_supports("sse4.2");
  return hw;
#else
  return false;
#endif
}

/// CRC-32C over `size` bytes at `data`. Dispatches to the SSE4.2 instruction
/// when the CPU supports it, else the portable table.
inline std::uint32_t crc32c(const std::uint8_t* data, std::size_t size) {
#ifdef LIBZARR_CRC32C_X86
  if (crc32c_uses_hardware()) {
    return crc32c_sse42(data, size);
  }
#endif
  return crc32c_table(data, size);
}

// ---- IEEE 754 binary16 (float16 has no native C++17 type) ------------------

/// double -> binary16 bits, round-to-nearest-even; overflows go to infinity.
inline std::uint16_t double_to_half_bits(double value) {
  const auto f = static_cast<float>(value);
  std::uint32_t bits = 0;
  std::memcpy(&bits, &f, 4);
  const std::uint32_t sign = (bits >> 16U) & 0x8000U;
  const std::uint32_t mantissa = bits & 0x007FFFFFU;
  const std::int32_t exponent = static_cast<std::int32_t>((bits >> 23U) & 0xFFU) - 127 + 15;
  if (exponent >= 0x1F) {  // overflow or inf/nan
    if (((bits >> 23U) & 0xFFU) == 0xFF && mantissa != 0) {
      return static_cast<std::uint16_t>(sign | 0x7E00U);  // pinned quiet NaN
    }
    return static_cast<std::uint16_t>(sign | 0x7C00U);  // infinity
  }
  if (exponent <= 0) {  // subnormal or zero
    if (exponent < -10) {
      return static_cast<std::uint16_t>(sign);
    }
    const std::uint32_t sub = (mantissa | 0x00800000U) >> static_cast<std::uint32_t>(14 - exponent);
    const std::uint32_t rounded =
        sub + ((sub & 0x1FFFU) > (0x1000U - ((sub >> 13U) & 1U)) ? 0x2000U : 0U);
    return static_cast<std::uint16_t>(sign | (rounded >> 13U));
  }
  std::uint32_t half = sign | (static_cast<std::uint32_t>(exponent) << 10U) | (mantissa >> 13U);
  // round to nearest even on the truncated 13 bits
  const std::uint32_t rest = mantissa & 0x1FFFU;
  if (rest > 0x1000U || (rest == 0x1000U && (half & 1U) != 0)) {
    ++half;  // may carry into the exponent, which correctly yields infinity
  }
  return static_cast<std::uint16_t>(half);
}

/// binary16 bits -> double (exact).
inline double half_bits_to_double(std::uint16_t half) {
  const std::uint32_t sign = (static_cast<std::uint32_t>(half) & 0x8000U) << 16U;
  const std::uint32_t exponent = (static_cast<std::uint32_t>(half) >> 10U) & 0x1FU;
  const std::uint32_t mantissa = static_cast<std::uint32_t>(half) & 0x3FFU;
  std::uint32_t bits = 0;
  if (exponent == 0x1F) {
    bits = sign | 0x7F800000U | (mantissa << 13U);
  } else if (exponent != 0) {
    bits = sign | ((exponent - 15 + 127) << 23U) | (mantissa << 13U);
  } else if (mantissa != 0) {  // subnormal: normalize
    std::uint32_t m = mantissa;
    std::int32_t e = -1;
    while ((m & 0x400U) == 0) {
      m <<= 1U;
      ++e;
    }
    bits = sign | (static_cast<std::uint32_t>(112 - e) << 23U) | ((m & 0x3FFU) << 13U);
  } else {
    bits = sign;  // zero
  }
  float f = 0;
  std::memcpy(&f, &bits, 4);
  return static_cast<double>(f);
}

/// numcodecs-style byte shuffle: element bytes are regrouped by byte
/// position (out[j*count + i] = in[i*es + j]); trailing bytes that do not
/// fill an element are copied unchanged. Its own inverse is unshuffle_bytes.
inline Bytes shuffle_bytes(const Bytes& src, std::uint32_t elementsize) {
  Bytes out(src.size());
  const std::size_t count = src.size() / elementsize;
  for (std::size_t i = 0; i < count; ++i) {
    for (std::size_t j = 0; j < elementsize; ++j) {
      out[j * count + i] = src[i * elementsize + j];
    }
  }
  std::memcpy(out.data() + count * elementsize, src.data() + count * elementsize,
              src.size() - count * elementsize);
  return out;
}

/// Inverse of shuffle_bytes.
inline Bytes unshuffle_bytes(const Bytes& src, std::uint32_t elementsize) {
  Bytes out(src.size());
  const std::size_t count = src.size() / elementsize;
  for (std::size_t i = 0; i < count; ++i) {
    for (std::size_t j = 0; j < elementsize; ++j) {
      out[i * elementsize + j] = src[j * count + i];
    }
  }
  std::memcpy(out.data() + count * elementsize, src.data() + count * elementsize,
              src.size() - count * elementsize);
  return out;
}

// ---- base64 (RFC 4648, with padding) — v2 fill_value for raw dtypes -------

inline std::string base64_encode(const std::uint8_t* data, std::size_t size) {
  constexpr std::string_view kAlphabet =
      "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
  std::string out;
  out.reserve(((size + 2) / 3) * 4);
  std::size_t i = 0;
  for (; i + 3 <= size; i += 3) {
    const std::uint32_t v = (static_cast<std::uint32_t>(data[i]) << 16U) |
                            (static_cast<std::uint32_t>(data[i + 1]) << 8U) |
                            static_cast<std::uint32_t>(data[i + 2]);
    out.push_back(kAlphabet[(v >> 18U) & 0x3FU]);
    out.push_back(kAlphabet[(v >> 12U) & 0x3FU]);
    out.push_back(kAlphabet[(v >> 6U) & 0x3FU]);
    out.push_back(kAlphabet[v & 0x3FU]);
  }
  if (i < size) {
    const bool two = (size - i) == 2;
    const std::uint32_t v = (static_cast<std::uint32_t>(data[i]) << 16U) |
                            (two ? static_cast<std::uint32_t>(data[i + 1]) << 8U : 0U);
    out.push_back(kAlphabet[(v >> 18U) & 0x3FU]);
    out.push_back(kAlphabet[(v >> 12U) & 0x3FU]);
    out.push_back(two ? kAlphabet[(v >> 6U) & 0x3FU] : '=');
    out.push_back('=');
  }
  return out;
}

inline std::uint32_t base64_value(char c, const char* what) {
  if (c >= 'A' && c <= 'Z') {
    return static_cast<std::uint32_t>(c - 'A');
  }
  if (c >= 'a' && c <= 'z') {
    return static_cast<std::uint32_t>(c - 'a') + 26U;
  }
  if (c >= '0' && c <= '9') {
    return static_cast<std::uint32_t>(c - '0') + 52U;
  }
  if (c == '+') {
    return 62U;
  }
  if (c == '/') {
    return 63U;
  }
  throw error(std::string(what) + ": invalid base64 character '" + std::string(1, c) + "'");
}

inline Bytes base64_decode(std::string_view text, const char* what) {
  const auto value_of = [what](char c) { return base64_value(c, what); };
  if (text.size() % 4 != 0) {
    throw error(std::string(what) + ": base64 length must be a multiple of 4");
  }
  Bytes out;
  out.reserve(text.size() / 4 * 3);
  for (std::size_t i = 0; i < text.size(); i += 4) {
    const bool last = i + 4 == text.size();
    std::size_t pad = 0;
    if (last && text[i + 3] == '=') {
      pad = text[i + 2] == '=' ? 2 : 1;
    } else if (last && text[i + 2] == '=') {
      throw error(std::string(what) + ": invalid base64 padding");
    }
    std::uint32_t v = (value_of(text[i]) << 18U) | (value_of(text[i + 1]) << 12U);
    if (pad < 2) {
      v |= value_of(text[i + 2]) << 6U;
    }
    if (pad < 1) {
      v |= value_of(text[i + 3]);
    }
    out.push_back(static_cast<std::uint8_t>((v >> 16U) & 0xFFU));
    if (pad < 2) {
      out.push_back(static_cast<std::uint8_t>((v >> 8U) & 0xFFU));
    }
    if (pad < 1) {
      out.push_back(static_cast<std::uint8_t>(v & 0xFFU));
    }
  }
  return out;
}

// ---- N-dimensional box math ------------------------------------------------

/// C-order strides in *bytes* for `shape` with `itemsize`-byte elements.
inline std::vector<std::uint64_t> c_strides_bytes(const std::vector<std::uint64_t>& shape,
                                                  std::uint32_t itemsize) {
  std::vector<std::uint64_t> strides(shape.size());
  std::uint64_t acc = itemsize;
  for (std::size_t d = shape.size(); d-- > 0;) {
    strides[d] = acc;
    acc *= shape[d];
  }
  return strides;
}

/// Copies box `box` (element counts per dimension) from position `src_origin`
/// of a C-order buffer of shape `src_shape` to position `dst_origin` of a
/// C-order buffer of shape `dst_shape`. Rank 0 copies a single element.
/// Bounds are the caller's responsibility (internal invariant).
inline void copy_box(const std::uint8_t* src, const std::vector<std::uint64_t>& src_shape,
                     const std::vector<std::uint64_t>& src_origin, std::uint8_t* dst,
                     const std::vector<std::uint64_t>& dst_shape,
                     const std::vector<std::uint64_t>& dst_origin,
                     const std::vector<std::uint64_t>& box, std::uint32_t itemsize) {
  const std::size_t rank = box.size();
  if (rank == 0) {
    std::memcpy(dst, src, itemsize);
    return;
  }
  for (const std::uint64_t extent : box) {
    if (extent == 0) {
      return;
    }
  }
  const std::vector<std::uint64_t> src_strides = c_strides_bytes(src_shape, itemsize);
  const std::vector<std::uint64_t> dst_strides = c_strides_bytes(dst_shape, itemsize);
  const std::size_t row_bytes = checked_size(box[rank - 1] * itemsize, "copy_box row");

  std::vector<std::uint64_t> index(rank, 0);  // index over `box`, last dim always 0
  while (true) {
    std::uint64_t src_off = 0;
    std::uint64_t dst_off = 0;
    for (std::size_t d = 0; d < rank; ++d) {
      src_off += (src_origin[d] + index[d]) * src_strides[d];
      dst_off += (dst_origin[d] + index[d]) * dst_strides[d];
    }
    std::memcpy(dst + dst_off, src + src_off, row_bytes);
    // odometer over dimensions [0, rank-1)
    std::size_t d = rank - 1;
    while (d-- > 0) {
      if (++index[d] < box[d]) {
        break;
      }
      index[d] = 0;
    }
    if (d == static_cast<std::size_t>(-1)) {
      return;
    }
  }
}

/// Gathers a C-order buffer of `shape` from `src`, reading the element for
/// C-index (i0, i1, ...) at byte offset sum(i_d * src_strides_bytes[d]).
/// Used to decode transposed (e.g. v2 order:"F") chunks.
inline void gather_strided(const std::uint8_t* src,
                           const std::vector<std::uint64_t>& src_strides_bytes, std::uint8_t* dst,
                           const std::vector<std::uint64_t>& shape, std::uint32_t itemsize) {
  const std::size_t rank = shape.size();
  if (rank == 0) {
    std::memcpy(dst, src, itemsize);
    return;
  }
  const std::uint64_t total = checked_product(shape, "gather_strided");
  if (total == 0) {
    return;
  }
  std::vector<std::uint64_t> index(rank, 0);
  std::uint8_t* out = dst;
  for (std::uint64_t n = 0; n < total; ++n) {
    std::uint64_t src_off = 0;
    for (std::size_t d = 0; d < rank; ++d) {
      src_off += index[d] * src_strides_bytes[d];
    }
    std::memcpy(out, src + src_off, itemsize);
    out += itemsize;
    std::size_t d = rank;
    while (d-- > 0) {
      if (++index[d] < shape[d]) {
        break;
      }
      index[d] = 0;
    }
  }
}

/// Fills `count` elements at `dst` with the `itemsize`-byte pattern `elem`,
/// or with zeros when `elem` is null.
inline void fill_elements(std::uint8_t* dst, std::uint64_t count, const std::uint8_t* elem,
                          std::uint32_t itemsize) {
  if (elem == nullptr) {
    std::memset(dst, 0, checked_size(count * itemsize, "fill"));
    return;
  }
  bool all_zero = true;
  for (std::uint32_t i = 0; i < itemsize; ++i) {
    all_zero = all_zero && elem[i] == 0;
  }
  if (all_zero) {
    std::memset(dst, 0, checked_size(count * itemsize, "fill"));
    return;
  }
  for (std::uint64_t i = 0; i < count; ++i) {
    std::memcpy(dst + i * itemsize, elem, itemsize);
  }
}

}  // namespace zarr::detail

#endif  // LIBZARR_DETAIL_COMMON_HPP
