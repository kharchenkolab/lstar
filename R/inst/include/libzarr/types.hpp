// SPDX-License-Identifier: MIT

#ifndef LIBZARR_TYPES_HPP
#define LIBZARR_TYPES_HPP

#include <cassert>
#include <cstdint>
#include <stdexcept>
#include <vector>

/// \file types.hpp
/// Fundamental types shared across libzarr: the error type, the byte buffer
/// alias, data types, and library version macros.

// Macros, not constants: consumers need these preprocessor-testable.
// NOLINTBEGIN(modernize-macro-to-enum,cppcoreguidelines-macro-to-enum)
/// Library major version.
#define LIBZARR_VERSION_MAJOR 0
/// Library minor version.
#define LIBZARR_VERSION_MINOR 1
/// Library patch version.
#define LIBZARR_VERSION_PATCH 0
// NOLINTEND(modernize-macro-to-enum,cppcoreguidelines-macro-to-enum)

namespace zarr {

/// Exception thrown for every failure reachable from user input or store
/// contents (malformed metadata, unknown codecs, out-of-range reads, ...).
/// Messages are precise and self-contained. Internal invariant violations use
/// assert() instead and never throw.
class error : public std::runtime_error {
 public:
  using std::runtime_error::runtime_error;
};

/// Owned byte buffer used throughout the value-based public API.
using Bytes = std::vector<std::uint8_t>;

/// Zarr storage format version.
enum class ZarrFormat : std::uint8_t {
  v2 = 2,
  v3 = 3,
};

/// Element type kinds. `raw` is a fixed-size byte string (v2 `|V<n>`, v3
/// `r<bits>`) whose size lives in DataType::itemsize.
enum class DType : std::uint8_t {
  boolean,
  int8,
  int16,
  int32,
  int64,
  uint8,
  uint16,
  uint32,
  uint64,
  float16,
  float32,
  float64,
  complex64,
  complex128,
  raw,
};

/// True for float16/float32/float64.
constexpr bool is_float(DType kind) {
  return kind == DType::float16 || kind == DType::float32 || kind == DType::float64;
}

/// True for int8..int64.
constexpr bool is_signed_int(DType kind) {
  return kind == DType::int8 || kind == DType::int16 || kind == DType::int32 ||
         kind == DType::int64;
}

/// True for uint8..uint64.
constexpr bool is_unsigned_int(DType kind) {
  return kind == DType::uint8 || kind == DType::uint16 || kind == DType::uint32 ||
         kind == DType::uint64;
}

/// True for complex64/complex128.
constexpr bool is_complex(DType kind) {
  return kind == DType::complex64 || kind == DType::complex128;
}

/// Element size in bytes fixed by the kind; 0 for DType::raw (whose size is
/// per-array).
constexpr std::uint32_t fixed_itemsize(DType kind) {
  switch (kind) {
    case DType::boolean:
    case DType::int8:
    case DType::uint8:
      return 1;
    case DType::int16:
    case DType::uint16:
    case DType::float16:
      return 2;
    case DType::int32:
    case DType::uint32:
    case DType::float32:
      return 4;
    case DType::int64:
    case DType::uint64:
    case DType::float64:
    case DType::complex64:
      return 8;
    case DType::complex128:
      return 16;
    case DType::raw:
      return 0;
  }
  return 0;  // unreachable; keeps compilers satisfied
}

/// A concrete element type: kind plus size (the size only varies for raw).
struct DataType {
  /// Element type kind.
  DType kind = DType::uint8;
  /// Element size in bytes.
  std::uint32_t itemsize = 1;

  /// A DataType of fixed-size `kind` (anything but raw).
  static constexpr DataType of(DType kind) {
    assert(kind != DType::raw);
    return DataType{kind, fixed_itemsize(kind)};
  }

  /// A raw byte-string type of `size` bytes (v2 `|V<size>`).
  static constexpr DataType raw_bytes(std::uint32_t size) { return DataType{DType::raw, size}; }

  /// Equal kind and size.
  friend constexpr bool operator==(DataType a, DataType b) {
    return a.kind == b.kind && a.itemsize == b.itemsize;
  }
  /// Differing kind or size.
  friend constexpr bool operator!=(DataType a, DataType b) { return !(a == b); }
};

}  // namespace zarr

#endif  // LIBZARR_TYPES_HPP
