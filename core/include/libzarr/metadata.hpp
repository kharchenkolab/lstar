// SPDX-License-Identifier: MIT

#ifndef LIBZARR_METADATA_HPP
#define LIBZARR_METADATA_HPP

#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <optional>
#include <string>
#include <vector>

#include <nlohmann/json.hpp>

#include "libzarr/detail/common.hpp"
#include "libzarr/types.hpp"

/// \file metadata.hpp
/// The normalized, version-independent array metadata (ArrayMeta) plus codec
/// descriptors and fill-value encoding. Version-specific parse/emit lives in
/// v2.hpp / v3.hpp, which lower into these types.

namespace zarr {

/// JSON type used throughout libzarr (vendored nlohmann/json by default;
/// see LIBZARR_EXTERNAL_JSON). Objects keep keys sorted, which makes every
/// serialization deterministic.
using json = nlohmann::json;

/// One codec in a chain, in v3 nomenclature: a name plus a JSON configuration.
/// v2 metadata is lowered into this form on read (compressor -> one
/// bytes->bytes codec; order:"F" -> a transpose codec; dtype byte order ->
/// the "bytes" codec's endian).
struct CodecSpec {
  /// Codec name ("transpose", "bytes", "gzip", "zlib", ...).
  std::string name;
  /// Codec-specific configuration.
  json configuration = json::object();
};

/// Convenience factory: gzip (RFC 1952) at `level` (0-9).
inline CodecSpec gzip(int level = 5) { return {"gzip", {{"level", level}}}; }

/// Convenience factory: zlib (RFC 1950) at `level` (0-9). Zarr v2 only.
inline CodecSpec zlib(int level = 5) { return {"zlib", {{"level", level}}}; }

/// Convenience factory: blosc (v3-style named shuffle: "noshuffle",
/// "shuffle" or "bitshuffle").
inline CodecSpec blosc(const std::string& cname = "lz4", int clevel = 5,
                       const std::string& shuffle = "shuffle") {
  return {"blosc", {{"cname", cname}, {"clevel", clevel}, {"shuffle", shuffle}}};
}

/// Convenience factory: zstd. Level 0 means zstd's default; `checksum`
/// appends a frame checksum.
inline CodecSpec zstd(int level = 0, bool checksum = false) {
  return {"zstd", {{"level", level}, {"checksum", checksum}}};
}

/// Chunk-key scheme. v2: indices joined by a separator ('.' default), rank 0
/// is "0". v3 "default": "c" prefix + separator-joined indices ('/' default),
/// rank 0 is "c". v3's "v2" encoding maps onto `v2`.
enum class ChunkKeyKind : std::uint8_t {
  v2,
  v3_default,
};

/// One level of sharding, outermost first (nested sharding stacks levels).
/// ArrayMeta::chunk_shape is always the innermost (true chunk) shape; each
/// level's shard shape is an integer multiple of the next level's.
struct ShardLevel {
  /// Shard (outer chunk) shape at this level.
  std::vector<std::uint64_t> shard_shape;
  /// Codecs of the shard index (fixed-size: `bytes` and optional `crc32c`).
  std::vector<CodecSpec> index_codecs;
  /// v3 sharding spec index_location: end (default) or start.
  bool index_at_end = true;
};

/// Normalized array metadata, shared by every format version. Version quirks
/// are resolved at parse time; everything downstream (codecs, chunk I/O)
/// consumes only this.
struct ArrayMeta {
  /// Format this array was read from / will be written as.
  ZarrFormat format = ZarrFormat::v2;
  /// Array shape; empty = 0-dimensional.
  std::vector<std::uint64_t> shape;
  /// Chunk shape, same rank as `shape`; chunks may exceed the array extent.
  std::vector<std::uint64_t> chunk_shape;
  /// Element type.
  DataType dtype;
  /// Fill value as one element in native byte order; std::nullopt when the
  /// source metadata had fill_value:null (legal in v2), which reads as zeros.
  std::optional<Bytes> fill;
  /// Chunk-key scheme (see ChunkKeyKind).
  ChunkKeyKind key_encoding = ChunkKeyKind::v2;
  /// Chunk-key separator ('.' or '/').
  char dimension_separator = '.';
  /// v3 dimension_names member, preserved verbatim (null when absent).
  json dimension_names;
  /// Codec chain of the (innermost) chunks, in v3 order: array->array*, one
  /// array->bytes ("bytes"), then bytes->bytes*. sharding_indexed never
  /// appears here — it is lowered into `shard_levels`.
  std::vector<CodecSpec> codecs;
  /// Sharding levels (empty = unsharded); see ShardLevel.
  std::vector<ShardLevel> shard_levels;
  /// User attributes (v2 .zattrs / v3 attributes).
  json attributes = json::object();

  /// Number of elements in the whole array (1 for rank 0).
  [[nodiscard]] std::uint64_t element_count() const {
    return detail::checked_product(shape, "array shape");
  }
  /// Number of elements in one (full) chunk (1 for rank 0).
  [[nodiscard]] std::uint64_t chunk_element_count() const {
    return detail::checked_product(chunk_shape, "chunk shape");
  }
  /// Chunk-grid extent per dimension.
  [[nodiscard]] std::vector<std::uint64_t> grid_shape() const {
    std::vector<std::uint64_t> grid(shape.size());
    for (std::size_t d = 0; d < shape.size(); ++d) {
      grid[d] = detail::ceil_div(shape[d], chunk_shape[d]);
    }
    return grid;
  }
};

/// Options for opening arrays and groups.
struct OpenOptions {
  /// v3 core requires rejecting unrecognized metadata members and null fill
  /// values; lenient mode (opt-in) ignores/defaults them instead. No effect
  /// on v2, which is read tolerantly by design.
  bool lenient = false;
};

/// Serializes JSON in libzarr's canonical form: 4-space indent, sorted keys
/// (nlohmann's storage order), UTF-8. Byte-stable across platforms.
inline Bytes canonical_json_bytes(const json& j) {
  const std::string text = j.dump(4);
  return {text.begin(), text.end()};
}

namespace detail {

/// Runs a metadata-parsing callable, converting any escaping nlohmann
/// exception (type_error from mis-typed members, out_of_range, ...) into
/// zarr::error: malformed input must never surface library-internal
/// exception types.
template <typename Fn>
auto guard_json(const std::string& ctx, const Fn& fn) -> decltype(fn()) {
  try {
    return fn();
  } catch (const json::exception& e) {
    throw error(ctx + ": malformed metadata (" + e.what() + ")");
  }
}

/// Reads a JSON value as uint64. Handles nlohmann's split integer storage:
/// parsed non-negative literals are number_unsigned, programmatically
/// constructed ints are number_integer.
inline std::uint64_t json_to_uint64(const json& v, const std::string& ctx) {
  if (v.is_number_unsigned()) {
    return v.get<std::uint64_t>();
  }
  if (v.is_number_integer()) {
    const auto i = v.get<std::int64_t>();
    if (i < 0) {
      throw error(ctx + ": expected a non-negative integer, got " + std::to_string(i));
    }
    return static_cast<std::uint64_t>(i);
  }
  throw error(ctx + ": expected a non-negative integer, got " + v.dump());
}

/// Reads an optional integer member, tolerating numeric strings: NCZarr
/// (libnetcdf 4.9.x) writes numbers like compressor levels and filter
/// element sizes as JSON strings ("1", "0").
inline std::int64_t lenient_int(const json& obj, const char* key, std::int64_t fallback,
                                const std::string& ctx) {
  const auto it = obj.find(key);
  if (it == obj.end()) {
    return fallback;
  }
  if (it->is_number_integer() || it->is_number_unsigned()) {
    return it->get<std::int64_t>();
  }
  if (it->is_string()) {
    const auto s = it->get<std::string>();
    if (!s.empty()) {
      char* end = nullptr;
      const long long v = std::strtoll(s.c_str(), &end, 10);
      if (end == s.c_str() + s.size()) {
        return static_cast<std::int64_t>(v);
      }
    }
  }
  throw error(ctx + ": '" + key + "' must be an integer, got " + it->dump());
}

template <typename T>
Bytes scalar_bytes(T value) {
  Bytes out(sizeof(T));
  std::memcpy(out.data(), &value, sizeof(T));
  return out;
}

/// Parses a JSON array of non-negative integers (shapes, chunk shapes).
inline std::vector<std::uint64_t> parse_extents(const json& v, const char* name,
                                                const std::string& ctx) {
  if (!v.is_array()) {
    throw error(ctx + ": '" + name + "' must be an array");
  }
  std::vector<std::uint64_t> out;
  out.reserve(v.size());
  for (const json& e : v) {
    out.push_back(json_to_uint64(e, ctx + ": " + name));
  }
  return out;
}

/// Pinned quiet-NaN bit patterns: the specs' "NaN" form carries no payload,
/// so we fix one for byte-stable output (f16 0x7e00, f32 0x7fc00000, f64
/// 0x7ff8000000000000).
inline Bytes quiet_nan_bytes(DType kind) {
  if (kind == DType::float16) {
    return scalar_bytes<std::uint16_t>(0x7e00U);
  }
  if (kind == DType::float32) {
    return scalar_bytes<std::uint32_t>(0x7fc00000U);
  }
  assert(kind == DType::float64);
  return scalar_bytes<std::uint64_t>(0x7ff8000000000000ULL);
}

/// +/- infinity, encoded per float kind.
inline Bytes infinity_bytes(DType kind, bool negative) {
  if (kind == DType::float16) {
    return scalar_bytes<std::uint16_t>(negative ? 0xfc00U : 0x7c00U);
  }
  if (kind == DType::float32) {
    const float inf = std::numeric_limits<float>::infinity();
    return scalar_bytes(negative ? -inf : inf);
  }
  assert(kind == DType::float64);
  const double inf = std::numeric_limits<double>::infinity();
  return scalar_bytes(negative ? -inf : inf);
}

inline Bytes fill_from_double(double value, DataType dt, const std::string& ctx);

/// Encodes a JSON integer fill for `dt`, range-checked with a precise error.
inline Bytes fill_from_int(std::int64_t value, DataType dt, const std::string& ctx) {
  const auto check = [&](std::int64_t lo, std::int64_t hi) {
    if (value < lo || value > hi) {
      throw error(ctx + ": fill_value " + std::to_string(value) + " out of range for dtype");
    }
  };
  switch (dt.kind) {
    case DType::boolean:
      check(0, 1);
      return scalar_bytes<std::uint8_t>(static_cast<std::uint8_t>(value));
    case DType::int8:
      check(std::numeric_limits<std::int8_t>::min(), std::numeric_limits<std::int8_t>::max());
      return scalar_bytes<std::int8_t>(static_cast<std::int8_t>(value));
    case DType::int16:
      check(std::numeric_limits<std::int16_t>::min(), std::numeric_limits<std::int16_t>::max());
      return scalar_bytes<std::int16_t>(static_cast<std::int16_t>(value));
    case DType::int32:
      check(std::numeric_limits<std::int32_t>::min(), std::numeric_limits<std::int32_t>::max());
      return scalar_bytes<std::int32_t>(static_cast<std::int32_t>(value));
    case DType::int64:
      return scalar_bytes<std::int64_t>(value);
    case DType::uint8:
    case DType::uint16:
    case DType::uint32:
    case DType::uint64: {
      if (value < 0) {
        throw error(ctx + ": fill_value " + std::to_string(value) +
                    " is negative for unsigned dtype");
      }
      const auto u = static_cast<std::uint64_t>(value);
      if (dt.kind == DType::uint8 && u > std::numeric_limits<std::uint8_t>::max()) {
        check(0, 255);
      }
      if (dt.kind == DType::uint16 && u > std::numeric_limits<std::uint16_t>::max()) {
        check(0, 65535);
      }
      if (dt.kind == DType::uint32 && u > std::numeric_limits<std::uint32_t>::max()) {
        check(0, 4294967295LL);
      }
      if (dt.kind == DType::uint8) {
        return scalar_bytes<std::uint8_t>(static_cast<std::uint8_t>(u));
      }
      if (dt.kind == DType::uint16) {
        return scalar_bytes<std::uint16_t>(static_cast<std::uint16_t>(u));
      }
      if (dt.kind == DType::uint32) {
        return scalar_bytes<std::uint32_t>(static_cast<std::uint32_t>(u));
      }
      return scalar_bytes<std::uint64_t>(u);
    }
    case DType::float32:
    case DType::float64:
      return fill_from_double(static_cast<double>(value), dt, ctx);
    default:
      throw error(ctx + ": numeric fill_value invalid for this dtype");
  }
}

/// Encodes a JSON unsigned fill (needed for uint64 values >= 2^63, which do
/// not fit int64_t — a known interop trap).
inline Bytes fill_from_uint(std::uint64_t value, DataType dt, const std::string& ctx) {
  if (value <= static_cast<std::uint64_t>(std::numeric_limits<std::int64_t>::max())) {
    return fill_from_int(static_cast<std::int64_t>(value), dt, ctx);
  }
  if (dt.kind != DType::uint64) {
    throw error(ctx + ": fill_value " + std::to_string(value) + " out of range for dtype");
  }
  return scalar_bytes<std::uint64_t>(value);
}

inline Bytes fill_from_double(double value, DataType dt, const std::string& ctx) {
  switch (dt.kind) {
    case DType::float16:
      return scalar_bytes<std::uint16_t>(double_to_half_bits(value));
    case DType::float32:
      return scalar_bytes<float>(static_cast<float>(value));
    case DType::float64:
      return scalar_bytes<double>(value);
    default:
      // Tolerance: integral fills occasionally arrive as JSON floats
      // (e.g. 1.0 for an int dtype); accept when exactly integral.
      if (std::nearbyint(value) == value &&
          value >= static_cast<double>(std::numeric_limits<std::int64_t>::min()) &&
          value <= static_cast<double>(std::numeric_limits<std::int64_t>::max())) {
        return fill_from_int(static_cast<std::int64_t>(value), dt, ctx);
      }
      throw error(ctx + ": non-integral fill_value for integer dtype");
  }
}

/// Emits a normalized fill (native-order element bytes) as JSON, dtype-directed.
/// Floats use the spec string forms for non-finite values; NaN/Infinity must
/// never reach nlohmann as doubles (it would serialize them as null).
inline json fill_to_json(const std::optional<Bytes>& fill, DataType dt) {
  if (!fill) {
    return nullptr;
  }
  const std::uint8_t* p = fill->data();
  const auto load = [&](auto probe) {
    decltype(probe) v;
    std::memcpy(&v, p, sizeof(v));
    return v;
  };
  switch (dt.kind) {
    case DType::boolean:
      return load(std::uint8_t{}) != 0;
    case DType::int8:
      return load(std::int8_t{});
    case DType::int16:
      return load(std::int16_t{});
    case DType::int32:
      return load(std::int32_t{});
    case DType::int64:
      return load(std::int64_t{});
    case DType::uint8:
      return load(std::uint8_t{});
    case DType::uint16:
      return load(std::uint16_t{});
    case DType::uint32:
      return load(std::uint32_t{});
    case DType::uint64:
      return load(std::uint64_t{});
    case DType::float16:
    case DType::float32:
    case DType::float64: {
      double v = 0;
      if (dt.kind == DType::float16) {
        v = half_bits_to_double(load(std::uint16_t{}));
      } else if (dt.kind == DType::float32) {
        v = static_cast<double>(load(float{}));
      } else {
        v = load(double{});
      }
      if (std::isnan(v)) {
        return "NaN";
      }
      if (std::isinf(v)) {
        return v > 0 ? "Infinity" : "-Infinity";
      }
      return v;
    }
    case DType::complex64:
    case DType::complex128: {
      const DType component = dt.kind == DType::complex64 ? DType::float32 : DType::float64;
      const std::uint32_t half = dt.itemsize / 2;
      return json::array({fill_to_json(Bytes(p, p + half), DataType::of(component)),
                          fill_to_json(Bytes(p + half, p + dt.itemsize), DataType::of(component))});
    }
    case DType::raw:
      return base64_encode(p, dt.itemsize);
    default:
      throw error("fill_value emission not implemented for this dtype");
  }
}

}  // namespace detail

}  // namespace zarr

#endif  // LIBZARR_METADATA_HPP
