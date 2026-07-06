// SPDX-License-Identifier: MIT

#ifndef LIBZARR_V3_HPP
#define LIBZARR_V3_HPP

#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <initializer_list>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

#include "libzarr/detail/common.hpp"
#include "libzarr/metadata.hpp"
#include "libzarr/types.hpp"
#include "libzarr/v2.hpp"

/// \file v3.hpp
/// Zarr v3 metadata (zarr.json): read-side parsing into normalized ArrayMeta.
/// Strict by the v3.1 core spec (unrecognized members are errors) with an
/// opt-in lenient mode; every deliberate read-tolerance cites its origin.

namespace zarr::v3 {

/// v3 metadata document name.
inline constexpr const char* kMetaKey = "zarr.json";

/// Store key of the metadata document for the node at `path` ("" = root).
inline std::string meta_key(const std::string& path) {
  return path.empty() ? kMetaKey : path + "/" + kMetaKey;
}

namespace detail_v3 {

/// v3 core: an implementation MUST fail on metadata members it does not
/// recognize, unless the member is an object carrying "must_understand":
/// false. Lenient mode (opt-in) ignores them instead.
inline void check_members(const json& j, std::initializer_list<std::string_view> known,
                          const std::string& ctx, bool lenient) {
  if (lenient) {
    return;
  }
  for (const auto& item : j.items()) {
    bool recognized = false;
    for (const std::string_view name : known) {
      recognized = recognized || item.key() == name;
    }
    if (recognized) {
      continue;
    }
    if (item.value().is_object() && !item.value().value("must_understand", true)) {
      continue;
    }
    throw error(ctx + ": unknown metadata member '" + item.key() +
                "' (v3 core requires rejecting unrecognized members; open in lenient mode to "
                "ignore it)");
  }
}

inline const json& require(const json& j, const char* name, const std::string& ctx) {
  const auto it = j.find(name);
  if (it == j.end()) {
    throw error(ctx + ": missing required member '" + name + "'");
  }
  return *it;
}

/// Value of one hex digit, or 99 when invalid.
inline std::uint32_t hex_digit(char c) {
  if (c >= '0' && c <= '9') {
    return static_cast<std::uint32_t>(c - '0');
  }
  if (c >= 'a' && c <= 'f') {
    return static_cast<std::uint32_t>(c - 'a') + 10U;
  }
  if (c >= 'A' && c <= 'F') {
    return static_cast<std::uint32_t>(c - 'A') + 10U;
  }
  return 99;
}

/// Parses "0x..."/"0b..." bit-pattern strings into `itemsize` bytes, given in
/// the written (big-endian numeral) order.
inline std::optional<Bytes> parse_bit_string(const std::string& s, std::uint32_t itemsize,
                                             const std::string& ctx) {
  const bool hex = detail::starts_with(s, "0x");
  const bool bin = detail::starts_with(s, "0b");
  if (!hex && !bin) {
    return std::nullopt;
  }
  const std::string_view digits = std::string_view(s).substr(2);
  const std::size_t per_byte = hex ? 2 : 8;
  const std::uint32_t base = hex ? 16U : 2U;
  if (digits.size() != itemsize * per_byte) {
    throw error(ctx + ": bit-pattern fill_value '" + s + "' must have " +
                std::to_string(itemsize * per_byte) + (hex ? " hex" : " binary") + " digits");
  }
  Bytes out(itemsize);
  for (std::uint32_t b = 0; b < itemsize; ++b) {
    std::uint32_t value = 0;
    for (std::size_t d = 0; d < per_byte; ++d) {
      const std::uint32_t digit = hex_digit(digits[b * per_byte + d]);
      if (digit >= base) {
        std::string msg = ctx;
        msg += ": invalid digit in fill_value '";
        msg += s;
        msg += "'";
        throw error(msg);
      }
      value = value * base + digit;
    }
    out[b] = static_cast<std::uint8_t>(value);
  }
  return out;
}

/// One float component per the v3 fill_value rules: number, "NaN",
/// "Infinity", "-Infinity", or a bit-pattern string.
inline Bytes float_fill(const json& v, DType kind, std::uint32_t itemsize, const std::string& ctx) {
  if (v.is_number()) {
    return detail::fill_from_double(v.get<double>(), DataType{kind, itemsize}, ctx);
  }
  if (v.is_string()) {
    const auto s = v.get<std::string>();
    if (s == "NaN") {
      return detail::quiet_nan_bytes(kind);
    }
    if (s == "Infinity") {
      return detail::infinity_bytes(kind, false);
    }
    if (s == "-Infinity") {
      return detail::infinity_bytes(kind, true);
    }
    if (auto bits = parse_bit_string(s, itemsize, ctx)) {
      // v3 core: the string is the numeral of the bit pattern (big-endian
      // digits); convert to native byte order.
      if (detail::host_is_little_endian()) {
        std::reverse(bits->begin(), bits->end());
      }
      return *std::move(bits);
    }
    throw error(ctx + ": cannot interpret float fill_value '" + s + "'");
  }
  throw error(ctx + ": cannot interpret float fill_value " + v.dump());
}

}  // namespace detail_v3

/// Parses a v3 data_type name. Extension (object) data types are rejected.
inline DataType parse_data_type(const json& v, const std::string& ctx) {
  if (!v.is_string()) {
    throw error(ctx + ": extension data_type objects are not supported");
  }
  const auto s = v.get<std::string>();
  struct NamedType {
    std::string_view name;
    DType kind;
  };
  static constexpr std::array<NamedType, 14> kNames{{{"bool", DType::boolean},
                                                     {"int8", DType::int8},
                                                     {"int16", DType::int16},
                                                     {"int32", DType::int32},
                                                     {"int64", DType::int64},
                                                     {"uint8", DType::uint8},
                                                     {"uint16", DType::uint16},
                                                     {"uint32", DType::uint32},
                                                     {"uint64", DType::uint64},
                                                     {"float16", DType::float16},
                                                     {"float32", DType::float32},
                                                     {"float64", DType::float64},
                                                     {"complex64", DType::complex64},
                                                     {"complex128", DType::complex128}}};
  for (const NamedType& named : kNames) {
    if (s == named.name) {
      return DataType::of(named.kind);
    }
  }
  if (s.size() > 1 && s[0] == 'r') {
    std::uint64_t bits = 0;
    for (std::size_t i = 1; i < s.size(); ++i) {
      if (s[i] < '0' || s[i] > '9' || bits > 0xFFFFFF) {
        bits = 0;
        break;
      }
      bits = bits * 10 + static_cast<std::uint64_t>(s[i] - '0');
    }
    if (bits == 0 || bits % 8 != 0) {
      throw error(ctx + ": raw data_type '" + s + "' must be r<bits> with bits a multiple of 8");
    }
    return DataType::raw_bytes(static_cast<std::uint32_t>(bits / 8));
  }
  throw error(ctx + ": unknown data_type '" + s + "'");
}

namespace detail_v3 {

/// Raw (r<bits>) fills: a bit-pattern string of exactly itemsize bytes (kept
/// in the written byte order), or an array of byte values.
inline Bytes raw_fill(const json& v, DataType dt, const std::string& ctx) {
  if (v.is_string()) {
    if (auto bits = parse_bit_string(v.get<std::string>(), dt.itemsize, ctx)) {
      return *std::move(bits);
    }
    throw error(ctx + ": raw fill_value string must be a 0x/0b bit pattern");
  }
  if (v.is_array() && v.size() == dt.itemsize) {
    Bytes out(dt.itemsize);
    for (std::uint32_t i = 0; i < dt.itemsize; ++i) {
      const std::uint64_t byte = detail::json_to_uint64(v[i], ctx + ": fill_value");
      if (byte > 0xFF) {
        throw error(ctx + ": raw fill_value bytes must be 0..255");
      }
      out[i] = static_cast<std::uint8_t>(byte);
    }
    return out;
  }
  throw error(ctx + ": raw fill_value must be a bit-pattern string or an array of " +
              std::to_string(dt.itemsize) + " byte values");
}

}  // namespace detail_v3

/// Parses a v3 fill_value, dtype-directed (v3 core "Fill value" table).
inline std::optional<Bytes> parse_fill(const json& v, DataType dt, const std::string& ctx,
                                       bool lenient) {
  if (v.is_null()) {
    // v3 requires a concrete fill_value; null appears from pre-final writers.
    if (lenient) {
      return std::nullopt;
    }
    throw error(ctx + ": v3 fill_value must not be null (open in lenient mode to read as zeros)");
  }
  switch (dt.kind) {
    case DType::boolean:
      if (!v.is_boolean()) {
        throw error(ctx + ": bool fill_value must be true or false");
      }
      return detail::scalar_bytes<std::uint8_t>(v.get<bool>() ? 1 : 0);
    case DType::float16:
    case DType::float32:
    case DType::float64:
      return detail_v3::float_fill(v, dt.kind, dt.itemsize, ctx);
    case DType::complex64:
    case DType::complex128: {
      // v3 core: complex fill_value is a [real, imaginary] two-element array.
      if (!v.is_array() || v.size() != 2) {
        throw error(ctx + ": complex fill_value must be a [re, im] array");
      }
      const DType component = dt.kind == DType::complex64 ? DType::float32 : DType::float64;
      const std::uint32_t half = dt.itemsize / 2;
      Bytes out = detail_v3::float_fill(v[0], component, half, ctx);
      const Bytes imag = detail_v3::float_fill(v[1], component, half, ctx);
      out.insert(out.end(), imag.begin(), imag.end());
      return out;
    }
    case DType::raw:
      return detail_v3::raw_fill(v, dt, ctx);
    default:
      if (v.is_number_unsigned()) {
        return detail::fill_from_uint(v.get<std::uint64_t>(), dt, ctx);
      }
      if (v.is_number_integer()) {
        return detail::fill_from_int(v.get<std::int64_t>(), dt, ctx);
      }
      if (v.is_number_float()) {
        return detail::fill_from_double(v.get<double>(), dt, ctx);
      }
      throw error(ctx + ": integer fill_value expected, got " + v.dump());
  }
}

namespace detail_v3 {

/// Normalizes one codecs-list entry, applying documented read-tolerances
/// for pre-final spellings.
inline CodecSpec parse_codec_entry(const json& c, std::size_t rank, const std::string& ctx) {
  CodecSpec spec;
  if (c.is_string()) {
    // Pre-final v3 writers emitted bare codec-name strings (read tolerance).
    spec.name = c.get<std::string>();
  } else if (c.is_object() && c.contains("name") && c["name"].is_string()) {
    spec.name = c["name"].get<std::string>();
    spec.configuration = c.value("configuration", json::object());
    if (!spec.configuration.is_object()) {
      throw error(ctx + ": codec '" + spec.name + "' configuration must be an object");
    }
  } else {
    throw error(ctx + ": each codec must be an object with a 'name'");
  }
  if (spec.name == "endian") {
    // zarr-python 2.x's experimental v3 wrote "endian" for what the final
    // spec names "bytes" (read tolerance).
    spec.name = "bytes";
  }
  if (spec.name == "transpose" && spec.configuration.value("order", json()).is_string()) {
    // 2022-draft transpose configs used "C"/"F" strings; the final spec
    // requires an explicit permutation (read tolerance).
    const auto order = spec.configuration["order"].get<std::string>();
    if (order != "C" && order != "F") {
      std::string msg = ctx;
      msg += ": transpose order '";
      msg += order;
      msg += R"(' is not "C", "F" or an array)";
      throw error(msg);
    }
    json perm = json::array();
    for (std::size_t d = 0; d < rank; ++d) {
      perm.push_back(order == "F" ? rank - 1 - d : d);
    }
    spec.configuration["order"] = perm;
  }
  return spec;
}

/// Lowers the v3 codecs member into CodecSpec form.
inline std::vector<CodecSpec> parse_codecs(const json& v, std::size_t rank,
                                           const std::string& ctx) {
  if (!v.is_array()) {
    throw error(ctx + ": 'codecs' must be an array");
  }
  std::vector<CodecSpec> out;
  out.reserve(v.size());
  for (const json& c : v) {
    out.push_back(parse_codec_entry(c, rank, ctx));
  }
  return out;
}

/// Fills meta.key_encoding / dimension_separator from chunk_key_encoding.
inline void parse_chunk_key_encoding(const json& cke, ArrayMeta& meta, const std::string& ctx) {
  if (!cke.is_object() || !cke.contains("name")) {
    throw error(ctx + ": chunk_key_encoding must be an object with a 'name'");
  }
  if (cke["name"] == "default") {
    meta.key_encoding = ChunkKeyKind::v3_default;
    meta.dimension_separator = '/';
  } else if (cke["name"] == "v2") {
    meta.key_encoding = ChunkKeyKind::v2;
    meta.dimension_separator = '.';
  } else {
    throw error(ctx + ": unknown chunk_key_encoding '" + cke["name"].dump() + "'");
  }
  const json config = cke.value("configuration", json::object());
  const json separator = config.value("separator", json());
  if (!separator.is_null()) {
    if (!separator.is_string() || (separator != "/" && separator != ".")) {
      throw error(ctx + R"(: chunk_key_encoding separator must be "/" or ".")");
    }
    meta.dimension_separator = separator.get<std::string>()[0];
  }
}

/// Validates and stores the optional dimension_names member.
inline void parse_dimension_names(const json& j, ArrayMeta& meta, const std::string& ctx) {
  const auto it = j.find("dimension_names");
  if (it == j.end()) {
    return;
  }
  if (!it->is_array() || it->size() != meta.shape.size()) {
    throw error(ctx + ": 'dimension_names' must be an array of rank length");
  }
  for (const json& name : *it) {
    if (!name.is_string() && !name.is_null()) {
      throw error(ctx + ": dimension names must be strings or null");
    }
  }
  meta.dimension_names = *it;
}

/// Lowers sharding_indexed codecs into ArrayMeta::shard_levels, recursively
/// (nested sharding stacks levels). After lowering, `codecs` holds only the
/// innermost chunk chain.
inline void lower_sharding(ArrayMeta& meta, const std::string& ctx) {
  while (!meta.codecs.empty() && meta.codecs.front().name == "sharding_indexed") {
    if (meta.codecs.size() != 1) {
      // Ranges into a shard must map 1:1 onto stored bytes; codecs wrapped
      // around the shard (outer transpose, whole-shard compression) break
      // that, so they are rejected rather than silently degraded.
      throw error(ctx +
                  ": sharding_indexed cannot be combined with other codecs at the same "
                  "level");
    }
    const json config = meta.codecs.front().configuration.is_object()
                            ? meta.codecs.front().configuration
                            : json::object();
    const std::vector<std::uint64_t> inner_shape = detail::parse_extents(
        require(config, "chunk_shape", ctx + ": sharding_indexed"), "chunk_shape", ctx);
    if (inner_shape.size() != meta.chunk_shape.size()) {
      throw error(ctx + ": sharding_indexed chunk_shape rank mismatch");
    }
    for (std::size_t d = 0; d < inner_shape.size(); ++d) {
      // v3 sharding spec: the inner chunk shape must evenly divide the shard.
      if (inner_shape[d] == 0 || meta.chunk_shape[d] % inner_shape[d] != 0) {
        throw error(ctx + ": sharding_indexed chunk_shape must evenly divide the shard shape");
      }
    }

    ShardLevel level;
    level.shard_shape = meta.chunk_shape;
    level.index_codecs =
        parse_codecs(require(config, "index_codecs", ctx + ": sharding_indexed"), 1, ctx);
    const json location = config.value("index_location", json("end"));
    if (location != "end" && location != "start") {
      throw error(ctx + R"(: index_location must be "end" or "start")");
    }
    level.index_at_end = location == "end";

    meta.shard_levels.push_back(std::move(level));
    meta.chunk_shape = inner_shape;
    meta.codecs =
        parse_codecs(require(config, "codecs", ctx + ": sharding_indexed"), meta.shape.size(), ctx);
  }
  for (const CodecSpec& codec : meta.codecs) {
    if (codec.name == "sharding_indexed") {
      throw error(ctx + ": sharding_indexed must be the sole codec of its level");
    }
  }
}

}  // namespace detail_v3

namespace detail_v3 {

inline ArrayMeta parse_array_meta_impl(const json& j, const std::string& ctx, bool lenient) {
  if (!j.is_object()) {
    throw error(ctx + ": expected a JSON object");
  }
  if (detail::json_to_uint64(detail_v3::require(j, "zarr_format", ctx), ctx + ": zarr_format") !=
      3) {
    throw error(ctx + ": zarr_format must be 3");
  }
  if (detail_v3::require(j, "node_type", ctx) != "array") {
    throw error(ctx + ": node_type must be 'array'");
  }
  detail_v3::check_members(
      j,
      {"zarr_format", "node_type", "shape", "data_type", "chunk_grid", "chunk_key_encoding",
       "fill_value", "codecs", "attributes", "dimension_names", "storage_transformers"},
      ctx, lenient);

  ArrayMeta meta;
  meta.format = ZarrFormat::v3;
  meta.shape = detail::parse_extents(detail_v3::require(j, "shape", ctx), "shape", ctx);
  meta.dtype = parse_data_type(detail_v3::require(j, "data_type", ctx), ctx);

  const json& grid = detail_v3::require(j, "chunk_grid", ctx);
  if (!grid.is_object() || grid.value("name", "") != std::string("regular")) {
    throw error(ctx + ": only the 'regular' chunk_grid is supported");
  }
  // Bind `grid_ctx` as a named lvalue: passing `ctx + "..."` directly trips
  // gcc's -Wdangling-reference on the reference initialization below.
  const std::string grid_ctx = ctx + ": chunk_grid";
  const json& grid_config = detail_v3::require(grid, "configuration", grid_ctx);
  meta.chunk_shape = detail::parse_extents(detail_v3::require(grid_config, "chunk_shape", grid_ctx),
                                           "chunk_shape", ctx);
  if (meta.chunk_shape.size() != meta.shape.size()) {
    throw error(ctx + ": chunk_shape rank " + std::to_string(meta.chunk_shape.size()) +
                " != shape rank " + std::to_string(meta.shape.size()));
  }
  for (const std::uint64_t c : meta.chunk_shape) {
    if (c == 0) {
      throw error(ctx + ": chunk extents must be positive");
    }
  }

  detail_v3::parse_chunk_key_encoding(detail_v3::require(j, "chunk_key_encoding", ctx), meta, ctx);

  meta.fill = parse_fill(detail_v3::require(j, "fill_value", ctx), meta.dtype, ctx + ": fill_value",
                         lenient);
  meta.codecs =
      detail_v3::parse_codecs(detail_v3::require(j, "codecs", ctx), meta.shape.size(), ctx);
  detail_v3::lower_sharding(meta, ctx);

  meta.attributes = j.value("attributes", json::object());
  if (!meta.attributes.is_object()) {
    throw error(ctx + ": 'attributes' must be an object");
  }
  detail_v3::parse_dimension_names(j, meta, ctx);

  const auto st_it = j.find("storage_transformers");
  if (st_it != j.end() && !(st_it->is_array() && st_it->empty())) {
    throw error(ctx + ": storage_transformers are not supported");
  }
  return meta;
}

}  // namespace detail_v3

/// Parses a v3 array zarr.json document into normalized ArrayMeta.
inline ArrayMeta parse_array_meta(const json& j, const std::string& ctx, bool lenient = false) {
  return detail::guard_json(ctx, [&] { return detail_v3::parse_array_meta_impl(j, ctx, lenient); });
}

/// Result of parsing a v3 group zarr.json.
struct GroupMeta {
  /// User attributes.
  json attributes = json::object();
  /// Inline consolidated metadata (node path -> zarr.json document), per the
  /// zarr-python convention (zarr-specs #309 — a convention, not yet an
  /// accepted spec); std::nullopt when absent.
  std::optional<json> consolidated;
};

/// Parses a v3 group zarr.json document.
inline GroupMeta parse_group_meta(const json& j, const std::string& ctx, bool lenient = false) {
  return detail::guard_json(ctx, [&]() -> GroupMeta {
    if (!j.is_object()) {
      throw error(ctx + ": expected a JSON object");
    }
    if (detail::json_to_uint64(detail_v3::require(j, "zarr_format", ctx), ctx + ": zarr_format") !=
        3) {
      throw error(ctx + ": zarr_format must be 3");
    }
    if (detail_v3::require(j, "node_type", ctx) != "group") {
      throw error(ctx + ": node_type must be 'group'");
    }
    detail_v3::check_members(j, {"zarr_format", "node_type", "attributes", "consolidated_metadata"},
                             ctx, lenient);

    GroupMeta meta;
    meta.attributes = j.value("attributes", json::object());
    const auto cons_it = j.find("consolidated_metadata");
    if (cons_it != j.end() && cons_it->is_object() && cons_it->contains("metadata") &&
        (*cons_it)["metadata"].is_object()) {
      meta.consolidated = (*cons_it)["metadata"];
    }
    return meta;
  });
}

/// v3 chunk key relative to the array (v3 core "chunk key encoding").
inline std::string chunk_key(const std::vector<std::uint64_t>& index, char separator) {
  std::string key = "c";  // rank 0: the key is exactly "c"
  for (const std::uint64_t i : index) {
    key += separator;
    key += std::to_string(i);
  }
  return key;
}

// ---- emission (canonical, deterministic) ------------------------------------

/// Emits the canonical v3 data_type name.
inline std::string emit_data_type(DataType dt) {
  switch (dt.kind) {
    case DType::boolean:
      return "bool";
    case DType::int8:
      return "int8";
    case DType::int16:
      return "int16";
    case DType::int32:
      return "int32";
    case DType::int64:
      return "int64";
    case DType::uint8:
      return "uint8";
    case DType::uint16:
      return "uint16";
    case DType::uint32:
      return "uint32";
    case DType::uint64:
      return "uint64";
    case DType::float16:
      return "float16";
    case DType::float32:
      return "float32";
    case DType::float64:
      return "float64";
    case DType::complex64:
      return "complex64";
    case DType::complex128:
      return "complex128";
    case DType::raw:
      return "r" + std::to_string(std::uint64_t{dt.itemsize} * 8);
  }
  throw error("v3 emission not implemented for this dtype");
}

namespace detail_v3 {

/// "0x..." form of native-order `bytes` (emitted as a big-endian numeral,
/// matching the parse direction).
inline std::string hex_bit_string(const Bytes& bytes, bool reverse_for_endianness) {
  constexpr std::string_view kDigits = "0123456789abcdef";
  std::string out = "0x";
  for (std::size_t i = 0; i < bytes.size(); ++i) {
    const std::size_t at =
        reverse_for_endianness && detail::host_is_little_endian() ? bytes.size() - 1 - i : i;
    out.push_back(kDigits[static_cast<std::size_t>(bytes[at]) >> 4U]);
    out.push_back(kDigits[static_cast<std::size_t>(bytes[at]) & 0x0FU]);
  }
  return out;
}

/// One float component as canonical v3 JSON. Non-finite values use the spec
/// strings; a NaN with a non-default payload must use the hex form
/// (v3 core: the hex form is the only NaN-payload representation).
inline json emit_float_fill(const std::uint8_t* data, DType kind, std::uint32_t width) {
  double v = 0;
  if (kind == DType::float16) {
    std::uint16_t bits = 0;
    std::memcpy(&bits, data, 2);
    v = detail::half_bits_to_double(bits);
  } else if (kind == DType::float32) {
    float f = 0;
    std::memcpy(&f, data, 4);
    v = static_cast<double>(f);
  } else {
    std::memcpy(&v, data, 8);
  }
  if (std::isnan(v)) {
    const Bytes bits(data, data + width);
    if (bits == detail::quiet_nan_bytes(kind)) {
      return "NaN";
    }
    return hex_bit_string(bits, /*reverse_for_endianness=*/true);
  }
  if (std::isinf(v)) {
    return v > 0 ? "Infinity" : "-Infinity";
  }
  return v;
}

}  // namespace detail_v3

/// Emits a normalized fill as canonical v3 JSON. A missing fill (legal only
/// on leniently-read metadata) is synthesized as zeros: v3 requires a
/// concrete fill_value.
inline json emit_fill(const std::optional<Bytes>& fill, DataType dt) {
  const Bytes zeros(dt.itemsize, 0);
  const Bytes& bytes = fill ? *fill : zeros;
  switch (dt.kind) {
    case DType::boolean:
      return bytes[0] != 0;
    case DType::float16:
    case DType::float32:
    case DType::float64:
      return detail_v3::emit_float_fill(bytes.data(), dt.kind, dt.itemsize);
    case DType::complex64:
    case DType::complex128: {
      const DType component = dt.kind == DType::complex64 ? DType::float32 : DType::float64;
      const std::uint32_t half = dt.itemsize / 2;
      return json::array({detail_v3::emit_float_fill(bytes.data(), component, half),
                          detail_v3::emit_float_fill(bytes.data() + half, component, half)});
    }
    case DType::raw:
      return detail_v3::hex_bit_string(bytes, /*reverse_for_endianness=*/false);
    default:
      // Integers reuse the version-independent emission (plain JSON numbers).
      return detail::fill_to_json(bytes, dt);
  }
}

namespace detail_v3 {

inline json emit_codec_list(const std::vector<CodecSpec>& codecs) {
  json out = json::array();
  for (const CodecSpec& codec : codecs) {
    if (codec.name == "shuffle") {
      throw error("the v2 shuffle filter cannot be represented in v3 metadata");
    }
    json c = {{"name", codec.name}};
    if (codec.configuration.is_object() && !codec.configuration.empty()) {
      c["configuration"] = codec.configuration;
    }
    out.push_back(std::move(c));
  }
  return out;
}

}  // namespace detail_v3

/// Emits canonical v3 array metadata. Deterministic: fixed member set (empty
/// attributes and absent dimension_names are omitted), sorted keys, stable
/// forms. Shard levels fold back into nested sharding_indexed codecs.
inline json emit_array_meta(const ArrayMeta& meta) {
  json j;
  j["zarr_format"] = 3;
  j["node_type"] = "array";
  j["shape"] = meta.shape;
  j["data_type"] = emit_data_type(meta.dtype);
  const std::vector<std::uint64_t>& grid_shape =
      meta.shard_levels.empty() ? meta.chunk_shape : meta.shard_levels.front().shard_shape;
  j["chunk_grid"] = {{"name", "regular"}, {"configuration", {{"chunk_shape", grid_shape}}}};
  j["chunk_key_encoding"] = {
      {"name", meta.key_encoding == ChunkKeyKind::v3_default ? "default" : "v2"},
      {"configuration", {{"separator", std::string(1, meta.dimension_separator)}}}};
  j["fill_value"] = emit_fill(meta.fill, meta.dtype);

  json codecs = detail_v3::emit_codec_list(meta.codecs);
  for (std::size_t i = meta.shard_levels.size(); i-- > 0;) {
    const ShardLevel& level = meta.shard_levels[i];
    const std::vector<std::uint64_t>& inner_shape =
        i + 1 < meta.shard_levels.size() ? meta.shard_levels[i + 1].shard_shape : meta.chunk_shape;
    codecs = json::array({{{"name", "sharding_indexed"},
                           {"configuration",
                            {{"chunk_shape", inner_shape},
                             {"codecs", std::move(codecs)},
                             {"index_codecs", detail_v3::emit_codec_list(level.index_codecs)},
                             {"index_location", level.index_at_end ? "end" : "start"}}}}});
  }
  j["codecs"] = std::move(codecs);
  if (meta.attributes.is_object() && !meta.attributes.empty()) {
    j["attributes"] = meta.attributes;
  }
  if (meta.dimension_names.is_array()) {
    j["dimension_names"] = meta.dimension_names;
  }
  return j;
}

/// Emits canonical v3 group metadata.
inline json emit_group_meta(const json& attributes) {
  json j;
  j["zarr_format"] = 3;
  j["node_type"] = "group";
  if (attributes.is_object() && !attributes.empty()) {
    j["attributes"] = attributes;
  }
  return j;
}

/// Builds (or rebuilds) the inline consolidated-metadata member of the root
/// zarr.json from every v3 document in the store. Opt-in and explicit: the
/// convention (zarr-specs #309) is not yet an accepted spec, so libzarr never
/// writes it unasked.
inline void consolidate(Store& store) {
  const auto root_bytes = store.read(kMetaKey);
  if (!root_bytes) {
    throw error("v3::consolidate: no zarr.json at the store root");
  }
  json root = v2::parse_json(*root_bytes, kMetaKey);
  json metadata = json::object();
  for (const std::string& key : store.list_prefix("")) {
    if (key == kMetaKey || !detail::ends_with(key, std::string("/") + kMetaKey)) {
      continue;
    }
    const std::string path = key.substr(0, key.size() - std::string(kMetaKey).size() - 1);
    if (const auto bytes = store.read(key)) {
      metadata[path] = v2::parse_json(*bytes, key);
    }
  }
  root["consolidated_metadata"] = {
      {"kind", "inline"}, {"must_understand", false}, {"metadata", std::move(metadata)}};
  store.write(kMetaKey, canonical_json_bytes(root));
}

}  // namespace zarr::v3

#endif  // LIBZARR_V3_HPP
