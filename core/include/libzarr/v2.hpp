// SPDX-License-Identifier: MIT

#ifndef LIBZARR_V2_HPP
#define LIBZARR_V2_HPP

#include <cstdint>
#include <cstdlib>
#include <limits>
#include <optional>
#include <string>
#include <vector>

#include "libzarr/detail/common.hpp"
#include "libzarr/metadata.hpp"
#include "libzarr/store.hpp"
#include "libzarr/types.hpp"

/// \file v2.hpp
/// Zarr v2 metadata: parsing (tolerant) and emission (canonical), chunk keys,
/// and consolidated metadata (.zmetadata). Parsing lowers v2 forms into the
/// normalized ArrayMeta: compressor -> a bytes->bytes codec, order:"F" -> a
/// transpose codec, dtype byte order -> the "bytes" codec's endian.

namespace zarr::v2 {

/// v2 array metadata document name.
inline constexpr const char* kArraySuffix = ".zarray";
/// v2 group metadata document name.
inline constexpr const char* kGroupSuffix = ".zgroup";
/// v2 attributes document name.
inline constexpr const char* kAttrsSuffix = ".zattrs";
/// v2 consolidated-metadata document (store root).
inline constexpr const char* kConsolidatedKey = ".zmetadata";

/// Store key of a metadata document for the node at `path` ("" = root).
[[nodiscard]] inline std::string meta_key(const std::string& path, const char* suffix) {
  return path.empty() ? suffix : path + "/" + suffix;
}

// ---- dtype ------------------------------------------------------------------

/// A parsed v2 dtype string: the element type plus its stored byte order.
struct ParsedDType {
  /// Element type.
  DataType dtype;
  /// Stored byte order (v2 '>' dtypes).
  bool big_endian = false;
};

namespace detail_v2 {

template <typename Fail>
DataType int_dtype(bool is_signed, std::uint64_t size, const Fail& fail) {
  DType kind{};
  if (size == 1) {
    kind = is_signed ? DType::int8 : DType::uint8;
  } else if (size == 2) {
    kind = is_signed ? DType::int16 : DType::uint16;
  } else if (size == 4) {
    kind = is_signed ? DType::int32 : DType::uint32;
  } else if (size == 8) {
    kind = is_signed ? DType::int64 : DType::uint64;
  } else {
    fail("integer size must be 1, 2, 4 or 8");
  }
  return DataType::of(kind);
}

/// Maps a validated (code, size) pair to a DataType; `fail` reports errors
/// with the full dtype context.
template <typename Fail>
DataType dtype_of_code(char code, std::uint64_t size, const Fail& fail) {
  switch (code) {
    case 'b':
      if (size != 1) {
        fail("bool must have size 1");
      }
      return DataType::of(DType::boolean);
    case 'i':
    case 'u':
      return int_dtype(code == 'i', size, fail);
    case 'f':
      if (size == 2) {
        return DataType::of(DType::float16);
      }
      if (size == 4) {
        return DataType::of(DType::float32);
      }
      if (size == 8) {
        return DataType::of(DType::float64);
      }
      fail("float size must be 2, 4 or 8");
      break;
    case 'V':
      if (size == 0) {
        fail("raw dtype must have a positive size");
      }
      return DataType::raw_bytes(static_cast<std::uint32_t>(size));
    case 'c':
      if (size == 8) {
        return DataType::of(DType::complex64);
      }
      if (size == 16) {
        return DataType::of(DType::complex128);
      }
      fail("complex size must be 8 or 16");
      break;
    case 'S':
    case 'U':
      fail("string dtypes are not supported");
      break;
    case 'm':
    case 'M':
      fail("datetime dtypes are not supported");
      break;
    default:
      fail("unknown type code");
  }
  return {};  // unreachable: fail() always throws
}

}  // namespace detail_v2

/// Parses a v2 dtype string (`[<>|][biufV]<size>`). Read-supported kinds only;
/// everything else fails with a precise error naming the dtype.
[[nodiscard]] inline ParsedDType parse_data_type(const std::string& text, const std::string& ctx) {
  const auto fail = [&](const char* why) { throw error(ctx + ": dtype '" + text + "': " + why); };
  if (text.size() < 3) {
    fail("expected the form '<f8', '|u1', ...");
  }
  const char order = text[0];
  if (order != '<' && order != '>' && order != '|') {
    fail("byte order must be '<', '>' or '|'");
  }
  std::uint64_t size = 0;
  for (std::size_t i = 2; i < text.size(); ++i) {
    if (text[i] < '0' || text[i] > '9' || size > 0xFFFFFF) {
      fail("malformed item size");
    }
    size = size * 10 + static_cast<std::uint64_t>(text[i] - '0');
  }
  const DataType dtype = detail_v2::dtype_of_code(text[1], size, fail);
  // Byte order only matters for multi-byte components; bool/raw never swap.
  const bool multi_byte = dtype.kind != DType::boolean && dtype.kind != DType::raw && size > 1;
  return {dtype, order == '>' && multi_byte};
}

/// Emits the canonical v2 dtype string. `big_endian` is preserved from an
/// opened array (we never *create* big-endian arrays, but re-emitting an
/// opened one must not lie about its chunk bytes).
[[nodiscard]] inline std::string emit_data_type(DataType dt, bool big_endian) {
  const char order_multi = big_endian ? '>' : '<';
  switch (dt.kind) {
    case DType::boolean:
      return "|b1";
    case DType::int8:
      return "|i1";
    case DType::uint8:
      return "|u1";
    case DType::int16:
      return std::string(1, order_multi) + "i2";
    case DType::int32:
      return std::string(1, order_multi) + "i4";
    case DType::int64:
      return std::string(1, order_multi) + "i8";
    case DType::uint16:
      return std::string(1, order_multi) + "u2";
    case DType::uint32:
      return std::string(1, order_multi) + "u4";
    case DType::uint64:
      return std::string(1, order_multi) + "u8";
    case DType::float16:
      return std::string(1, order_multi) + "f2";
    case DType::float32:
      return std::string(1, order_multi) + "f4";
    case DType::float64:
      return std::string(1, order_multi) + "f8";
    case DType::complex64:
      return std::string(1, order_multi) + "c8";
    case DType::complex128:
      return std::string(1, order_multi) + "c16";
    case DType::raw:
      return "|V" + std::to_string(dt.itemsize);
    default:
      throw error("v2 emission not implemented for this dtype");
  }
}

// ---- fill_value ---------------------------------------------------------------

/// Parses a v2 fill_value, dtype-directed. Tolerant per our READ policy; each
/// tolerance cites its origin.
[[nodiscard]] inline std::optional<Bytes> parse_fill(const json& v, DataType dt,
                                                     const std::string& ctx);

namespace detail_v2 {

inline std::optional<Bytes> non_finite_fill(const std::string& s, DType kind) {
  // v2 spec: "NaN", "Infinity" and "-Infinity" are the sanctioned string
  // encodings of non-finite floats. "+Infinity" appears in the wild (it is
  // the v3 spelling); accept it on read.
  if (s == "NaN") {
    return detail::quiet_nan_bytes(kind);
  }
  if (s == "Infinity" || s == "+Infinity") {
    return detail::infinity_bytes(kind, false);
  }
  if (s == "-Infinity") {
    return detail::infinity_bytes(kind, true);
  }
  return std::nullopt;
}

/// One complex component: a number or a non-finite string (zarr-python's v2
/// complex fill form mirrors v3's [re, im] pairs).
inline Bytes complex_component_fill(const json& v, DType kind, const std::string& ctx) {
  if (v.is_number()) {
    return detail::fill_from_double(v.get<double>(), DataType::of(kind), ctx);
  }
  if (v.is_string()) {
    if (auto fill = non_finite_fill(v.get<std::string>(), kind)) {
      return *std::move(fill);
    }
  }
  throw error(ctx + ": cannot interpret complex fill_value component " + v.dump());
}

/// GDAL's Zarr driver emits numeric fills as JSON strings (read tolerance);
/// returns nullopt when `s` is not fully numeric.
inline std::optional<Bytes> numeric_string_fill(const std::string& s, DataType dt,
                                                const std::string& ctx) {
  if (s.empty()) {
    return std::nullopt;
  }
  const char* begin = s.c_str();
  char* end = nullptr;
  if (is_float(dt.kind)) {
    const double d = std::strtod(begin, &end);
    if (end == begin + s.size()) {
      return detail::fill_from_double(d, dt, ctx);
    }
  } else if (is_signed_int(dt.kind) || dt.kind == DType::boolean) {
    const long long i = std::strtoll(begin, &end, 10);
    if (end == begin + s.size()) {
      return detail::fill_from_int(static_cast<std::int64_t>(i), dt, ctx);
    }
  } else if (is_unsigned_int(dt.kind) && s[0] != '-') {
    const unsigned long long u = std::strtoull(begin, &end, 10);
    if (end == begin + s.size()) {
      return detail::fill_from_uint(static_cast<std::uint64_t>(u), dt, ctx);
    }
  }
  return std::nullopt;
}

inline Bytes string_fill(const std::string& s, DataType dt, const std::string& ctx) {
  if (is_float(dt.kind)) {
    if (auto fill = non_finite_fill(s, dt.kind)) {
      return *std::move(fill);
    }
  }
  if (dt.kind == DType::raw) {
    // v2 spec: raw ("V") fill values are base64-encoded.
    Bytes decoded = detail::base64_decode(s, (ctx + ": fill_value").c_str());
    if (decoded.size() != dt.itemsize) {
      throw error(ctx + ": base64 fill_value decodes to " + std::to_string(decoded.size()) +
                  " bytes, dtype needs " + std::to_string(dt.itemsize));
    }
    return decoded;
  }
  if (auto fill = numeric_string_fill(s, dt, ctx)) {
    return *std::move(fill);
  }
  throw error(ctx + ": cannot interpret fill_value '" + s + "' for this dtype");
}

}  // namespace detail_v2

[[nodiscard]] inline std::optional<Bytes> parse_fill(const json& v, DataType dt,
                                                     const std::string& ctx) {
  if (v.is_null()) {
    return std::nullopt;  // v2 spec: null = fill value undefined (reads as zeros)
  }
  if (v.is_array()) {
    if (is_complex(dt.kind) && v.size() == 2) {
      // zarr-python encodes v2 complex fills as [re, im], like v3.
      const DType component = dt.kind == DType::complex64 ? DType::float32 : DType::float64;
      Bytes out = detail_v2::complex_component_fill(v[0], component, ctx);
      const Bytes imag = detail_v2::complex_component_fill(v[1], component, ctx);
      out.insert(out.end(), imag.begin(), imag.end());
      return out;
    }
    // NCZarr 4.8.0 wraps fill_value in a 1-element array (read tolerance).
    if (v.size() == 1) {
      return parse_fill(v[0], dt, ctx);
    }
    throw error(ctx + ": fill_value must be a scalar, got an array of " + std::to_string(v.size()));
  }
  if (v.is_boolean()) {
    if (dt.kind != DType::boolean) {
      throw error(ctx + ": boolean fill_value for non-bool dtype");
    }
    return detail::scalar_bytes<std::uint8_t>(v.get<bool>() ? 1 : 0);
  }
  if (v.is_string()) {
    return detail_v2::string_fill(v.get<std::string>(), dt, ctx);
  }
  if (v.is_number_unsigned()) {
    return detail::fill_from_uint(v.get<std::uint64_t>(), dt, ctx);
  }
  if (v.is_number_integer()) {
    return detail::fill_from_int(v.get<std::int64_t>(), dt, ctx);
  }
  if (v.is_number_float()) {
    return detail::fill_from_double(v.get<double>(), dt, ctx);
  }
  throw error(ctx + ": unsupported fill_value " + v.dump());
}

// ---- .zarray ----------------------------------------------------------------

/// Parses a .zarray document into normalized ArrayMeta (without attributes,
/// which live in .zattrs). Unknown members are ignored: v2 predates v3's
/// must-understand rule and extra keys are common in the wild.
namespace detail_v2 {

inline char parse_separator(const json& j, const std::string& ctx) {
  const auto it = j.find("dimension_separator");
  if (it == j.end()) {
    return '.';
  }
  if (!it->is_string() || (it->get<std::string>() != "." && it->get<std::string>() != "/")) {
    throw error(ctx + R"(: 'dimension_separator' must be "." or "/")");
  }
  return it->get<std::string>()[0];
}

/// Lowers the v2 filters member into bytes->bytes codec specs (applied
/// before the compressor, per numcodecs). Only the shuffle filter is
/// supported: netCDF's NCZarr writes it by default alongside zlib.
inline std::vector<CodecSpec> parse_filters(const json& j, const std::string& ctx) {
  std::vector<CodecSpec> out;
  const auto it = j.find("filters");
  if (it == j.end() || it->is_null()) {
    return out;
  }
  // Tolerance: no-filters is canonically null, but [] appears in the wild.
  if (!it->is_array()) {
    throw error(ctx + ": 'filters' must be null or an array");
  }
  for (const json& f : *it) {
    if (!f.is_object() || !f.contains("id") || !f["id"].is_string()) {
      throw error(ctx + ": each filter must be an object with an 'id'");
    }
    const auto id = f["id"].get<std::string>();
    if (id == "shuffle") {
      // NCZarr 4.9.x writes elementsize "0" (a string, and zero) meaning
      // "the dtype's item size"; the resolved size is filled in at codec
      // resolution.
      const std::int64_t elementsize = detail::lenient_int(f, "elementsize", 0, ctx);
      out.push_back(CodecSpec{"shuffle", {{"elementsize", elementsize}}});
    } else {
      std::string msg = ctx;
      msg += ": unsupported v2 filter '";
      msg += id;
      msg += "'";
      throw error(msg);
    }
  }
  return out;
}

/// Lowers the v2 compressor member into 0 or 1 bytes->bytes codec specs.
inline std::optional<CodecSpec> parse_compressor(const json& j, const std::string& ctx) {
  const auto it = j.find("compressor");
  if (it == j.end() || it->is_null()) {
    return std::nullopt;
  }
  if (!it->is_object() || !it->contains("id") || !(*it)["id"].is_string()) {
    throw error(ctx + ": 'compressor' must be null or an object with an 'id'");
  }
  const auto id = (*it)["id"].get<std::string>();
  if (id == "zlib" || id == "gzip") {
    // numcodecs defaults level to 1 when absent.
    const std::int64_t level = detail::lenient_int(*it, "level", 1, ctx);
    if (level < 0 || level > 9) {
      throw error(ctx + ": compressor level must be in 0..9");
    }
    return CodecSpec{id, {{"level", level}}};
  }
  if (id == "blosc") {
    // numcodecs Blosc: numeric shuffle (-1 = auto), no typesize member (the
    // dtype's itemsize applies). Values are validated at codec resolution.
    // Evaluated before the braced list: a .value() throw during json
    // initializer-list construction leaks json_ref temporaries (fuzz+LSan).
    const std::string cname = it->value("cname", "lz4");
    const std::int64_t clevel = detail::lenient_int(*it, "clevel", 5, ctx);
    const json shuffle = it->value("shuffle", json(1));
    const std::int64_t blocksize = detail::lenient_int(*it, "blocksize", 0, ctx);
    return CodecSpec{
        "blosc",
        {{"cname", cname}, {"clevel", clevel}, {"shuffle", shuffle}, {"blocksize", blocksize}}};
  }
  if (id == "zstd") {
    // numcodecs Zstd (zarr-python 3's default for v2-format arrays).
    const std::int64_t level = detail::lenient_int(*it, "level", 0, ctx);
    return CodecSpec{"zstd", {{"level", level}}};
  }
  throw error(ctx + ": unsupported v2 compressor '" + id + "'");
}

}  // namespace detail_v2

namespace detail_v2 {

inline ArrayMeta parse_array_meta_impl(const json& j, const std::string& ctx) {
  if (!j.is_object()) {
    throw error(ctx + ": expected a JSON object");
  }
  const auto require = [&](const char* name) -> const json& {
    const auto it = j.find(name);
    if (it == j.end()) {
      throw error(ctx + ": missing required member '" + name + "'");
    }
    return *it;
  };
  if (detail::json_to_uint64(require("zarr_format"), ctx + ": zarr_format") != 2) {
    throw error(ctx + ": zarr_format must be 2");
  }

  ArrayMeta meta;
  meta.format = ZarrFormat::v2;
  meta.shape = detail::parse_extents(require("shape"), "shape", ctx);
  meta.chunk_shape = detail::parse_extents(require("chunks"), "chunks", ctx);
  if (meta.chunk_shape.size() != meta.shape.size()) {
    throw error(ctx + ": 'chunks' must be an array of the same rank as 'shape'");
  }
  for (const std::uint64_t c : meta.chunk_shape) {
    if (c == 0) {
      throw error(ctx + ": chunk extents must be positive");
    }
  }

  const json& dtype = require("dtype");
  if (!dtype.is_string()) {
    throw error(ctx + ": structured dtypes are not supported");
  }
  const ParsedDType parsed = parse_data_type(dtype.get<std::string>(), ctx);
  meta.dtype = parsed.dtype;

  const json& order = require("order");
  if (!order.is_string() || (order != "C" && order != "F")) {
    throw error(ctx + R"(: 'order' must be "C" or "F")");
  }

  // Tolerance: fill_value/compressor are required members, but minimal
  // writers omit them; missing reads as null.
  const auto fill_it = j.find("fill_value");
  meta.fill = fill_it == j.end() ? std::nullopt : parse_fill(*fill_it, meta.dtype, ctx);

  meta.dimension_separator = parse_separator(j, ctx);

  // Lowering into the normalized codec chain: array->array, bytes, then
  // filters (numcodecs applies them before the compressor), then compressor.
  if (order == "F" && meta.shape.size() >= 2) {
    json perm = json::array();
    for (std::size_t d = meta.shape.size(); d-- > 0;) {
      perm.push_back(d);
    }
    meta.codecs.push_back({"transpose", {{"order", perm}}});
  }
  meta.codecs.push_back({"bytes", {{"endian", parsed.big_endian ? "big" : "little"}}});
  for (CodecSpec& filter : parse_filters(j, ctx)) {
    meta.codecs.push_back(std::move(filter));
  }
  if (auto compressor = parse_compressor(j, ctx)) {
    meta.codecs.push_back(*std::move(compressor));
  }
  return meta;
}

/// The canonical v2 filters member for a lowered codec chain.
inline json emit_filters(const std::vector<CodecSpec>& codecs) {
  json filters = json::array();
  for (const CodecSpec& codec : codecs) {
    if (codec.name == "shuffle") {
      filters.push_back(
          {{"id", "shuffle"},
           {"elementsize", codec.configuration.value("elementsize", std::int64_t{0})}});
    }
  }
  return filters.empty() ? json(nullptr) : filters;
}

}  // namespace detail_v2

/// Parses a .zarray document into normalized ArrayMeta (without attributes,
/// which live in .zattrs). Unknown members are ignored: v2 predates v3's
/// must-understand rule and extra keys are common in the wild.
[[nodiscard]] inline ArrayMeta parse_array_meta(const json& j, const std::string& ctx) {
  return detail::guard_json(ctx, [&] { return detail_v2::parse_array_meta_impl(j, ctx); });
}

/// Emits canonical .zarray JSON. Deterministic: sorted keys, stable forms;
/// dimension_separator appears only when '/' (matching common practice).
[[nodiscard]] inline json emit_array_meta(const ArrayMeta& meta) {
  json j;
  j["zarr_format"] = 2;
  j["shape"] = meta.shape;
  j["chunks"] = meta.chunk_shape;
  j["filters"] = detail_v2::emit_filters(meta.codecs);
  j["fill_value"] = detail::fill_to_json(meta.fill, meta.dtype);
  if (meta.dimension_separator == '/') {
    j["dimension_separator"] = "/";
  }

  bool big_endian = false;
  bool f_order = false;
  j["compressor"] = nullptr;
  for (const CodecSpec& codec : meta.codecs) {
    if (codec.name == "bytes") {
      big_endian = codec.configuration.value("endian", "little") == std::string("big");
    } else if (codec.name == "transpose") {
      // v2 can only express the full reversal (order:"F").
      const json& perm = codec.configuration.at("order");
      for (std::size_t i = 0; i < perm.size(); ++i) {
        if (perm[i].get<std::uint64_t>() != perm.size() - 1 - i) {
          throw error(R"(v2 cannot represent this transpose; only order:"F" (full reversal))");
        }
      }
      f_order = true;
    } else if (codec.name == "gzip" || codec.name == "zlib") {
      j["compressor"] = {{"id", codec.name},
                         {"level", codec.configuration.value("level", std::int64_t{1})}};
    } else if (codec.name == "blosc") {
      // Canonical numcodecs form: numeric shuffle.
      const json shuffle = codec.configuration.value("shuffle", json(1));
      std::int64_t shuffle_num = 1;
      if (shuffle.is_number_integer()) {
        shuffle_num = shuffle.get<std::int64_t>();
      } else if (shuffle == "noshuffle") {
        shuffle_num = 0;
      } else if (shuffle == "bitshuffle") {
        shuffle_num = 2;
      }
      j["compressor"] = {{"id", "blosc"},
                         {"cname", codec.configuration.value("cname", "lz4")},
                         {"clevel", codec.configuration.value("clevel", std::int64_t{5})},
                         {"shuffle", shuffle_num},
                         {"blocksize", codec.configuration.value("blocksize", std::int64_t{0})}};
    } else if (codec.name == "zstd") {
      j["compressor"] = {{"id", "zstd"},
                         {"level", codec.configuration.value("level", std::int64_t{0})}};
    } else if (codec.name == "shuffle") {
      // already emitted into the filters member
    } else {
      throw error("v2 cannot represent codec '" + codec.name + "'");
    }
  }
  j["dtype"] = emit_data_type(meta.dtype, big_endian);
  j["order"] = f_order ? "F" : "C";
  return j;
}

// ---- groups / chunk keys ------------------------------------------------------

/// The (only) content of a v2 group document.
[[nodiscard]] inline json emit_group_meta() { return json{{"zarr_format", 2}}; }

/// Validates a .zgroup document.
inline void check_group_meta(const json& j, const std::string& ctx) {
  detail::guard_json(ctx, [&] {
    if (!j.is_object() || j.find("zarr_format") == j.end() ||
        detail::json_to_uint64(j.at("zarr_format"), ctx) != 2) {
      throw error(ctx + ": not a v2 group (zarr_format must be 2)");
    }
  });
}

/// v2 chunk key relative to the array: indices joined by the separator;
/// rank-0 arrays use the fixed key "0".
[[nodiscard]] inline std::string chunk_key(const std::vector<std::uint64_t>& index,
                                           char separator) {
  if (index.empty()) {
    return "0";
  }
  std::string key = std::to_string(index[0]);
  for (std::size_t d = 1; d < index.size(); ++d) {
    key += separator;
    key += std::to_string(index[d]);
  }
  return key;
}

// ---- consolidated metadata (.zmetadata) ---------------------------------------

/// Writes a metadata document and, when a consolidated .zmetadata exists at
/// the store root, keeps it in sync (decision: consolidated metadata is
/// maintained by default so it can never go stale through libzarr writes).
inline void write_meta_key(Store& store, const std::string& key, const json& value) {
  store.write(key, canonical_json_bytes(value));
  if (auto existing = store.read(kConsolidatedKey)) {
    json c = detail::parse_json(*existing, kConsolidatedKey);
    c["metadata"][key] = value;
    store.write(kConsolidatedKey, canonical_json_bytes(c));
  }
}

/// Removes a metadata document, keeping .zmetadata in sync (see write_meta_key).
inline void erase_meta_key(Store& store, const std::string& key) {
  store.erase(key);
  if (auto existing = store.read(kConsolidatedKey)) {
    json c = detail::parse_json(*existing, kConsolidatedKey);
    auto meta_it = c.find("metadata");
    if (meta_it != c.end()) {
      meta_it->erase(key);
    }
    store.write(kConsolidatedKey, canonical_json_bytes(c));
  }
}

/// Builds (or rebuilds) .zmetadata from every v2 metadata document in the
/// store, in zarr-python's consolidated format 1.
inline void consolidate(Store& store) {
  json metadata = json::object();
  for (const std::string& key : store.list_prefix("")) {
    const std::string_view k = key;
    const auto leaf_is = [&](const char* suffix) {
      return k == suffix || detail::ends_with(k, std::string("/") + suffix);
    };
    if (leaf_is(kArraySuffix) || leaf_is(kGroupSuffix) || leaf_is(kAttrsSuffix)) {
      const auto bytes = store.read(key);
      if (bytes) {
        metadata[key] = detail::parse_json(*bytes, key);
      }
    }
  }
  const json c = {{"metadata", metadata}, {"zarr_consolidated_format", 1}};
  store.write(kConsolidatedKey, canonical_json_bytes(c));
}

/// Loads the consolidated metadata map if present and well-formed.
[[nodiscard]] inline std::optional<json> read_consolidated(Store& store) {
  const auto bytes = store.read(kConsolidatedKey);
  if (!bytes) {
    return std::nullopt;
  }
  json c = detail::parse_json(*bytes, kConsolidatedKey);
  if (!c.is_object() || c.value("zarr_consolidated_format", std::int64_t{0}) != 1 ||
      !c.contains("metadata") || !c["metadata"].is_object()) {
    throw error(std::string(kConsolidatedKey) + ": unrecognized consolidated metadata format");
  }
  return c["metadata"];
}

}  // namespace zarr::v2

#endif  // LIBZARR_V2_HPP
