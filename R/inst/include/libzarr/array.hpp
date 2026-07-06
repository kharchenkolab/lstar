// SPDX-License-Identifier: MIT

#ifndef LIBZARR_ARRAY_HPP
#define LIBZARR_ARRAY_HPP

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <memory>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include "libzarr/codecs.hpp"
#include "libzarr/detail/common.hpp"
#include "libzarr/metadata.hpp"
#include "libzarr/sharding.hpp"
#include "libzarr/store.hpp"
#include "libzarr/types.hpp"
#include "libzarr/v2.hpp"
#include "libzarr/v3.hpp"

/// \file array.hpp
/// The Array API: create/open, whole-array and per-chunk read/write, and
/// byte-range sub-chunk reads. Buffers in, buffers out; all data is native
/// byte order, C layout.

namespace zarr {

namespace detail {

/// Validates a node path: "" (root) or '/'-separated non-empty segments that
/// do not start with '.' (which would collide with metadata documents).
inline void validate_path(const std::string& path) {
  if (path.empty()) {
    return;
  }
  if (path.front() == '/' || path.back() == '/') {
    throw error("node path must not start or end with '/': '" + path + "'");
  }
  std::size_t start = 0;
  while (start <= path.size()) {
    const std::size_t slash = path.find('/', start);
    const std::size_t end = slash == std::string::npos ? path.size() : slash;
    if (end == start) {
      throw error("node path has an empty segment: '" + path + "'");
    }
    if (path[start] == '.') {
      throw error("node path segments must not start with '.': '" + path + "'");
    }
    if (slash == std::string::npos) {
      break;
    }
    start = slash + 1;
  }
}

/// Advances a C-order odometer over `extents`; returns false after the last
/// index. Rank 0 iterates exactly once.
inline bool next_index(std::vector<std::uint64_t>& index,
                       const std::vector<std::uint64_t>& extents) {
  std::size_t d = extents.size();
  while (d-- > 0) {
    if (++index[d] < extents[d]) {
      return true;
    }
    index[d] = 0;
  }
  return false;
}

}  // namespace detail

/// Parameters for Array::create.
struct ArraySpec {
  /// Storage format version to write.
  ZarrFormat format = ZarrFormat::v2;
  /// Array shape; empty = 0-dimensional.
  std::vector<std::uint64_t> shape;
  /// Chunk shape, same rank as `shape`, extents >= 1.
  std::vector<std::uint64_t> chunks;
  /// Element type.
  DataType dtype;
  /// bytes->bytes codecs (v2: at most one of zarr::gzip / zarr::zlib;
  /// v3 additionally zarr::CodecSpec{"blosc", ...} / {"crc32c", {}}).
  std::vector<CodecSpec> codecs;
  /// Fill value as one native-order element; defaults to zeros.
  std::optional<Bytes> fill;
  /// Initial user attributes.
  json attributes = json::object();
  /// v2 chunk-key separator; '.' is canonical, '/' supported. (v3 arrays are
  /// created with the "default" encoding and '/' regardless.)
  char dimension_separator = '.';
  /// v3 dimension_names (array of strings/null, rank length), or null.
  json dimension_names;
  /// v3 sharding: shard (outer chunk) shape; each extent must be a multiple
  /// of the corresponding `chunks` extent. Empty = unsharded. `codecs` apply
  /// to the inner chunks; the index gets `bytes` + `crc32c`.
  std::vector<std::uint64_t> shards;
};

/// A Zarr array bound to a Store. Value-semantics handle: cheap to move,
/// holds a shared reference to the store.
class Array {
 public:
  /// Creates an array at `path` ("" = the store root) in spec.format,
  /// writing canonical metadata. Fails if the spec is invalid; overwrites
  /// existing metadata.
  static Array create(std::shared_ptr<Store> store, const std::string& path,
                      const ArraySpec& spec) {
    if (!store) {
      throw error("Array::create: null store");
    }
    detail::validate_path(path);
    const bool v3 = spec.format == ZarrFormat::v3;
    const std::string ctx = v3 ? v3::meta_key(path) : v2::meta_key(path, v2::kArraySuffix);
    if (spec.chunks.size() != spec.shape.size()) {
      throw error(ctx + ": chunks rank " + std::to_string(spec.chunks.size()) + " != shape rank " +
                  std::to_string(spec.shape.size()));
    }
    for (const std::uint64_t c : spec.chunks) {
      if (c == 0) {
        throw error(ctx + ": chunk extents must be positive");
      }
    }
    if (spec.dimension_separator != '.' && spec.dimension_separator != '/') {
      throw error(ctx + ": dimension_separator must be '.' or '/'");
    }

    ArrayMeta meta;
    meta.format = spec.format;
    meta.shape = spec.shape;
    meta.chunk_shape = spec.chunks;
    meta.dtype = spec.dtype;
    meta.attributes = spec.attributes;
    apply_format_members(spec, meta, ctx);
    if (spec.fill) {
      if (spec.fill->size() != spec.dtype.itemsize) {
        throw error(ctx + ": fill is " + std::to_string(spec.fill->size()) +
                    " bytes, dtype needs " + std::to_string(spec.dtype.itemsize));
      }
      meta.fill = spec.fill;
    } else {
      meta.fill = Bytes(spec.dtype.itemsize, 0);  // canonical default: zeros
    }
    meta.codecs.push_back({"bytes", {{"endian", "little"}}});
    for (const CodecSpec& codec : spec.codecs) {
      meta.codecs.push_back(codec);
    }

    Array array(std::move(store), path, std::move(meta));  // resolves + validates codecs
    if (v3) {
      array.store_->write(v3::meta_key(path),
                          canonical_json_bytes(v3::emit_array_meta(array.meta_)));
    } else {
      v2::write_meta_key(*array.store_, array.meta_store_key(), v2::emit_array_meta(array.meta_));
      array.write_attributes();
    }
    return array;
  }

  /// Opens an existing array at `path`, probing v3 (zarr.json) first, then
  /// v2 (.zarray). `consolidated` is a pre-fetched store-key -> document map
  /// (supplied by Group when consolidated metadata is present).
  static Array open(std::shared_ptr<Store> store, const std::string& path, OpenOptions options = {},
                    const std::shared_ptr<const json>& consolidated = nullptr) {
    if (!store) {
      throw error("Array::open: null store");
    }
    detail::validate_path(path);
    const auto read_doc = [&](const std::string& key) -> std::optional<json> {
      if (consolidated) {
        const auto it = consolidated->find(key);
        if (it == consolidated->end()) {
          return std::nullopt;
        }
        return *it;
      }
      const auto bytes = store->read(key);
      if (!bytes) {
        return std::nullopt;
      }
      return v2::parse_json(*bytes, key);
    };

    // Probe order: zarr.json first, so v3 opens cost one round-trip.
    const std::string v3_key = v3::meta_key(path);
    if (const auto doc = read_doc(v3_key)) {
      if (doc->is_object() && doc->value("node_type", "") == std::string("group")) {
        throw error("'" + path + "' is a group, not an array");
      }
      return {std::move(store), path, v3::parse_array_meta(*doc, v3_key, options.lenient)};
    }

    const std::string meta_key = v2::meta_key(path, v2::kArraySuffix);
    auto doc = read_doc(meta_key);
    if (!doc) {
      if (store->exists(v2::meta_key(path, v2::kGroupSuffix))) {
        throw error("'" + path + "' is a group, not an array");
      }
      throw error("no array at '" + path + "' (neither " + v3_key + " nor " + meta_key + " found)");
    }
    ArrayMeta meta = v2::parse_array_meta(*doc, meta_key);
    if (const auto attrs = read_doc(v2::meta_key(path, v2::kAttrsSuffix))) {
      meta.attributes = *attrs;
    }
    return {std::move(store), path, std::move(meta)};
  }

  /// Normalized metadata (shape, chunks, dtype, codecs, attributes).
  [[nodiscard]] const ArrayMeta& meta() const { return meta_; }
  /// Node path within the store ("" = root).
  [[nodiscard]] const std::string& path() const { return path_; }
  /// Chunk-grid extent per dimension.
  [[nodiscard]] std::vector<std::uint64_t> grid_shape() const { return meta_.grid_shape(); }
  /// Whole-array size in bytes (elements x itemsize).
  [[nodiscard]] std::uint64_t nbytes() const {
    return meta_.element_count() * meta_.dtype.itemsize;
  }
  /// One full chunk's size in bytes.
  [[nodiscard]] std::uint64_t chunk_nbytes() const { return pipeline_.decoded_chunk_bytes(); }

  /// Reads one chunk as a full native/C-order buffer. Missing chunks read as
  /// fill; edge chunks come back full-sized (fill-padded), per the format.
  [[nodiscard]] Bytes read_chunk(const std::vector<std::uint64_t>& index) const {
    auto stored = chunk_store_->read(chunk_store_key(index));
    if (!stored) {
      return filled_chunk();
    }
    return pipeline_.decode(std::move(*stored));
  }

  /// Writes one full chunk (native byte order, C layout, exactly
  /// chunk_nbytes() bytes — edge chunks included, fill-padded).
  void write_chunk(const std::vector<std::uint64_t>& index, const void* data, std::size_t size) {
    Bytes chunk(static_cast<const std::uint8_t*>(data),
                static_cast<const std::uint8_t*>(data) + size);
    chunk_store_->write(chunk_store_key(index), pipeline_.encode(std::move(chunk)));
    chunk_store_->flush();
  }

  /// Reads `element_count` consecutive elements of one chunk starting at
  /// linear element `element_offset` (in stored C order), using a store
  /// byte-range read — a single hosted object can be range-served. Requires
  /// an uncompressed, untransposed chunk layout
  /// (CodecPipeline::supports_partial_read).
  [[nodiscard]] Bytes read_chunk_range(const std::vector<std::uint64_t>& index,
                                       std::uint64_t element_offset,
                                       std::uint64_t element_count) const {
    if (!pipeline_.supports_partial_read()) {
      throw error("byte-range chunk reads need an uncompressed, untransposed layout");
    }
    const std::uint32_t itemsize = meta_.dtype.itemsize;
    const std::uint64_t chunk_elements = meta_.chunk_element_count();
    if (element_count > chunk_elements || element_offset > chunk_elements - element_count) {
      throw error("chunk range [" + std::to_string(element_offset) + ", +" +
                  std::to_string(element_count) + ") exceeds " + std::to_string(chunk_elements) +
                  " elements");
    }
    auto stored = chunk_store_->read_range(
        chunk_store_key(index),
        ByteRange::slice(element_offset * itemsize, element_count * itemsize));
    if (!stored) {
      Bytes out(detail::checked_size(element_count * itemsize, "chunk range"));
      detail::fill_elements(out.data(), element_count, meta_.fill ? meta_.fill->data() : nullptr,
                            itemsize);
      return out;
    }
    return pipeline_.decode_range(std::move(*stored));
  }

  /// Reads the whole array into `dst` (native order, C layout); `size` must
  /// equal nbytes().
  void read(void* dst, std::size_t size) const {
    read_region(std::vector<std::uint64_t>(meta_.shape.size(), 0), meta_.shape, dst, size);
  }

  /// Writes the whole array from `src` (native order, C layout); `size` must
  /// equal nbytes(). Edge chunks are fill-padded, per the format.
  void write(const void* src, std::size_t size) {
    write_region(std::vector<std::uint64_t>(meta_.shape.size(), 0), meta_.shape, src, size);
  }

  /// Reads the hyperslab of `shape` starting at `origin` (array coordinates)
  /// into `dst` as a C-order native buffer of exactly `size` bytes. Missing
  /// chunks read as fill.
  void read_region(const std::vector<std::uint64_t>& origin,
                   const std::vector<std::uint64_t>& shape, void* dst, std::size_t size) const {
    validate_region(origin, shape, size, "read_region");
    if (size == 0) {
      return;
    }
    auto* out = static_cast<std::uint8_t*>(dst);
    for_each_region_chunk(origin, shape, [&](const RegionChunk& rc) {
      const Bytes chunk = read_chunk(rc.index);
      detail::copy_box(chunk.data(), meta_.chunk_shape, rc.origin_in_chunk, out, shape,
                       rc.origin_in_region, rc.box, meta_.dtype.itemsize);
    });
  }

  /// Writes the hyperslab of `shape` starting at `origin` from `src` (a
  /// C-order native buffer of exactly `size` bytes). Partially covered
  /// chunks are read-modify-written, preserving their other elements;
  /// chunks the region covers entirely are rebuilt without a read.
  void write_region(const std::vector<std::uint64_t>& origin,
                    const std::vector<std::uint64_t>& shape, const void* src, std::size_t size) {
    validate_region(origin, shape, size, "write_region");
    if (size == 0) {
      return;
    }
    const auto* in = static_cast<const std::uint8_t*>(src);
    for_each_region_chunk(origin, shape, [&](const RegionChunk& rc) {
      Bytes chunk = rc.covered ? filled_chunk() : read_chunk(rc.index);
      detail::copy_box(in, shape, rc.origin_in_region, chunk.data(), meta_.chunk_shape,
                       rc.origin_in_chunk, rc.box, meta_.dtype.itemsize);
      chunk_store_->write(chunk_store_key(rc.index), pipeline_.encode(std::move(chunk)));
    });
    chunk_store_->flush();
  }

  /// User attributes (.zattrs).
  [[nodiscard]] const json& attributes() const { return meta_.attributes; }

  /// Replaces the user attributes and persists them. For v3, the stored
  /// zarr.json is patched in place, preserving any extension members.
  void set_attributes(json attributes) {
    meta_.attributes = std::move(attributes);
    if (meta_.format == ZarrFormat::v3) {
      const std::string key = v3::meta_key(path_);
      const auto bytes = store_->read(key);
      if (!bytes) {
        throw error(key + ": metadata disappeared");
      }
      json doc = v2::parse_json(*bytes, key);
      if (meta_.attributes.empty()) {
        doc.erase("attributes");
      } else {
        doc["attributes"] = meta_.attributes;
      }
      store_->write(key, canonical_json_bytes(doc));
      return;
    }
    write_attributes();
  }

  /// Store key of the chunk at `index` (bounds-checked against the grid).
  [[nodiscard]] std::string chunk_store_key(const std::vector<std::uint64_t>& index) const {
    const auto grid = meta_.grid_shape();
    if (index.size() != grid.size()) {
      throw error("chunk index rank " + std::to_string(index.size()) + " != array rank " +
                  std::to_string(grid.size()));
    }
    for (std::size_t d = 0; d < grid.size(); ++d) {
      if (index[d] >= grid[d]) {
        throw error("chunk index " + std::to_string(index[d]) + " out of range for dimension " +
                    std::to_string(d) + " (grid extent " + std::to_string(grid[d]) + ")");
      }
    }
    const std::string relative = meta_.key_encoding == ChunkKeyKind::v3_default
                                     ? v3::chunk_key(index, meta_.dimension_separator)
                                     : v2::chunk_key(index, meta_.dimension_separator);
    return path_.empty() ? relative : path_ + "/" + relative;
  }

 private:
  /// Applies the format-specific ArraySpec members (chunk-key encoding,
  /// dimension_names, shards) with their validation.
  static void apply_format_members(const ArraySpec& spec, ArrayMeta& meta, const std::string& ctx) {
    if (spec.format != ZarrFormat::v3) {
      meta.dimension_separator = spec.dimension_separator;
      if (!spec.dimension_names.is_null()) {
        throw error(ctx + ": dimension_names is a v3 feature");
      }
      if (!spec.shards.empty()) {
        throw error(ctx + ": sharding is a v3 feature");
      }
      return;
    }
    // Canonical v3 creation: the "default" chunk-key encoding with '/'.
    meta.key_encoding = ChunkKeyKind::v3_default;
    meta.dimension_separator = '/';
    if (spec.dimension_names.is_array()) {
      if (spec.dimension_names.size() != spec.shape.size()) {
        throw error(ctx + ": dimension_names must have rank length");
      }
      meta.dimension_names = spec.dimension_names;
    }
    if (spec.shards.empty()) {
      return;
    }
    if (spec.shards.size() != spec.chunks.size()) {
      throw error(ctx + ": shards rank must match chunks rank");
    }
    for (std::size_t d = 0; d < spec.shards.size(); ++d) {
      // v3 sharding spec: chunks must evenly divide the shard.
      if (spec.shards[d] == 0 || spec.shards[d] % spec.chunks[d] != 0) {
        throw error(ctx + ": each shard extent must be a positive multiple of the chunk extent");
      }
    }
    ShardLevel level;
    level.shard_shape = spec.shards;
    level.index_codecs = {{"bytes", {{"endian", "little"}}}, {"crc32c", {}}};
    meta.shard_levels.push_back(std::move(level));
  }

  Array(std::shared_ptr<Store> store, std::string path, ArrayMeta meta)
      : store_(std::move(store)),
        path_(std::move(path)),
        meta_(std::move(meta)),
        pipeline_(CodecPipeline::resolve(meta_)),
        chunk_store_(wrap_shards(store_, meta_, path_)) {
    // Materializing a chunk must be possible on this platform (wasm32!).
    detail::checked_size(pipeline_.decoded_chunk_bytes(), "chunk");
  }

  /// Builds the chunk-I/O store: the raw store for plain arrays, or a chain
  /// of ShardStore adapters (outermost level first) for sharded ones.
  /// Metadata I/O always uses the raw store.
  static std::shared_ptr<Store> wrap_shards(std::shared_ptr<Store> store, const ArrayMeta& meta,
                                            const std::string& path) {
    std::shared_ptr<Store> chunks = std::move(store);
    const std::string prefix = path.empty() ? "" : path + "/";
    for (std::size_t i = 0; i < meta.shard_levels.size(); ++i) {
      const ShardLevel& level = meta.shard_levels[i];
      const std::vector<std::uint64_t>& inner_shape = i + 1 < meta.shard_levels.size()
                                                          ? meta.shard_levels[i + 1].shard_shape
                                                          : meta.chunk_shape;
      ShardParams params;
      params.chunk_prefix = prefix;
      params.key_encoding = meta.key_encoding;
      params.separator = meta.dimension_separator;
      params.index_codecs = level.index_codecs;
      params.index_at_end = level.index_at_end;
      params.per_shard.resize(inner_shape.size());
      params.inner_grid.resize(inner_shape.size());
      for (std::size_t d = 0; d < inner_shape.size(); ++d) {
        params.per_shard[d] = level.shard_shape[d] / inner_shape[d];
        params.inner_grid[d] = detail::ceil_div(meta.shape[d], inner_shape[d]);
      }
      chunks = std::make_shared<ShardStore>(std::move(chunks), std::move(params));
    }
    return chunks;
  }

  [[nodiscard]] std::string meta_store_key() const { return v2::meta_key(path_, v2::kArraySuffix); }

  void write_attributes() {
    const std::string key = v2::meta_key(path_, v2::kAttrsSuffix);
    if (meta_.attributes.empty()) {
      v2::erase_meta_key(*store_, key);  // canonical: no empty .zattrs documents
    } else {
      v2::write_meta_key(*store_, key, meta_.attributes);
    }
  }

  [[nodiscard]] Bytes filled_chunk() const {
    Bytes chunk(detail::checked_size(pipeline_.decoded_chunk_bytes(), "chunk"));
    detail::fill_elements(chunk.data(), meta_.chunk_element_count(),
                          meta_.fill ? meta_.fill->data() : nullptr, meta_.dtype.itemsize);
    return chunk;
  }

  /// The intersection of one chunk with a requested region.
  struct RegionChunk {
    std::vector<std::uint64_t> index;             ///< chunk-grid index
    std::vector<std::uint64_t> origin_in_chunk;   ///< intersection start, chunk coords
    std::vector<std::uint64_t> origin_in_region;  ///< intersection start, region coords
    std::vector<std::uint64_t> box;               ///< intersection extents
    /// True when the region covers the chunk's whole in-array portion (so a
    /// write need not read the existing chunk first).
    bool covered = true;
  };

  void validate_region(const std::vector<std::uint64_t>& origin,
                       const std::vector<std::uint64_t>& shape, std::size_t size,
                       const char* what) const {
    const std::size_t rank = meta_.shape.size();
    if (origin.size() != rank || shape.size() != rank) {
      throw error(std::string(what) + ": origin/shape rank must be " + std::to_string(rank));
    }
    for (std::size_t d = 0; d < rank; ++d) {
      if (shape[d] > meta_.shape[d] || origin[d] > meta_.shape[d] - shape[d]) {
        throw error(std::string(what) + ": region [" + std::to_string(origin[d]) + ", " +
                    std::to_string(origin[d]) + "+" + std::to_string(shape[d]) +
                    ") exceeds dimension " + std::to_string(d) + " (extent " +
                    std::to_string(meta_.shape[d]) + ")");
      }
    }
    const std::uint64_t bytes = detail::checked_product(shape, what) * meta_.dtype.itemsize;
    if (size != detail::checked_size(bytes, what)) {
      throw error(std::string(what) + ": buffer is " + std::to_string(size) +
                  " bytes, region needs " + std::to_string(bytes));
    }
  }

  /// Invokes `fn(RegionChunk)` for every chunk intersecting the (non-empty,
  /// validated) region, in C order.
  template <typename Fn>
  void for_each_region_chunk(const std::vector<std::uint64_t>& origin,
                             const std::vector<std::uint64_t>& shape, const Fn& fn) const {
    const std::size_t rank = meta_.shape.size();
    std::vector<std::uint64_t> first(rank, 0);
    std::vector<std::uint64_t> last(rank, 0);
    for (std::size_t d = 0; d < rank; ++d) {
      first[d] = origin[d] / meta_.chunk_shape[d];
      last[d] = (origin[d] + shape[d] - 1) / meta_.chunk_shape[d];
    }

    RegionChunk rc;
    rc.index = first;
    rc.origin_in_chunk.assign(rank, 0);
    rc.origin_in_region.assign(rank, 0);
    rc.box.assign(rank, 0);
    while (true) {
      rc.covered = true;
      for (std::size_t d = 0; d < rank; ++d) {
        const std::uint64_t chunk_start = rc.index[d] * meta_.chunk_shape[d];
        const std::uint64_t valid_end =
            std::min(chunk_start + meta_.chunk_shape[d], meta_.shape[d]);
        const std::uint64_t begin = std::max(chunk_start, origin[d]);
        const std::uint64_t end = std::min(valid_end, origin[d] + shape[d]);
        rc.origin_in_chunk[d] = begin - chunk_start;
        rc.origin_in_region[d] = begin - origin[d];
        rc.box[d] = end - begin;
        rc.covered = rc.covered && begin == chunk_start && end == valid_end;
      }
      fn(rc);
      // odometer over [first, last]
      std::size_t d = rank;
      bool advanced = false;
      while (d-- > 0) {
        if (rc.index[d] < last[d]) {
          ++rc.index[d];
          advanced = true;
          break;
        }
        rc.index[d] = first[d];
      }
      if (!advanced) {
        return;
      }
    }
  }

  std::shared_ptr<Store> store_;
  std::string path_;
  ArrayMeta meta_;
  CodecPipeline pipeline_;
  std::shared_ptr<Store> chunk_store_;
};

}  // namespace zarr

#endif  // LIBZARR_ARRAY_HPP
