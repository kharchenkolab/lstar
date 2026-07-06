// SPDX-License-Identifier: MIT

#ifndef LIBZARR_SHARDING_HPP
#define LIBZARR_SHARDING_HPP

#include <cstdint>
#include <limits>
#include <memory>
#include <optional>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "libzarr/codecs.hpp"
#include "libzarr/detail/common.hpp"
#include "libzarr/metadata.hpp"
#include "libzarr/store.hpp"
#include "libzarr/types.hpp"
#include "libzarr/v2.hpp"
#include "libzarr/v3.hpp"

/// \file sharding.hpp
/// The sharding_indexed codec, modeled as a Store adapter rather than a
/// codec: a shard is simply the outer chunk's stored object, and ShardStore
/// maps inner-chunk keys onto byte ranges of it via the trailing (or
/// leading) index. The array machinery above stays completely unaware, and
/// nested sharding falls out as ShardStore-wrapping-ShardStore.

namespace zarr {

namespace detail_shard {

inline constexpr std::uint64_t kSentinel = std::numeric_limits<std::uint64_t>::max();

struct IndexEntry {
  std::uint64_t offset = kSentinel;
  std::uint64_t nbytes = kSentinel;
  [[nodiscard]] bool missing() const { return offset == kSentinel && nbytes == kSentinel; }
};

/// Validates index_codecs (sharding spec: the encoded index must have a
/// fixed size, so only `bytes` and `crc32c` qualify) and returns the encoded
/// byte count for `entry_count` [offset, nbytes] pairs.
inline std::uint64_t index_encoded_size(const std::vector<CodecSpec>& index_codecs,
                                        std::uint64_t entry_count, const std::string& ctx) {
  std::uint64_t size = entry_count * 2 * 8;
  bool have_bytes = false;
  for (const CodecSpec& codec : index_codecs) {
    if (codec.name == "bytes") {
      have_bytes = true;
    } else if (codec.name == "crc32c") {
      size += 4;
    } else {
      throw error(ctx +
                  ": index_codecs may only contain 'bytes' and 'crc32c' (the index must "
                  "have a fixed encoded size), got '" +
                  codec.name + "'");
    }
  }
  if (!have_bytes) {
    throw error(ctx + ": index_codecs must contain the 'bytes' codec");
  }
  return size;
}

}  // namespace detail_shard

/// How a ShardStore maps inner-chunk keys onto shards.
struct ShardParams {
  /// Store-key prefix of the array's chunks ("" for a root array, else
  /// "path/").
  std::string chunk_prefix;
  /// Chunk-key scheme shared by the inner keys and the shard keys.
  ChunkKeyKind key_encoding = ChunkKeyKind::v3_default;
  /// Separator of both key layers.
  char separator = '/';
  /// Inner chunks per shard, per dimension (outer shape / inner shape).
  std::vector<std::uint64_t> per_shard;
  /// Total inner-chunk grid extents (for key validation).
  std::vector<std::uint64_t> inner_grid;
  /// Codecs of the shard index (v3 sharding spec: fixed-size only).
  std::vector<CodecSpec> index_codecs;
  /// v3 sharding spec: index_location "end" (default) or "start".
  bool index_at_end = true;

  /// The params for shard `level` of a (possibly nested) sharded array — the exact computation the
  /// array machinery uses to wrap its chunk store, factored out so external byte-range readers can
  /// build the SAME mapping from metadata alone. `prefix` is the store-key prefix of the array's
  /// chunks ("" yields leaf-relative keys). Inner shape is the next level's shard shape, or the
  /// chunk shape at the innermost level.
  [[nodiscard]] static ShardParams for_level(const ArrayMeta& meta, std::size_t level,
                                             const std::string& prefix) {
    const std::vector<std::uint64_t>& inner_shape =
        level + 1 < meta.shard_levels.size() ? meta.shard_levels[level + 1].shard_shape
                                             : meta.chunk_shape;
    const ShardLevel& lvl = meta.shard_levels[level];
    ShardParams params;
    params.chunk_prefix = prefix;
    params.key_encoding = meta.key_encoding;
    params.separator = meta.dimension_separator;
    params.index_codecs = lvl.index_codecs;
    params.index_at_end = lvl.index_at_end;
    params.per_shard.resize(inner_shape.size());
    params.inner_grid.resize(inner_shape.size());
    for (std::size_t d = 0; d < inner_shape.size(); ++d) {
      params.per_shard[d] = lvl.shard_shape[d] / inner_shape[d];
      params.inner_grid[d] = detail::ceil_div(meta.shape[d], inner_shape[d]);
    }
    return params;
  }
};

/// Store adapter presenting the inner chunks of a sharded array as ordinary
/// keys. Reads cost one index fetch per shard (cached, one suffix/prefix
/// range request) plus one range request per inner chunk. Writes assemble
/// whole shards in memory — the peak footprint is one encoded shard — and
/// are completed by flush() (called automatically when writes move to
/// another shard, and by Array after every write operation).
///
/// Internal adapter: listing is not supported.
class ShardStore final : public Store {
 public:
  /// Binds inner-chunk keys of `params` onto shard objects in `source`.
  ShardStore(std::shared_ptr<Store> source, ShardParams params)
      : source_(std::move(source)),
        params_(std::move(params)),
        entry_count_(detail::checked_product(params_.per_shard, "shard grid")),
        index_size_(detail_shard::index_encoded_size(params_.index_codecs, entry_count_,
                                                     "sharding_indexed")),
        index_pipeline_(CodecPipeline::resolve(index_meta())) {
    if (!source_) {
      throw error("ShardStore: null store");
    }
  }

  ~ShardStore() override = default;
  ShardStore(const ShardStore&) = delete;
  ShardStore& operator=(const ShardStore&) = delete;
  ShardStore(ShardStore&&) = delete;
  ShardStore& operator=(ShardStore&&) = delete;

  [[nodiscard]] std::optional<Bytes> read(std::string_view key) override {
    return read_range(key, ByteRange::full());
  }

  [[nodiscard]] std::optional<Bytes> read_range(std::string_view key, ByteRange range) override {
    const Location loc = locate(key);
    if (assembly_ && assembly_->shard_key == loc.shard_key) {
      // Serve pending writes so read-modify-write sequences stay coherent.
      const auto& entry = assembly_->entries[loc.intra];
      if (!entry) {
        return std::nullopt;
      }
      return slice_value(*entry, range, key);
    }
    const auto* index = load_index(loc.shard_key);
    if (index == nullptr || (*index)[loc.intra].missing()) {
      return std::nullopt;
    }
    const detail_shard::IndexEntry entry = (*index)[loc.intra];
    std::uint64_t begin = 0;
    std::uint64_t count = entry.nbytes;
    resolve_range(range, entry.nbytes, begin, count, key);
    return source_->read_range(loc.shard_key, ByteRange::slice(entry.offset + begin, count));
  }

  [[nodiscard]] std::optional<std::uint64_t> size(std::string_view key) override {
    const Location loc = locate(key);
    if (assembly_ && assembly_->shard_key == loc.shard_key) {
      const auto& entry = assembly_->entries[loc.intra];
      if (!entry) {
        return std::nullopt;
      }
      return entry->size();
    }
    const auto* index = load_index(loc.shard_key);
    if (index == nullptr || (*index)[loc.intra].missing()) {
      return std::nullopt;
    }
    return (*index)[loc.intra].nbytes;
  }

  [[nodiscard]] bool exists(std::string_view key) override { return size(key).has_value(); }

  // --- Byte-range support for external readers that own their fetch loop ---------------------------
  // The two-phase resolve behind read_range, exposed so a reader (e.g. the WASM/JS viewer) can fetch
  // only a shard's index and then only the wanted chunk's bytes, instead of the whole shard object.

  /// Where an inner chunk lives: the store key of its owning shard, the chunk's C-order slot within
  /// that shard, and the shard index's on-disk layout (encoded size + end/start) so the reader knows
  /// which bytes of the shard to fetch for the index. Pure — no store read.
  struct ChunkPlacement {
    std::string shard_key;
    std::uint64_t intra = 0;
    std::uint64_t index_size = 0;
    bool index_at_end = true;
  };
  /// The chunk's bytes within its shard: [offset, nbytes), or `missing` when the chunk is a fill
  /// (all-default) chunk absent from the shard.
  struct ChunkExtent {
    std::uint64_t offset = 0;
    std::uint64_t nbytes = 0;
    bool missing = false;
  };

  [[nodiscard]] ChunkPlacement place(std::string_view inner_key) const {
    const Location loc = locate(inner_key);
    return {loc.shard_key, static_cast<std::uint64_t>(loc.intra), index_size_, params_.index_at_end};
  }

  /// Decode a fetched shard index (the `index_size` bytes at the index location) and return the
  /// slot's extent within the shard. The caller supplies the bytes; no store read happens here.
  [[nodiscard]] ChunkExtent extent(const Bytes& index_bytes, std::uint64_t intra) const {
    if (intra >= entry_count_) {
      throw error("ShardStore::extent: slot out of range");
    }
    const Bytes decoded = index_pipeline_.decode(Bytes(index_bytes));
    detail_shard::IndexEntry e;
    std::memcpy(&e.offset, decoded.data() + intra * 16, 8);
    std::memcpy(&e.nbytes, decoded.data() + intra * 16 + 8, 8);
    if (e.missing()) {
      return {0, 0, true};
    }
    return {e.offset, e.nbytes, false};
  }

  void write(std::string_view key, Bytes value) override {
    put(locate(key), std::optional<Bytes>(std::move(value)));
  }

  void erase(std::string_view key) override { put(locate(key), std::nullopt); }

  /// Writes out the shard being assembled, then flushes the source (which
  /// matters when the source is itself a ShardStore — nested sharding).
  void flush() override {
    flush_assembly();
    source_->flush();
  }

  [[nodiscard]] std::vector<std::string> list_prefix(std::string_view /*prefix*/) override {
    throw error("ShardStore is an internal adapter; listing is not supported");
  }
  [[nodiscard]] DirListing list_dir(std::string_view /*prefix*/) override {
    throw error("ShardStore is an internal adapter; listing is not supported");
  }

 private:
  struct Location {
    std::string shard_key;
    std::size_t intra = 0;  // C-order position within the shard
  };

  struct Assembly {
    std::string shard_key;
    std::vector<std::optional<Bytes>> entries;
    bool dirty = false;
  };

  /// Metadata describing the index as a decodable "array": uint64 pairs.
  [[nodiscard]] ArrayMeta index_meta() const {
    ArrayMeta meta;
    meta.shape = {entry_count_ * 2};
    meta.chunk_shape = {entry_count_ * 2};
    meta.dtype = DataType::of(DType::uint64);
    meta.codecs = params_.index_codecs;
    return meta;
  }

  static void resolve_range(ByteRange range, std::uint64_t size, std::uint64_t& begin,
                            std::uint64_t& count, std::string_view key) {
    if (range.kind == ByteRange::Kind::slice) {
      if (range.length > size || range.offset > size - range.length) {
        throw error("read_range: slice out of bounds for inner chunk '" + std::string(key) + "' (" +
                    std::to_string(size) + " bytes)");
      }
      begin = range.offset;
      count = range.length;
    } else if (range.kind == ByteRange::Kind::suffix) {
      if (range.length > size) {
        throw error("read_range: suffix out of bounds for inner chunk '" + std::string(key) +
                    "' (" + std::to_string(size) + " bytes)");
      }
      begin = size - range.length;
      count = range.length;
    }
  }

  [[nodiscard]] static std::optional<Bytes> slice_value(const Bytes& value, ByteRange range,
                                                        std::string_view key) {
    std::uint64_t begin = 0;
    std::uint64_t count = value.size();
    resolve_range(range, value.size(), begin, count, key);
    const auto first = value.begin() + static_cast<std::ptrdiff_t>(begin);
    return Bytes(first, first + static_cast<std::ptrdiff_t>(count));
  }

  /// Parses an inner-chunk key back into grid indices and splits it into the
  /// owning shard's key plus the C-order position inside that shard.
  [[nodiscard]] Location locate(std::string_view key) const {
    const std::size_t rank = params_.per_shard.size();
    std::string_view rest = key;
    if (!detail::starts_with(rest, params_.chunk_prefix)) {
      throw error("ShardStore: key '" + std::string(key) + "' is outside the array's chunks");
    }
    rest = rest.substr(params_.chunk_prefix.size());
    if (params_.key_encoding == ChunkKeyKind::v3_default) {
      if (rest.empty() || rest[0] != 'c' || (rest.size() > 1 && rest[1] != params_.separator)) {
        throw error("ShardStore: malformed chunk key '" + std::string(key) + "'");
      }
      rest = rest.size() > 1 ? rest.substr(2) : rest.substr(1);
    }
    std::vector<std::uint64_t> index(rank, 0);
    std::size_t pos = 0;
    for (std::size_t d = 0; d < rank; ++d) {
      std::uint64_t value = 0;
      const std::size_t start = pos;
      while (pos < rest.size() && rest[pos] >= '0' && rest[pos] <= '9') {
        value = value * 10 + static_cast<std::uint64_t>(rest[pos] - '0');
        ++pos;
      }
      if (pos == start || value >= params_.inner_grid[d]) {
        throw error("ShardStore: malformed chunk key '" + std::string(key) + "'");
      }
      index[d] = value;
      if (d + 1 < rank) {
        if (pos >= rest.size() || rest[pos] != params_.separator) {
          throw error("ShardStore: malformed chunk key '" + std::string(key) + "'");
        }
        ++pos;
      }
    }
    if (pos != rest.size()) {
      throw error("ShardStore: malformed chunk key '" + std::string(key) + "'");
    }

    Location loc;
    std::vector<std::uint64_t> outer(rank, 0);
    std::uint64_t intra = 0;
    for (std::size_t d = 0; d < rank; ++d) {
      outer[d] = index[d] / params_.per_shard[d];
      intra = intra * params_.per_shard[d] + index[d] % params_.per_shard[d];
    }
    loc.intra = detail::checked_size(intra, "shard entry");
    const std::string relative = params_.key_encoding == ChunkKeyKind::v3_default
                                     ? v3::chunk_key(outer, params_.separator)
                                     : v2::chunk_key(outer, params_.separator);
    loc.shard_key = params_.chunk_prefix + relative;
    return loc;
  }

  /// Fetches and decodes a shard's index with one range request; nullptr if
  /// the whole shard is absent. Decoded indices are kept in a small LRU.
  [[nodiscard]] const std::vector<detail_shard::IndexEntry>* load_index(
      const std::string& shard_key) {
    for (std::size_t i = 0; i < cache_.size(); ++i) {
      if (cache_[i].first == shard_key) {
        if (i != 0) {
          std::rotate(cache_.begin(), cache_.begin() + static_cast<std::ptrdiff_t>(i),
                      cache_.begin() + static_cast<std::ptrdiff_t>(i) + 1);
        }
        return &cache_.front().second;
      }
    }
    auto stored =
        source_->read_range(shard_key, params_.index_at_end ? ByteRange::suffix(index_size_)
                                                            : ByteRange::slice(0, index_size_));
    if (!stored) {
      return nullptr;
    }
    const Bytes decoded = index_pipeline_.decode(std::move(*stored));
    std::vector<detail_shard::IndexEntry> entries(static_cast<std::size_t>(entry_count_));
    for (std::size_t i = 0; i < entries.size(); ++i) {
      std::memcpy(&entries[i].offset, decoded.data() + i * 16, 8);
      std::memcpy(&entries[i].nbytes, decoded.data() + i * 16 + 8, 8);
      const bool sentinel_mismatch = (entries[i].offset == detail_shard::kSentinel) !=
                                     (entries[i].nbytes == detail_shard::kSentinel);
      const bool overflow =
          !entries[i].missing() && entries[i].nbytes > detail_shard::kSentinel - entries[i].offset;
      if (sentinel_mismatch || overflow) {
        throw error(shard_key + ": corrupt shard index");
      }
    }
    cache_.insert(cache_.begin(), {shard_key, std::move(entries)});
    if (cache_.size() > kCacheCapacity) {
      cache_.pop_back();
    }
    return &cache_.front().second;
  }

  void put(const Location& loc, std::optional<Bytes> value) {
    if (assembly_ && assembly_->shard_key != loc.shard_key) {
      flush_assembly();  // early flush: writes moved on to another shard
    }
    if (!assembly_) {
      assembly_ = load_assembly(loc.shard_key);
    }
    assembly_->entries[loc.intra] = std::move(value);
    assembly_->dirty = true;
  }

  /// Read-modify-write: one full read of the existing shard seeds the
  /// assembly; a missing shard seeds it empty.
  [[nodiscard]] Assembly load_assembly(const std::string& shard_key) {
    Assembly assembly;
    assembly.shard_key = shard_key;
    assembly.entries.assign(static_cast<std::size_t>(entry_count_), std::nullopt);
    const auto stored = source_->read(shard_key);
    if (stored) {
      const auto* index = load_index(shard_key);
      assert(index != nullptr);
      for (std::size_t i = 0; i < index->size(); ++i) {
        const detail_shard::IndexEntry entry = (*index)[i];
        if (entry.missing()) {
          continue;
        }
        if (entry.offset > stored->size() || entry.nbytes > stored->size() - entry.offset) {
          throw error(shard_key + ": shard index points outside the shard");
        }
        const auto first = stored->begin() + static_cast<std::ptrdiff_t>(entry.offset);
        assembly.entries[i] = Bytes(first, first + static_cast<std::ptrdiff_t>(entry.nbytes));
      }
    }
    return assembly;
  }

  void flush_assembly() {
    if (!assembly_) {
      return;
    }
    Assembly assembly = *std::move(assembly_);
    assembly_.reset();
    if (!assembly.dirty) {
      return;
    }
    drop_cached(assembly.shard_key);
    bool any = false;
    for (const auto& entry : assembly.entries) {
      any = any || entry.has_value();
    }
    if (!any) {
      source_->erase(assembly.shard_key);  // all-fill shards are not stored
      return;
    }

    Bytes body;
    std::vector<std::uint64_t> raw_index(assembly.entries.size() * 2, detail_shard::kSentinel);
    const std::uint64_t base = params_.index_at_end ? 0 : index_size_;
    for (std::size_t i = 0; i < assembly.entries.size(); ++i) {
      const std::optional<Bytes>& entry = assembly.entries[i];
      if (!entry) {
        continue;
      }
      raw_index[i * 2] = base + body.size();
      raw_index[i * 2 + 1] = entry->size();
      body.insert(body.end(), entry->begin(), entry->end());
    }
    Bytes index_bytes(raw_index.size() * 8);
    std::memcpy(index_bytes.data(), raw_index.data(), index_bytes.size());
    index_bytes = index_pipeline_.encode(std::move(index_bytes));

    Bytes shard;
    shard.reserve(body.size() + index_bytes.size());
    if (params_.index_at_end) {
      shard = std::move(body);
      shard.insert(shard.end(), index_bytes.begin(), index_bytes.end());
    } else {
      shard = std::move(index_bytes);
      shard.insert(shard.end(), body.begin(), body.end());
    }
    source_->write(assembly.shard_key, std::move(shard));
  }

  void drop_cached(const std::string& shard_key) {
    for (std::size_t i = 0; i < cache_.size(); ++i) {
      if (cache_[i].first == shard_key) {
        cache_.erase(cache_.begin() + static_cast<std::ptrdiff_t>(i));
        return;
      }
    }
  }

  static constexpr std::size_t kCacheCapacity = 16;

  std::shared_ptr<Store> source_;
  ShardParams params_;
  std::uint64_t entry_count_;
  std::uint64_t index_size_;
  CodecPipeline index_pipeline_;
  std::vector<std::pair<std::string, std::vector<detail_shard::IndexEntry>>> cache_;
  std::optional<Assembly> assembly_;
};

}  // namespace zarr

#endif  // LIBZARR_SHARDING_HPP
