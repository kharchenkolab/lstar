// SPDX-License-Identifier: MIT

#ifndef LIBZARR_STORE_HPP
#define LIBZARR_STORE_HPP

#include <cstddef>
#include <cstdint>
#include <map>
#include <optional>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "libzarr/types.hpp"

/// \file store.hpp
/// The key->bytes Store abstraction all libzarr I/O goes through, and the
/// built-in MemoryStore backend.
///
/// Keys are '/'-separated UTF-8 paths with no leading or trailing '/' and no
/// empty segments (e.g. "group/array/c/0/0"). A prefix is either "" (the
/// root) or ends with '/'.

namespace zarr {

/// Byte-range request for Store::read_range.
struct ByteRange {
  /// Which part of the value to read.
  enum class Kind : std::uint8_t {
    full,    ///< the whole value
    slice,   ///< `length` bytes starting at `offset`
    suffix,  ///< the final `length` bytes
  };

  /// Which part of the value to read.
  Kind kind = Kind::full;
  /// Start of the range; used by Kind::slice only.
  std::uint64_t offset = 0;
  /// Number of bytes; used by Kind::slice and Kind::suffix.
  std::uint64_t length = 0;

  /// The whole value.
  [[nodiscard]] static constexpr ByteRange full() { return ByteRange{}; }

  /// `length` bytes starting at byte `offset`. A range that lies outside the
  /// value is an error, not a truncation.
  [[nodiscard]] static constexpr ByteRange slice(std::uint64_t offset, std::uint64_t length) {
    return ByteRange{Kind::slice, offset, length};
  }

  /// The final `length` bytes — fetches a trailing shard index in one round
  /// trip without knowing the value size. `length` greater than the value
  /// size is an error.
  [[nodiscard]] static constexpr ByteRange suffix(std::uint64_t length) {
    return ByteRange{Kind::suffix, 0, length};
  }
};

/// One entry of a batched read (Store::read_many): a key and the byte range
/// to fetch from it (the whole value by default).
struct ReadRequest {
  /// Key to read; must outlive the read_many call (as with any Store key).
  std::string_view key;
  /// Range within the value.
  ByteRange range = ByteRange::full();
};

/// Immediate children of a prefix, as returned by Store::list_dir.
struct DirListing {
  /// Child keys, relative to the queried prefix, sorted.
  std::vector<std::string> keys;
  /// Child prefixes ("directories"), relative, without trailing '/', sorted.
  std::vector<std::string> prefixes;
};

/// Abstract key->bytes store. All array and metadata logic in libzarr is
/// written against this interface; backends supply bytes from memory, files,
/// ZIP archives, HTTP range requests, ...
///
/// The core is single-threaded by design; implementations are not required to
/// be thread-safe.
class Store {
 public:
  Store() = default;
  Store(const Store&) = delete;
  Store& operator=(const Store&) = delete;
  Store(Store&&) = delete;
  Store& operator=(Store&&) = delete;
  virtual ~Store() = default;

  /// Full value at `key`, or std::nullopt if the key is absent.
  [[nodiscard]] virtual std::optional<Bytes> read(std::string_view key) = 0;

  /// Part of the value at `key`: std::nullopt if the key is absent, throws
  /// zarr::error if the range lies outside the value. The default
  /// implementation reads the full value and slices it; backends with native
  /// range reads (files, HTTP) should override.
  [[nodiscard]] virtual std::optional<Bytes> read_range(std::string_view key, ByteRange range);

  /// Size in bytes of the value at `key`, or std::nullopt if absent. The
  /// default implementation reads the whole value; backends with cheap stat
  /// (files: fstat, HTTP: HEAD) should override.
  [[nodiscard]] virtual std::optional<std::uint64_t> size(std::string_view key) {
    const auto value = read(key);
    if (!value) {
      return std::nullopt;
    }
    return value->size();
  }

  /// Reads several ranges in one call, order-preserving: result[i] is the
  /// value for requests[i], or std::nullopt if that key is absent. The
  /// default implementation loops over read_range; backends whose I/O has
  /// per-request latency (HTTP, object stores) should override to issue the
  /// requests concurrently or as coalesced multi-range reads, cutting round
  /// trips. Still synchronous: it returns only once every range is resolved.
  [[nodiscard]] virtual std::vector<std::optional<Bytes>> read_many(
      const std::vector<ReadRequest>& requests) {
    std::vector<std::optional<Bytes>> out;
    out.reserve(requests.size());
    for (const ReadRequest& req : requests) {
      out.push_back(read_range(req.key, req.range));
    }
    return out;
  }

  /// Create or replace the value at `key`.
  virtual void write(std::string_view key, Bytes value) = 0;

  /// True if `key` holds a value.
  [[nodiscard]] virtual bool exists(std::string_view key) = 0;

  /// Remove `key`; removing an absent key is a no-op.
  virtual void erase(std::string_view key) = 0;

  /// Completes any buffered writes. Plain stores write through and need
  /// nothing; adapters that assemble objects (sharding) override.
  virtual void flush() {}

  /// All keys starting with `prefix` ("" or ending in '/'), sorted.
  [[nodiscard]] virtual std::vector<std::string> list_prefix(std::string_view prefix) = 0;

  /// Immediate children under `prefix` ("" or ending in '/').
  [[nodiscard]] virtual DirListing list_dir(std::string_view prefix) = 0;
};

/// In-memory Store backed by a sorted map: the default store for tests, WASM
/// builds, and assembling stores to be packed into archives.
class MemoryStore final : public Store {
 public:
  [[nodiscard]] std::optional<Bytes> read(std::string_view key) override {
    const auto it = map_.find(key);
    if (it == map_.end()) {
      return std::nullopt;
    }
    return it->second;
  }

  void write(std::string_view key, Bytes value) override {
    map_.insert_or_assign(std::string(key), std::move(value));
  }

  [[nodiscard]] std::optional<std::uint64_t> size(std::string_view key) override {
    const auto it = map_.find(key);
    if (it == map_.end()) {
      return std::nullopt;
    }
    return it->second.size();
  }

  [[nodiscard]] bool exists(std::string_view key) override { return map_.find(key) != map_.end(); }

  void erase(std::string_view key) override {
    const auto it = map_.find(key);
    if (it != map_.end()) {
      map_.erase(it);
    }
  }

  [[nodiscard]] std::vector<std::string> list_prefix(std::string_view prefix) override {
    check_prefix(prefix);
    std::vector<std::string> out;
    for (auto it = map_.lower_bound(prefix); it != map_.end() && starts_with(it->first, prefix);
         ++it) {
      out.push_back(it->first);
    }
    return out;
  }

  [[nodiscard]] DirListing list_dir(std::string_view prefix) override {
    check_prefix(prefix);
    DirListing out;
    for (auto it = map_.lower_bound(prefix); it != map_.end() && starts_with(it->first, prefix);
         ++it) {
      const auto rest = std::string_view(it->first).substr(prefix.size());
      const auto slash = rest.find('/');
      if (slash == std::string_view::npos) {
        out.keys.emplace_back(rest);
      } else {
        // Keys under one child prefix are contiguous in a sorted map, so
        // comparing against the last emitted prefix deduplicates.
        const auto child = rest.substr(0, slash);
        if (out.prefixes.empty() || out.prefixes.back() != child) {
          out.prefixes.emplace_back(child);
        }
      }
    }
    return out;
  }

  /// Number of keys held. (Named `key_count`, not `size`, to avoid colliding
  /// with `Store::size(key)`, which returns a value's byte length.)
  [[nodiscard]] std::size_t key_count() const { return map_.size(); }

 private:
  static bool starts_with(std::string_view text, std::string_view prefix) {
    return text.size() >= prefix.size() && text.compare(0, prefix.size(), prefix) == 0;
  }

  static void check_prefix(std::string_view prefix) {
    if (!prefix.empty() && prefix.back() != '/') {
      throw error("store prefix must be empty or end with '/', got '" + std::string(prefix) + "'");
    }
  }

  std::map<std::string, Bytes, std::less<>> map_;
};

inline std::optional<Bytes> Store::read_range(std::string_view key, ByteRange range) {
  auto value = read(key);
  if (!value) {
    return std::nullopt;
  }
  const std::uint64_t size = value->size();
  std::uint64_t begin = 0;
  std::uint64_t count = 0;
  switch (range.kind) {
    case ByteRange::Kind::full:
      return value;
    case ByteRange::Kind::slice:
      if (range.length > size || range.offset > size - range.length) {
        throw error("read_range: slice at offset " + std::to_string(range.offset) + " of length " +
                    std::to_string(range.length) + " out of bounds for \"" + std::string(key) +
                    "\" (" + std::to_string(size) + " bytes)");
      }
      begin = range.offset;
      count = range.length;
      break;
    case ByteRange::Kind::suffix:
      if (range.length > size) {
        throw error("read_range: suffix of length " + std::to_string(range.length) +
                    " out of bounds for \"" + std::string(key) + "\" (" + std::to_string(size) +
                    " bytes)");
      }
      begin = size - range.length;
      count = range.length;
      break;
  }
  const auto first = value->begin() + static_cast<std::ptrdiff_t>(begin);
  return Bytes(first, first + static_cast<std::ptrdiff_t>(count));
}

}  // namespace zarr

#endif  // LIBZARR_STORE_HPP
