// SPDX-License-Identifier: MIT

#ifndef LIBZARR_ZIP_HPP
#define LIBZARR_ZIP_HPP

#include <algorithm>
#include <array>
#include <cstdint>
#include <map>
#include <memory>
#include <optional>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "libzarr/detail/common.hpp"
#include "libzarr/store.hpp"
#include "libzarr/types.hpp"

/// \file zip.hpp
/// Single-file archives: a whole store packed into one ZIP with STORED
/// (uncompressed) entries only, so every entry — and any byte range inside
/// it — stays byte-range-readable through the archive (chunk codecs still
/// apply; the zip layer never re-compresses). ZIP64-aware. ZipStore is a
/// read-only Store view over an archive held in any other Store, using
/// nothing but read_range.

namespace zarr {

namespace detail_zip {

// APPNOTE 4.3.x structure signatures and fixed sizes.
inline constexpr std::uint32_t kLocalSig = 0x04034b50;
inline constexpr std::uint32_t kCentralSig = 0x02014b50;
inline constexpr std::uint32_t kEocdSig = 0x06054b50;
inline constexpr std::uint32_t kEocd64Sig = 0x06064b50;
inline constexpr std::uint32_t kEocd64LocatorSig = 0x07064b50;
inline constexpr std::size_t kEocdSize = 22;
inline constexpr std::size_t kEocd64Size = 56;
inline constexpr std::size_t kEocd64LocatorSize = 20;
inline constexpr std::size_t kLocalHeaderSize = 30;
inline constexpr std::size_t kCentralHeaderSize = 46;
inline constexpr std::uint32_t kMax32 = 0xFFFFFFFFU;
inline constexpr std::uint16_t kMax16 = 0xFFFFU;

inline std::uint16_t rd16(const std::uint8_t* p) {
  return static_cast<std::uint16_t>(p[0] | (static_cast<std::uint16_t>(p[1]) << 8U));
}
inline std::uint32_t rd32(const std::uint8_t* p) {
  return static_cast<std::uint32_t>(p[0]) | (static_cast<std::uint32_t>(p[1]) << 8U) |
         (static_cast<std::uint32_t>(p[2]) << 16U) | (static_cast<std::uint32_t>(p[3]) << 24U);
}
inline std::uint64_t rd64(const std::uint8_t* p) {
  return static_cast<std::uint64_t>(rd32(p)) | (static_cast<std::uint64_t>(rd32(p + 4)) << 32U);
}
inline void wr16(Bytes& out, std::uint16_t v) {
  out.push_back(static_cast<std::uint8_t>(v & 0xFFU));
  out.push_back(static_cast<std::uint8_t>(v >> 8U));
}
inline void wr32(Bytes& out, std::uint32_t v) {
  for (int i = 0; i < 4; ++i) {
    out.push_back(static_cast<std::uint8_t>(v >> (8 * i)));
  }
}
inline void wr64(Bytes& out, std::uint64_t v) {
  for (int i = 0; i < 8; ++i) {
    out.push_back(static_cast<std::uint8_t>(v >> (8 * i)));
  }
}

/// CRC-32 (IEEE 802.3, as required by the zip format). Implemented here so
/// zero-dependency builds do not need zlib.
inline std::uint32_t crc32(const std::uint8_t* data, std::size_t size) {
  static const std::array<std::uint32_t, 256> table = [] {
    std::array<std::uint32_t, 256> t{};
    for (std::uint32_t i = 0; i < 256; ++i) {
      std::uint32_t c = i;
      for (int k = 0; k < 8; ++k) {
        c = ((c & 1U) != 0) ? 0xEDB88320U ^ (c >> 1U) : c >> 1U;
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

struct Entry {
  std::uint64_t header_offset = 0;  ///< of the local file header
  std::uint64_t size = 0;           ///< uncompressed == compressed (STORED)
  std::uint16_t method = 0;
  /// Start of the entry's bytes; resolved lazily from the local header
  /// (whose name/extra lengths may differ from the central directory's).
  std::optional<std::uint64_t> data_offset;
};

struct CentralDirectory {
  std::uint64_t count = 0;
  std::uint64_t size = 0;
  std::uint64_t offset = 0;
};

/// An entry being written by zip_pack.
struct PackEntry {
  std::string name;
  std::uint32_t crc = 0;
  std::uint64_t size = 0;
  std::uint64_t header_offset = 0;
  bool wide = false;  ///< sizes need ZIP64
  bool far = false;   ///< header offset needs ZIP64
};

// APPNOTE 4.3.7. Determinism: zero timestamps (DOS epoch 1980-01-01), no
// flags, no comments; version 4.5 only where ZIP64 structures are present.
inline void append_local_header(Bytes& out, const PackEntry& e) {
  wr32(out, kLocalSig);
  wr16(out, e.wide ? 45 : 20);
  wr16(out, 0);       // flags
  wr16(out, 0);       // method: STORED
  wr16(out, 0);       // mod time
  wr16(out, 0x0021);  // mod date: 1980-01-01
  wr32(out, e.crc);
  const auto size32 = e.wide ? kMax32 : static_cast<std::uint32_t>(e.size);
  wr32(out, size32);  // compressed == uncompressed for STORED
  wr32(out, size32);
  wr16(out, static_cast<std::uint16_t>(e.name.size()));
  wr16(out, e.wide ? 20 : 0);  // extra length
  out.insert(out.end(), e.name.begin(), e.name.end());
  if (e.wide) {
    wr16(out, 0x0001);  // ZIP64 extra: usize + csize
    wr16(out, 16);
    wr64(out, e.size);
    wr64(out, e.size);
  }
}

// APPNOTE 4.3.12.
inline void append_central_header(Bytes& cen, const PackEntry& e) {
  const bool any64 = e.wide || e.far;
  const auto size32 = e.wide ? kMax32 : static_cast<std::uint32_t>(e.size);
  wr32(cen, kCentralSig);
  wr16(cen, any64 ? 45 : 20);  // version made by
  wr16(cen, any64 ? 45 : 20);  // version needed
  wr16(cen, 0);                // flags
  wr16(cen, 0);                // method
  wr16(cen, 0);                // time
  wr16(cen, 0x0021);           // date
  wr32(cen, e.crc);
  wr32(cen, size32);
  wr32(cen, size32);
  wr16(cen, static_cast<std::uint16_t>(e.name.size()));
  wr16(cen, any64 ? static_cast<std::uint16_t>(4 + (e.wide ? 16 : 0) + (e.far ? 8 : 0)) : 0);
  wr16(cen, 0);  // comment length
  wr16(cen, 0);  // disk number
  wr16(cen, 0);  // internal attributes
  wr32(cen, 0);  // external attributes
  wr32(cen, e.far ? kMax32 : static_cast<std::uint32_t>(e.header_offset));
  cen.insert(cen.end(), e.name.begin(), e.name.end());
  if (any64) {
    wr16(cen, 0x0001);
    wr16(cen, static_cast<std::uint16_t>((e.wide ? 16 : 0) + (e.far ? 8 : 0)));
    if (e.wide) {
      wr64(cen, e.size);
      wr64(cen, e.size);
    }
    if (e.far) {
      wr64(cen, e.header_offset);
    }
  }
}

// APPNOTE 4.3.14-16: EOCD64 + locator when any value overflows, then EOCD.
inline void append_end_records(Bytes& out, std::uint64_t count, std::uint64_t cen_size,
                               std::uint64_t cen_offset, bool force_zip64) {
  const bool eocd64 = force_zip64 || count >= kMax16 || cen_size >= kMax32 || cen_offset >= kMax32;
  if (eocd64) {
    const std::uint64_t eocd64_offset = out.size();
    wr32(out, kEocd64Sig);
    wr64(out, kEocd64Size - 12);  // size of the remainder
    wr16(out, 45);                // version made by
    wr16(out, 45);                // version needed
    wr32(out, 0);                 // this disk
    wr32(out, 0);                 // central directory disk
    wr64(out, count);
    wr64(out, count);
    wr64(out, cen_size);
    wr64(out, cen_offset);
    wr32(out, kEocd64LocatorSig);
    wr32(out, 0);  // disk with the EOCD64
    wr64(out, eocd64_offset);
    wr32(out, 1);  // total disks
  }
  wr32(out, kEocdSig);
  wr16(out, 0);  // this disk
  wr16(out, 0);  // central directory disk
  const auto count16 = eocd64 ? kMax16 : static_cast<std::uint16_t>(count);
  wr16(out, count16);
  wr16(out, count16);
  wr32(out, eocd64 ? kMax32 : static_cast<std::uint32_t>(cen_size));
  wr32(out, eocd64 ? kMax32 : static_cast<std::uint32_t>(cen_offset));
  wr16(out, 0);  // comment length
}

}  // namespace detail_zip

/// Read-only Store view over a STORED-entry ZIP archive that itself lives at
/// `archive_key` inside another Store. All access goes through read_range, so
/// a remote archive is never downloaded whole; entry reads cost one range
/// request (plus one 30-byte header read the first time an entry is touched).
class ZipStore final : public Store {
 public:
  /// Opens the archive at `archive_key` in `source` (parses its central
  /// directory with one size probe plus one suffix range read).
  ZipStore(std::shared_ptr<Store> source, std::string archive_key)
      : source_(std::move(source)), key_(std::move(archive_key)) {
    if (!source_) {
      throw error("ZipStore: null store");
    }
    parse_directory();
  }

  [[nodiscard]] std::optional<Bytes> read(std::string_view key) override {
    return read_range(key, ByteRange::full());
  }

  [[nodiscard]] std::optional<Bytes> read_range(std::string_view key, ByteRange range) override {
    const auto it = entries_.find(key);
    if (it == entries_.end()) {
      return std::nullopt;
    }
    detail_zip::Entry& entry = it->second;
    if (entry.method != 0) {
      // Scope guard: only STORED entries stay byte-range-readable.
      throw error(context() + ": entry '" + std::string(key) + "' uses compression method " +
                  std::to_string(entry.method) + "; only STORED entries are supported");
    }
    std::uint64_t begin = 0;
    std::uint64_t count = entry.size;
    if (range.kind == ByteRange::Kind::slice) {
      if (range.length > entry.size || range.offset > entry.size - range.length) {
        throw error(context() + ": range out of bounds for entry '" + std::string(key) + "' (" +
                    std::to_string(entry.size) + " bytes)");
      }
      begin = range.offset;
      count = range.length;
    } else if (range.kind == ByteRange::Kind::suffix) {
      if (range.length > entry.size) {
        throw error(context() + ": suffix out of bounds for entry '" + std::string(key) + "' (" +
                    std::to_string(entry.size) + " bytes)");
      }
      begin = entry.size - range.length;
      count = range.length;
    }
    return must_read(data_offset(key, entry) + begin, count);
  }

  [[nodiscard]] std::optional<std::uint64_t> size(std::string_view key) override {
    const auto it = entries_.find(key);
    if (it == entries_.end()) {
      return std::nullopt;
    }
    return it->second.size;
  }

  [[nodiscard]] bool exists(std::string_view key) override {
    return entries_.find(key) != entries_.end();
  }

  void write(std::string_view /*key*/, Bytes /*value*/) override {
    throw error("ZipStore is read-only");
  }
  void erase(std::string_view /*key*/) override { throw error("ZipStore is read-only"); }

  [[nodiscard]] std::vector<std::string> list_prefix(std::string_view prefix) override {
    check_prefix(prefix);
    std::vector<std::string> out;
    for (auto it = entries_.lower_bound(prefix);
         it != entries_.end() && detail::starts_with(it->first, prefix); ++it) {
      out.push_back(it->first);
    }
    return out;
  }

  [[nodiscard]] DirListing list_dir(std::string_view prefix) override {
    check_prefix(prefix);
    DirListing out;
    for (auto it = entries_.lower_bound(prefix);
         it != entries_.end() && detail::starts_with(it->first, prefix); ++it) {
      const auto rest = std::string_view(it->first).substr(prefix.size());
      const auto slash = rest.find('/');
      if (slash == std::string_view::npos) {
        out.keys.emplace_back(rest);
      } else {
        const auto child = rest.substr(0, slash);
        if (out.prefixes.empty() || out.prefixes.back() != child) {
          out.prefixes.emplace_back(child);
        }
      }
    }
    return out;
  }

  /// Number of entries in the archive.
  [[nodiscard]] std::size_t entry_count() const { return entries_.size(); }

 private:
  [[nodiscard]] std::string context() const { return key_.empty() ? "zip archive" : key_; }

  static void check_prefix(std::string_view prefix) {
    if (!prefix.empty() && prefix.back() != '/') {
      throw error("store prefix must be empty or end with '/', got '" + std::string(prefix) + "'");
    }
  }

  [[nodiscard]] Bytes must_read(std::uint64_t offset, std::uint64_t length) {
    auto bytes = source_->read_range(key_, ByteRange::slice(offset, length));
    if (!bytes) {
      throw error(context() + ": archive disappeared mid-read");
    }
    return *std::move(bytes);
  }

  void parse_directory() {
    namespace z = detail_zip;
    const auto archive_size = source_->size(key_);
    if (!archive_size) {
      throw error(context() + ": not found");
    }
    if (*archive_size < z::kEocdSize) {
      throw error(context() + ": too small to be a zip archive");
    }
    // The EOCD sits within the final 22 + 65535 bytes (max comment); one
    // suffix read also covers a possible ZIP64 locator just before it.
    const std::uint64_t tail_len =
        std::min<std::uint64_t>(*archive_size, z::kEocdSize + z::kMax16 + z::kEocd64LocatorSize);
    auto tail_bytes = source_->read_range(key_, ByteRange::suffix(tail_len));
    if (!tail_bytes) {
      throw error(context() + ": archive disappeared mid-read");
    }
    const Bytes tail = *std::move(tail_bytes);
    const std::uint64_t tail_start = *archive_size - tail_len;

    const std::size_t eocd = find_eocd(tail);
    const detail_zip::CentralDirectory dir = locate_directory(tail, tail_start, eocd);
    const Bytes cen = must_read(dir.offset, dir.size);
    std::size_t pos = 0;
    for (std::uint64_t n = 0; n < dir.count; ++n) {
      parse_entry(cen, pos);
    }
  }

  /// Finds the EOCD: scan back from the end; the comment must reach the end.
  [[nodiscard]] std::size_t find_eocd(const Bytes& tail) const {
    namespace z = detail_zip;
    for (std::size_t i = tail.size() - z::kEocdSize + 1; i-- > 0;) {
      if (z::rd32(tail.data() + i) == z::kEocdSig &&
          i + z::kEocdSize + z::rd16(tail.data() + i + 20) == tail.size()) {
        return i;
      }
    }
    throw error(context() + ": end-of-central-directory record not found (not a zip archive?)");
  }

  /// Reads the (possibly ZIP64) central-directory location from the EOCD.
  [[nodiscard]] detail_zip::CentralDirectory locate_directory(const Bytes& tail,
                                                              std::uint64_t tail_start,
                                                              std::size_t eocd) {
    namespace z = detail_zip;
    if (z::rd16(tail.data() + eocd + 4) != 0 || z::rd16(tail.data() + eocd + 6) != 0) {
      throw error(context() + ": multi-disk archives are not supported");
    }
    z::CentralDirectory dir;
    dir.count = z::rd16(tail.data() + eocd + 10);
    dir.size = z::rd32(tail.data() + eocd + 12);
    dir.offset = z::rd32(tail.data() + eocd + 16);
    if (dir.count != z::kMax16 && dir.size != z::kMax32 && dir.offset != z::kMax32) {
      return dir;
    }
    // ZIP64: the locator sits immediately before the EOCD.
    if (eocd < z::kEocd64LocatorSize ||
        z::rd32(tail.data() + eocd - z::kEocd64LocatorSize) != z::kEocd64LocatorSig) {
      throw error(context() + ": ZIP64 locator not found");
    }
    const std::uint64_t eocd64_offset = z::rd64(tail.data() + eocd - z::kEocd64LocatorSize + 8);
    Bytes eocd64;
    if (eocd64_offset >= tail_start) {
      const auto local = static_cast<std::size_t>(eocd64_offset - tail_start);
      // APPNOTE 4.3.14: the ZIP64 EOCD must lie fully within the archive. A
      // crafted locator can point past the tail we hold; assigning past end()
      // is an out-of-bounds read (fuzz-found SEGV).
      if (local > tail.size() || z::kEocd64Size > tail.size() - local) {
        throw error(context() + ": ZIP64 end-of-central-directory record out of range");
      }
      eocd64.assign(tail.begin() + static_cast<std::ptrdiff_t>(local),
                    tail.begin() + static_cast<std::ptrdiff_t>(local + z::kEocd64Size));
    } else {
      eocd64 = must_read(eocd64_offset, z::kEocd64Size);
    }
    // must_read returns fewer bytes when the offset runs past EOF; the fixed
    // field offsets below (rd32/rd64) require the whole record.
    if (eocd64.size() < z::kEocd64Size) {
      throw error(context() + ": ZIP64 end-of-central-directory record truncated");
    }
    if (z::rd32(eocd64.data()) != z::kEocd64Sig) {
      throw error(context() + ": bad ZIP64 end-of-central-directory record");
    }
    dir.count = z::rd64(eocd64.data() + 32);
    dir.size = z::rd64(eocd64.data() + 40);
    dir.offset = z::rd64(eocd64.data() + 48);
    return dir;
  }

  /// Parses one central-directory record at `pos` (advanced past it).
  void parse_entry(const Bytes& cen, std::size_t& pos) {
    namespace z = detail_zip;
    if (pos + z::kCentralHeaderSize > cen.size() || z::rd32(cen.data() + pos) != z::kCentralSig) {
      throw error(context() + ": corrupt central directory");
    }
    const std::uint8_t* h = cen.data() + pos;
    z::Entry entry;
    entry.method = z::rd16(h + 10);
    std::uint64_t csize = z::rd32(h + 20);
    entry.size = z::rd32(h + 24);
    const std::size_t name_len = z::rd16(h + 28);
    const std::size_t extra_len = z::rd16(h + 30);
    const std::size_t comment_len = z::rd16(h + 32);
    entry.header_offset = z::rd32(h + 42);
    if (pos + z::kCentralHeaderSize + name_len + extra_len + comment_len > cen.size()) {
      throw error(context() + ": corrupt central directory");
    }
    std::string name(reinterpret_cast<const char*>(h + z::kCentralHeaderSize), name_len);
    apply_zip64_extra(h + z::kCentralHeaderSize + name_len, extra_len, entry, csize, name);

    if (entry.method == 0 && csize != entry.size) {
      throw error(context() + ": STORED entry '" + name + "' has mismatched sizes");
    }
    if (!name.empty() && name.back() != '/') {  // skip directory placeholders
      entries_.insert_or_assign(std::move(name), entry);
    }
    pos += z::kCentralHeaderSize + name_len + extra_len + comment_len;
  }

  /// ZIP64 extra field (id 0x0001): 64-bit values, present only for the
  /// 32-bit fields that hold the 0xFFFFFFFF marker, in the fixed order
  /// usize, csize, offset.
  void apply_zip64_extra(const std::uint8_t* extra, std::size_t extra_len, detail_zip::Entry& entry,
                         std::uint64_t& csize, const std::string& name) const {
    namespace z = detail_zip;
    std::size_t epos = 0;
    while (epos + 4 <= extra_len) {
      const std::uint16_t id = z::rd16(extra + epos);
      const std::uint16_t len = z::rd16(extra + epos + 2);
      if (epos + 4 + len > extra_len) {
        throw error(context() + ": corrupt extra field in '" + name + "'");
      }
      if (id == 0x0001) {
        const std::uint8_t* f = extra + epos + 4;
        std::size_t fpos = 0;
        const auto take64 = [&](std::uint64_t& value) {
          if (fpos + 8 > len) {
            throw error(context() + ": truncated ZIP64 extra field in '" + name + "'");
          }
          value = z::rd64(f + fpos);
          fpos += 8;
        };
        if (entry.size == z::kMax32) {
          take64(entry.size);
        }
        if (csize == z::kMax32) {
          take64(csize);
        }
        if (entry.header_offset == z::kMax32) {
          take64(entry.header_offset);
        }
      }
      epos += std::size_t{4} + len;
    }
  }

  /// Resolves where an entry's bytes start: the central directory does not
  /// record it, and the local header's name/extra lengths can differ from
  /// the central ones.
  std::uint64_t data_offset(std::string_view key, detail_zip::Entry& entry) {
    namespace z = detail_zip;
    if (!entry.data_offset) {
      const Bytes lfh = must_read(entry.header_offset, z::kLocalHeaderSize);
      if (z::rd32(lfh.data()) != z::kLocalSig) {
        throw error(context() + ": corrupt local header for entry '" + std::string(key) + "'");
      }
      entry.data_offset = entry.header_offset + z::kLocalHeaderSize + z::rd16(lfh.data() + 26) +
                          z::rd16(lfh.data() + 28);
    }
    return *entry.data_offset;
  }

  std::shared_ptr<Store> source_;
  std::string key_;
  std::map<std::string, detail_zip::Entry, std::less<>> entries_;
};

namespace detail_zip {

/// Implementation of zip_pack. `force_zip64` forces the ZIP64 structures even
/// when no value overflows 32 bits — used by tests to exercise them on small
/// archives; not part of the public API.
inline void zip_pack_impl(Store& source, Store& dest, const std::string& dest_key,
                          const std::string& prefix, bool force_zip64) {
  namespace z = detail_zip;
  Bytes out;
  Bytes cen;
  std::uint64_t count = 0;

  for (const std::string& key : source.list_prefix(prefix)) {
    const auto value = source.read(key);
    if (!value) {
      continue;  // key vanished between list and read
    }
    z::PackEntry entry;
    entry.name = key.substr(prefix.size());
    if (entry.name.size() > z::kMax16) {
      throw error("zip_pack: entry name too long: '" + entry.name + "'");
    }
    entry.crc = z::crc32(value->data(), value->size());
    entry.size = value->size();
    entry.header_offset = out.size();
    entry.wide = force_zip64 || entry.size >= z::kMax32;
    entry.far = force_zip64 || entry.header_offset >= z::kMax32;
    z::append_local_header(out, entry);
    out.insert(out.end(), value->begin(), value->end());
    z::append_central_header(cen, entry);
    ++count;
  }

  const std::uint64_t cen_offset = out.size();
  out.insert(out.end(), cen.begin(), cen.end());
  z::append_end_records(out, count, cen.size(), cen_offset, force_zip64);
  dest.write(dest_key, std::move(out));
}

}  // namespace detail_zip

/// Packs every key of `source` under `prefix` into a STORED-entry ZIP written
/// at `dest_key` in `dest`. Deterministic byte-for-byte: sorted entries, zero
/// timestamps (DOS epoch), no comments.
inline void zip_pack(Store& source, Store& dest, const std::string& dest_key,
                     const std::string& prefix = "") {
  detail_zip::zip_pack_impl(source, dest, dest_key, prefix, /*force_zip64=*/false);
}

}  // namespace zarr

#endif  // LIBZARR_ZIP_HPP
