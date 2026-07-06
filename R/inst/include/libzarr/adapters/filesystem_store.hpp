// SPDX-License-Identifier: MIT

#ifndef LIBZARR_ADAPTERS_FILESYSTEM_STORE_HPP
#define LIBZARR_ADAPTERS_FILESYSTEM_STORE_HPP

#include <algorithm>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <ios>
#include <optional>
#include <string>
#include <string_view>
#include <system_error>
#include <utility>
#include <vector>

#include "libzarr/detail/common.hpp"
#include "libzarr/store.hpp"
#include "libzarr/types.hpp"

/// \file filesystem_store.hpp
/// Directory-backed Store. This is an *adapter*: it is the only libzarr
/// header that touches the filesystem, and WASM builds simply omit it — the
/// core never includes it.

namespace zarr {

/// Store mapping keys to files under a root directory.
class FilesystemStore final : public Store {
 public:
  /// Binds to `root`, creating the directory when `create` is true.
  explicit FilesystemStore(std::filesystem::path root, bool create = true)
      : root_(std::move(root)) {
    if (create) {
      std::filesystem::create_directories(root_);
    }
  }

  [[nodiscard]] std::optional<Bytes> read(std::string_view key) override {
    std::ifstream in(key_path(key), std::ios::binary);
    if (!in) {
      return std::nullopt;
    }
    Bytes out;
    in.seekg(0, std::ios::end);
    const std::streamoff size = in.tellg();
    in.seekg(0, std::ios::beg);
    out.resize(static_cast<std::size_t>(size));
    in.read(reinterpret_cast<char*>(out.data()), size);
    if (!in) {
      throw error("failed to read '" + std::string(key) + "'");
    }
    return out;
  }

  [[nodiscard]] std::optional<Bytes> read_range(std::string_view key, ByteRange range) override {
    if (range.kind == ByteRange::Kind::full) {
      return read(key);
    }
    std::ifstream in(key_path(key), std::ios::binary);
    if (!in) {
      return std::nullopt;
    }
    in.seekg(0, std::ios::end);
    const auto size = static_cast<std::uint64_t>(in.tellg());
    std::uint64_t begin = 0;
    if (range.kind == ByteRange::Kind::slice) {
      if (range.length > size || range.offset > size - range.length) {
        throw error("read_range: slice at offset " + std::to_string(range.offset) + " of length " +
                    std::to_string(range.length) + " out of bounds for '" + std::string(key) +
                    "' (" + std::to_string(size) + " bytes)");
      }
      begin = range.offset;
    } else {  // suffix
      if (range.length > size) {
        throw error("read_range: suffix of length " + std::to_string(range.length) +
                    " out of bounds for '" + std::string(key) + "' (" + std::to_string(size) +
                    " bytes)");
      }
      begin = size - range.length;
    }
    Bytes out(detail::checked_size(range.length, "read_range"));
    in.seekg(static_cast<std::streamoff>(begin), std::ios::beg);
    in.read(reinterpret_cast<char*>(out.data()), static_cast<std::streamsize>(out.size()));
    if (!in) {
      throw error("failed to read range from '" + std::string(key) + "'");
    }
    return out;
  }

  void write(std::string_view key, Bytes value) override {
    const std::filesystem::path path = key_path(key);
    std::filesystem::create_directories(path.parent_path());
    std::ofstream out(path, std::ios::binary | std::ios::trunc);
    if (!out) {
      throw error("cannot open '" + std::string(key) + "' for writing");
    }
    out.write(reinterpret_cast<const char*>(value.data()),
              static_cast<std::streamsize>(value.size()));
    if (!out) {
      throw error("failed to write '" + std::string(key) + "'");
    }
  }

  [[nodiscard]] std::optional<std::uint64_t> size(std::string_view key) override {
    std::error_code ec;
    const auto bytes = std::filesystem::file_size(key_path(key), ec);
    if (ec) {
      return std::nullopt;
    }
    return bytes;
  }

  [[nodiscard]] bool exists(std::string_view key) override {
    return std::filesystem::is_regular_file(key_path(key));
  }

  void erase(std::string_view key) override {
    std::error_code ec;
    std::filesystem::remove(key_path(key), ec);  // absent key is a no-op
  }

  [[nodiscard]] std::vector<std::string> list_prefix(std::string_view prefix) override {
    check_prefix(prefix);
    std::vector<std::string> out;
    if (!std::filesystem::is_directory(root_)) {
      return out;
    }
    for (const auto& entry : std::filesystem::recursive_directory_iterator(root_)) {
      if (!entry.is_regular_file()) {
        continue;
      }
      const std::string key = std::filesystem::relative(entry.path(), root_).generic_string();
      if (detail::starts_with(key, prefix)) {
        out.push_back(key);
      }
    }
    std::sort(out.begin(), out.end());
    return out;
  }

  [[nodiscard]] DirListing list_dir(std::string_view prefix) override {
    check_prefix(prefix);
    DirListing out;
    std::filesystem::path dir = root_;
    if (!prefix.empty()) {
      dir /= std::filesystem::path(prefix.substr(0, prefix.size() - 1));
    }
    if (!std::filesystem::is_directory(dir)) {
      return out;
    }
    for (const auto& entry : std::filesystem::directory_iterator(dir)) {
      const std::string name = entry.path().filename().generic_string();
      if (entry.is_regular_file()) {
        out.keys.push_back(name);
      } else if (entry.is_directory()) {
        out.prefixes.push_back(name);
      }
    }
    std::sort(out.keys.begin(), out.keys.end());
    std::sort(out.prefixes.begin(), out.prefixes.end());
    return out;
  }

 private:
  static void check_prefix(std::string_view prefix) {
    if (!prefix.empty() && prefix.back() != '/') {
      throw error("store prefix must be empty or end with '/', got '" + std::string(prefix) + "'");
    }
  }

  [[nodiscard]] std::filesystem::path key_path(std::string_view key) const {
    if (key.empty()) {
      throw error("store key must not be empty");
    }
    std::filesystem::path path = root_;
    std::size_t start = 0;
    while (start <= key.size()) {
      const std::size_t slash = key.find('/', start);
      const std::size_t end = slash == std::string_view::npos ? key.size() : slash;
      const std::string_view segment = key.substr(start, end - start);
      // Reject anything that could escape the root directory.
      if (segment.empty() || segment == "." || segment == "..") {
        throw error("invalid store key '" + std::string(key) + "'");
      }
      path /= segment;
      if (slash == std::string_view::npos) {
        break;
      }
      start = slash + 1;
    }
    return path;
  }

  std::filesystem::path root_;
};

}  // namespace zarr

#endif  // LIBZARR_ADAPTERS_FILESYSTEM_STORE_HPP
