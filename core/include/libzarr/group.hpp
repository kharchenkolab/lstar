// SPDX-License-Identifier: MIT

#ifndef LIBZARR_GROUP_HPP
#define LIBZARR_GROUP_HPP

#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "libzarr/array.hpp"
#include "libzarr/metadata.hpp"
#include "libzarr/store.hpp"
#include "libzarr/types.hpp"
#include "libzarr/v2.hpp"
#include "libzarr/v3.hpp"

/// \file group.hpp
/// The Group API: hierarchy create/open/traverse. Groups are value-semantics
/// handles sharing the store; opening the root through consolidated metadata
/// (.zmetadata) is automatic when present, so remote stores pay one metadata
/// round-trip.

namespace zarr {

/// Immediate children of a group, by kind.
struct GroupChildren {
  /// Names of child arrays, sorted.
  std::vector<std::string> arrays;
  /// Names of child groups, sorted.
  std::vector<std::string> groups;
};

/// A Zarr group bound to a Store.
class Group {
 public:
  /// Creates a group at `path` ("" = store root) in the requested format,
  /// including any missing ancestor groups (writers that skip intermediate
  /// group documents are a known interop hazard). Existing group documents
  /// along the chain are left untouched.
  static Group create(std::shared_ptr<Store> store, const std::string& path = "",
                      ZarrFormat format = ZarrFormat::v2) {
    if (!store) {
      throw error("Group::create: null store");
    }
    detail::validate_path(path);
    write_group_chain(*store, path, format);
    return {std::move(store), path, json::object(), nullptr, format, OpenOptions{}};
  }

  /// Opens the group at `path`, probing v3 (zarr.json) first, then v2.
  /// Consolidated metadata — v2 .zmetadata at the store root, or the v3
  /// inline convention — is used when present and shared with every node
  /// opened through this group.
  static Group open(std::shared_ptr<Store> store, const std::string& path = "",
                    OpenOptions options = {}) {
    if (!store) {
      throw error("Group::open: null store");
    }
    detail::validate_path(path);

    // v3 probe.
    const std::string v3_key = v3::meta_key(path);
    if (const auto bytes = store->read(v3_key)) {
      const json doc = v2::parse_json(*bytes, v3_key);
      if (doc.is_object() && doc.value("node_type", "") == std::string("array")) {
        throw error("'" + path + "' is an array, not a group");
      }
      v3::GroupMeta meta = v3::parse_group_meta(doc, v3_key, options.lenient);
      std::shared_ptr<const json> consolidated;
      if (meta.consolidated) {
        // The inline map is node-path -> document; rekey by store key so
        // child opens look up uniformly for both formats.
        json by_key = json::object();
        for (const auto& item : meta.consolidated->items()) {
          const std::string child = path.empty() ? item.key() : path + "/" + item.key();
          by_key[v3::meta_key(child)] = item.value();
        }
        by_key[v3_key] = doc;
        consolidated = std::make_shared<const json>(std::move(by_key));
      }
      return {std::move(store),        path,           std::move(meta.attributes),
              std::move(consolidated), ZarrFormat::v3, options};
    }

    // v2 path, reading through .zmetadata at the store root when present.
    std::shared_ptr<const json> consolidated;
    if (const auto c = v2::read_consolidated(*store)) {
      consolidated = std::make_shared<const json>(*c);
    }
    return open_with(std::move(store), path, std::move(consolidated), options);
  }

  /// Node path within the store ("" = root).
  [[nodiscard]] const std::string& path() const { return path_; }

  /// User attributes (.zattrs).
  [[nodiscard]] const json& attributes() const { return attributes_; }

  /// Replaces the user attributes and persists them. For v3, the stored
  /// zarr.json is patched in place, preserving any extension members.
  void set_attributes(json attributes) {
    attributes_ = std::move(attributes);
    if (format_ == ZarrFormat::v3) {
      const std::string key = v3::meta_key(path_);
      const auto bytes = store_->read(key);
      if (!bytes) {
        throw error(key + ": metadata disappeared");
      }
      json doc = v2::parse_json(*bytes, key);
      if (attributes_.empty()) {
        doc.erase("attributes");
      } else {
        doc["attributes"] = attributes_;
      }
      store_->write(key, canonical_json_bytes(doc));
      return;
    }
    const std::string key = v2::meta_key(path_, v2::kAttrsSuffix);
    if (attributes_.empty()) {
      v2::erase_meta_key(*store_, key);  // canonical: no empty .zattrs documents
    } else {
      v2::write_meta_key(*store_, key, attributes_);
    }
  }

  /// Creates a child (possibly nested, e.g. "a/b") group in this group's
  /// format.
  Group create_group(const std::string& name) { return create(store_, child_path(name), format_); }

  /// Creates a child (possibly nested) array, writing ancestor groups first.
  /// The group's format governs; spec.format is ignored.
  Array create_array(const std::string& name, ArraySpec spec) {
    spec.format = format_;
    const std::string target = child_path(name);
    const std::size_t slash = target.rfind('/');
    if (slash != std::string::npos) {
      write_group_chain(*store_, target.substr(0, slash), format_);
    }
    return Array::create(store_, target, spec);
  }

  /// Opens a child (possibly nested) group.
  [[nodiscard]] Group open_group(const std::string& name) const {
    const std::string target = child_path(name);
    if (format_ == ZarrFormat::v3) {
      return open_v3_child_group(target);
    }
    return open_with(store_, target, consolidated_, options_);
  }

  /// Opens a child (possibly nested) array.
  [[nodiscard]] Array open_array(const std::string& name) const {
    return Array::open(store_, child_path(name), options_, consolidated_);
  }

  /// Lists immediate children, classified by their metadata documents.
  [[nodiscard]] GroupChildren children() const {
    GroupChildren out;
    const std::string prefix = path_.empty() ? "" : path_ + "/";
    for (const std::string& child : store_->list_dir(prefix).prefixes) {
      const std::string child_full = prefix + child;
      if (store_->exists(child_full + "/" + v2::kArraySuffix)) {
        out.arrays.push_back(child);
      } else if (store_->exists(child_full + "/" + v2::kGroupSuffix)) {
        out.groups.push_back(child);
      } else if (const auto doc = read_doc(v3::meta_key(child_full))) {
        // v3: the node kind lives inside zarr.json.
        const std::string node_type = doc->is_object() ? doc->value("node_type", "") : "";
        if (node_type == "array") {
          out.arrays.push_back(child);
        } else if (node_type == "group") {
          out.groups.push_back(child);
        }
      }
      // Other directories are not Zarr nodes; ignore them.
    }
    return out;
  }

 private:
  Group(std::shared_ptr<Store> store, std::string path, json attributes,
        std::shared_ptr<const json> consolidated, ZarrFormat format = ZarrFormat::v2,
        OpenOptions options = {})
      : store_(std::move(store)),
        path_(std::move(path)),
        attributes_(std::move(attributes)),
        consolidated_(std::move(consolidated)),
        format_(format),
        options_(options) {}

  [[nodiscard]] std::optional<json> read_doc(const std::string& key) const {
    if (consolidated_) {
      const auto it = consolidated_->find(key);
      if (it == consolidated_->end()) {
        return std::nullopt;
      }
      return *it;
    }
    const auto bytes = store_->read(key);
    if (!bytes) {
      return std::nullopt;
    }
    return v2::parse_json(*bytes, key);
  }

  [[nodiscard]] Group open_v3_child_group(const std::string& target) const {
    const std::string key = v3::meta_key(target);
    const auto doc = read_doc(key);
    if (!doc) {
      throw error("no group at '" + target + "' (" + key + " not found)");
    }
    if (doc->is_object() && doc->value("node_type", "") == std::string("array")) {
      throw error("'" + target + "' is an array, not a group");
    }
    v3::GroupMeta meta = v3::parse_group_meta(*doc, key, options_.lenient);
    // Children of this group keep reading through the root's consolidated map.
    return {store_, target, std::move(meta.attributes), consolidated_, ZarrFormat::v3, options_};
  }

  static Group open_with(std::shared_ptr<Store> store, const std::string& path,
                         std::shared_ptr<const json> consolidated, OpenOptions options) {
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

    const std::string group_key = v2::meta_key(path, v2::kGroupSuffix);
    const auto doc = read_doc(group_key);
    if (!doc) {
      if (store->exists(v2::meta_key(path, v2::kArraySuffix))) {
        throw error("'" + path + "' is an array, not a group");
      }
      throw error("no group at '" + path + "' (neither " + v3::meta_key(path) + " nor " +
                  group_key + " found)");
    }
    v2::check_group_meta(*doc, group_key);
    json attributes = json::object();
    if (const auto attrs = read_doc(v2::meta_key(path, v2::kAttrsSuffix))) {
      attributes = *attrs;
    }
    return {std::move(store),        path,           std::move(attributes),
            std::move(consolidated), ZarrFormat::v2, options};
  }

  /// Writes group metadata at `path` and every missing ancestor.
  static void write_group_chain(Store& store, const std::string& path, ZarrFormat format) {
    std::vector<std::string> chain;
    chain.emplace_back("");
    std::size_t start = 0;
    while (!path.empty()) {
      const std::size_t slash = path.find('/', start);
      if (slash == std::string::npos) {
        chain.push_back(path);
        break;
      }
      chain.push_back(path.substr(0, slash));
      start = slash + 1;
    }
    for (const std::string& node : chain) {
      if (store.exists(v2::meta_key(node, v2::kArraySuffix))) {
        throw error("'" + node + "' is an array; cannot create a group inside it");
      }
      if (format == ZarrFormat::v3) {
        const std::string key = v3::meta_key(node);
        if (const auto bytes = store.read(key)) {
          const json doc = v2::parse_json(*bytes, key);
          if (doc.is_object() && doc.value("node_type", "") == std::string("array")) {
            throw error("'" + node + "' is an array; cannot create a group inside it");
          }
          continue;  // existing group document: leave it (and its attributes) alone
        }
        store.write(key, canonical_json_bytes(v3::emit_group_meta(json::object())));
      } else {
        const std::string key = v2::meta_key(node, v2::kGroupSuffix);
        if (!store.exists(key)) {
          v2::write_meta_key(store, key, v2::group_meta_json());
        }
      }
    }
  }

  [[nodiscard]] std::string child_path(const std::string& name) const {
    detail::validate_path(name);
    if (name.empty()) {
      throw error("child name must not be empty");
    }
    return path_.empty() ? name : path_ + "/" + name;
  }

  std::shared_ptr<Store> store_;
  std::string path_;
  json attributes_;
  std::shared_ptr<const json> consolidated_;
  ZarrFormat format_ = ZarrFormat::v2;
  OpenOptions options_;
};

}  // namespace zarr

#endif  // LIBZARR_GROUP_HPP
