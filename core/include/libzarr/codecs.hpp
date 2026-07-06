// SPDX-License-Identifier: MIT

#ifndef LIBZARR_CODECS_HPP
#define LIBZARR_CODECS_HPP

#include <cstdint>
#include <limits>
#include <optional>
#include <string>
#include <vector>

#include "libzarr/detail/common.hpp"
#include "libzarr/metadata.hpp"
#include "libzarr/types.hpp"

#ifdef LIBZARR_HAS_ZLIB
#include "libzarr/codecs_gzip.hpp"
#endif
#ifdef LIBZARR_HAS_BLOSC
#include "libzarr/codecs_blosc.hpp"
#endif
#ifdef LIBZARR_HAS_ZSTD
#include "libzarr/codecs_zstd.hpp"
#endif

/// \file codecs.hpp
/// The codec pipeline: a chain of CodecSpec resolved once per array into an
/// executable encode/decode plan. The chain is partitioned per the v3 model —
/// array->array codecs, exactly one array->bytes codec ("bytes"), then
/// bytes->bytes codecs — and no-op stages are elided at resolve time.

namespace zarr {

/// A codec chain resolved against one array's dtype and chunk shape.
/// encode() maps a native-order, C-layout full chunk to stored bytes;
/// decode() is the exact inverse. Both are whole-buffer and value-based.
class CodecPipeline {
 public:
  /// Validates `meta.codecs` (partition order, exactly one "bytes" stage,
  /// known names, well-formed configurations) and builds the plan. Codecs
  /// compiled out of this build (e.g. gzip without LIBZARR_HAS_ZLIB) fail
  /// here with a precise error, never at link time.
  [[nodiscard]] static CodecPipeline resolve(const ArrayMeta& meta) {
    CodecPipeline p;
    p.chunk_shape_ = meta.chunk_shape;
    p.itemsize_ = meta.dtype.itemsize;
    const std::uint64_t count = meta.chunk_element_count();
    if (meta.dtype.itemsize != 0 &&
        count > std::numeric_limits<std::uint64_t>::max() / meta.dtype.itemsize) {
      throw error("chunk byte size overflows uint64");
    }
    p.chunk_bytes_ = count * meta.dtype.itemsize;
    // v3 core: byte order applies per complex component (two floats).
    p.swap_width_ = is_complex(meta.dtype.kind) ? meta.dtype.itemsize / 2 : meta.dtype.itemsize;

    // v3 core: array->array*, then exactly one array->bytes, then bytes->bytes*.
    bool past_bytes = false;
    bool have_bytes = false;
    for (const CodecSpec& codec : meta.codecs) {
      if (codec.name == "transpose") {
        if (past_bytes) {
          throw error("codec 'transpose' must precede the 'bytes' codec");
        }
        p.set_transpose(codec);
      } else if (codec.name == "bytes") {
        if (have_bytes) {
          throw error("codec chain has more than one 'bytes' codec");
        }
        have_bytes = true;
        past_bytes = true;
        p.set_byte_order(codec);
      } else if (codec.name == "sharding_indexed") {
        // Sharding is not executed as a codec: metadata parsing lowers it
        // into ArrayMeta::shard_levels and Array wraps the store instead.
        throw error("codec 'sharding_indexed' must be lowered into shard levels, not resolved");
      } else {
        if (!past_bytes) {
          throw error("codec '" + codec.name + "' must follow the 'bytes' codec");
        }
        p.add_byte_stage(codec, meta);
      }
    }
    if (!have_bytes) {
      throw error("codec chain is missing the 'bytes' (array->bytes) codec");
    }
    p.compute_expected_sizes();
    return p;
  }

  /// Size of a decoded full chunk in bytes.
  [[nodiscard]] std::uint64_t decoded_chunk_bytes() const { return chunk_bytes_; }

  /// True when stored bytes equal decoded bytes element-for-element (no
  /// transpose, no byteswap, no bytes->bytes stage) — the precondition for
  /// byte-range sub-chunk reads without post-processing.
  [[nodiscard]] bool is_identity() const {
    return !transpose_order_ && !byteswap_ && byte_stages_.empty();
  }

  /// True when the stored byte at a given linear element offset is
  /// independent of the rest of the chunk (no bytes->bytes stage) —
  /// byte-range reads work, possibly with a byteswap. Transposed layouts do
  /// not qualify.
  [[nodiscard]] bool supports_partial_read() const {
    return !transpose_order_ && byte_stages_.empty();
  }

  /// Encodes a native, C-order full chunk. The buffer must be exactly
  /// decoded_chunk_bytes() long.
  [[nodiscard]] Bytes encode(Bytes chunk) const {
    if (chunk.size() != chunk_bytes_) {
      throw error("encode: chunk buffer is " + std::to_string(chunk.size()) + " bytes, expected " +
                  std::to_string(chunk_bytes_));
    }
    if (transpose_order_) {
      // Write support for transposed layouts (v2 order:"F") is deliberately
      // absent: we emit canonical C-order arrays only.
      throw error("writing to a transposed (order:'F') array is not supported");
    }
    if (byteswap_) {
      detail::byteswap_inplace(chunk.data(), chunk.size() / swap_width_, swap_width_);
    }
    for (const ByteStage& stage : byte_stages_) {
      chunk = encode_stage(stage, std::move(chunk));
    }
    return chunk;
  }

  /// Post-processes a partial (byte-range) chunk read; valid only when
  /// supports_partial_read(). The sole possible transform is a byteswap.
  [[nodiscard]] Bytes decode_range(Bytes raw) const {
    assert(supports_partial_read());
    if (byteswap_) {
      detail::byteswap_inplace(raw.data(), raw.size() / swap_width_, swap_width_);
    }
    return raw;
  }

  /// Decodes stored chunk bytes to a native, C-order full chunk of exactly
  /// decoded_chunk_bytes() bytes; anything inconsistent is a zarr::error.
  [[nodiscard]] Bytes decode(Bytes stored) const {
    for (std::size_t i = byte_stages_.size(); i-- > 0;) {
      stored = decode_stage(byte_stages_[i], std::move(stored), decode_expected_[i]);
    }
    if (stored.size() != chunk_bytes_) {
      throw error("decode: chunk is " + std::to_string(stored.size()) + " bytes, expected " +
                  std::to_string(chunk_bytes_));
    }
    if (byteswap_) {
      detail::byteswap_inplace(stored.data(), stored.size() / swap_width_, swap_width_);
    }
    if (transpose_order_) {
      Bytes out(stored.size());
      detail::gather_strided(stored.data(), gather_strides_, out.data(), chunk_shape_, itemsize_);
      return out;
    }
    return stored;
  }

 private:
  struct ByteStage {
    enum class Kind : std::uint8_t {
      deflate,  // gzip/zlib framing
      crc32c,
      blosc,
      zstd,
      shuffle,  // v2 numcodecs shuffle filter
    };
    Kind kind = Kind::deflate;
    // deflate
    int level = 5;
    bool gzip_framing = true;
    // blosc (encode-side parameters; decode is self-describing)
    std::string blosc_cname = "lz4";
    int blosc_clevel = 5;
    int blosc_shuffle = 1;
    std::uint32_t blosc_typesize = 1;
    std::uint64_t blosc_blocksize = 0;
    // zstd
    int zstd_level = 0;
    bool zstd_checksum = false;
    // shuffle
    std::uint32_t shuffle_elementsize = 1;
  };

  [[nodiscard]] static Bytes encode_stage(const ByteStage& stage, Bytes data) {
    switch (stage.kind) {
      case ByteStage::Kind::deflate:
#ifdef LIBZARR_HAS_ZLIB
        return detail::deflate_bytes(data, stage.level, stage.gzip_framing, "encode");
#else
        throw error("codec requires zlib but LIBZARR_HAS_ZLIB is not defined");
#endif
      case ByteStage::Kind::crc32c: {
        // v3 crc32c codec: little-endian CRC-32C of the payload, appended.
        const std::uint32_t checksum = detail::crc32c(data.data(), data.size());
        for (int i = 0; i < 4; ++i) {
          data.push_back(static_cast<std::uint8_t>(checksum >> (8 * i)));
        }
        return data;
      }
      case ByteStage::Kind::blosc:
#ifdef LIBZARR_HAS_BLOSC
      {
        detail::BloscParams params;
        params.cname = stage.blosc_cname;
        params.clevel = stage.blosc_clevel;
        params.shuffle = stage.blosc_shuffle;
        params.typesize = stage.blosc_typesize;
        params.blocksize = detail::checked_size(stage.blosc_blocksize, "blosc blocksize");
        return detail::blosc_compress_bytes(data, params, "encode");
      }
#else
        throw error("codec requires blosc but LIBZARR_HAS_BLOSC is not defined");
#endif
      case ByteStage::Kind::zstd:
#ifdef LIBZARR_HAS_ZSTD
        return detail::zstd_compress_bytes(data, stage.zstd_level, stage.zstd_checksum, "encode");
#else
        throw error("codec requires zstd but LIBZARR_HAS_ZSTD is not defined");
#endif
      case ByteStage::Kind::shuffle:
        return detail::shuffle_bytes(data, stage.shuffle_elementsize);
    }
    return data;  // unreachable
  }

  [[nodiscard]] static Bytes decode_stage(const ByteStage& stage, Bytes data,
                                          [[maybe_unused]] std::optional<std::uint64_t> expected) {
    switch (stage.kind) {
      case ByteStage::Kind::deflate:
#ifdef LIBZARR_HAS_ZLIB
        return detail::inflate_bytes(data, expected, "decode");
#else
        throw error("codec requires zlib but LIBZARR_HAS_ZLIB is not defined");
#endif
      case ByteStage::Kind::crc32c: {
        if (data.size() < 4) {
          throw error("decode: crc32c codec needs at least 4 bytes");
        }
        const std::size_t payload = data.size() - 4;
        std::uint32_t stored_crc = 0;
        for (int i = 3; i >= 0; --i) {
          stored_crc = (stored_crc << 8U) | data[payload + static_cast<std::size_t>(i)];
        }
        if (detail::crc32c(data.data(), payload) != stored_crc) {
          throw error("decode: crc32c checksum mismatch (corrupt chunk)");
        }
        data.resize(payload);
        return data;
      }
      case ByteStage::Kind::blosc:
#ifdef LIBZARR_HAS_BLOSC
        return detail::blosc_decompress_bytes(data, expected, "decode");
#else
        throw error("codec requires blosc but LIBZARR_HAS_BLOSC is not defined");
#endif
      case ByteStage::Kind::zstd:
#ifdef LIBZARR_HAS_ZSTD
        return detail::zstd_decompress_bytes(data, expected, "decode");
#else
        throw error("codec requires zstd but LIBZARR_HAS_ZSTD is not defined");
#endif
      case ByteStage::Kind::shuffle:
        return detail::unshuffle_bytes(data, stage.shuffle_elementsize);
    }
    return data;  // unreachable
  }

  void set_transpose(const CodecSpec& codec) {
    const auto it = codec.configuration.find("order");
    if (it == codec.configuration.end() || !it->is_array()) {
      throw error("codec 'transpose' requires an 'order' array");
    }
    const std::size_t rank = chunk_shape_.size();
    std::vector<std::uint32_t> order;
    std::vector<bool> seen(rank, false);
    for (const json& v : *it) {
      const std::uint64_t dim = detail::json_to_uint64(v, "codec 'transpose': 'order'");
      if (dim >= rank) {
        throw error("codec 'transpose': 'order' must be a permutation of 0.." +
                    std::to_string(rank == 0 ? 0 : rank - 1));
      }
      const auto d = static_cast<std::uint32_t>(dim);
      if (seen[d]) {
        throw error("codec 'transpose': repeated dimension " + std::to_string(d) + " in 'order'");
      }
      seen[d] = true;
      order.push_back(d);
    }
    if (order.size() != rank) {
      throw error("codec 'transpose': 'order' has " + std::to_string(order.size()) +
                  " entries for a rank-" + std::to_string(rank) + " array");
    }
    bool identity = true;
    for (std::size_t i = 0; i < rank; ++i) {
      identity = identity && order[i] == i;
    }
    if (identity) {
      return;  // no-op elision
    }
    transpose_order_ = order;
    // Stored (encoded) dimension i holds source dimension order[i]; build the
    // per-source-dimension byte strides used to gather back to C order.
    std::vector<std::uint64_t> stored_shape(rank);
    for (std::size_t i = 0; i < rank; ++i) {
      stored_shape[i] = chunk_shape_[order[i]];
    }
    const std::vector<std::uint64_t> stored_strides =
        detail::c_strides_bytes(stored_shape, itemsize_);
    gather_strides_.assign(rank, 0);
    for (std::size_t i = 0; i < rank; ++i) {
      gather_strides_[order[i]] = stored_strides[i];
    }
  }

  void set_byte_order(const CodecSpec& codec) {
    std::string endian = "little";
    const auto it = codec.configuration.find("endian");
    if (it != codec.configuration.end()) {
      if (!it->is_string()) {
        throw error("codec 'bytes': 'endian' must be a string");
      }
      endian = it->get<std::string>();
    }
    if (endian != "little" && endian != "big") {
      throw error("codec 'bytes': unknown endian '" + endian + "'");
    }
    const bool stored_little = endian == "little";
    byteswap_ = swap_width_ > 1 && stored_little != detail::host_is_little_endian();
  }

  void add_byte_stage(const CodecSpec& codec, const ArrayMeta& meta) {
    if (codec.name == "gzip" || codec.name == "zlib") {
      add_deflate(codec);
    } else if (codec.name == "crc32c") {
      ByteStage stage;
      stage.kind = ByteStage::Kind::crc32c;
      byte_stages_.push_back(stage);
    } else if (codec.name == "blosc") {
      add_blosc(codec, meta);
    } else if (codec.name == "zstd") {
      add_zstd(codec);
    } else if (codec.name == "shuffle") {
      add_shuffle(codec, meta);
    } else {
      // sharding_indexed never reaches here: resolve()'s dispatch intercepts
      // it (it is lowered into shard levels at parse time, not resolved).
      throw error("unknown codec '" + codec.name + "'");
    }
  }

  void add_deflate(const CodecSpec& codec) {
    ByteStage stage;
    stage.kind = ByteStage::Kind::deflate;
    stage.gzip_framing = codec.name == "gzip";
    const auto it = codec.configuration.find("level");
    if (it != codec.configuration.end()) {
      if (!it->is_number_integer() || it->get<std::int64_t>() < 0 || it->get<std::int64_t>() > 9) {
        throw error("codec '" + codec.name + "': 'level' must be an integer in 0..9");
      }
      stage.level = it->get<int>();
    }
#ifndef LIBZARR_HAS_ZLIB
    throw error("codec '" + codec.name +
                "' is not built into this libzarr (compile with LIBZARR_HAS_ZLIB and link zlib)");
#endif
    byte_stages_.push_back(stage);
  }

  void add_zstd(const CodecSpec& codec) {
    ByteStage stage;
    stage.kind = ByteStage::Kind::zstd;
    const json config = codec.configuration.is_object() ? codec.configuration : json::object();
    const json level = config.value("level", json(std::int64_t{0}));
    if (!level.is_number_integer()) {
      throw error("codec 'zstd': 'level' must be an integer");
    }
    stage.zstd_level = level.get<int>();
    const json checksum = config.value("checksum", json(false));
    if (!checksum.is_boolean()) {
      throw error("codec 'zstd': 'checksum' must be a boolean");
    }
    stage.zstd_checksum = checksum.get<bool>();
#ifndef LIBZARR_HAS_ZSTD
    throw error(
        "codec 'zstd' is not built into this libzarr (compile with LIBZARR_HAS_ZSTD and link "
        "zstd)");
#endif
    byte_stages_.push_back(stage);
  }

  void add_shuffle(const CodecSpec& codec, const ArrayMeta& meta) {
    ByteStage stage;
    stage.kind = ByteStage::Kind::shuffle;
    const json config = codec.configuration.is_object() ? codec.configuration : json::object();
    const std::int64_t elementsize = config.value("elementsize", std::int64_t{0});
    if (elementsize < 0 || elementsize > 0xFFFF) {
      throw error("filter 'shuffle': invalid elementsize " + std::to_string(elementsize));
    }
    // NCZarr writes elementsize 0 for "the dtype's item size"; the actual
    // stored bytes are shuffled with the item size.
    stage.shuffle_elementsize =
        elementsize == 0 ? meta.dtype.itemsize : static_cast<std::uint32_t>(elementsize);
    if (stage.shuffle_elementsize == 0) {
      throw error("filter 'shuffle': element size cannot be zero");
    }
    byte_stages_.push_back(stage);
  }

  void add_blosc(const CodecSpec& codec, const ArrayMeta& meta) {
    ByteStage stage;
    stage.kind = ByteStage::Kind::blosc;
    // A default-constructed CodecSpec has a *null* configuration, on which
    // json::value() throws; normalize to an empty object.
    const json config = codec.configuration.is_object() ? codec.configuration : json::object();
    stage.blosc_cname = config.value("cname", "lz4");
    const std::int64_t clevel = config.value("clevel", std::int64_t{5});
    if (clevel < 0 || clevel > 9) {
      throw error("codec 'blosc': 'clevel' must be in 0..9");
    }
    stage.blosc_clevel = static_cast<int>(clevel);
    const json shuffle = config.value("shuffle", json("shuffle"));
    if (shuffle == "noshuffle") {
      stage.blosc_shuffle = 0;
    } else if (shuffle == "shuffle") {
      stage.blosc_shuffle = 1;
    } else if (shuffle == "bitshuffle") {
      stage.blosc_shuffle = 2;
    } else if (shuffle.is_number_integer() && shuffle.get<std::int64_t>() >= -1 &&
               shuffle.get<std::int64_t>() <= 2) {
      // v2 numcodecs uses numeric shuffle; -1 = automatic (bitshuffle for
      // 1-byte types, else byte shuffle) — read tolerance.
      const auto n = shuffle.get<std::int64_t>();
      if (n == -1) {
        stage.blosc_shuffle = meta.dtype.itemsize == 1 ? 2 : 1;
      } else {
        stage.blosc_shuffle = static_cast<int>(n);
      }
    } else {
      throw error("codec 'blosc': unknown shuffle " + shuffle.dump());
    }
    stage.blosc_typesize =
        static_cast<std::uint32_t>(config.value("typesize", std::int64_t{meta.dtype.itemsize}));
    stage.blosc_blocksize = static_cast<std::uint64_t>(config.value("blocksize", std::int64_t{0}));
    bool known_cname = false;
    for (const char* name : {"blosclz", "lz4", "lz4hc", "snappy", "zlib", "zstd"}) {
      known_cname = known_cname || stage.blosc_cname == name;
    }
    if (!known_cname) {
      throw error("codec 'blosc': unknown cname '" + stage.blosc_cname + "'");
    }
#ifndef LIBZARR_HAS_BLOSC
    throw error(
        "codec 'blosc' is not built into this libzarr (compile with LIBZARR_HAS_BLOSC and link "
        "c-blosc)");
#endif
    byte_stages_.push_back(stage);
  }

  /// Precomputes, per stage, the byte size its decode must produce: known for
  /// the segment of size-preserving/size-shifting stages nearest the array
  /// (crc32c adds exactly 4), unknowable outside the first compressor.
  void compute_expected_sizes() {
    decode_expected_.assign(byte_stages_.size(), std::nullopt);
    std::optional<std::uint64_t> size = chunk_bytes_;
    for (std::size_t i = 0; i < byte_stages_.size(); ++i) {
      decode_expected_[i] = size;  // decode of stage i must yield its encode input
      if (!size) {
        continue;
      }
      if (byte_stages_[i].kind == ByteStage::Kind::crc32c) {
        size = *size + 4;
      } else if (byte_stages_[i].kind == ByteStage::Kind::shuffle) {
        // size-preserving
      } else {
        size = std::nullopt;  // compressed size is unknowable
      }
    }
  }

  std::vector<std::uint64_t> chunk_shape_;
  std::uint32_t itemsize_ = 1;
  std::uint64_t chunk_bytes_ = 0;
  std::uint32_t swap_width_ = 1;
  bool byteswap_ = false;
  std::optional<std::vector<std::uint32_t>> transpose_order_;
  std::vector<std::uint64_t> gather_strides_;
  std::vector<ByteStage> byte_stages_;
  std::vector<std::optional<std::uint64_t>> decode_expected_;
};

}  // namespace zarr

#endif  // LIBZARR_CODECS_HPP
