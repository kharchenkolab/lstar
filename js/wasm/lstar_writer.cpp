// lstar._writer — Emscripten/embind binding for the WRITE-side pure functions from libzarr: chunk codec
// ENCODE (gzip/zstd) and shard-object assembly (shard::pack). The write twin of lstar_reader.cpp's
// shardLocate/shardEntry: JS owns chunking + store writes + metadata; the drift-prone bytes (codec
// encode, and the shard index + crc32c) are libzarr compiled to WASM, so a JS-written store is
// byte-identical to what the C++/Python writers produce. Built with FULL zstd (LIBZARR_HAS_ZSTD, NO
// DECODE_ONLY) so encode covers zstd -- separate from the decode-only reader module (lstar_io).
#include <cstdint>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include <emscripten/bind.h>
#include <emscripten/val.h>

#include "lstar/lstar.hpp"

using emscripten::val;
using zarr::Bytes;

static Bytes to_bytes(const val& v) { return emscripten::convertJSArrayToNumberVector<std::uint8_t>(v); }
static val to_u8(const Bytes& b) {  // copy into a JS-owned Uint8Array (a memory view would dangle)
    return val::global("Uint8Array").new_(val(emscripten::typed_memory_view(b.size(), b.data())));
}

// Encode one chunk's raw little-endian bytes through a v3 codec chain [bytes, (gzip|zstd)?]. The chunk is
// treated as opaque uint8 (its bytes are already the on-disk LE layout), so the `bytes` codec is identity
// and only the compressor transforms it -- and it's the SAME libzarr codec code the reader decodes with,
// so the output is byte-consistent across surfaces. compressor: "none" | "gzip" | "zstd".
static val encodeChunk(const val& rawJS, const std::string& compressor, int level) {
    Bytes raw = to_bytes(rawJS);
    const std::uint64_t n = raw.empty() ? 1 : raw.size();
    zarr::ArrayMeta meta;
    meta.shape = {n};
    meta.chunk_shape = {n};
    meta.dtype = zarr::DataType::of(zarr::DType::uint8);
    meta.codecs.push_back(zarr::CodecSpec{"bytes", {{"endian", "little"}}});
    if (compressor == "gzip") meta.codecs.push_back(zarr::CodecSpec{"gzip", {{"level", level}}});
    else if (compressor == "zstd") meta.codecs.push_back(zarr::CodecSpec{"zstd", {{"level", level}}});
    else if (compressor != "none") throw zarr::error("encodeChunk: unknown compressor '" + compressor + "'");
    return to_u8(zarr::CodecPipeline::resolve(meta).encode(std::move(raw)));
}

// Assemble one shard object from its slot-ordered, already-ENCODED inner chunks (a null/undefined entry is
// a fill/absent slot). Binds zarr::shard::pack via a minimal 1-D meta whose single shard level has
// `entries.length` slots and the standard index codecs ([bytes, crc32c]). pack is dimension-agnostic given
// the slot order, so the bytes are identical to the real multi-dim array's shard object. Returns the shard
// bytes (empty when every slot is absent -> JS skips writing that shard key). `indexAtEnd` must match what
// the array's metadata declares (lstar emits index_location:"end").
static val packShard(const val& entriesJS, bool indexAtEnd) {
    const unsigned n = entriesJS["length"].as<unsigned>();
    std::vector<std::optional<Bytes>> entries(n);
    for (unsigned i = 0; i < n; ++i) {
        const val e = entriesJS[i];
        if (!e.isNull() && !e.isUndefined()) entries[i] = to_bytes(e);
    }
    const std::uint64_t slots = n ? n : 1;
    zarr::ArrayMeta meta;
    meta.shape = {slots};
    meta.chunk_shape = {1};
    meta.dtype = zarr::DataType::of(zarr::DType::uint8);
    zarr::ShardLevel lvl;
    lvl.shard_shape = {slots};
    lvl.index_codecs = {zarr::CodecSpec{"bytes", {{"endian", "little"}}}, zarr::CodecSpec{"crc32c", {}}};
    lvl.index_at_end = indexAtEnd;
    meta.shard_levels.push_back(std::move(lvl));
    return to_u8(zarr::shard::pack(meta, entries, 0));
}

EMSCRIPTEN_BINDINGS(lstar_writer) {
    emscripten::function("encodeChunk", &encodeChunk);
    emscripten::function("packShard", &packShard);
}
