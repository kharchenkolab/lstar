// lstar._io — Emscripten/embind binding over libstar's Zarr I/O (the libzarr-backed reader).
//
// This is the SAME C++ core that reads L* stores for the R (cpp11) and Python (pybind11) packages,
// compiled to WebAssembly. It retires the parallel zarrita.js reimplementation: the browser/Node viewer
// reads a store through one recipe (libzarr), not a second one. The libzarr Store is SYNCHRONOUS, so the
// JS host supplies bytes through a synchronous callback -- in Node that callback reads the filesystem; in
// the browser the JS layer prefetches the keys `keysFor()` names into a cache the callback reads from
// (prefetch-then-sync-decode; no Asyncify). Reads both v2 (.zarray/.zgroup) and v3 (zarr.json).
#include <cstdint>
#include <memory>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

#include <emscripten/bind.h>
#include <emscripten/val.h>

#include "lstar/lstar.hpp"

using emscripten::val;
using zarr::ByteRange;
using zarr::Bytes;

static Bytes to_bytes(const val& v) { return emscripten::convertJSArrayToNumberVector<std::uint8_t>(v); }
static val to_u8(const Bytes& b) {  // copy into a JS-owned Uint8Array (a memory view would dangle)
    return val::global("Uint8Array").new_(val(emscripten::typed_memory_view(b.size(), b.data())));
}

// A read-only Store whose read()/read_range() delegate to a synchronous JS callback: read(key) returns a
// Uint8Array, or null/undefined for a missing object. Listing is unused (the reader navigates by name).
struct JsStore : zarr::Store {
    val jread;
    explicit JsStore(val r) : jread(std::move(r)) {}
    std::optional<Bytes> read(std::string_view k) override {
        val v = jread(std::string(k));
        if (v.isNull() || v.isUndefined()) return std::nullopt;
        return to_bytes(v);
    }
    std::optional<Bytes> read_range(std::string_view k, ByteRange r) override {
        auto whole = read(k);                          // whole-object read + slice (host callback is key-granular)
        if (!whole) return std::nullopt;
        const std::uint64_t n = whole->size();
        std::uint64_t off = 0, len = n;
        if (r.kind == ByteRange::Kind::slice) { off = r.offset; len = std::min<std::uint64_t>(r.length, n - std::min(off, n)); }
        else if (r.kind == ByteRange::Kind::suffix) { len = std::min<std::uint64_t>(r.length, n); off = n - len; }
        return Bytes(whole->begin() + off, whole->begin() + off + len);
    }
    std::optional<std::uint64_t> size(std::string_view k) override {
        auto b = read(k); return b ? std::optional<std::uint64_t>(b->size()) : std::nullopt;
    }
    bool exists(std::string_view k) override { return read(k).has_value(); }
    void write(std::string_view, Bytes) override { throw zarr::error("JsStore is read-only"); }
    void erase(std::string_view) override { throw zarr::error("JsStore is read-only"); }
    std::vector<std::string> list_prefix(std::string_view) override { throw zarr::error("JsStore: no listing"); }
    zarr::DirListing list_dir(std::string_view) override { throw zarr::error("JsStore: no listing"); }
};

// enumerate every chunk key of an array's shape/chunk grid, in the array's own format + separator, using
// libzarr's own chunk_key -- so JS learns exactly what to prefetch without reimplementing key math.
static void enumerate_chunks(const zarr::ArrayMeta& m, std::vector<std::string>& out) {
    // The store objects to prefetch are the STORE-key grid, which for a sharded array is the OUTERMOST
    // shard grid -- NOT the inner-chunk grid. ArrayMeta::chunk_shape is always the innermost (true) chunk
    // shape; a sharded array packs many inner chunks into each shard object. Enumerating inner keys would
    // name c/0..c/N when only the M<=N shard objects exist -> N-M spurious 404s per array (decode is
    // correctness-safe -- it reads only the shards -- but the over-fetch defeats sharding's fewer-objects
    // goal). Use the outermost shard_shape so keysFor names exactly the objects that exist.
    const std::vector<std::uint64_t>& store_chunk =
        m.shard_levels.empty() ? m.chunk_shape : m.shard_levels.front().shard_shape;
    std::vector<std::uint64_t> grid(m.shape.size());
    std::uint64_t total = 1;
    for (size_t d = 0; d < m.shape.size(); ++d) {
        grid[d] = store_chunk[d] ? (m.shape[d] + store_chunk[d] - 1) / store_chunk[d] : 0;
        total *= grid[d];
    }
    if (m.shape.empty()) { out.push_back(m.format == zarr::ZarrFormat::v3 ? "c" : "0"); return; }  // 0-d
    std::vector<std::uint64_t> idx(m.shape.size(), 0);
    for (std::uint64_t c = 0; c < total; ++c) {
        out.push_back(m.format == zarr::ZarrFormat::v3 ? zarr::v3::chunk_key(idx, m.dimension_separator)
                                                       : zarr::v2::chunk_key(idx, m.dimension_separator));
        for (size_t d = idx.size(); d-- > 0;) { if (++idx[d] < grid[d]) break; idx[d] = 0; }
    }
}

// The reader lstar's JS layer drives: one Store, many name-addressed reads. Mirrors lstar.hpp's own
// per-node open (Group::open / Array::open both probe v3 then v2), so it reads exactly what R/Python do.
class Reader {
 public:
    explicit Reader(val read_cb) : store_(std::make_shared<JsStore>(std::move(read_cb))) {}

    // group user-attributes (the L* manifest lives at root under "lstar"); "" = root.
    std::string groupAttrs(const std::string& path) {
        return zarr::Group::open(store_, path).attributes().dump();
    }
    std::string rootAttrs() { return groupAttrs(""); }

    // an array's shape from its metadata only (no chunk reads) -> number[]. Used for axis lengths.
    val shape(const std::string& path) {
        zarr::Array a = zarr::Array::open(store_, path);
        val out = val::array();
        for (size_t i = 0; i < a.meta().shape.size(); ++i) out.set(i, (double)a.meta().shape[i]);
        return out;
    }

    // Pure-sync descriptor the JS byte-range fast path needs, format-agnostically (v2 or v3): element
    // size, chunk shape, and whether a chunk stores raw contiguous element bytes (no bytes->bytes
    // compressor), so a sub-range of the array maps 1:1 onto a store byte range. JS owns the fetch; this
    // owns the Zarr metadata interpretation.
    val arrayInfo(const std::string& path) {
        zarr::Array a = zarr::Array::open(store_, path);
        bool raw = true;                                // raw iff no compression codec in the chain
        for (const auto& c : a.meta().codecs)
            if (c.name == "gzip" || c.name == "zlib" || c.name == "blosc" || c.name == "zstd") raw = false;
        val out = val::object();
        out.set("dtype", zarr::v2::emit_data_type(a.meta().dtype, /*big_endian=*/false));
        out.set("itemsize", (double)a.meta().dtype.itemsize);
        val cs = val::array();
        for (size_t i = 0; i < a.meta().chunk_shape.size(); ++i) cs.set(i, (double)a.meta().chunk_shape[i]);
        out.set("chunkShape", cs);
        out.set("uncompressed", raw);
        // A sharded array's inner chunks live INSIDE shard objects, so the plain chunk-key byte-range
        // fast path does not apply; the reader falls back to a whole-array read (correct) unless the JS
        // shard-resolve path is used. (sharding_indexed is lowered into shard_levels, not codecs.)
        out.set("sharded", !a.meta().shard_levels.empty());
        return out;
    }

    // The store key of one chunk, in the array's own chunk-key encoding — v2 "0" / "i.j", v3 default
    // "c/0" / "c/i/j". The drift-prone bit: computed by libzarr, never reimplemented in JS.
    std::string chunkKey(const std::string& path, val idx_js) {
        zarr::Array a = zarr::Array::open(store_, path);
        std::vector<int> t = emscripten::convertJSArrayToNumberVector<int>(idx_js);   // JS Numbers (chunk idx < 2^31)
        std::vector<std::uint64_t> idx(t.begin(), t.end());
        return a.meta().format == zarr::ZarrFormat::v3
                   ? zarr::v3::chunk_key(idx, a.meta().dimension_separator)
                   : zarr::v2::chunk_key(idx, a.meta().dimension_separator);
    }

    // Sharded byte-range resolve (v3 sharding), phase 1: where inner chunk `idx` lives. Returns the
    // leaf shard-object key (JS prepends the array path, exactly like chunkKey), the chunk's slot in
    // that shard, and the shard index's layout (encoded size + end/start) so JS knows which bytes of
    // the shard to fetch for the index. Metadata-only; no chunk read. libzarr owns the shard math.
    val shardLocate(const std::string& path, val idx_js) {
        zarr::Array a = zarr::Array::open(store_, path);
        std::vector<int> t = emscripten::convertJSArrayToNumberVector<int>(idx_js);
        std::vector<std::uint64_t> idx(t.begin(), t.end());
        // libzarr's public shard-resolution façade (v1.0). path="" -> the shard key is leaf-relative
        // (JS prepends the array path, same convention as chunkKey). Pure: no store read.
        const auto p = zarr::shard::place(a.meta(), /*path=*/"", idx, /*level=*/0);
        val out = val::object();
        out.set("shardKey", p.shard_key);
        out.set("intra", (double)p.slot);
        out.set("indexSize", (double)p.index_size);
        out.set("indexAtEnd", p.index_at_end);
        return out;
    }

    // Sharded byte-range resolve, phase 2: decode the shard index bytes JS fetched (per shardLocate)
    // and return inner chunk `intra`'s [offset, nbytes) within the shard, or missing == a fill chunk.
    // JS then range-reads exactly the chunk's bytes off the shard object. No store read here.
    val shardEntry(const std::string& path, val index_js, double intra) {
        zarr::Array a = zarr::Array::open(store_, path);
        const Bytes idx = to_bytes(index_js);
        const auto e = zarr::shard::extent(a.meta(), idx, (std::uint64_t)intra, /*level=*/0);
        val out = val::object();
        out.set("offset", (double)e.offset);
        out.set("nbytes", (double)e.nbytes);
        out.set("missing", e.missing);
        return out;
    }

    // a whole array -> {dtype: "<f4"|..., shape: number[], bytes: Uint8Array (C-order, native)}.
    val array(const std::string& path) {
        zarr::Array a = zarr::Array::open(store_, path);
        Bytes b(a.nbytes());
        if (!b.empty()) a.read(b.data(), b.size());
        val out = val::object();
        out.set("dtype", zarr::v2::emit_data_type(a.meta().dtype, /*big_endian=*/false));
        val shape = val::array();
        for (size_t i = 0; i < a.meta().shape.size(); ++i) shape.set(i, (double)a.meta().shape[i]);
        out.set("shape", shape);
        out.set("bytes", to_u8(b));
        return out;
    }

    // the store keys required to read array(path): its metadata document + every chunk key. JS prefetches
    // these (async) into the callback's cache before calling array(); the metadata comes from expanding the
    // store's consolidated metadata, the chunks are listed here in libzarr's own key scheme.
    val keysFor(const std::string& path) {
        zarr::Array a = zarr::Array::open(store_, path);
        std::vector<std::string> keys;
        keys.push_back(path + (a.meta().format == zarr::ZarrFormat::v3 ? "/zarr.json" : "/.zarray"));
        std::vector<std::string> chunks;
        enumerate_chunks(a.meta(), chunks);
        for (auto& c : chunks) keys.push_back(path + "/" + c);
        val out = val::array();
        for (size_t i = 0; i < keys.size(); ++i) out.set(i, keys[i]);
        return out;
    }

    std::string version() { return "lstar-io 0.1.0"; }

 private:
    std::shared_ptr<JsStore> store_;
};

EMSCRIPTEN_BINDINGS(lstar_io) {
    emscripten::class_<Reader>("Reader")
        .constructor<val>()
        .function("rootAttrs", &Reader::rootAttrs)
        .function("groupAttrs", &Reader::groupAttrs)
        .function("shape", &Reader::shape)
        .function("arrayInfo", &Reader::arrayInfo)
        .function("chunkKey", &Reader::chunkKey)
        .function("shardLocate", &Reader::shardLocate)
        .function("shardEntry", &Reader::shardEntry)
        .function("array", &Reader::array)
        .function("keysFor", &Reader::keysFor)
        .function("version", &Reader::version);
}
