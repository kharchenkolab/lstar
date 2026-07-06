// lstar — C++ core (libstar): the L* model + a Zarr v2 reader/writer.
//
// Reads multi-chunk arrays over an arbitrary chunk grid, little-endian, with optional gzip/zlib
// chunk compression (when built with zlib; define LSTAR_HAVE_ZLIB). Edge chunks are full-size and
// fill-padded per the v2 spec; missing chunks read as fill_value 0. The writer emits single-chunk
// arrays plus a consolidated .zmetadata. Blosc and Zarr v3 sharding are still planned (see
// misc/plan1.md). Written fresh (MIT); pagoda2's IO is GPL-3 and was used only as a schema reference.
#pragma once

#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

// nlohmann/json instantiates std::char_traits<unsigned char> (via its
// std::basic_string<unsigned char> output adapters). libc++ on Xcode 26.5+
// deprecates char_traits<T> for non-standard T, which -Werror-style CI turns
// into a build failure. Silence it just for this third-party header; our own
// code names no such instantiation. (GCC-form pragma is honored by Clang too.)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#include "nlohmann/json.hpp"
#pragma GCC diagnostic pop

// libzarr: the header-only Zarr v2+v3 core that backs lstar's array I/O (read_array/read_array_range
// today; the writer + consolidated orchestration next). Shares this file's vendored nlohmann/json.
#include <libzarr/libzarr.hpp>
#include <libzarr/adapters/filesystem_store.hpp>

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef LSTAR_HAVE_ZLIB
#include <zlib.h>
#endif

namespace lstar {

namespace fs = std::filesystem;
using json = nlohmann::json;

// ---------------------------------------------------------------- dtypes ----

// Size in bytes of a zarr numeric dtype string ("<f8", "<i4", "|u1", ...).
inline size_t dtype_size(const std::string& z) {
    char c = z.empty() ? '0' : z.back();
    if (c >= '1' && c <= '9') return static_cast<size_t>(c - '0');
    throw std::runtime_error("unsupported zarr dtype: " + z);
}

// A raw N-D numeric array: zarr dtype + shape + bytes (C order, little-endian).
struct NdArray {
    std::string dtype;
    std::vector<int64_t> shape;
    std::vector<uint8_t> bytes;

    int64_t nelem() const {
        int64_t n = 1;
        for (auto s : shape) n *= s;
        return n;
    }
    template <class T> const T* as() const { return reinterpret_cast<const T*>(bytes.data()); }
    template <class T> T* as() { return reinterpret_cast<T*>(bytes.data()); }
};

// Integer index arrays (sparse indptr/indices) come in either width -- scipy emits int32 for
// small matrices and int64 once a dimension or nnz exceeds 2^31. Normalize to int64 so callers
// don't have to branch on dtype.
inline std::vector<int64_t> as_i64(const NdArray& a) {
    const int64_t n = a.nelem();
    std::vector<int64_t> out(static_cast<size_t>(n));
    if (a.dtype == "<i8" || a.dtype == "<u8") {
        const int64_t* p = a.as<int64_t>();
        for (int64_t i = 0; i < n; ++i) out[i] = p[i];
    } else if (a.dtype == "<i4" || a.dtype == "<u4") {
        const int32_t* p = a.as<int32_t>();
        for (int64_t i = 0; i < n; ++i) out[i] = p[i];
    } else {
        throw std::runtime_error("as_i64: non-integer dtype " + a.dtype);
    }
    return out;
}

// Numeric values to double regardless of stored width -- measures are commonly float32 (<f4) in
// AnnData/Seurat as well as float64 (<f8). The reduction/transpose primitives take double, so
// normalize first.
inline std::vector<double> as_f64(const NdArray& a) {
    const int64_t n = a.nelem();
    std::vector<double> out(static_cast<size_t>(n));
    if (a.dtype == "<f8") {
        const double* p = a.as<double>();
        for (int64_t i = 0; i < n; ++i) out[i] = p[i];
    } else if (a.dtype == "<f4") {
        const float* p = a.as<float>();
        for (int64_t i = 0; i < n; ++i) out[i] = p[i];
    } else {
        throw std::runtime_error("as_f64: non-float dtype " + a.dtype);
    }
    return out;
}

// ---------------------------------------------------------------- model -----

struct Axis {
    std::string name;
    std::string origin = "observed";
    std::string role;                 // "" == none
    std::string induced_by;           // "" == none; the field this axis was induced from (factor/coordinate)
    std::vector<std::string> labels;
    json provenance = json::object();
};

struct Field {
    std::string name;
    std::string role;
    std::string encoding;             // dense | utf8 | csr | csc | coo
    std::string state;                // "" == none
    std::string subtype;              // "" == none
    std::vector<std::string> span;
    json provenance = json::object();
    json uncertainty = json(nullptr);   // optional per-value uncertainty metadata (round-trips verbatim)
    bool directed = false, has_directed = false;
    bool weighted = false, has_weighted = false;

    // value payload, by encoding:
    NdArray dense;                    // "dense"
    std::vector<std::string> strings; // "utf8"
    NdArray data, indices, indptr;    // "csr"/"csc"
    std::vector<int64_t> shape;       // sparse 2-D shape
    NdArray codes;                    // "categorical": int codes (-1 = missing)
    std::vector<std::string> categories;  // "categorical": inline category labels
    bool ordered = false, has_ordered = false;
    NdArray mask;                     // optional uint8 validity mask, 1 == missing (nullable Int/bool/string)
    bool has_mask = false;
    std::string coverage = "full";    // "full" or "partial" (-> index keys the covered subset of index_axis)
    NdArray index;                    // partial coverage: int64 positions into index_axis (one per value row)
    bool has_index = false;
    std::string index_axis;           // which span axis `index` keys into
};

// Lossless passthrough (AnnData `uns`, Seurat `@misc`). The core round-trips it *verbatim* and never
// interprets it: `attrs` is the group's whole `lstar` attribute dict (kind + `tree` JSON string + the
// `arrays` manifest), preserved as-is; `leaves` are the dense/utf8 arrays the manifest names. Only the
// originating profile (Python) walks the tree to rebuild the live object.
struct AuxLeaf {
    std::string id, kind;             // kind: "dense" | "utf8"
    NdArray dense;
    std::vector<std::string> strings;
};
struct Aux {
    std::string ns;                   // namespace, e.g. "anndata.uns"
    json attrs;                       // the group's lstar attrs (tree string + arrays manifest), verbatim
    std::vector<AuxLeaf> leaves;
};

struct Dataset {
    std::string kind = "sample";
    std::string spec_version = "0.1";
    std::vector<std::string> profiles, dropped;
    std::vector<Axis> axes;
    std::vector<Field> fields;
    std::vector<Aux> aux;

    Axis* axis(const std::string& n) {
        for (auto& a : axes) if (a.name == n) return &a;
        return nullptr;
    }
    Field* field(const std::string& n) {
        for (auto& f : fields) if (f.name == n) return &f;
        return nullptr;
    }
};

// ------------------------------------------------------------ low-level IO --

inline std::string read_text(const fs::path& p) {
    std::ifstream f(p, std::ios::binary);
    if (!f) throw std::runtime_error("cannot open " + p.string());
    std::ostringstream ss;
    ss << f.rdbuf();
    return ss.str();
}
inline std::vector<uint8_t> read_bytes(const fs::path& p) {
    // Bulk read: size via ate, then one read() call. (istreambuf_iterator reads byte-by-byte through
    // the streambuf -- ~10-50x slower; it dominated chunk-decode time on large stores.)
    std::ifstream f(p, std::ios::binary | std::ios::ate);
    if (!f) throw std::runtime_error("cannot open " + p.string());
    const std::streamsize n = f.tellg();
    std::vector<uint8_t> buf(n > 0 ? static_cast<size_t>(n) : 0);
    if (n > 0) {
        f.seekg(0);
        if (!f.read(reinterpret_cast<char*>(buf.data()), n))
            throw std::runtime_error("short read " + p.string());
    }
    return buf;
}
inline void write_text(const fs::path& p, const std::string& s) {
    std::ofstream f(p, std::ios::binary);
    f.write(s.data(), static_cast<std::streamsize>(s.size()));
}
inline void write_bytes(const fs::path& p, const uint8_t* b, size_t n) {
    std::ofstream f(p, std::ios::binary);
    f.write(reinterpret_cast<const char*>(b), static_cast<std::streamsize>(n));
}
inline json read_json(const fs::path& p) { return json::parse(read_text(p)); }

// A group's user attributes, format-agnostically: libzarr's Group::open probes v3 (zarr.json's
// "attributes") before v2 (.zattrs), so lstar reads either format's stores without branching. Opening a
// fresh FilesystemStore rooted at the group dir mirrors read_array's per-node pattern; the root group
// reads through .zmetadata/inline-consolidated when present, subgroups read their local metadata.
inline json read_group_attrs(const fs::path& dir) {
    auto store = std::make_shared<zarr::FilesystemStore>(dir, /*create=*/false);
    return zarr::Group::open(store, "").attributes();
}

// chunk key for an all-zero chunk index: "0" (1-D), "0.0" (2-D), ...
inline std::string zero_chunk_key(size_t ndim) {
    std::string k = "0";
    for (size_t i = 1; i < ndim; ++i) k += ".0";
    return k;
}

// ------------------------------------------------------------- codecs -------

// Inflate a zlib- or gzip-wrapped DEFLATE stream (the numcodecs "zlib"/"gzip" codecs). zlib's
// windowBits=47 (32 + 15) auto-detects either header.
inline std::vector<uint8_t> inflate_stream(const uint8_t* src, size_t n, size_t hint) {
#ifdef LSTAR_HAVE_ZLIB
    std::vector<uint8_t> out;
    out.reserve(hint ? hint : n * 3);
    z_stream zs;
    std::memset(&zs, 0, sizeof(zs));
    if (inflateInit2(&zs, 47) != Z_OK) throw std::runtime_error("inflateInit2 failed");
    zs.next_in = const_cast<Bytef*>(src);
    zs.avail_in = static_cast<uInt>(n);
    std::vector<uint8_t> buf(1 << 20);
    int ret;
    do {
        zs.next_out = buf.data();
        zs.avail_out = static_cast<uInt>(buf.size());
        ret = inflate(&zs, Z_NO_FLUSH);
        if (ret != Z_OK && ret != Z_STREAM_END && ret != Z_BUF_ERROR) {
            inflateEnd(&zs);
            throw std::runtime_error("inflate failed (" + std::to_string(ret) + ")");
        }
        out.insert(out.end(), buf.data(), buf.data() + (buf.size() - zs.avail_out));
        if (ret == Z_BUF_ERROR && zs.avail_in == 0) break;
    } while (ret != Z_STREAM_END);
    inflateEnd(&zs);
    return out;
#else
    (void)src; (void)n; (void)hint;
    throw std::runtime_error("compressed chunk needs zlib; build libstar with LSTAR_HAVE_ZLIB");
#endif
}

inline std::vector<uint8_t> decode_chunk(const json& compressor, std::vector<uint8_t> raw,
                                         size_t hint) {
    if (compressor.is_null()) return raw;
    std::string id = compressor.value("id", std::string());
    if (id == "gzip" || id == "zlib")
        return inflate_stream(raw.data(), raw.size(), hint);
    throw std::runtime_error("unsupported compressor id: " + id);
}

// Deflate to a gzip- (windowBits 15+16) or zlib- (15) wrapped stream -- the encode counterpart of
// inflate_stream, so the C++/R writer can emit the numcodecs "gzip"/"zlib" chunk codecs the reader
// (and zarr-python) already understand.
inline std::vector<uint8_t> deflate_stream(const uint8_t* src, size_t n, int level, bool gzip) {
#ifdef LSTAR_HAVE_ZLIB
    z_stream zs;
    std::memset(&zs, 0, sizeof(zs));
    if (deflateInit2(&zs, level, Z_DEFLATED, gzip ? (15 + 16) : 15, 8, Z_DEFAULT_STRATEGY) != Z_OK)
        throw std::runtime_error("deflateInit2 failed");
    zs.next_in = const_cast<Bytef*>(src);
    zs.avail_in = static_cast<uInt>(n);
    std::vector<uint8_t> out;
    out.reserve(n / 2 + 64);
    std::vector<uint8_t> buf(1 << 20);
    int ret;
    do {
        zs.next_out = buf.data();
        zs.avail_out = static_cast<uInt>(buf.size());
        ret = deflate(&zs, Z_FINISH);
        out.insert(out.end(), buf.data(), buf.data() + (buf.size() - zs.avail_out));
    } while (ret != Z_STREAM_END);
    deflateEnd(&zs);
    return out;
#else
    (void)src; (void)n; (void)level; (void)gzip;
    throw std::runtime_error("compression needs zlib; build libstar with LSTAR_HAVE_ZLIB");
#endif
}

inline std::vector<uint8_t> encode_chunk(const json& compressor, std::vector<uint8_t> raw) {
    if (compressor.is_null()) return raw;
    std::string id = compressor.value("id", std::string());
    int level = compressor.value("level", 5);
    if (id == "gzip") return deflate_stream(raw.data(), raw.size(), level, true);
    if (id == "zlib") return deflate_stream(raw.data(), raw.size(), level, false);
    throw std::runtime_error("unsupported compressor id for write: " + id);
}

// ------------------------------------------------------------ array IO ------

// Read a (possibly chunked, possibly compressed) zarr v2 array into a contiguous C-order buffer.
// Edge chunks are stored full-size and fill-padded; missing chunk files default to fill_value 0.
inline NdArray read_array(const fs::path& dir) {
    // A single zarr array lives at `dir` (its .zarray + chunk files); libzarr reads it as the root of a
    // filesystem store. Chunk decode, fill-padding, gzip, and multi-dim assembly are libzarr's job now.
    auto store = std::make_shared<zarr::FilesystemStore>(dir, false);
    zarr::Array za = zarr::Array::open(store, "");
    NdArray a;
    a.dtype = zarr::v2::emit_dtype(za.meta().dtype, /*big_endian=*/false);  // little-endian, C-order (lstar always)
    a.shape.assign(za.meta().shape.begin(), za.meta().shape.end());          // uint64 -> int64
    a.bytes.resize(za.nbytes());
    if (!a.bytes.empty()) za.read(a.bytes.data(), a.bytes.size());
    return a;
}

// Read elements [lo, hi) of a 1-D zarr array, decoding only the chunks that overlap the range.
// This is the bounded-memory primitive behind the blocked reader: a column block of a CSC measure
// touches only its slice of `data`/`indices`, so a per-gene reduction never reads the whole matrix.
// (When the array was written as a single chunk -- the unchunked default -- this still decodes that
// one chunk, so bounded streaming needs a chunked store, e.g. one written with chunk_elems set.)
inline NdArray read_array_range(const fs::path& dir, int64_t lo, int64_t hi) {
    // Elements [lo, hi) of a 1-D array. libzarr's read_region reads only the overlapping chunks (decoding
    // each), padding missing chunks with fill 0 -- the same bounded read the hand-rolled loop did.
    auto store = std::make_shared<zarr::FilesystemStore>(dir, false);
    zarr::Array za = zarr::Array::open(store, "");
    if (za.meta().shape.size() != 1) throw std::runtime_error("read_array_range: 1-D only: " + dir.string());
    const int64_t n = static_cast<int64_t>(za.meta().shape[0]);
    lo = std::max<int64_t>(0, lo);
    hi = std::min<int64_t>(n, hi);
    const int64_t len = std::max<int64_t>(0, hi - lo);
    NdArray a;
    a.dtype = zarr::v2::emit_dtype(za.meta().dtype, /*big_endian=*/false);
    a.shape = {len};
    const size_t dsz = za.meta().dtype.itemsize;
    a.bytes.assign(static_cast<size_t>(len) * dsz, 0);
    if (len == 0) return a;
    za.read_region({static_cast<std::uint64_t>(lo)}, {static_cast<std::uint64_t>(len)}, a.bytes.data(), a.bytes.size());
    return a;
}

// Chunk shape splitting the first axis so each chunk holds ~chunk_elems elements (0 -> single chunk,
// the portable default). Mirrors the Python writer's `_chunks_for`.
inline std::vector<int64_t> chunk_shape_for(const std::vector<int64_t>& shape, int64_t chunk_elems) {
    if (chunk_elems <= 0 || shape.empty() || shape[0] == 0) return shape;
    int64_t inner = 1;
    for (size_t i = 1; i < shape.size(); ++i) inner *= shape[i];
    int64_t rows = std::max<int64_t>(1, chunk_elems / std::max<int64_t>(1, inner));
    if (rows >= shape[0]) return shape;
    std::vector<int64_t> cs = shape;
    cs[0] = rows;
    return cs;
}

// Shard (outer-chunk) shape packing ~shard_elems into one shard object -- a v3 hosting optimization:
// many inner chunks in one object (fewer HTTP requests), byte-range-readable via the shard index. The
// shard first-axis extent is a POSITIVE MULTIPLE of the chunk extent (v3 spec), capped so one shard
// covers at most the whole array. Returns {} (unsharded) when there's nothing to gain (single chunk).
inline std::vector<int64_t> shard_shape_for(const std::vector<int64_t>& shape,
                                            const std::vector<int64_t>& chunks, int64_t shard_elems) {
    if (shard_elems <= 0 || shape.empty() || chunks.empty() || chunks[0] <= 0) return {};
    const int64_t chunk_rows = chunks[0];
    if (chunk_rows >= shape[0]) return {};                     // single chunk -> nothing to shard
    int64_t inner = 1;
    for (size_t i = 1; i < shape.size(); ++i) inner *= shape[i];
    const int64_t num_chunks = (shape[0] + chunk_rows - 1) / chunk_rows;
    const int64_t target_rows = shard_elems / std::max<int64_t>(1, inner);
    const int64_t k = std::max<int64_t>(1, std::min<int64_t>(num_chunks, target_rows / chunk_rows));
    std::vector<int64_t> ss = chunks;                          // inner dims: shard extent == chunk extent (a multiple)
    ss[0] = k * chunk_rows;
    return ss;
}

// Write a zarr v2 array, optionally chunked (along the first axis) and/or compressed. chunk_elems<=0
// and compressor=null reproduce the original single-chunk uncompressed output exactly (so existing
// stores are byte-identical). Because only the first (outermost, slowest) axis is chunked, each chunk
// is a contiguous byte range of the C-order buffer -- no scatter; edge chunks are full-size,
// fill-padded per the v2 spec.
inline void write_array(const fs::path& dir, const NdArray& a,
                        int64_t chunk_elems = 0, const json& compressor = json(nullptr),
                        zarr::ZarrFormat fmt = zarr::ZarrFormat::v2, int64_t shard_elems = 0) {
    // libzarr writes the array metadata + chunks in `fmt` (v2 .zarray, or v3 zarr.json).
    // chunk_shape_for keeps lstar's first-axis-only chunking and single-chunk default; fill 0 matches
    // lstar's layout. For v3, libzarr prepends the required `bytes` codec and uses the 'default' chunk-key
    // encoding ('/'), so the same compressor spec + dtype carry over. shard_elems>0 (v3 only) packs inner
    // chunks into shard objects. Canonical metadata formatting differs from the old hand-rolled dump, but
    // readers are value-based (conformance checks values, not store bytes).
    auto store = std::make_shared<zarr::FilesystemStore>(dir, /*create=*/true);
    zarr::ArraySpec spec;
    spec.format = fmt;
    spec.shape.assign(a.shape.begin(), a.shape.end());
    const std::vector<int64_t> chunks = chunk_shape_for(a.shape, chunk_elems);
    spec.chunks.assign(chunks.begin(), chunks.end());
    if (fmt == zarr::ZarrFormat::v3 && shard_elems > 0) {      // v3 sharding: shard extent = multiple of chunk
        const std::vector<int64_t> shards = shard_shape_for(a.shape, chunks, shard_elems);
        if (!shards.empty()) spec.shards.assign(shards.begin(), shards.end());
    }
    spec.dtype = zarr::v2::parse_dtype(a.dtype, dir.string()).dtype;
    spec.dimension_separator = '.';
    if (!compressor.is_null()) {
        const std::string id = compressor.value("id", std::string());
        if (id != "gzip" && id != "zlib") throw std::runtime_error("unsupported compressor id for write: " + id);
        spec.codecs.push_back(zarr::CodecSpec{id, {{"level", compressor.value("level", 5)}}});
    }
    zarr::Array za = zarr::Array::create(store, "", spec);
    if (!a.bytes.empty()) za.write(a.bytes.data(), a.bytes.size());
}

// ---------------------------------------------------------- string codec ----

inline std::vector<std::string> read_strings(const fs::path& gdir, const std::string& name) {
    NdArray data = read_array(gdir / name);
    NdArray offs = read_array(gdir / (name + "_offsets"));
    const int64_t* off = offs.as<int64_t>();
    const char* buf = reinterpret_cast<const char*>(data.bytes.data());
    int64_t n = offs.nelem() - 1;
    std::vector<std::string> out;
    out.reserve(static_cast<size_t>(n));
    for (int64_t i = 0; i < n; ++i) out.emplace_back(buf + off[i], buf + off[i + 1]);
    return out;
}

inline void write_strings(const fs::path& gdir, const std::string& name,
                          const std::vector<std::string>& strs,
                          int64_t chunk_elems = 0, const json& compressor = json(nullptr),
                          zarr::ZarrFormat fmt = zarr::ZarrFormat::v2, int64_t shard_elems = 0) {
    std::vector<int64_t> off(strs.size() + 1, 0);
    std::string buf;
    for (size_t i = 0; i < strs.size(); ++i) {
        buf += strs[i];
        off[i + 1] = off[i] + static_cast<int64_t>(strs[i].size());
    }
    NdArray data;
    data.dtype = "|u1";
    data.shape = {static_cast<int64_t>(buf.size())};
    data.bytes.assign(buf.begin(), buf.end());
    NdArray offsets;
    offsets.dtype = "<i8";
    offsets.shape = {static_cast<int64_t>(off.size())};
    offsets.bytes.resize(off.size() * 8);
    std::memcpy(offsets.bytes.data(), off.data(), off.size() * 8);
    write_array(gdir / name, data, chunk_elems, compressor, fmt, shard_elems);
    write_array(gdir / (name + "_offsets"), offsets, chunk_elems, compressor, fmt, shard_elems);
}

// ------------------------------------------------------------ group IO ------

inline void write_group(const fs::path& dir, const json& attrs,
                        zarr::ZarrFormat fmt = zarr::ZarrFormat::v2) {
    fs::create_directories(dir);
    if (fmt == zarr::ZarrFormat::v2) {                    // v2: hand-write .zgroup + .zattrs (proven default path)
        write_text(dir / ".zgroup", json{{"zarr_format", 2}}.dump());
        write_text(dir / ".zattrs", attrs.dump());
    } else {                                              // v3: libzarr emits zarr.json (node_type + attributes)
        auto store = std::make_shared<zarr::FilesystemStore>(dir, /*create=*/true);
        zarr::Group g = zarr::Group::create(store, "", zarr::ZarrFormat::v3);
        g.set_attributes(attrs);
    }
}

// Walk a finished store and emit consolidated metadata (so a reader gets the tree in one read instead
// of many small stats): v2 -> a `.zmetadata` document; v3 -> the inline convention in root zarr.json.
inline void consolidate_metadata(const fs::path& root, zarr::ZarrFormat fmt = zarr::ZarrFormat::v2) {
    // libzarr emits deterministically, parent-before-child. (The old fs-walk order was a latent hazard
    // for strict consolidated readers -- the class of bug we hit on the JS writer side.)
    zarr::FilesystemStore store(root, /*create=*/false);
    if (fmt == zarr::ZarrFormat::v2) zarr::v2::consolidate(store);
    else zarr::v3::consolidate(store);
}

inline std::string opt_str(const json& m, const char* k) {
    return (m.contains(k) && !m[k].is_null()) ? m[k].get<std::string>() : std::string();
}

// ------------------------------------------------------------ STORED zip ----
//
// A `.lstar.zarr.zip` is a normal store packed into ONE file with every entry STORED (no deflate):
// zarr chunks are already codec-compressed, so re-deflating wastes CPU for ~no gain, and -- the
// load-bearing reason -- only a STORED entry stays byte-range-readable inside the archive, which is
// the point of a hosted single file. C++/R read a `.zip` by extracting its STORED entries to a temp
// dir (a seek + copy, NO decompression -- the STORED win) and running the normal reader; they write
// one by writing a normal store and packing it. (Reading a chunk's byte range *without* extracting --
// the remote-hosted case -- is handled on the JS surface, where HTTP Range into the zip is the point.)

namespace zipfmt {
inline uint16_t rd16(const uint8_t* p) { return (uint16_t)(p[0] | (p[1] << 8)); }
inline uint32_t rd32(const uint8_t* p) {
    return (uint32_t)(p[0] | (p[1] << 8) | (p[2] << 16) | ((uint32_t)p[3] << 24));
}
inline uint64_t rd64(const uint8_t* p) { uint64_t v = 0; for (int i = 0; i < 8; i++) v |= (uint64_t)p[i] << (8 * i); return v; }
inline void wr16(std::vector<uint8_t>& o, uint16_t v) { o.push_back(v & 0xff); o.push_back((v >> 8) & 0xff); }
inline void wr32(std::vector<uint8_t>& o, uint32_t v) { for (int i = 0; i < 4; i++) o.push_back((v >> (8 * i)) & 0xff); }
inline void wr64(std::vector<uint8_t>& o, uint64_t v) { for (int i = 0; i < 8; i++) o.push_back((v >> (8 * i)) & 0xff); }

inline uint32_t crc32_bytes(const uint8_t* data, size_t n) {  // standard CRC-32 (poly 0xEDB88320), no zlib dep
    static const std::array<uint32_t, 256> table = [] {
        std::array<uint32_t, 256> a{};
        for (uint32_t i = 0; i < 256; i++) {
            uint32_t c = i;
            for (int k = 0; k < 8; k++) c = (c & 1) ? (0xEDB88320u ^ (c >> 1)) : (c >> 1);
            a[i] = c;
        }
        return a;
    }();
    uint32_t c = 0xFFFFFFFFu;
    for (size_t i = 0; i < n; i++) c = table[(c ^ data[i]) & 0xff] ^ (c >> 8);
    return c ^ 0xFFFFFFFFu;
}

inline std::vector<uint8_t> read_range(std::istream& f, uint64_t off, uint64_t n) {
    f.seekg((std::streamoff)off);
    std::vector<uint8_t> buf(n);
    if (n && !f.read((char*)buf.data(), (std::streamsize)n)) throw std::runtime_error("zip: short read");
    return buf;
}
}  // namespace zipfmt

struct ZipIndexEntry { uint64_t data_off; uint64_t size; };

// Parse a zip's central directory (ZIP64-aware) into name -> (data offset, size). STORED-only: a
// DEFLATE-compressed entry throws a clear, actionable error (a hosted single-file store must be STORED
// so its chunks stay byte-range-readable).
inline std::map<std::string, ZipIndexEntry> zip_read_index(std::istream& f, uint64_t fsize,
                                                           const std::string& forwhat) {
    using namespace zipfmt;
    const uint64_t taillen = std::min<uint64_t>(fsize, 65557);
    std::vector<uint8_t> tail = read_range(f, fsize - taillen, taillen);
    // locate the End Of Central Directory record (sig 0x06054b50), scanning backward
    int64_t e = -1;
    for (int64_t i = (int64_t)tail.size() - 22; i >= 0; --i)
        if (rd32(&tail[i]) == 0x06054b50u) { e = i; break; }
    if (e < 0) throw std::runtime_error(forwhat + ": not a zip (no end-of-central-directory record)");
    uint64_t n_entries = rd16(&tail[e + 10]);
    uint64_t cd_size = rd32(&tail[e + 12]);
    uint64_t cd_off = rd32(&tail[e + 16]);
    if (n_entries == 0xFFFFu || cd_size == 0xFFFFFFFFu || cd_off == 0xFFFFFFFFu) {  // ZIP64
        int64_t loc = e - 20;
        if (loc < 0 || rd32(&tail[loc]) != 0x07064b50u)
            throw std::runtime_error(forwhat + ": ZIP64 end-of-central-directory locator missing");
        uint64_t z64off = rd64(&tail[loc + 8]);
        std::vector<uint8_t> z = read_range(f, z64off, 56);
        if (rd32(&z[0]) != 0x06064b50u) throw std::runtime_error(forwhat + ": bad ZIP64 EOCD record");
        n_entries = rd64(&z[32]);
        cd_size = rd64(&z[40]);
        cd_off = rd64(&z[48]);
    }
    std::vector<uint8_t> cd = read_range(f, cd_off, cd_size);
    std::map<std::string, ZipIndexEntry> idx;
    std::vector<std::string> deflated;
    uint64_t p = 0;
    for (uint64_t i = 0; i < n_entries; i++) {
        if (p + 46 > cd.size() || rd32(&cd[p]) != 0x02014b50u)
            throw std::runtime_error(forwhat + ": corrupt central directory");
        uint16_t method = rd16(&cd[p + 10]);
        uint64_t usize = rd32(&cd[p + 24]);
        uint16_t nlen = rd16(&cd[p + 28]);
        uint16_t elen = rd16(&cd[p + 30]);
        uint16_t clen = rd16(&cd[p + 32]);
        uint64_t lho = rd32(&cd[p + 42]);
        std::string name((const char*)&cd[p + 46], nlen);
        // ZIP64 extra: fill in whichever of usize/lho were 0xFFFFFFFF, in spec order
        uint64_t ep = p + 46 + nlen, eend = ep + elen;
        while (ep + 4 <= eend) {
            uint16_t hid = rd16(&cd[ep]), hsz = rd16(&cd[ep + 2]);
            uint64_t q = ep + 4;
            if (hid == 0x0001) {
                if (usize == 0xFFFFFFFFu && q + 8 <= eend) { usize = rd64(&cd[q]); q += 8; }
                if (rd32(&cd[p + 20]) == 0xFFFFFFFFu && q + 8 <= eend) { q += 8; }  // compressed (==usize for STORED)
                if (lho == 0xFFFFFFFFu && q + 8 <= eend) { lho = rd64(&cd[q]); q += 8; }
            }
            ep += 4 + hsz;
        }
        if (method != 0) deflated.push_back(name);
        // data starts after the LOCAL header, whose name/extra lengths can differ from the central's
        std::vector<uint8_t> lh = zipfmt::read_range(f, lho, 30);
        uint64_t data_off = lho + 30 + rd16(&lh[26]) + rd16(&lh[28]);
        idx[name] = ZipIndexEntry{data_off, usize};
        p += 46 + nlen + elen + clen;
    }
    if (!deflated.empty())
        throw std::runtime_error(
            forwhat + ": this .lstar.zarr.zip is DEFLATE-compressed (" + std::to_string(deflated.size()) +
            " entries, e.g. '" + deflated[0] + "') -- a hosted single-file store must be written STORED so "
            "its chunks stay byte-range-readable. Repack it STORED (lstar convert, or `zip -0 -r`).");
    return idx;
}

inline fs::path unique_temp_dir(const std::string& tag) {
    fs::path base = fs::temp_directory_path();
    static std::atomic<uint64_t> ctr{0};
    for (int i = 0; i < 100000; i++) {
        uint64_t salt = ((uintptr_t)&ctr) ^ (ctr.fetch_add(1) * 0x9E3779B97F4A7C15ull);
        fs::path p = base / (tag + std::to_string(salt));
        std::error_code ec;
        if (fs::create_directory(p, ec)) return p;
    }
    throw std::runtime_error("cannot create a unique temp directory under " + base.string());
}

// Extract every STORED entry of `zippath` into `outdir` (a seek + copy per entry; no decompression).
inline void zip_extract_to_dir(const fs::path& zippath, const fs::path& outdir) {
    std::ifstream f(zippath, std::ios::binary | std::ios::ate);
    if (!f) throw std::runtime_error("cannot open " + zippath.string());
    uint64_t fsize = (uint64_t)f.tellg();
    auto idx = zip_read_index(f, fsize, zippath.string());
    for (auto& kv : idx) {
        if (kv.first.empty() || kv.first.back() == '/') continue;   // skip directory entries
        fs::path out = outdir / kv.first;
        fs::create_directories(out.parent_path());
        auto bytes = zipfmt::read_range(f, kv.second.data_off, kv.second.size);
        write_bytes(out, bytes.data(), bytes.size());
    }
}

// Pack a directory store into ONE file with every entry STORED (ZIP64-aware for >4GB / >65535 entries).
inline void pack_stored_zip(const fs::path& srcdir, const fs::path& zippath) {
    using namespace zipfmt;
    std::vector<std::pair<std::string, fs::path>> files;
    for (auto& en : fs::recursive_directory_iterator(srcdir)) {
        if (!en.is_regular_file()) continue;
        files.push_back({fs::relative(en.path(), srcdir).generic_string(), en.path()});
    }
    auto is_meta = [](const std::string& a) {                       // basename starts ".z" (.zarray/.zattrs/.z*)
        size_t s = a.rfind('/'); return a.compare(s == std::string::npos ? 0 : s + 1, 2, ".z") == 0;
    };
    std::sort(files.begin(), files.end(), [&](const auto& a, const auto& b) {  // metadata first, then by name
        bool am = is_meta(a.first), bm = is_meta(b.first);
        return am != bm ? am : a.first < b.first;
    });
    std::ofstream out(zippath, std::ios::binary);
    if (!out) throw std::runtime_error("cannot write " + zippath.string());
    struct CD { std::string name; uint32_t crc; uint64_t size, off; };
    std::vector<CD> cds;
    uint64_t offset = 0;
    for (auto& fp : files) {
        std::vector<uint8_t> data = read_bytes(fp.second);
        uint32_t crc = crc32_bytes(data.data(), data.size());
        uint64_t sz = data.size();
        bool z64 = sz >= 0xFFFFFFFFull;
        std::vector<uint8_t> extra;
        if (z64) { wr16(extra, 0x0001); wr16(extra, 16); wr64(extra, sz); wr64(extra, sz); }
        std::vector<uint8_t> lh;
        wr32(lh, 0x04034b50); wr16(lh, z64 ? 45 : 20); wr16(lh, 0); wr16(lh, 0);
        wr16(lh, 0); wr16(lh, 0x21); wr32(lh, crc);
        wr32(lh, z64 ? 0xFFFFFFFFu : (uint32_t)sz); wr32(lh, z64 ? 0xFFFFFFFFu : (uint32_t)sz);
        wr16(lh, (uint16_t)fp.first.size()); wr16(lh, (uint16_t)extra.size());
        lh.insert(lh.end(), fp.first.begin(), fp.first.end());
        lh.insert(lh.end(), extra.begin(), extra.end());
        out.write((char*)lh.data(), (std::streamsize)lh.size());
        out.write((char*)data.data(), (std::streamsize)data.size());
        cds.push_back({fp.first, crc, sz, offset});
        offset += lh.size() + data.size();
    }
    uint64_t cd_start = offset;
    std::vector<uint8_t> cd;
    for (auto& c : cds) {
        bool zsz = c.size >= 0xFFFFFFFFull, zoff = c.off >= 0xFFFFFFFFull, z64 = zsz || zoff;
        std::vector<uint8_t> extra;
        if (z64) {
            std::vector<uint8_t> body;
            if (zsz) { wr64(body, c.size); wr64(body, c.size); }
            if (zoff) { wr64(body, c.off); }
            wr16(extra, 0x0001); wr16(extra, (uint16_t)body.size());
            extra.insert(extra.end(), body.begin(), body.end());
        }
        wr32(cd, 0x02014b50); wr16(cd, z64 ? 45 : 20); wr16(cd, z64 ? 45 : 20); wr16(cd, 0);
        wr16(cd, 0); wr16(cd, 0); wr16(cd, 0x21); wr32(cd, c.crc);
        wr32(cd, zsz ? 0xFFFFFFFFu : (uint32_t)c.size); wr32(cd, zsz ? 0xFFFFFFFFu : (uint32_t)c.size);
        wr16(cd, (uint16_t)c.name.size()); wr16(cd, (uint16_t)extra.size());
        wr16(cd, 0); wr16(cd, 0); wr16(cd, 0); wr32(cd, 0);
        wr32(cd, zoff ? 0xFFFFFFFFu : (uint32_t)c.off);
        cd.insert(cd.end(), c.name.begin(), c.name.end());
        cd.insert(cd.end(), extra.begin(), extra.end());
    }
    out.write((char*)cd.data(), (std::streamsize)cd.size());
    uint64_t cd_size = cd.size(), nrec = cds.size();
    bool need_z64 = nrec >= 0xFFFFu || cd_start >= 0xFFFFFFFFull || cd_size >= 0xFFFFFFFFull;
    if (need_z64) {
        uint64_t z64eocd = cd_start + cd_size;
        std::vector<uint8_t> z;
        wr32(z, 0x06064b50); wr64(z, 44); wr16(z, 45); wr16(z, 45); wr32(z, 0); wr32(z, 0);
        wr64(z, nrec); wr64(z, nrec); wr64(z, cd_size); wr64(z, cd_start);
        wr32(z, 0x07064b50); wr32(z, 0); wr64(z, z64eocd); wr32(z, 1);
        out.write((char*)z.data(), (std::streamsize)z.size());
    }
    std::vector<uint8_t> eo;
    wr32(eo, 0x06054b50); wr16(eo, 0); wr16(eo, 0);
    wr16(eo, nrec >= 0xFFFFu ? 0xFFFFu : (uint16_t)nrec);
    wr16(eo, nrec >= 0xFFFFu ? 0xFFFFu : (uint16_t)nrec);
    wr32(eo, cd_size >= 0xFFFFFFFFull ? 0xFFFFFFFFu : (uint32_t)cd_size);
    wr32(eo, cd_start >= 0xFFFFFFFFull ? 0xFFFFFFFFu : (uint32_t)cd_start);
    wr16(eo, 0);
    out.write((char*)eo.data(), (std::streamsize)eo.size());
}

// ------------------------------------------------------------ dataset IO ----

inline Dataset read(const fs::path& root) {
    if (root.extension() == ".zip") {                 // single-file .lstar.zarr.zip: extract, read, clean up
        fs::path tmp = unique_temp_dir("lstar_unzip_");
        try {
            zip_extract_to_dir(root, tmp);
            Dataset ds = read(tmp);
            fs::remove_all(tmp);
            return ds;
        } catch (...) { std::error_code ec; fs::remove_all(tmp, ec); throw; }
    }
    json rmeta = read_group_attrs(root)["lstar"];
    Dataset ds;
    ds.kind = rmeta.value("kind", std::string("sample"));
    ds.spec_version = rmeta.value("spec_version", std::string("0.1"));
    ds.profiles = rmeta.value("profiles", std::vector<std::string>{});
    ds.dropped = rmeta.value("dropped", std::vector<std::string>{});

    for (auto& an : rmeta["axes"]) {
        std::string name = an.get<std::string>();
        fs::path g = root / "axes" / name;
        json m = read_group_attrs(g)["lstar"];
        Axis a;
        a.name = name;
        a.origin = m.value("origin", std::string("observed"));
        a.role = opt_str(m, "role");
        a.induced_by = opt_str(m, "induced_by");
        if (m.contains("provenance") && !m["provenance"].is_null()) a.provenance = m["provenance"];
        a.labels = read_strings(g, "labels");
        ds.axes.push_back(std::move(a));
    }

    for (auto& fn : rmeta["fields"]) {
        std::string name = fn.get<std::string>();
        fs::path g = root / "fields" / name;
        json m = read_group_attrs(g)["lstar"];
        Field f;
        f.name = name;
        f.role = opt_str(m, "role");
        f.encoding = opt_str(m, "encoding");
        f.state = opt_str(m, "state");
        f.subtype = opt_str(m, "subtype");
        f.coverage = m.contains("coverage") && !m["coverage"].is_null() ? m["coverage"].get<std::string>() : "full";
        f.index_axis = opt_str(m, "index_axis");
        if (m.contains("span") && !m["span"].is_null())
            f.span = m["span"].get<std::vector<std::string>>();
        if (m.contains("provenance") && !m["provenance"].is_null()) f.provenance = m["provenance"];
        if (m.contains("uncertainty") && !m["uncertainty"].is_null()) f.uncertainty = m["uncertainty"];
        if (m.contains("directed") && !m["directed"].is_null()) {
            f.directed = m["directed"].get<bool>();
            f.has_directed = true;
        }
        if (m.contains("weighted") && !m["weighted"].is_null()) {
            f.weighted = m["weighted"].get<bool>();
            f.has_weighted = true;
        }
        if (f.encoding == "csr" || f.encoding == "csc") {
            f.data = read_array(g / "data");
            f.indices = read_array(g / "indices");
            f.indptr = read_array(g / "indptr");
            f.shape = m["shape"].get<std::vector<int64_t>>();
        } else if (f.encoding == "utf8") {
            f.strings = read_strings(g, "values");
        } else if (f.encoding == "categorical") {
            f.codes = read_array(g / "codes");
            if (fs::exists(g / "categories")) f.categories = read_strings(g, "categories");  // inline
            f.ordered = m.value("ordered", false);
            f.has_ordered = true;
        } else {  // dense
            f.dense = read_array(g / "values");
        }
        if (m.contains("nullable") && !m["nullable"].is_null() && m["nullable"].get<bool>()
            && fs::exists(g / "mask")) {                               // optional validity mask
            f.mask = read_array(g / "mask");
            f.has_mask = true;
        }
        if (fs::exists(g / "index")) {                                 // partial coverage index
            f.index = read_array(g / "index");
            f.has_index = true;
        }
        ds.fields.push_back(std::move(f));
    }

    if (rmeta.contains("passthrough") && !rmeta["passthrough"].is_null()) {  // verbatim passthrough subtree
        for (auto& an : rmeta["passthrough"]) {
            std::string ns = an.get<std::string>();
            fs::path g = root / "passthrough" / ns;
            Aux ax;
            ax.ns = ns;
            ax.attrs = read_group_attrs(g)["lstar"];
            if (ax.attrs.contains("arrays") && !ax.attrs["arrays"].is_null()) {
                for (auto& a : ax.attrs["arrays"]) {
                    AuxLeaf leaf;
                    leaf.id = a["id"].get<std::string>();
                    leaf.kind = a["kind"].get<std::string>();
                    if (leaf.kind == "utf8") leaf.strings = read_strings(g, leaf.id);
                    else leaf.dense = read_array(g / leaf.id);
                    ax.leaves.push_back(std::move(leaf));
                }
            }
            ds.aux.push_back(std::move(ax));
        }
    }
    return ds;
}

inline void write(const Dataset& ds, const fs::path& root,
                  int64_t chunk_elems = 0, const json& compressor = json(nullptr),
                  zarr::ZarrFormat fmt = zarr::ZarrFormat::v2, int64_t shard_elems = 0) {
    if (root.extension() == ".zip") {                 // single-file .lstar.zarr.zip: write a dir, pack STORED
        fs::path tmp = unique_temp_dir("lstar_zip_");
        try {
            write(ds, tmp, chunk_elems, compressor, fmt, shard_elems);  // recurse into the directory path (writes + consolidates)
            pack_stored_zip(tmp, root);
            fs::remove_all(tmp);
            return;
        } catch (...) { std::error_code ec; fs::remove_all(tmp, ec); throw; }
    }
    if (fs::exists(root)) fs::remove_all(root);

    std::vector<std::string> axnames, fnames;
    for (auto& a : ds.axes) axnames.push_back(a.name);
    for (auto& f : ds.fields) fnames.push_back(f.name);

    json rl;
    rl["spec_version"] = ds.spec_version;
    rl["kind"] = ds.kind;
    rl["profiles"] = ds.profiles;
    rl["dropped"] = ds.dropped;
    rl["axes"] = axnames;
    rl["fields"] = fnames;
    std::vector<std::string> auxnames;
    for (auto& a : ds.aux) auxnames.push_back(a.ns);
    rl["passthrough"] = auxnames;
    write_group(root, json{{"lstar", rl}}, fmt);
    write_group(root / "axes", json::object(), fmt);
    write_group(root / "fields", json::object(), fmt);
    write_group(root / "models", json::object(), fmt);

    for (auto& a : ds.axes) {
        fs::path g = root / "axes" / a.name;
        json al;
        al["kind"] = "axis";
        al["origin"] = a.origin;
        al["role"] = a.role.empty() ? json(nullptr) : json(a.role);
        al["induced_by"] = a.induced_by.empty() ? json(nullptr) : json(a.induced_by);
        al["provenance"] = a.provenance;
        write_group(g, json{{"lstar", al}}, fmt);
        write_strings(g, "labels", a.labels, chunk_elems, compressor, fmt, shard_elems);
    }

    for (auto& f : ds.fields) {
        fs::path g = root / "fields" / f.name;
        json fl;
        fl["kind"] = "field";
        fl["role"] = f.role;
        fl["encoding"] = f.encoding;
        fl["span"] = f.span;
        fl["state"] = f.state.empty() ? json(nullptr) : json(f.state);
        fl["subtype"] = f.subtype.empty() ? json(nullptr) : json(f.subtype);
        fl["coverage"] = f.has_index ? json("partial") : json(f.coverage.empty() ? "full" : f.coverage);
        fl["index_axis"] = f.has_index && !f.index_axis.empty() ? json(f.index_axis) : json(nullptr);
        fl["uncertainty"] = f.uncertainty;   // preserve verbatim (was hardcoded null -> silent data loss)
        fl["provenance"] = f.provenance;
        fl["directed"] = f.has_directed ? json(f.directed) : json(nullptr);
        fl["weighted"] = f.has_weighted ? json(f.weighted) : json(nullptr);
        fl["nullable"] = f.has_mask ? json(true) : json(nullptr);
        if (f.encoding == "csr" || f.encoding == "csc") fl["shape"] = f.shape;
        if (f.encoding == "utf8") fl["shape"] = std::vector<int64_t>{(int64_t)f.strings.size()};
        if (f.encoding == "categorical") {
            fl["shape"] = std::vector<int64_t>{f.codes.nelem()};
            fl["ordered"] = f.has_ordered ? json(f.ordered) : json(false);
        }
        if (!(f.encoding == "csr" || f.encoding == "csc" || f.encoding == "utf8" || f.encoding == "categorical"))
            fl["shape"] = f.dense.shape;   // dense: shape in the manifest too (parity; a reader shouldn't need values/.zarray)
        write_group(g, json{{"lstar", fl}}, fmt);
        if (f.encoding == "csr" || f.encoding == "csc") {
            write_array(g / "data", f.data, chunk_elems, compressor, fmt, shard_elems);
            write_array(g / "indices", f.indices, chunk_elems, compressor, fmt, shard_elems);
            write_array(g / "indptr", f.indptr, chunk_elems, compressor, fmt, shard_elems);
        } else if (f.encoding == "utf8") {
            write_strings(g, "values", f.strings, chunk_elems, compressor, fmt, shard_elems);
        } else if (f.encoding == "categorical") {
            write_array(g / "codes", f.codes, chunk_elems, compressor, fmt, shard_elems);
            write_strings(g, "categories", f.categories, chunk_elems, compressor, fmt, shard_elems);  // inline (P1)
        } else {
            write_array(g / "values", f.dense, chunk_elems, compressor, fmt, shard_elems);
        }
        if (f.has_mask) write_array(g / "mask", f.mask, chunk_elems, compressor, fmt, shard_elems);
        if (f.has_index) write_array(g / "index", f.index, chunk_elems, compressor, fmt, shard_elems);
    }

    if (!ds.aux.empty()) {                                            // verbatim passthrough subtree
        write_group(root / "passthrough", json::object(), fmt);
        for (auto& ax : ds.aux) {
            fs::path g = root / "passthrough" / ax.ns;
            write_group(g, json{{"lstar", ax.attrs}}, fmt);
            for (auto& leaf : ax.leaves) {
                if (leaf.kind == "utf8") write_strings(g, leaf.id, leaf.strings, chunk_elems, compressor, fmt, shard_elems);
                else write_array(g / leaf.id, leaf.dense, chunk_elems, compressor, fmt, shard_elems);
            }
        }
    }
    consolidate_metadata(root, fmt);
}

// ---------------------------------------------------- translation primitives --
//
// Performance-critical kernels (the reason for a C++ core). These parallelize over the
// outer (column) axis with OpenMP -- the embarrassingly-parallel pattern that conversions,
// variance modeling, and per-group summaries need. Thread count is an explicit argument
// (caller's policy), matching the convention proven in pagoda2 (n_threads<=0 -> default).

struct ColStats {
    std::vector<double> mean, var;
    std::vector<int64_t> nnz;
};

// CSC <-> CSR conversion (an exact, O(nnz) storage transpose). A cells x genes measure stored CSC
// (compressed by gene) becomes the same matrix stored CSR (compressed by cell) -- the orientation
// flip the Seurat/SCE profiles need (genes x cells assays) without densifying. Run it the other
// direction to flip back; it is its own inverse on the (rows<->cols) labels.
template <class T, class Idx = int64_t>
struct CsxArrays {
    std::vector<T> data;                  // preserves the input value dtype (no widening)
    std::vector<Idx> indices;             // index width follows the caller: int64 by default, int32 for
    std::vector<Idx> indptr;              // the memory-lean WASM/browser prep (halves the per-nonzero term)
    int64_t nrows = 0, ncols = 0;
};

// `data` is templated on the stored value dtype so a float32 measure transposes as float32 (no widening).
// `Idx` is the index width: L* stores hold int32 indices, so the default int64 means every existing caller
// (Python/R/core, which normalize to int64 via as_i64()) is unchanged by template *deduction* from its
// int64 pointers, while the WASM prep passes int32 pointers to keep indices at native width -- the indices
// are the DOMINANT per-nonzero term at scale (two int64 nnz arrays, in + out), so halving them is the bulk
// of the browser memory win. Precondition for the int32 instantiation: nnz (and ncols/nrows) all < 2^31 --
// always true when the caller's int32 indptr can even represent nnz.
template <class T, class Idx = int64_t>
inline CsxArrays<T, Idx> csc_to_csr(const T* data, const Idx* indices, const Idx* indptr,
                                    int64_t nrows, int64_t ncols) {
    const int64_t nnz = static_cast<int64_t>(indptr[ncols]);
    CsxArrays<T, Idx> out;
    out.nrows = nrows;
    out.ncols = ncols;
    out.data.resize(static_cast<size_t>(nnz));
    out.indices.resize(static_cast<size_t>(nnz));
    out.indptr.assign(static_cast<size_t>(nrows) + 1, 0);
    for (int64_t k = 0; k < nnz; ++k) out.indptr[indices[k] + 1]++;     // count per row
    for (int64_t i = 0; i < nrows; ++i) out.indptr[i + 1] += out.indptr[i];
    std::vector<Idx> next(out.indptr.begin(), out.indptr.end() - 1);
    for (int64_t j = 0; j < ncols; ++j) {                               // scatter by row
        for (int64_t k = static_cast<int64_t>(indptr[j]); k < static_cast<int64_t>(indptr[j + 1]); ++k) {
            int64_t dst = static_cast<int64_t>(next[indices[k]]++);
            out.data[static_cast<size_t>(dst)] = data[k];
            out.indices[static_cast<size_t>(dst)] = static_cast<Idx>(j);
        }
    }
    return out;
}

// CSR -> CSC storage transpose (orientation flip), value dtype preserved. A CSR (nrows x ncols) is
// byte-identical to a CSC of the transpose (ncols x nrows), so this is csc_to_csr with the dims swapped.
// Lets a binding normalize CSR counts to the CSC the viewer/DE kernels expect (Python/R get this from
// scipy/Matrix; JS-WASM has no such library, so it calls this).
template <class T, class Idx = int64_t>
inline CsxArrays<T, Idx> csr_to_csc(const T* data, const Idx* indices, const Idx* indptr,
                                    int64_t nrows, int64_t ncols) {
    return csc_to_csr<T, Idx>(data, indices, indptr, ncols, nrows);
}

// Zero-aware per-column mean/variance of a CSC matrix with `nrows` rows and `ncols` columns.
// `data` are the nnz values; `indptr` is length ncols+1. Implicit zeros are accounted for.
//
// When `lognorm` is set, each nonzero is transformed on the fly via log1p(value) before the
// moments are accumulated. This mirrors pagoda2's lazy normalized-view pattern: the dense
// log-normalized matrix is never materialized; the per-gene dispersion statistics are computed
// straight off the raw CSC store in one streaming pass. log1p(0)==0, so the implicit zeros stay
// implicit and the zero-aware variance correction is unchanged. The extra transcendental per
// nonzero makes the kernel compute-bound, which is where the column-parallel OpenMP schedule pays
// off (the plain-sum variant above is memory-bandwidth-bound and barely scales).
// `indptr` is int64 (length ncols+1) -- normalize with as_i64() when reading from a store.
// `data` is templated on the stored value dtype (float or double): a float32 measure (common in
// AnnData/Seurat) is read in place -- no widening copy -- while sums accumulate in double, so the
// kernel is memory-lean yet float64-accurate. n_threads is the explicit threading policy from the
// caller: 1 = serial, N = N threads, <=0 = the OpenMP default. Results are thread-count invariant.
// Optional per-cell depth normalization: when `depth` is given (length nrows, keyed by row index via
// `indices`), each nonzero is normalized to `x * depthScale / depth[row]` *before* the optional log1p
// -- i.e. the value of pagoda2's "plain" analysis view, computed straight off the raw store. With
// `population=true` the variance divides by `nrows` (E[X^2]-E[X]^2 over all cells, pagoda2's
// convention) rather than the sample `nrows-1`. Defaults (no depth, sample variance) are unchanged.
template <class T>
inline ColStats csc_col_mean_var(const T* data, const int64_t* indptr,
                                 int64_t ncols, int64_t nrows, int n_threads = 0,
                                 bool lognorm = false,
                                 const int64_t* indices = nullptr,
                                 const double* depth = nullptr, double depthScale = 1.0,
                                 bool population = false) {
    ColStats s;
    s.mean.assign((size_t)ncols, 0.0);
    s.var.assign((size_t)ncols, 0.0);
    s.nnz.assign((size_t)ncols, 0);
#ifdef _OPENMP
    if (n_threads > 0) omp_set_num_threads(n_threads);
    #pragma omp parallel for schedule(static) if (n_threads != 1)
#endif
    for (int64_t j = 0; j < ncols; ++j) {
        int64_t a = indptr[j], b = indptr[j + 1];
        double sum = 0.0;
        for (int64_t k = a; k < b; ++k) {
            double x = static_cast<double>(data[k]);     // widen the scalar, not the array
            if (depth) x = x * depthScale / depth[indices[k]];
            sum += lognorm ? std::log1p(x) : x;
        }
        double m = sum / (double)nrows;
        double ss = 0.0;
        for (int64_t k = a; k < b; ++k) {
            double x = static_cast<double>(data[k]);
            if (depth) x = x * depthScale / depth[indices[k]];
            double v = lognorm ? std::log1p(x) : x;
            double d = v - m;
            ss += d * d;
        }
        ss += m * m * (double)(nrows - (b - a));  // contribution of the implicit zeros
        const double denom = population ? (double)nrows : (double)(nrows - 1);
        s.mean[(size_t)j] = m;
        s.var[(size_t)j] = (denom > 0.0) ? ss / denom : 0.0;
        s.nnz[(size_t)j] = b - a;
    }
    return s;
}

// Dispatch csc_col_mean_var on a block's stored value dtype, so a float32 measure is reduced in
// place (no widening copy) -- the memory-lean path. `indptr` is the block-local pointer (int64).
inline ColStats col_mean_var_dispatch(const NdArray& data, const int64_t* indptr,
                                      int64_t ncols, int64_t nrows, int n_threads, bool lognorm,
                                      const int64_t* indices = nullptr, const double* depth = nullptr,
                                      double depthScale = 1.0, bool population = false) {
    const std::string& dt = data.dtype;
    if (dt == "<f8") return csc_col_mean_var(data.as<double>(),  indptr, ncols, nrows, n_threads, lognorm, indices, depth, depthScale, population);
    if (dt == "<f4") return csc_col_mean_var(data.as<float>(),   indptr, ncols, nrows, n_threads, lognorm, indices, depth, depthScale, population);
    if (dt == "<i4") return csc_col_mean_var(data.as<int32_t>(), indptr, ncols, nrows, n_threads, lognorm, indices, depth, depthScale, population);
    if (dt == "<i8") return csc_col_mean_var(data.as<int64_t>(), indptr, ncols, nrows, n_threads, lognorm, indices, depth, depthScale, population);
    throw std::runtime_error("col_mean_var: unsupported data dtype " + dt);
}

// Read a contiguous gene (column) range [g_lo, g_hi) of a CSC measure field as its own CSC arrays,
// touching only the `data`/`indices` chunks that overlap the range (via read_array_range). This is
// the general bounded block-read primitive a consumer drives to build out-of-core ops over an L*
// store *without* each reduction having to live in libstar -- e.g. pagoda2 reads gene blocks and
// applies its own view kernels. `data`/`indices` keep their stored dtype; `indptr` is rebased local.
struct CscBlock {
    NdArray data, indices;            // the block's nonzeros (stored dtype preserved)
    std::vector<int64_t> indptr;      // length (g_hi-g_lo+1), local (starts at 0)
    int64_t nrows = 0, ncols = 0;     // nrows = #cells; ncols = g_hi-g_lo (genes in the block)
};
inline CscBlock read_csc_block(const fs::path& field_group, int64_t g_lo, int64_t g_hi) {
    json m = read_group_attrs(field_group)["lstar"];
    if (opt_str(m, "encoding") != "csc")
        throw std::runtime_error("read_csc_block: field is not CSC (gene-major)");
    auto shape = m["shape"].get<std::vector<int64_t>>();
    const int64_t nrows = shape[0], ngenes = shape[1];
    g_lo = std::max<int64_t>(0, g_lo);
    g_hi = std::min<int64_t>(ngenes, g_hi);
    if (g_hi < g_lo) g_hi = g_lo;
    std::vector<int64_t> ip = as_i64(read_array(field_group / "indptr"));   // whole (small)
    const int64_t lo = ip[g_lo], hi = ip[g_hi];
    CscBlock b;
    b.nrows = nrows;
    b.ncols = g_hi - g_lo;
    b.data = read_array_range(field_group / "data", lo, hi);                // only overlapping chunks
    b.indices = read_array_range(field_group / "indices", lo, hi);
    b.indptr.resize((size_t)(b.ncols + 1));
    for (int64_t j = 0; j <= b.ncols; ++j) b.indptr[(size_t)j] = ip[g_lo + j] - lo;
    return b;
}

// A 1-D chunked zarr array read one element-range at a time, caching the last decoded chunk. Because
// a gather issues ranges in ascending order, each chunk is decoded at most once across all ranges --
// the key to an efficient scattered-column gather (vs decoding a chunk once per touched column).
struct ChunkReader {
    zarr::Array arr;
    int64_t cs = 1, n = 0;
    size_t dsz = 1;
    std::string dtype;
    int64_t cached = -1;
    std::vector<uint8_t> buf;                              // the decoded chunk (cs*dsz bytes)

    explicit ChunkReader(const fs::path& d)
        : arr(zarr::Array::open(std::make_shared<zarr::FilesystemStore>(d, false), "")) {
        cs = static_cast<int64_t>(arr.meta().chunk_shape.at(0));
        n = static_cast<int64_t>(arr.meta().shape.at(0));
        dsz = arr.meta().dtype.itemsize;
        dtype = zarr::v2::emit_dtype(arr.meta().dtype, /*big_endian=*/false);
    }
    const uint8_t* chunk(int64_t ci) {
        if (ci != cached) {                                // libzarr decodes + fill-pads the whole chunk
            buf = arr.read_chunk({static_cast<std::uint64_t>(ci)});
            cached = ci;
        }
        return buf.data();
    }
    // copy source elements [lo, hi) to `out` starting at element `dst` (out is dsz-byte elements)
    void copy_range(uint8_t* out, int64_t dst, int64_t lo, int64_t hi) {
        int64_t k = 0;
        for (int64_t p = lo; p < hi; ) {
            int64_t ci = p / cs;
            int64_t cend = std::min<int64_t>(hi, (ci + 1) * cs);
            const uint8_t* cb = chunk(ci);
            std::memcpy(out + static_cast<size_t>(dst + k) * dsz,
                        cb + static_cast<size_t>(p - ci * cs) * dsz,
                        static_cast<size_t>(cend - p) * dsz);
            k += cend - p;
            p = cend;
        }
    }
};

// Gather an arbitrary set of gene (column) indices of a CSC measure, decoding each touched data/index
// chunk at most once (the efficient scattered-subset read; `read_csc_block`-per-run re-decodes chunks
// for scattered columns). `cols` must be sorted ascending and unique (the R wrapper enforces this and
// restores the caller's order). Returns the gathered columns as their own CSC arrays.
inline CscBlock read_csc_cols(const fs::path& field_group, const std::vector<int64_t>& cols) {
    json m = read_group_attrs(field_group)["lstar"];
    if (opt_str(m, "encoding") != "csc")
        throw std::runtime_error("read_csc_cols: field is not CSC (gene-major)");
    auto shape = m["shape"].get<std::vector<int64_t>>();
    std::vector<int64_t> ip = as_i64(read_array(field_group / "indptr"));
    CscBlock b;
    b.nrows = shape[0];
    b.ncols = static_cast<int64_t>(cols.size());
    b.indptr.assign(static_cast<size_t>(b.ncols + 1), 0);
    for (int64_t j = 0; j < b.ncols; ++j)
        b.indptr[(size_t)(j + 1)] = b.indptr[(size_t)j] + (ip[cols[(size_t)j] + 1] - ip[cols[(size_t)j]]);
    const int64_t total = b.indptr.back();
    ChunkReader dr(field_group / "data"), ir(field_group / "indices");
    b.data.dtype = dr.dtype;
    b.indices.dtype = ir.dtype;
    b.data.shape = {total};
    b.indices.shape = {total};
    b.data.bytes.assign(static_cast<size_t>(total) * dr.dsz, 0);
    b.indices.bytes.assign(static_cast<size_t>(total) * ir.dsz, 0);
    for (int64_t j = 0; j < b.ncols; ++j) {
        int64_t g = cols[(size_t)j], lo = ip[g], hi = ip[g + 1], dst = b.indptr[(size_t)j];
        dr.copy_range(b.data.bytes.data(), dst, lo, hi);
        ir.copy_range(b.indices.bytes.data(), dst, lo, hi);
    }
    return b;
}

// Zero-aware per-column mean/variance of a CSC measure read straight from a store *block-by-block*,
// so the whole matrix never lands in memory -- the C++/R counterpart of Python's `stream_col_stats`.
// `field_group` is the field's directory (…/fields/<name>). `indptr` (small) is read whole; for each
// column block only that block's slice of `data` is read (via read_array_range, touching just the
// overlapping chunks) and reduced. Blocking by column is *exact* (each column's nonzeros sit fully
// in one block), so per-column results need no cross-block merge. Bounded memory requires a chunked
// `data` array (a store written with chunk_elems set, e.g. by streamed conversion). n_threads is the
// usual policy (1=serial, N=N threads, <=0=OpenMP default), applied within each block's reduction.
// `depth` (optional, length nrows) + `depthScale` + `population` add pagoda2's "plain" view in the
// same streamed pass: each nonzero becomes log1p(x*depthScale/depth[row]) and the variance divides by
// nrows. When depth is given the block's `indices` are read too (to recover each nonzero's row/depth),
// so a normalized pass moves ~2x the bytes of a raw one -- still bounded, still one streaming pass.
inline ColStats stream_csc_col_mean_var(const fs::path& field_group, int64_t block = 4096,
                                        int n_threads = 0, bool lognorm = false,
                                        const std::vector<double>* depth = nullptr,
                                        double depthScale = 1.0, bool population = false) {
    json m = read_group_attrs(field_group)["lstar"];
    if (opt_str(m, "encoding") != "csc")
        throw std::runtime_error("stream_csc_col_mean_var: field is not CSC (need gene-major)");
    auto shape = m["shape"].get<std::vector<int64_t>>();
    const int64_t nrows = shape[0], ncols = shape[1];
    std::vector<int64_t> indptr = as_i64(read_array(field_group / "indptr"));
    if (block <= 0) block = 4096;
    const double* depthp = (depth && !depth->empty()) ? depth->data() : nullptr;
    if (depthp && (int64_t)depth->size() != nrows)
        throw std::runtime_error("stream_csc_col_mean_var: depth length must equal nrows");
    ColStats out;
    out.mean.assign((size_t)ncols, 0.0);
    out.var.assign((size_t)ncols, 0.0);
    out.nnz.assign((size_t)ncols, 0);
    for (int64_t a = 0; a < ncols; a += block) {
        const int64_t b = std::min(a + block, ncols);
        const int64_t lo = indptr[a], hi = indptr[b];
        NdArray dblk = read_array_range(field_group / "data", lo, hi);     // only this block's nonzeros
        std::vector<int64_t> idxblk;
        const int64_t* idxp = nullptr;
        if (depthp) { idxblk = as_i64(read_array_range(field_group / "indices", lo, hi)); idxp = idxblk.data(); }
        std::vector<int64_t> liptr((size_t)(b - a + 1));
        for (int64_t j = 0; j <= b - a; ++j) liptr[(size_t)j] = indptr[a + j] - lo;
        ColStats s = col_mean_var_dispatch(dblk, liptr.data(), b - a, nrows, n_threads, lognorm,
                                           idxp, depthp, depthScale, population);
        std::copy(s.mean.begin(), s.mean.end(), out.mean.begin() + a);
        std::copy(s.var.begin(), s.var.end(), out.var.begin() + a);
        std::copy(s.nnz.begin(), s.nnz.end(), out.nnz.begin() + a);
    }
    return out;
}

// Per-group sufficient statistics over a CSC measure (cells x genes, gene-major): for each gene
// column, accumulate sum / sumsq / n_expr per cell-group. `group_of_cell` (length nrows) maps each
// cell to a group in [0,ngroups) (or <0 to skip). With `lognorm`, nonzeros are log1p'd on the fly.
// Returns flat (ngroups x ncols) arrays. This is the exporter workhorse for the viewer profile's
// cluster stats + grouped heatmaps + cluster markers — computed straight off the raw store, no
// dense matrix. Parallel over genes (columns): each column writes distinct (g*ncols+j) slots, so
// there is no cross-thread race. Accumulates in double; templated on the stored dtype (memory-lean).
struct GroupStats { std::vector<double> sum, sumsq, n_expr; int64_t ngroups = 0, ngenes = 0; };
template <class T>
inline GroupStats csc_col_sum_by_group(const T* data, const int64_t* indptr, const int64_t* indices,
                                       int64_t /*nrows*/, int64_t ncols, const int* group_of_cell,
                                       int ngroups, bool lognorm = false, int n_threads = 0) {
    GroupStats s; s.ngroups = ngroups; s.ngenes = ncols;
    const size_t sz = (size_t)ngroups * (size_t)ncols;
    s.sum.assign(sz, 0.0); s.sumsq.assign(sz, 0.0); s.n_expr.assign(sz, 0.0);
#ifdef _OPENMP
    if (n_threads > 0) omp_set_num_threads(n_threads);
    #pragma omp parallel for schedule(static) if (n_threads != 1)
#endif
    for (int64_t j = 0; j < ncols; ++j) {
        for (int64_t k = indptr[j]; k < indptr[j + 1]; ++k) {
            int g = group_of_cell[indices[k]];
            if (g < 0 || g >= ngroups) continue;
            double v = static_cast<double>(data[k]);
            if (lognorm) v = std::log1p(v);
            size_t idx = (size_t)g * (size_t)ncols + (size_t)j;
            s.sum[idx] += v; s.sumsq[idx] += v * v; s.n_expr[idx] += 1.0;
        }
    }
    return s;
}

// ---- viewer@0.1 recipe kernels (docs/format.md "The viewer profile") --------------------------
// One implementation, bound to R / Python / WASM, so a precomputed navigator equals the viewer's
// on-the-fly value ("prepped == live"). See the profile spec for field definitions and orientation.

// 1-vs-rest marker table from per-(group,gene) sufficient stats (the `sum`/`n_expr` of
// csc_col_sum_by_group computed over log1p). Inputs S, NE are GROUP-major (g*ngenes + gene); `nper`
// = group sizes (length ngroups); `ncells` = total cells. Outputs are GENE-major (gene*ngroups + g)
// -- the spec's ng x K orientation: lfc[gene,g] = mean_log_in_g - mean_log_in_rest;
// padj[gene,g] = clip(exp(-|lfc|*sqrt(nexpr_g + 1)), 1e-12, 1). A fast surrogate, not a calibrated test.
struct Markers { std::vector<double> lfc, padj; int64_t ngenes = 0; int ngroups = 0; };
inline Markers markers_one_vs_rest(const double* S, const double* NE, const int64_t* nper,
                                   int ngroups, int64_t ngenes, int64_t ncells) {
    Markers m; m.ngenes = ngenes; m.ngroups = ngroups;
    m.lfc.assign((size_t)ngenes * ngroups, 0.0);
    m.padj.assign((size_t)ngenes * ngroups, 1.0);
    std::vector<double> grand((size_t)ngenes, 0.0);
    for (int g = 0; g < ngroups; ++g)
        for (int64_t j = 0; j < ngenes; ++j) grand[j] += S[(size_t)g * ngenes + j];
    for (int g = 0; g < ngroups; ++g) {
        double ng1 = std::max<double>(nper[g], 1.0);
        double nr  = std::max<double>((double)ncells - (double)nper[g], 1.0);
        for (int64_t j = 0; j < ngenes; ++j) {
            double sg = S[(size_t)g * ngenes + j];
            double d  = sg / ng1 - (grand[j] - sg) / nr;
            size_t o  = (size_t)j * ngroups + g;
            m.lfc[o]  = d;
            double p  = std::exp(-std::fabs(d * std::sqrt(NE[(size_t)g * ngenes + j] + 1.0)));
            m.padj[o] = p < 1e-12 ? 1e-12 : (p > 1.0 ? 1.0 : p);
        }
    }
    return m;
}

// Lanczos log-gamma + Numerical-Recipes continued-fraction for the regularized incomplete beta, in
// log space (the upper F-tail underflows to 0 for strongly overdispersed genes). reg_incbeta_log
// returns log I_x(a,b).
inline double od_gammln_(double xx) {
    static const double cof[6] = {76.18009172947146, -86.50532032941677, 24.01409824083091,
        -1.231739572450155, 0.1208650973866179e-2, -0.5395239384953e-5};
    double x = xx, y = xx, tmp = x + 5.5; tmp -= (x + 0.5) * std::log(tmp);
    double ser = 1.000000000190015;
    for (int j = 0; j < 6; j++) { y += 1.0; ser += cof[j] / y; }
    return -tmp + std::log(2.5066282746310005 * ser / x);
}
inline double od_betacf_(double a, double b, double x) {
    const int MAXIT = 300; const double EPS = 3e-12, FPMIN = 1e-300;
    double qab = a + b, qap = a + 1.0, qam = a - 1.0;
    double c = 1.0, d = 1.0 - qab * x / qap; if (std::fabs(d) < FPMIN) d = FPMIN; d = 1.0 / d; double h = d;
    for (int mm = 1; mm <= MAXIT; mm++) {
        int m2 = 2 * mm; double aa = mm * (b - mm) * x / ((qam + m2) * (a + m2));
        d = 1.0 + aa * d; if (std::fabs(d) < FPMIN) d = FPMIN; c = 1.0 + aa / c; if (std::fabs(c) < FPMIN) c = FPMIN;
        d = 1.0 / d; h *= d * c;
        aa = -(a + mm) * (qab + mm) * x / ((a + m2) * (qap + m2));
        d = 1.0 + aa * d; if (std::fabs(d) < FPMIN) d = FPMIN; c = 1.0 + aa / c; if (std::fabs(c) < FPMIN) c = FPMIN;
        d = 1.0 / d; double del = d * c; h *= del;
        if (std::fabs(del - 1.0) < EPS) break;
    }
    return h;
}
inline double reg_incbeta_log(double a, double b, double x) {
    if (x <= 0.0) return -700.0;                          // I_0(a,b) = 0
    if (x >= 1.0) return 0.0;                              // I_1(a,b) = 1
    double bt_log = od_gammln_(a + b) - od_gammln_(a) - od_gammln_(b)
                    + a * std::log(x) + b * std::log1p(-x);
    if (x < (a + 1.0) / (a + b + 2.0)) return bt_log + std::log(od_betacf_(a, b, x) / a);   // small-x: stable
    double v = 1.0 - std::exp(bt_log) * od_betacf_(b, a, 1.0 - x) / b;
    return std::log(v > 1e-300 ? v : 1e-300);
}

// tricube local-linear LOWESS, evaluated at `nanchor` evenly-spaced anchors and linearly interpolated
// back to each xs (constant edge extrapolation). A faithful port of the viewer's prep lowess.
inline std::vector<double> lowess_predict(const std::vector<double>& xs, const std::vector<double>& ys,
                                          double span = 0.3, int nanchor = 200) {
    const int64_t n = (int64_t)xs.size();
    std::vector<double> out((size_t)n, 0.0);
    if (n < 3) { double my = 0.0; for (double v : ys) my += v; my = n ? my / n : 0.0;
        for (auto& o : out) o = my; return out; }
    std::vector<int64_t> ord((size_t)n); for (int64_t i = 0; i < n; i++) ord[(size_t)i] = i;
    std::sort(ord.begin(), ord.end(), [&](int64_t a, int64_t b){ return xs[(size_t)a] < xs[(size_t)b]; });
    std::vector<double> sx((size_t)n), sy((size_t)n);
    for (int64_t i = 0; i < n; i++) { sx[(size_t)i] = xs[(size_t)ord[(size_t)i]]; sy[(size_t)i] = ys[(size_t)ord[(size_t)i]]; }
    const int64_t win = std::max<int64_t>(2, (int64_t)(span * (double)n));
    std::vector<double> ax((size_t)nanchor), ay((size_t)nanchor);
    for (int a = 0; a < nanchor; a++) {
        double x0 = sx[0] + (sx[(size_t)n - 1] - sx[0]) * (double)a / (double)(nanchor - 1);
        int64_t lo = (int64_t)(std::lower_bound(sx.begin(), sx.end(), x0) - sx.begin());
        int64_t l = std::max<int64_t>(0, lo - (win >> 1)); int64_t r = std::min<int64_t>(n, l + win); l = std::max<int64_t>(0, r - win);
        double maxd = 1e-9; for (int64_t i = l; i < r; i++) maxd = std::max(maxd, std::fabs(sx[(size_t)i] - x0));
        double sw = 0, swx = 0, swy = 0, swxx = 0, swxy = 0;
        for (int64_t i = l; i < r; i++) { double dd = std::fabs(sx[(size_t)i] - x0) / maxd, w = std::pow(1.0 - dd * dd * dd, 3);
            sw += w; swx += w * sx[(size_t)i]; swy += w * sy[(size_t)i]; swxx += w * sx[(size_t)i] * sx[(size_t)i]; swxy += w * sx[(size_t)i] * sy[(size_t)i]; }
        double den = sw * swxx - swx * swx;
        ay[(size_t)a] = (std::fabs(den) < 1e-12) ? swy / sw
                        : ((swy - (sw * swxy - swx * swy) / den * swx) / sw + (sw * swxy - swx * swy) / den * x0);
        ax[(size_t)a] = x0;
    }
    for (int64_t i = 0; i < n; i++) {                     // np.interp: constant-extrapolated at the edges
        double x = xs[(size_t)i];
        if (x <= ax[0]) { out[(size_t)i] = ay[0]; continue; }
        if (x >= ax[(size_t)nanchor - 1]) { out[(size_t)i] = ay[(size_t)nanchor - 1]; continue; }
        int64_t k = (int64_t)(std::lower_bound(ax.begin(), ax.end(), x) - ax.begin());
        out[(size_t)i] = ay[(size_t)k - 1] + (ay[(size_t)k] - ay[(size_t)k - 1]) * (x - ax[(size_t)k - 1]) / (ax[(size_t)k] - ax[(size_t)k - 1]);
    }
    return out;
}

// Per-gene overdispersion score (pagoda2 adjustVariance, D4): residual of log(var) about a lowess fit
// of log(var)~log(mean) over log1p(counts), scored by the upper-tail variance-ratio F-test:
// od = -log P(F > exp(residual); df1=df2=nobs), nobs = #expressing cells. Genes with nobs<3 or
// mean/var<=0 score 0. `mean`/`var`/`nobs` are per-gene (length ngenes) -- e.g. from csc_col_mean_var.
inline std::vector<double> overdispersion(const double* mean, const double* var, const int64_t* nobs,
                                          int64_t ngenes, double span = 0.3, int nanchor = 200) {
    std::vector<double> od((size_t)ngenes, 0.0);
    std::vector<int64_t> ok; std::vector<double> xs, ys;
    for (int64_t j = 0; j < ngenes; ++j)
        if (nobs[j] >= 3 && mean[j] > 0.0 && var[j] > 0.0) {
            ok.push_back(j); xs.push_back(std::log(mean[j])); ys.push_back(std::log(var[j]));
        }
    if ((int64_t)ok.size() <= 10) return od;
    std::vector<double> trend = lowess_predict(xs, ys, span, nanchor);
    for (size_t k = 0; k < ok.size(); ++k) {
        double f = std::exp(ys[k] - trend[k]);            // variance ratio
        double a = (double)nobs[ok[(size_t)k]] / 2.0;     // df/2
        double od_k = -reg_incbeta_log(a, a, 1.0 / (1.0 + f));   // -log P(F>f; nobs,nobs)
        od[(size_t)ok[(size_t)k]] = std::isfinite(od_k) ? od_k : 0.0;
    }
    return od;
}

// Canonical xy -> Hilbert index on an N x N curve (N a power of two); the reflection uses N-1.
inline int64_t hilbert_xy2d(int64_t N, int64_t x, int64_t y) {
    int64_t d = 0;
    for (int64_t s = N >> 1; s > 0; s >>= 1) {
        int64_t rx = (x & s) > 0 ? 1 : 0, ry = (y & s) > 0 ? 1 : 0;
        d += s * s * ((3 * rx) ^ ry);
        if (ry == 0) { if (rx == 1) { x = N - 1 - x; y = N - 1 - y; } int64_t t = x; x = y; y = t; }
    }
    return d;
}
// Per-cell Hilbert index over an N x N grid of the min-max-scaled 2-D embedding (`emb` is row-major
// ncells x 2). Used as the secondary (within-cluster spatial-locality) key for cell ordering.
inline std::vector<int64_t> hilbert_index(const double* emb, int64_t ncells, int64_t N = 1024) {
    std::vector<int64_t> out((size_t)ncells, 0);
    double xmin = 1e300, xmax = -1e300, ymin = 1e300, ymax = -1e300;
    for (int64_t i = 0; i < ncells; i++) { double x = emb[2 * i], y = emb[2 * i + 1];
        xmin = std::min(xmin, x); xmax = std::max(xmax, x); ymin = std::min(ymin, y); ymax = std::max(ymax, y); }
    double xr = (xmax - xmin); if (xr == 0) xr = 1.0; double yr = (ymax - ymin); if (yr == 0) yr = 1.0;
    for (int64_t i = 0; i < ncells; i++) {
        int64_t gx = (int64_t)std::floor((emb[2 * i] - xmin) / xr * (double)(N - 1));     gx = std::min<int64_t>(N - 1, std::max<int64_t>(0, gx));
        int64_t gy = (int64_t)std::floor((emb[2 * i + 1] - ymin) / yr * (double)(N - 1)); gy = std::min<int64_t>(N - 1, std::max<int64_t>(0, gy));
        out[(size_t)i] = hilbert_xy2d(N, gx, gy);
    }
    return out;
}
// pos_of[cell] = physical row, after a stable sort by (primary cluster code, optional secondary key,
// cell index). `secondary` (e.g. hilbert_index) may be null for a cluster-only order. The matrix is
// then physically reordered so each cluster is byte-contiguous (the spec's `counts_cellmajor_order`).
inline std::vector<int64_t> cell_order_pos(const int* primary_code, const int64_t* secondary,
                                           int64_t ncells) {
    std::vector<int64_t> perm((size_t)ncells); for (int64_t i = 0; i < ncells; i++) perm[(size_t)i] = i;
    std::stable_sort(perm.begin(), perm.end(), [&](int64_t a, int64_t b) {
        if (primary_code[a] != primary_code[b]) return primary_code[a] < primary_code[b];
        if (secondary && secondary[a] != secondary[b]) return secondary[a] < secondary[b];
        return a < b;
    });
    std::vector<int64_t> pos((size_t)ncells);
    for (int64_t p = 0; p < ncells; p++) pos[(size_t)perm[(size_t)p]] = p;
    return pos;
}
// Canonical viewer cell order -- the SINGLE source of truth every binding (Python/R/WASM) calls, so all
// surfaces emit a byte-identical `counts_cellmajor_order`. pos_of[cell] = physical row, stable-sorted by
// (cluster `primary_code`, then Hilbert index of the embedding when `emb2d` is given, else cell index).
// `emb2d` is row-major ncells x 2 (or null for a cluster-only order when no embedding is available).
inline std::vector<int64_t> viewer_cell_order(const int* primary_code, const double* emb2d,
                                              int64_t ncells, int64_t grid = 1024) {
    if (emb2d) {
        std::vector<int64_t> hil = hilbert_index(emb2d, ncells, grid);
        return cell_order_pos(primary_code, hil.data(), ncells);
    }
    return cell_order_pos(primary_code, nullptr, ncells);
}

// Streaming per-(group, gene) SUM of a CSC measure read straight from a store, with the same optional
// "plain" view (depth-normalize + log1p) as stream_csc_col_mean_var -- the fused counterpart of
// pagoda2's colSumByFacView (pseudobulk / marker sums), one threaded streamed pass, no R marshalling.
// `group_of_cell` (length nrows) maps each cell to a bucket in [0, ngroups) (out of range -> skipped;
// pass NA cells as bucket 0 to mirror pagoda2's `<NA>` row). Output is flat row-major (g*ncols + gene).
template <class T>
inline void sum_by_group_block(const T* data, const int64_t* lindptr, int64_t bcols, int64_t col0,
                               const int64_t* indices, int64_t ncols_total, const int* group_of_cell,
                               int ngroups, bool lognorm, const double* depth, double depthScale,
                               double* out, int n_threads) {
#ifdef _OPENMP
    if (n_threads > 0) omp_set_num_threads(n_threads);
    #pragma omp parallel for schedule(static) if (n_threads != 1)
#endif
    for (int64_t j = 0; j < bcols; ++j) {            // each j writes distinct out columns -> race-free
        for (int64_t k = lindptr[j]; k < lindptr[j + 1]; ++k) {
            const int64_t row = indices[k];
            const int g = group_of_cell[row];
            if (g < 0 || g >= ngroups) continue;
            double v = static_cast<double>(data[k]);
            if (depth) v = v * depthScale / depth[row];
            if (lognorm) v = std::log1p(v);
            out[(size_t)g * (size_t)ncols_total + (size_t)(col0 + j)] += v;
        }
    }
}
inline void sum_by_group_block_dispatch(const NdArray& data, const int64_t* lindptr, int64_t bcols,
        int64_t col0, const int64_t* indices, int64_t ncols_total, const int* group_of_cell,
        int ngroups, bool lognorm, const double* depth, double depthScale, double* out, int n_threads) {
    const std::string& dt = data.dtype;
    if (dt == "<f8") sum_by_group_block(data.as<double>(),  lindptr, bcols, col0, indices, ncols_total, group_of_cell, ngroups, lognorm, depth, depthScale, out, n_threads);
    else if (dt == "<f4") sum_by_group_block(data.as<float>(),   lindptr, bcols, col0, indices, ncols_total, group_of_cell, ngroups, lognorm, depth, depthScale, out, n_threads);
    else if (dt == "<i4") sum_by_group_block(data.as<int32_t>(), lindptr, bcols, col0, indices, ncols_total, group_of_cell, ngroups, lognorm, depth, depthScale, out, n_threads);
    else if (dt == "<i8") sum_by_group_block(data.as<int64_t>(), lindptr, bcols, col0, indices, ncols_total, group_of_cell, ngroups, lognorm, depth, depthScale, out, n_threads);
    else throw std::runtime_error("sum_by_group: unsupported data dtype " + dt);
}
inline std::vector<double> stream_csc_col_sum_by_group(const fs::path& field_group,
        const std::vector<int>& group_of_cell, int ngroups, bool lognorm,
        const std::vector<double>* depth, double depthScale, int64_t block, int n_threads) {
    json m = read_group_attrs(field_group)["lstar"];
    if (opt_str(m, "encoding") != "csc")
        throw std::runtime_error("stream_csc_col_sum_by_group: field is not CSC (need gene-major)");
    auto shape = m["shape"].get<std::vector<int64_t>>();
    const int64_t nrows = shape[0], ncols = shape[1];
    if ((int64_t)group_of_cell.size() != nrows)
        throw std::runtime_error("stream_csc_col_sum_by_group: group length must equal nrows");
    const double* depthp = (depth && !depth->empty()) ? depth->data() : nullptr;
    if (depthp && (int64_t)depth->size() != nrows)
        throw std::runtime_error("stream_csc_col_sum_by_group: depth length must equal nrows");
    std::vector<int64_t> indptr = as_i64(read_array(field_group / "indptr"));
    if (block <= 0) block = 4096;
    std::vector<double> out((size_t)ngroups * (size_t)ncols, 0.0);
    for (int64_t a = 0; a < ncols; a += block) {
        const int64_t b = std::min(a + block, ncols);
        const int64_t lo = indptr[a], hi = indptr[b];
        NdArray dblk = read_array_range(field_group / "data", lo, hi);
        std::vector<int64_t> idxblk = as_i64(read_array_range(field_group / "indices", lo, hi));
        std::vector<int64_t> liptr((size_t)(b - a + 1));
        for (int64_t j = 0; j <= b - a; ++j) liptr[(size_t)j] = indptr[a + j] - lo;
        sum_by_group_block_dispatch(dblk, liptr.data(), b - a, a, idxblk.data(), ncols,
                                    group_of_cell.data(), ngroups, lognorm, depthp, depthScale,
                                    out.data(), n_threads);
    }
    return out;
}

// Subsample differential-expression ranker over a CSR submatrix (sampled cells x genes, cell-major):
// `membership` (length nrows) labels each row 0 (group A), 1 (group B), or <0 (skip). `indices` are
// gene ids, `indptr` is length nrows+1. Returns per-gene log1p group means + log-fold-change
// (lfc = meanA - meanB); the caller ranks by |lfc|. This is the constant-cost selection-DE kernel:
// fed a few hundred rows of the cell-major DE panel it gives ranking-grade results without a full
// matrix pass. (AUC/Wilcoxon is a later refinement; fold-change ranks the dominant signal cheaply.)
struct DERank { std::vector<double> meanA, meanB, lfc; int64_t nA = 0, nB = 0; };
template <class T>
inline DERank subsample_de_rank(const T* data, const int64_t* indptr, const int64_t* indices,
                                int64_t nrows, int64_t ngenes, const int* membership, bool lognorm = true) {
    DERank r;
    std::vector<double> sumA((size_t)ngenes, 0.0), sumB((size_t)ngenes, 0.0);
    int64_t nA = 0, nB = 0;
    for (int64_t row = 0; row < nrows; ++row) {
        int m = membership[row];
        if (m < 0) continue;
        if (m == 0) ++nA; else ++nB;
        for (int64_t k = indptr[row]; k < indptr[row + 1]; ++k) {
            double v = static_cast<double>(data[k]);
            if (lognorm) v = std::log1p(v);
            (m == 0 ? sumA : sumB)[(size_t)indices[k]] += v;
        }
    }
    r.nA = nA; r.nB = nB;
    r.meanA.resize((size_t)ngenes); r.meanB.resize((size_t)ngenes); r.lfc.resize((size_t)ngenes);
    const double invA = 1.0 / (double)std::max<int64_t>(nA, 1), invB = 1.0 / (double)std::max<int64_t>(nB, 1);
    for (int64_t g = 0; g < ngenes; ++g) {
        double ma = sumA[(size_t)g] * invA, mb = sumB[(size_t)g] * invB;
        r.meanA[(size_t)g] = ma; r.meanB[(size_t)g] = mb; r.lfc[(size_t)g] = ma - mb;
    }
    return r;
}

}  // namespace lstar
