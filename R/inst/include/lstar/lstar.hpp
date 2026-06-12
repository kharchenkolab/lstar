// lstar — C++ core (libstar): the L* model + a Zarr v2 reader/writer.
//
// Reads multi-chunk arrays over an arbitrary chunk grid, little-endian, with optional gzip/zlib
// chunk compression (when built with zlib; define LSTAR_HAVE_ZLIB). Edge chunks are full-size and
// fill-padded per the v2 spec; missing chunks read as fill_value 0. The writer emits single-chunk
// arrays plus a consolidated .zmetadata. Blosc and Zarr v3 sharding are still planned (see
// misc/plan1.md). Written fresh (MIT); pagoda2's IO is GPL-3 and was used only as a schema reference.
#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "nlohmann/json.hpp"

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
    json za = read_json(dir / ".zarray");
    if (za.contains("order") && za["order"].is_string() && za["order"].get<std::string>() != "C")
        throw std::runtime_error("F-order array unsupported: " + dir.string());
    NdArray a;
    a.dtype = za["dtype"].get<std::string>();
    a.shape = za["shape"].get<std::vector<int64_t>>();
    auto chunks = za["chunks"].get<std::vector<int64_t>>();
    json compressor = za.contains("compressor") ? za["compressor"] : json(nullptr);
    const size_t ndim = a.shape.size();
    const size_t dsz = dtype_size(a.dtype);
    a.bytes.assign(static_cast<size_t>(a.nelem()) * dsz, 0);  // fill_value 0
    if (ndim == 0 || a.nelem() == 0) return a;

    std::vector<int64_t> grid(ndim), ostride(ndim, 1), cstride(ndim, 1);
    int64_t nchunks = 1, chunk_elems = 1;
    for (size_t i = 0; i < ndim; ++i) {
        grid[i] = (a.shape[i] + chunks[i] - 1) / chunks[i];
        nchunks *= grid[i];
        chunk_elems *= chunks[i];
    }
    for (int i = static_cast<int>(ndim) - 2; i >= 0; --i) {
        ostride[i] = ostride[i + 1] * a.shape[i + 1];
        cstride[i] = cstride[i + 1] * chunks[i + 1];
    }

    std::vector<int64_t> cc(ndim, 0);
    for (int64_t ci = 0; ci < nchunks; ++ci) {
        std::string key = std::to_string(cc[0]);
        for (size_t i = 1; i < ndim; ++i) key += "." + std::to_string(cc[i]);
        fs::path cf = dir / key;
        if (fs::exists(cf)) {
            std::vector<uint8_t> raw =
                decode_chunk(compressor, read_bytes(cf), static_cast<size_t>(chunk_elems) * dsz);
            raw.resize(static_cast<size_t>(chunk_elems) * dsz);  // pad a short final chunk

            std::vector<int64_t> base(ndim), ext(ndim);
            for (size_t i = 0; i < ndim; ++i) {
                base[i] = cc[i] * chunks[i];
                ext[i] = std::min<int64_t>(chunks[i], a.shape[i] - base[i]);
            }
            const int64_t runlen = ext[ndim - 1];
            int64_t nruns = 1;
            for (size_t i = 0; i + 1 < ndim; ++i) nruns *= ext[i];
            std::vector<int64_t> idx(ndim, 0);
            for (int64_t r = 0; r < nruns; ++r) {
                int64_t src = 0, dst = base[ndim - 1];
                for (size_t i = 0; i + 1 < ndim; ++i) {
                    src += idx[i] * cstride[i];
                    dst += (base[i] + idx[i]) * ostride[i];
                }
                std::memcpy(a.bytes.data() + static_cast<size_t>(dst) * dsz,
                            raw.data() + static_cast<size_t>(src) * dsz,
                            static_cast<size_t>(runlen) * dsz);
                for (int i = static_cast<int>(ndim) - 2; i >= 0; --i) {
                    if (++idx[i] < ext[i]) break;
                    idx[i] = 0;
                }
            }
        }
        for (int i = static_cast<int>(ndim) - 1; i >= 0; --i) {
            if (++cc[i] < grid[i]) break;
            cc[i] = 0;
        }
    }
    return a;
}

// Read elements [lo, hi) of a 1-D zarr array, decoding only the chunks that overlap the range.
// This is the bounded-memory primitive behind the blocked reader: a column block of a CSC measure
// touches only its slice of `data`/`indices`, so a per-gene reduction never reads the whole matrix.
// (When the array was written as a single chunk -- the unchunked default -- this still decodes that
// one chunk, so bounded streaming needs a chunked store, e.g. one written with chunk_elems set.)
inline NdArray read_array_range(const fs::path& dir, int64_t lo, int64_t hi) {
    json za = read_json(dir / ".zarray");
    if (za.contains("order") && za["order"].is_string() && za["order"].get<std::string>() != "C")
        throw std::runtime_error("F-order array unsupported: " + dir.string());
    auto shape = za["shape"].get<std::vector<int64_t>>();
    if (shape.size() != 1) throw std::runtime_error("read_array_range: 1-D only: " + dir.string());
    const int64_t n = shape[0];
    const int64_t cs = za["chunks"].get<std::vector<int64_t>>()[0];
    json compressor = za.contains("compressor") ? za["compressor"] : json(nullptr);
    NdArray a;
    a.dtype = za["dtype"].get<std::string>();
    const size_t dsz = dtype_size(a.dtype);
    lo = std::max<int64_t>(0, lo);
    hi = std::min<int64_t>(n, hi);
    const int64_t len = std::max<int64_t>(0, hi - lo);
    a.shape = {len};
    a.bytes.assign(static_cast<size_t>(len) * dsz, 0);   // fill_value 0 for any missing chunk
    if (len == 0) return a;
    for (int64_t ci = lo / cs; ci <= (hi - 1) / cs; ++ci) {
        const int64_t cstart = ci * cs;
        const int64_t s = std::max(lo, cstart), e = std::min(hi, cstart + cs);  // overlap [s,e)
        fs::path cf = dir / std::to_string(ci);
        if (e <= s || !fs::exists(cf)) continue;
        std::vector<uint8_t> raw =
            decode_chunk(compressor, read_bytes(cf), static_cast<size_t>(cs) * dsz);
        raw.resize(static_cast<size_t>(cs) * dsz);       // pad a short final chunk
        std::memcpy(a.bytes.data() + static_cast<size_t>(s - lo) * dsz,
                    raw.data() + static_cast<size_t>(s - cstart) * dsz,
                    static_cast<size_t>(e - s) * dsz);
    }
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

// Write a zarr v2 array, optionally chunked (along the first axis) and/or compressed. chunk_elems<=0
// and compressor=null reproduce the original single-chunk uncompressed output exactly (so existing
// stores are byte-identical). Because only the first (outermost, slowest) axis is chunked, each chunk
// is a contiguous byte range of the C-order buffer -- no scatter; edge chunks are full-size,
// fill-padded per the v2 spec.
inline void write_array(const fs::path& dir, const NdArray& a,
                        int64_t chunk_elems = 0, const json& compressor = json(nullptr)) {
    fs::create_directories(dir);
    std::vector<int64_t> chunks = chunk_shape_for(a.shape, chunk_elems);
    json za;
    za["zarr_format"] = 2;
    za["shape"] = a.shape;
    za["chunks"] = chunks;
    za["dtype"] = a.dtype;
    za["compressor"] = compressor;
    za["fill_value"] = 0;
    za["order"] = "C";
    za["filters"] = nullptr;
    write_text(dir / ".zarray", za.dump());
    write_text(dir / ".zattrs", json::object().dump());

    const size_t ndim = a.shape.size();
    const size_t dsz = dtype_size(a.dtype);
    if (ndim == 0) {                                          // scalar -> single chunk "0"
        write_bytes(dir / "0", a.bytes.data(), a.bytes.size());
        return;
    }
    int64_t inner = 1;
    for (size_t i = 1; i < ndim; ++i) inner *= a.shape[i];
    const int64_t crows = chunks[0];
    const int64_t nchunks0 = (a.shape[0] + crows - 1) / crows;
    const size_t chunk_bytes = static_cast<size_t>(crows) * static_cast<size_t>(inner) * dsz;
    for (int64_t ci = 0; ci < nchunks0; ++ci) {
        const int64_t r0 = ci * crows;
        const int64_t rows = std::min<int64_t>(crows, a.shape[0] - r0);
        std::vector<uint8_t> chunk(chunk_bytes, 0);          // full-size, fill_value 0 padded
        std::memcpy(chunk.data(), a.bytes.data() + static_cast<size_t>(r0) * inner * dsz,
                    static_cast<size_t>(rows) * inner * dsz);
        std::vector<uint8_t> enc = encode_chunk(compressor, std::move(chunk));
        std::string key = std::to_string(ci);
        for (size_t i = 1; i < ndim; ++i) key += ".0";       // only the first axis is chunked
        write_bytes(dir / key, enc.data(), enc.size());
    }
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
                          int64_t chunk_elems = 0, const json& compressor = json(nullptr)) {
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
    write_array(gdir / name, data, chunk_elems, compressor);
    write_array(gdir / (name + "_offsets"), offsets, chunk_elems, compressor);
}

// ------------------------------------------------------------ group IO ------

inline void write_group(const fs::path& dir, const json& attrs) {
    fs::create_directories(dir);
    write_text(dir / ".zgroup", json{{"zarr_format", 2}}.dump());
    write_text(dir / ".zattrs", attrs.dump());
}

// Walk a finished store and emit a consolidated .zmetadata (so zarr.open_consolidated works and
// the metadata is one read instead of many small stats).
inline void consolidate_metadata(const fs::path& root) {
    json meta = json::object();
    for (auto& e : fs::recursive_directory_iterator(root)) {
        if (!e.is_regular_file()) continue;
        std::string fn = e.path().filename().string();
        if (fn == ".zgroup" || fn == ".zattrs" || fn == ".zarray")
            meta[fs::relative(e.path(), root).generic_string()] = read_json(e.path());
    }
    json out;
    out["zarr_consolidated_format"] = 1;
    out["metadata"] = meta;
    write_text(root / ".zmetadata", out.dump());
}

inline std::string opt_str(const json& m, const char* k) {
    return (m.contains(k) && !m[k].is_null()) ? m[k].get<std::string>() : std::string();
}

// ------------------------------------------------------------ dataset IO ----

inline Dataset read(const fs::path& root) {
    json rmeta = read_json(root / ".zattrs")["lstar"];
    Dataset ds;
    ds.kind = rmeta.value("kind", std::string("sample"));
    ds.spec_version = rmeta.value("spec_version", std::string("0.1"));
    ds.profiles = rmeta.value("profiles", std::vector<std::string>{});
    ds.dropped = rmeta.value("dropped", std::vector<std::string>{});

    for (auto& an : rmeta["axes"]) {
        std::string name = an.get<std::string>();
        fs::path g = root / "axes" / name;
        json m = read_json(g / ".zattrs")["lstar"];
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
        json m = read_json(g / ".zattrs")["lstar"];
        Field f;
        f.name = name;
        f.role = opt_str(m, "role");
        f.encoding = opt_str(m, "encoding");
        f.state = opt_str(m, "state");
        f.subtype = opt_str(m, "subtype");
        if (m.contains("span") && !m["span"].is_null())
            f.span = m["span"].get<std::vector<std::string>>();
        if (m.contains("provenance") && !m["provenance"].is_null()) f.provenance = m["provenance"];
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
        ds.fields.push_back(std::move(f));
    }

    if (rmeta.contains("aux") && !rmeta["aux"].is_null()) {            // verbatim passthrough subtree
        for (auto& an : rmeta["aux"]) {
            std::string ns = an.get<std::string>();
            fs::path g = root / "aux" / ns;
            Aux ax;
            ax.ns = ns;
            ax.attrs = read_json(g / ".zattrs")["lstar"];
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
                  int64_t chunk_elems = 0, const json& compressor = json(nullptr)) {
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
    rl["aux"] = auxnames;
    write_group(root, json{{"lstar", rl}});
    write_group(root / "axes", json::object());
    write_group(root / "fields", json::object());
    write_group(root / "models", json::object());

    for (auto& a : ds.axes) {
        fs::path g = root / "axes" / a.name;
        json al;
        al["kind"] = "axis";
        al["origin"] = a.origin;
        al["role"] = a.role.empty() ? json(nullptr) : json(a.role);
        al["induced_by"] = a.induced_by.empty() ? json(nullptr) : json(a.induced_by);
        al["provenance"] = a.provenance;
        write_group(g, json{{"lstar", al}});
        write_strings(g, "labels", a.labels, chunk_elems, compressor);
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
        fl["coverage"] = "full";
        fl["uncertainty"] = nullptr;
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
        write_group(g, json{{"lstar", fl}});
        if (f.encoding == "csr" || f.encoding == "csc") {
            write_array(g / "data", f.data, chunk_elems, compressor);
            write_array(g / "indices", f.indices, chunk_elems, compressor);
            write_array(g / "indptr", f.indptr, chunk_elems, compressor);
        } else if (f.encoding == "utf8") {
            write_strings(g, "values", f.strings, chunk_elems, compressor);
        } else if (f.encoding == "categorical") {
            write_array(g / "codes", f.codes, chunk_elems, compressor);
            write_strings(g, "categories", f.categories, chunk_elems, compressor);  // inline (P1)
        } else {
            write_array(g / "values", f.dense, chunk_elems, compressor);
        }
        if (f.has_mask) write_array(g / "mask", f.mask, chunk_elems, compressor);
    }

    if (!ds.aux.empty()) {                                            // verbatim passthrough subtree
        write_group(root / "aux", json::object());
        for (auto& ax : ds.aux) {
            fs::path g = root / "aux" / ax.ns;
            write_group(g, json{{"lstar", ax.attrs}});
            for (auto& leaf : ax.leaves) {
                if (leaf.kind == "utf8") write_strings(g, leaf.id, leaf.strings, chunk_elems, compressor);
                else write_array(g / leaf.id, leaf.dense, chunk_elems, compressor);
            }
        }
    }
    consolidate_metadata(root);
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
template <class T>
struct CsxArrays {
    std::vector<T> data;                  // preserves the input value dtype (no widening)
    std::vector<int64_t> indices;
    std::vector<int64_t> indptr;
    int64_t nrows = 0, ncols = 0;
};

// indptr (length ncols+1) and indices are int64 -- normalize with as_i64() before calling. `data`
// is templated on the stored value dtype so a float32 measure transposes as float32 (memory-lean).
template <class T>
inline CsxArrays<T> csc_to_csr(const T* data, const int64_t* indices, const int64_t* indptr,
                               int64_t nrows, int64_t ncols) {
    const int64_t nnz = indptr[ncols];
    CsxArrays<T> out;
    out.nrows = nrows;
    out.ncols = ncols;
    out.data.resize(static_cast<size_t>(nnz));
    out.indices.resize(static_cast<size_t>(nnz));
    out.indptr.assign(static_cast<size_t>(nrows) + 1, 0);
    for (int64_t k = 0; k < nnz; ++k) out.indptr[indices[k] + 1]++;     // count per row
    for (int64_t i = 0; i < nrows; ++i) out.indptr[i + 1] += out.indptr[i];
    std::vector<int64_t> next(out.indptr.begin(), out.indptr.end() - 1);
    for (int64_t j = 0; j < ncols; ++j) {                               // scatter by row
        for (int64_t k = indptr[j]; k < indptr[j + 1]; ++k) {
            int64_t dst = next[indices[k]]++;
            out.data[static_cast<size_t>(dst)] = data[k];
            out.indices[static_cast<size_t>(dst)] = j;
        }
    }
    return out;
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
    json m = read_json(field_group / ".zattrs")["lstar"];
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
    fs::path dir;
    int64_t cs = 1, n = 0;
    size_t dsz = 1;
    json compressor;
    int64_t cached = -1;
    std::vector<uint8_t> buf;                              // the decoded chunk (cs*dsz bytes)

    explicit ChunkReader(const fs::path& d) : dir(d) {
        json za = read_json(d / ".zarray");
        cs = za["chunks"].get<std::vector<int64_t>>()[0];
        n = za["shape"].get<std::vector<int64_t>>()[0];
        dsz = dtype_size(za["dtype"].get<std::string>());
        compressor = za.contains("compressor") ? za["compressor"] : json(nullptr);
    }
    const uint8_t* chunk(int64_t ci) {
        if (ci != cached) {
            fs::path cf = dir / std::to_string(ci);
            if (fs::exists(cf)) {
                std::vector<uint8_t> raw =
                    decode_chunk(compressor, read_bytes(cf), static_cast<size_t>(cs) * dsz);
                raw.resize(static_cast<size_t>(cs) * dsz);
                buf.swap(raw);
            } else {
                buf.assign(static_cast<size_t>(cs) * dsz, 0);   // missing chunk -> fill 0
            }
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
    json m = read_json(field_group / ".zattrs")["lstar"];
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
    b.data.dtype = read_json(field_group / "data" / ".zarray")["dtype"].get<std::string>();
    b.indices.dtype = read_json(field_group / "indices" / ".zarray")["dtype"].get<std::string>();
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
    json m = read_json(field_group / ".zattrs")["lstar"];
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
    json m = read_json(field_group / ".zattrs")["lstar"];
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
