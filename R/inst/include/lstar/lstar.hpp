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
};

struct Dataset {
    std::string kind = "sample";
    std::string spec_version = "0.1";
    std::vector<std::string> profiles, dropped;
    std::vector<Axis> axes;
    std::vector<Field> fields;

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
    std::ifstream f(p, std::ios::binary);
    if (!f) throw std::runtime_error("cannot open " + p.string());
    return std::vector<uint8_t>(std::istreambuf_iterator<char>(f), std::istreambuf_iterator<char>());
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

inline void write_array(const fs::path& dir, const NdArray& a) {
    fs::create_directories(dir);
    json za;
    za["zarr_format"] = 2;
    za["shape"] = a.shape;
    za["chunks"] = a.shape;
    za["dtype"] = a.dtype;
    za["compressor"] = nullptr;
    za["fill_value"] = 0;
    za["order"] = "C";
    za["filters"] = nullptr;
    write_text(dir / ".zarray", za.dump());
    write_text(dir / ".zattrs", json::object().dump());
    write_bytes(dir / zero_chunk_key(a.shape.size()), a.bytes.data(), a.bytes.size());
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
                          const std::vector<std::string>& strs) {
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
    write_array(gdir / name, data);
    write_array(gdir / (name + "_offsets"), offsets);
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
        } else {  // dense
            f.dense = read_array(g / "values");
        }
        ds.fields.push_back(std::move(f));
    }
    return ds;
}

inline void write(const Dataset& ds, const fs::path& root) {
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
        al["provenance"] = a.provenance;
        write_group(g, json{{"lstar", al}});
        write_strings(g, "labels", a.labels);
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
        if (f.encoding == "csr" || f.encoding == "csc") fl["shape"] = f.shape;
        if (f.encoding == "utf8") fl["shape"] = std::vector<int64_t>{(int64_t)f.strings.size()};
        write_group(g, json{{"lstar", fl}});
        if (f.encoding == "csr" || f.encoding == "csc") {
            write_array(g / "data", f.data);
            write_array(g / "indices", f.indices);
            write_array(g / "indptr", f.indptr);
        } else if (f.encoding == "utf8") {
            write_strings(g, "values", f.strings);
        } else {
            write_array(g / "values", f.dense);
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
template <class T>
inline ColStats csc_col_mean_var(const T* data, const int64_t* indptr,
                                 int64_t ncols, int64_t nrows, int n_threads = 0,
                                 bool lognorm = false) {
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
            sum += lognorm ? std::log1p(x) : x;
        }
        double m = sum / (double)nrows;
        double ss = 0.0;
        for (int64_t k = a; k < b; ++k) {
            double x = static_cast<double>(data[k]);
            double v = lognorm ? std::log1p(x) : x;
            double d = v - m;
            ss += d * d;
        }
        ss += m * m * (double)(nrows - (b - a));  // contribution of the implicit zeros
        s.mean[(size_t)j] = m;
        s.var[(size_t)j] = (nrows > 1) ? ss / (double)(nrows - 1) : 0.0;
        s.nnz[(size_t)j] = b - a;
    }
    return s;
}

}  // namespace lstar
