// cpp11 bridge between R and the libstar C++ core.
// Returns/accepts L* datasets as R lists; the R layer assembles Matrix/Seurat objects.
#include <cpp11.hpp>

#include <cmath>
#include <cstdint>
#include <string>
#include <vector>

#include "lstar/lstar.hpp"

using namespace cpp11;

namespace {

writable::doubles nd_doubles(const lstar::NdArray& a) {
  int64_t n = a.nelem();
  writable::doubles out(n);
  const std::string& dt = a.dtype;
  if (dt == "<f8") { auto p = a.as<double>();   for (R_xlen_t i = 0; i < (R_xlen_t)n; ++i) out[i] = p[i]; }
  else if (dt == "<f4") { auto p = a.as<float>();    for (R_xlen_t i = 0; i < (R_xlen_t)n; ++i) out[i] = p[i]; }
  else if (dt == "<i4") { auto p = a.as<int32_t>();  for (R_xlen_t i = 0; i < (R_xlen_t)n; ++i) out[i] = p[i]; }
  else if (dt == "<i8") { auto p = a.as<int64_t>();  for (R_xlen_t i = 0; i < (R_xlen_t)n; ++i) out[i] = (double)p[i]; }
  else if (dt == "|u1" || dt == "<u1") { auto p = a.as<uint8_t>(); for (R_xlen_t i = 0; i < (R_xlen_t)n; ++i) out[i] = p[i]; }
  else if (dt == "<u4") { auto p = a.as<uint32_t>(); for (R_xlen_t i = 0; i < (R_xlen_t)n; ++i) out[i] = p[i]; }
  else if (dt == "|b1" || dt == "<b1") { auto p = a.as<uint8_t>(); for (R_xlen_t i = 0; i < (R_xlen_t)n; ++i) out[i] = p[i] ? 1.0 : 0.0; }
  else throw std::runtime_error("nd_doubles: unsupported dtype " + dt);
  return out;
}

writable::integers nd_integers(const lstar::NdArray& a) {
  int64_t n = a.nelem();
  writable::integers out(n);
  const std::string& dt = a.dtype;
  if (dt == "<i4") { auto p = a.as<int32_t>(); for (R_xlen_t i = 0; i < (R_xlen_t)n; ++i) out[i] = p[i]; }
  else if (dt == "<i8") { auto p = a.as<int64_t>(); for (R_xlen_t i = 0; i < (R_xlen_t)n; ++i) out[i] = (int)p[i]; }
  else if (dt == "<u4") { auto p = a.as<uint32_t>(); for (R_xlen_t i = 0; i < (R_xlen_t)n; ++i) out[i] = (int)p[i]; }
  else if (dt == "|u1" || dt == "<u1") { auto p = a.as<uint8_t>(); for (R_xlen_t i = 0; i < (R_xlen_t)n; ++i) out[i] = (int)p[i]; }
  else throw std::runtime_error("nd_integers: unsupported dtype " + dt);
  return out;
}

writable::strings to_strings(const std::vector<std::string>& v) {
  writable::strings out(v.size());
  for (size_t i = 0; i < v.size(); ++i) out[i] = v[i];
  return out;
}

// Aux (passthrough) dense leaves round-trip as raw bytes + dtype + shape, so *any* dtype is preserved
// exactly across the R boundary (no f8 widening), the same way the store holds them.
writable::raws nd_to_raws(const lstar::NdArray& a) {
  writable::raws out(a.bytes.size());
  for (size_t i = 0; i < a.bytes.size(); ++i) out[i] = (Rbyte)a.bytes[i];
  return out;
}

writable::integers to_ints(const std::vector<int64_t>& v) {
  writable::integers out(v.size());
  for (size_t i = 0; i < v.size(); ++i) out[i] = (int)v[i];
  return out;
}

// Decode an R numeric vector (C-order) into a libstar NdArray of a given zarr dtype.
lstar::NdArray nd_from_doubles(const doubles& v, std::vector<int64_t> shape, const std::string& dtype) {
  lstar::NdArray a;
  a.dtype = dtype;
  a.shape = shape;
  size_t n = (size_t)v.size();
  size_t dsz = lstar::dtype_size(dtype);
  a.bytes.resize(n * dsz);
  if (dtype == "<f8") { auto p = a.as<double>();  for (size_t i = 0; i < n; ++i) p[i] = v[i]; }
  else if (dtype == "<f4") { auto p = a.as<float>();   for (size_t i = 0; i < n; ++i) p[i] = (float)v[i]; }
  else if (dtype == "<i4") { auto p = a.as<int32_t>(); for (size_t i = 0; i < n; ++i) p[i] = (int32_t)v[i]; }
  else if (dtype == "<i8") { auto p = a.as<int64_t>(); for (size_t i = 0; i < n; ++i) p[i] = (int64_t)v[i]; }
  else throw std::runtime_error("nd_from_doubles: unsupported dtype " + dtype);
  return a;
}

// Pick a compact on-disk dtype for a raw-counts data vector. Used ONLY for fields declared
// state=="raw": integer-valued counts become i4 (or i8 if out of int32 range) instead of 8-byte
// floats. We gate on the declared semantics, not on the values, so a genuinely floating layer
// (lognorm/scaled, or a float measure that merely happens to be integer-valued) is never narrowed
// — that keeps dtype predictable across round-trips through outside formats (AnnData/Seurat/SCE).
// A raw layer that is itself non-integer (e.g. corrected counts) also stays f8.
static std::string pick_data_dtype(const doubles& v) {
  bool all_int = true;
  double mn = 0.0, mx = 0.0;
  const R_xlen_t n = v.size();
  for (R_xlen_t i = 0; i < n; ++i) {
    const double x = v[i];
    if (!std::isfinite(x) || x != std::trunc(x)) { all_int = false; break; }
    if (x < mn) mn = x;
    if (x > mx) mx = x;
  }
  if (!all_int) return "<f8";
  if (mn >= -2147483648.0 && mx <= 2147483647.0) return "<i4";
  return "<i8";
}

static writable::doubles vec_to_dbl(const std::vector<double>& v) {
  writable::doubles o((R_xlen_t)v.size());
  for (R_xlen_t i = 0; i < (R_xlen_t)v.size(); ++i) o[i] = v[(size_t)i];
  return o;
}

}  // namespace

// Per-group sufficient stats over a CSC measure (cells x genes). `group`: per-cell group in
// [0,ngroups) or <0 to skip. Returns sum/sumsq/n_expr flat (ngroups x ngenes). The same libstar
// kernel the Python and WASM bindings call ("one kernel, every runtime").
[[cpp11::register]]
list lstar_cpp_col_sum_by_group(doubles data, integers indptr, integers indices,
                                int nrows, int ncols, integers group, int ngroups, bool lognorm) {
  std::vector<double> dv((size_t)data.size()); for (R_xlen_t i = 0; i < data.size(); ++i) dv[(size_t)i] = data[i];
  std::vector<int64_t> ip((size_t)indptr.size()); for (R_xlen_t i = 0; i < indptr.size(); ++i) ip[(size_t)i] = indptr[i];
  std::vector<int64_t> ix((size_t)indices.size()); for (R_xlen_t i = 0; i < indices.size(); ++i) ix[(size_t)i] = indices[i];
  std::vector<int> grp((size_t)group.size()); for (R_xlen_t i = 0; i < group.size(); ++i) grp[(size_t)i] = group[i];
  auto s = lstar::csc_col_sum_by_group(dv.data(), ip.data(), ix.data(), nrows, ncols, grp.data(), ngroups, lognorm, 0);
  return writable::list({"sum"_nm = vec_to_dbl(s.sum), "sumsq"_nm = vec_to_dbl(s.sumsq),
                         "n_expr"_nm = vec_to_dbl(s.n_expr), "ngroups"_nm = ngroups, "ngenes"_nm = (int)s.ngenes});
}

// Subsample DE ranker over a CSR submatrix (sampled cells x genes). membership: 0=A, 1=B, <0=skip.
[[cpp11::register]]
list lstar_cpp_subsample_de_rank(doubles data, integers indptr, integers indices,
                                 int nrows, int ngenes, integers membership, bool lognorm) {
  std::vector<double> dv((size_t)data.size()); for (R_xlen_t i = 0; i < data.size(); ++i) dv[(size_t)i] = data[i];
  std::vector<int64_t> ip((size_t)indptr.size()); for (R_xlen_t i = 0; i < indptr.size(); ++i) ip[(size_t)i] = indptr[i];
  std::vector<int64_t> ix((size_t)indices.size()); for (R_xlen_t i = 0; i < indices.size(); ++i) ix[(size_t)i] = indices[i];
  std::vector<int> mem((size_t)membership.size()); for (R_xlen_t i = 0; i < membership.size(); ++i) mem[(size_t)i] = membership[i];
  auto r = lstar::subsample_de_rank(dv.data(), ip.data(), ix.data(), nrows, ngenes, mem.data(), lognorm);
  return writable::list({"meanA"_nm = vec_to_dbl(r.meanA), "meanB"_nm = vec_to_dbl(r.meanB),
                         "lfc"_nm = vec_to_dbl(r.lfc), "nA"_nm = (int)r.nA, "nB"_nm = (int)r.nB});
}

// Zero-aware per-gene mean/var/nnz of a CSC measure in a store, read block-by-block so the whole
// matrix never lands in memory (the bounded-memory reduction; same libstar kernel as Python/WASM).
[[cpp11::register]]
list lstar_cpp_stream_col_stats(std::string path, std::string field, int block,
                                int n_threads, bool lognorm, doubles depth,
                                double depthScale, bool population) {
  std::vector<double> dv;                          // empty -> no depth normalization
  if (depth.size() > 0) {
    dv.resize((size_t)depth.size());
    for (R_xlen_t i = 0; i < depth.size(); ++i) dv[(size_t)i] = depth[i];
  }
  const std::vector<double>* dp = dv.empty() ? nullptr : &dv;
  lstar::ColStats s = lstar::stream_csc_col_mean_var(
      path + "/fields/" + field, (int64_t)block, n_threads, lognorm, dp, depthScale, population);
  return writable::list({"mean"_nm = vec_to_dbl(s.mean), "var"_nm = vec_to_dbl(s.var),
                         "nnz"_nm = to_ints(s.nnz)});
}

// Streaming per-(group, gene) SUM of a CSC measure (the fused pseudobulk / colSumByFacView), with the
// optional depth-normalized log1p view. Returns the (ngroups x ncols) sums flat, row-major (g*ncols+j).
[[cpp11::register]]
doubles lstar_cpp_stream_col_sum_by_group(std::string path, std::string field, integers group,
                                          int ngroups, bool lognorm, doubles depth, double depthScale,
                                          int block, int n_threads) {
  std::vector<int> g((size_t)group.size());
  for (R_xlen_t i = 0; i < group.size(); ++i) g[(size_t)i] = group[i];
  std::vector<double> dv;
  if (depth.size() > 0) { dv.resize((size_t)depth.size()); for (R_xlen_t i = 0; i < depth.size(); ++i) dv[(size_t)i] = depth[i]; }
  const std::vector<double>* dp = dv.empty() ? nullptr : &dv;
  std::vector<double> out = lstar::stream_csc_col_sum_by_group(
      path + "/fields/" + field, g, ngroups, lognorm, dp, depthScale, (int64_t)block, n_threads);
  writable::doubles o((R_xlen_t)out.size());
  for (size_t i = 0; i < out.size(); ++i) o[(R_xlen_t)i] = out[i];
  return o;
}

// Read a contiguous gene (column) range [g_lo, g_hi) of a CSC measure as CSC arrays, touching only
// the overlapping chunks. The general bounded block-read primitive (R assembles the dgCMatrix).
[[cpp11::register]]
list lstar_cpp_read_csc_block(std::string path, std::string field, int g_lo, int g_hi) {
  lstar::CscBlock b = lstar::read_csc_block(path + "/fields/" + field, (int64_t)g_lo, (int64_t)g_hi);
  return writable::list({"data"_nm = nd_doubles(b.data), "indices"_nm = nd_integers(b.indices),
                         "indptr"_nm = to_ints(b.indptr), "nrows"_nm = (int)b.nrows,
                         "ncols"_nm = (int)b.ncols});
}

// Gather an arbitrary (sorted, unique, 0-based) set of gene columns of a CSC measure, decoding each
// touched chunk at most once -- the efficient scattered-subset read.
[[cpp11::register]]
list lstar_cpp_read_csc_cols(std::string path, std::string field, integers cols) {
  std::vector<int64_t> cv((size_t)cols.size());
  for (R_xlen_t i = 0; i < cols.size(); ++i) cv[(size_t)i] = (int64_t)cols[i];
  lstar::CscBlock b = lstar::read_csc_cols(path + "/fields/" + field, cv);
  return writable::list({"data"_nm = nd_doubles(b.data), "indices"_nm = nd_integers(b.indices),
                         "indptr"_nm = to_ints(b.indptr), "nrows"_nm = (int)b.nrows,
                         "ncols"_nm = (int)b.ncols});
}

[[cpp11::register]]
list lstar_cpp_read(std::string path) {
  lstar::Dataset ds = lstar::read(path);

  writable::list axes(ds.axes.size());
  writable::strings axnames(ds.axes.size());
  for (size_t k = 0; k < ds.axes.size(); ++k) {
    const auto& a = ds.axes[k];
    axes[k] = writable::list({"labels"_nm = to_strings(a.labels),
                              "origin"_nm = a.origin, "role"_nm = a.role,
                              "induced_by"_nm = a.induced_by});
    axnames[k] = a.name;
  }
  axes.attr("names") = axnames;

  writable::list fields(ds.fields.size());
  writable::strings fnames(ds.fields.size());
  for (size_t k = 0; k < ds.fields.size(); ++k) {
    const auto& f = ds.fields[k];
    writable::list fl;
    if (f.encoding == "csc" || f.encoding == "csr") {
      fl = writable::list({"role"_nm = f.role, "span"_nm = to_strings(f.span),
                           "encoding"_nm = f.encoding, "state"_nm = f.state, "subtype"_nm = f.subtype,
                           "data"_nm = nd_doubles(f.data), "indices"_nm = nd_integers(f.indices),
                           "indptr"_nm = nd_integers(f.indptr), "shape"_nm = to_ints(f.shape)});
    } else if (f.encoding == "utf8") {
      fl = writable::list({"role"_nm = f.role, "span"_nm = to_strings(f.span),
                           "encoding"_nm = f.encoding, "state"_nm = f.state, "subtype"_nm = f.subtype,
                           "strings"_nm = to_strings(f.strings)});
    } else if (f.encoding == "categorical") {
      fl = writable::list({"role"_nm = f.role, "span"_nm = to_strings(f.span),
                           "encoding"_nm = f.encoding, "state"_nm = f.state, "subtype"_nm = f.subtype,
                           "codes"_nm = nd_integers(f.codes), "categories"_nm = to_strings(f.categories),
                           "ordered"_nm = f.ordered});
    } else {
      fl = writable::list({"role"_nm = f.role, "span"_nm = to_strings(f.span),
                           "encoding"_nm = f.encoding, "state"_nm = f.state, "subtype"_nm = f.subtype,
                           "dense"_nm = nd_doubles(f.dense), "shape"_nm = to_ints(f.dense.shape)});
    }
    if (f.has_mask) fl.push_back("mask"_nm = nd_integers(f.mask));   // 1 == missing (nullable)
    if (f.has_index) {                                              // partial coverage over index_axis
      fl.push_back("index"_nm = nd_integers(f.index));
      fl.push_back("index_axis"_nm = f.index_axis);
      fl.push_back("coverage"_nm = std::string("partial"));
    }
    if (f.provenance.is_object() && !f.provenance.empty())          // provenance as an opaque JSON string
      fl.push_back("provenance"_nm = f.provenance.dump());          // (preserves arbitrary nesting verbatim)
    fields[k] = fl;
    fnames[k] = f.name;
  }
  fields.attr("names") = fnames;

  // Aux passthrough: attrs (tree string + manifest) carried as a JSON string; each leaf as raw bytes
  // (dense) or a character vector (utf8). The R layer round-trips it verbatim, never interpreting it.
  writable::list aux(ds.aux.size());
  writable::strings auxnames(ds.aux.size());
  for (size_t k = 0; k < ds.aux.size(); ++k) {
    const auto& a = ds.aux[k];
    writable::list leaves(a.leaves.size());
    writable::strings lnames(a.leaves.size());
    for (size_t j = 0; j < a.leaves.size(); ++j) {
      const auto& lf = a.leaves[j];
      if (lf.kind == "utf8")
        leaves[j] = writable::list({"kind"_nm = lf.kind, "strings"_nm = to_strings(lf.strings)});
      else
        leaves[j] = writable::list({"kind"_nm = lf.kind, "dtype"_nm = lf.dense.dtype,
                                    "shape"_nm = to_ints(lf.dense.shape), "bytes"_nm = nd_to_raws(lf.dense)});
      lnames[j] = lf.id;
    }
    leaves.attr("names") = lnames;
    aux[k] = writable::list({"attrs"_nm = a.attrs.dump(), "leaves"_nm = leaves});
    auxnames[k] = a.ns;
  }
  aux.attr("names") = auxnames;

  return writable::list({"kind"_nm = ds.kind, "spec_version"_nm = ds.spec_version,
                         "profiles"_nm = to_strings(ds.profiles),
                         "dropped"_nm = to_strings(ds.dropped),
                         "axes"_nm = axes, "fields"_nm = fields, "aux"_nm = aux});
}

// Write: R list -> libstar Dataset -> Zarr store. The R layer disassembles Matrix objects
// into (data, indices, indptr) / dense vectors before calling this.
[[cpp11::register]]
void lstar_cpp_write(list ds, std::string path, int chunk_elems = 0,
                     std::string compression = "", int level = 5) {
  lstar::json compressor = nullptr;                 // "" -> uncompressed; else numcodecs gzip/zlib codec
  if (compression == "gzip" || compression == "zlib")
    compressor = lstar::json{{"id", compression}, {"level", level}};
  else if (!compression.empty())
    throw std::runtime_error("unsupported compression: " + compression + " (use 'gzip', 'zlib', or '')");
  lstar::Dataset out;
  out.kind = as_cpp<std::string>(ds["kind"]);
  out.spec_version = as_cpp<std::string>(ds["spec_version"]);
  out.profiles = as_cpp<std::vector<std::string>>(ds["profiles"]);
  out.dropped = as_cpp<std::vector<std::string>>(ds["dropped"]);

  list axes = ds["axes"];
  strings axnames = axes.names();
  for (R_xlen_t i = 0; i < axes.size(); ++i) {
    list a = axes[i];
    lstar::Axis ax;
    ax.name = axnames[i];
    ax.labels = as_cpp<std::vector<std::string>>(a["labels"]);
    ax.origin = as_cpp<std::string>(a["origin"]);
    ax.role = as_cpp<std::string>(a["role"]);
    {                                                   // induced_by is optional (older payloads omit it)
      strings ns = a.names();
      for (R_xlen_t j = 0; j < ns.size(); ++j)
        if (std::string(ns[j]) == "induced_by") { ax.induced_by = as_cpp<std::string>(a["induced_by"]); break; }
    }
    out.axes.push_back(std::move(ax));
  }

  list fields = ds["fields"];
  strings fnames = fields.names();
  for (R_xlen_t i = 0; i < fields.size(); ++i) {
    list f = fields[i];
    lstar::Field fl;
    fl.name = fnames[i];
    fl.role = as_cpp<std::string>(f["role"]);
    fl.span = as_cpp<std::vector<std::string>>(f["span"]);
    fl.encoding = as_cpp<std::string>(f["encoding"]);
    fl.state = as_cpp<std::string>(f["state"]);
    fl.subtype = as_cpp<std::string>(f["subtype"]);
    {                                                     // provenance: an opaque JSON string -> json
      strings ns = f.names();
      for (R_xlen_t j = 0; j < ns.size(); ++j) if (std::string(ns[j]) == "provenance") {
        std::string pj = as_cpp<std::string>(f["provenance"]);
        if (!pj.empty()) {
          auto p = decltype(fl.provenance)::parse(pj, nullptr, false);  // lenient: discarded on bad input
          if (p.is_object()) fl.provenance = p;
        }
        break;
      }
    }
    if (fl.encoding == "csc" || fl.encoding == "csr") {
      integers shp = f["shape"];
      for (R_xlen_t j = 0; j < shp.size(); ++j) fl.shape.push_back((int64_t)shp[j]);
      doubles dat = f["data"], ind = f["indices"], ptr = f["indptr"];
      const std::string ddt = (fl.state == "raw") ? pick_data_dtype(dat) : std::string("<f8");
      fl.data    = nd_from_doubles(dat, {(int64_t)dat.size()}, ddt);
      fl.indices = nd_from_doubles(ind, {(int64_t)ind.size()}, "<i4");
      fl.indptr  = nd_from_doubles(ptr, {(int64_t)ptr.size()}, "<i4");
    } else if (fl.encoding == "utf8") {
      fl.strings = as_cpp<std::vector<std::string>>(f["strings"]);
    } else if (fl.encoding == "categorical") {
      integers codes = f["codes"];                        // 0-based int, -1 = missing
      fl.codes.dtype = "<i4";
      fl.codes.shape = {(int64_t)codes.size()};
      fl.codes.bytes.resize((size_t)codes.size() * 4);
      { auto p = fl.codes.as<int32_t>(); for (R_xlen_t j = 0; j < codes.size(); ++j) p[j] = codes[j]; }
      fl.categories = as_cpp<std::vector<std::string>>(f["categories"]);
      fl.ordered = as_cpp<bool>(f["ordered"]);
      fl.has_ordered = true;
    } else {
      integers shp = f["shape"];
      std::vector<int64_t> sh;
      for (R_xlen_t j = 0; j < shp.size(); ++j) sh.push_back((int64_t)shp[j]);
      doubles dn = f["dense"];
      fl.dense = nd_from_doubles(dn, sh, "<f8");          // dense stays predictable f8 (round-trip-safe)
    }
    {                                                     // optional validity mask (1 == missing)
      strings ns = f.names();
      for (R_xlen_t j = 0; j < ns.size(); ++j) if (std::string(ns[j]) == "mask") {
        integers mk = f["mask"];
        fl.mask.dtype = "|u1"; fl.mask.shape = {(int64_t)mk.size()};
        fl.mask.bytes.resize((size_t)mk.size());
        auto p = fl.mask.as<uint8_t>();
        for (R_xlen_t i = 0; i < mk.size(); ++i) p[i] = (mk[i] != NA_INTEGER && mk[i] != 0) ? 1 : 0;
        fl.has_mask = true;
        break;
      }
    }
    {                                                     // partial coverage: int index into index_axis
      strings ns = f.names();
      bool has_ix = false, has_iax = false;
      for (R_xlen_t j = 0; j < ns.size(); ++j) {
        std::string nm(ns[j]);
        if (nm == "index") has_ix = true;
        else if (nm == "index_axis") has_iax = true;
      }
      if (has_ix) {
        integers ix = f["index"];
        fl.index.dtype = "<i8"; fl.index.shape = {(int64_t)ix.size()};
        fl.index.bytes.resize((size_t)ix.size() * 8);
        auto p = fl.index.as<int64_t>();
        for (R_xlen_t i = 0; i < ix.size(); ++i) p[i] = (int64_t)ix[i];
        fl.has_index = true;
        fl.coverage = "partial";
        if (has_iax) fl.index_axis = as_cpp<std::string>(f["index_axis"]);
      }
    }
    out.fields.push_back(std::move(fl));
  }

  {                                                       // aux passthrough (optional; older lists omit it)
    strings dn = ds.names();
    bool has_aux = false;
    for (R_xlen_t j = 0; j < dn.size(); ++j) if (std::string(dn[j]) == "aux") { has_aux = true; break; }
    if (has_aux) {
      list aux = ds["aux"];
      strings auxnames = aux.names();
      for (R_xlen_t i = 0; i < aux.size(); ++i) {
        list a = aux[i];
        lstar::Aux ax;
        ax.ns = auxnames[i];
        ax.attrs = lstar::json::parse(as_cpp<std::string>(a["attrs"]));
        list leaves = a["leaves"];
        strings lnames = leaves.names();
        for (R_xlen_t j = 0; j < leaves.size(); ++j) {
          list lf = leaves[j];
          lstar::AuxLeaf leaf;
          leaf.id = lnames[j];
          leaf.kind = as_cpp<std::string>(lf["kind"]);
          if (leaf.kind == "utf8") {
            leaf.strings = as_cpp<std::vector<std::string>>(lf["strings"]);
          } else {
            leaf.dense.dtype = as_cpp<std::string>(lf["dtype"]);
            integers sh = lf["shape"];
            for (R_xlen_t s = 0; s < sh.size(); ++s) leaf.dense.shape.push_back((int64_t)sh[s]);
            raws bv = lf["bytes"];
            leaf.dense.bytes.resize((size_t)bv.size());
            for (R_xlen_t b = 0; b < bv.size(); ++b) leaf.dense.bytes[b] = (uint8_t)bv[b];
          }
          ax.leaves.push_back(std::move(leaf));
        }
        out.aux.push_back(std::move(ax));
      }
    }
  }
  lstar::write(out, path, (int64_t)chunk_elems, compressor);
}
