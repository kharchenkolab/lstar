// cpp11 bridge between R and the libstar C++ core.
// Returns/accepts L* datasets as R lists; the R layer assembles Matrix/Seurat objects.
#include <cpp11.hpp>

#include <string>
#include <vector>

#include "lstar/lstar.hpp"

using namespace cpp11;

namespace {

writable::doubles nd_doubles(const lstar::NdArray& a) {
  int64_t n = a.nelem();
  writable::doubles out(n);
  const std::string& dt = a.dtype;
  if (dt == "<f8") { auto p = a.as<double>();   for (int64_t i = 0; i < n; ++i) out[i] = p[i]; }
  else if (dt == "<f4") { auto p = a.as<float>();    for (int64_t i = 0; i < n; ++i) out[i] = p[i]; }
  else if (dt == "<i4") { auto p = a.as<int32_t>();  for (int64_t i = 0; i < n; ++i) out[i] = p[i]; }
  else if (dt == "<i8") { auto p = a.as<int64_t>();  for (int64_t i = 0; i < n; ++i) out[i] = (double)p[i]; }
  else if (dt == "|u1" || dt == "<u1") { auto p = a.as<uint8_t>(); for (int64_t i = 0; i < n; ++i) out[i] = p[i]; }
  else if (dt == "<u4") { auto p = a.as<uint32_t>(); for (int64_t i = 0; i < n; ++i) out[i] = p[i]; }
  else if (dt == "|b1" || dt == "<b1") { auto p = a.as<uint8_t>(); for (int64_t i = 0; i < n; ++i) out[i] = p[i] ? 1.0 : 0.0; }
  else throw std::runtime_error("nd_doubles: unsupported dtype " + dt);
  return out;
}

writable::integers nd_integers(const lstar::NdArray& a) {
  int64_t n = a.nelem();
  writable::integers out(n);
  const std::string& dt = a.dtype;
  if (dt == "<i4") { auto p = a.as<int32_t>(); for (int64_t i = 0; i < n; ++i) out[i] = p[i]; }
  else if (dt == "<i8") { auto p = a.as<int64_t>(); for (int64_t i = 0; i < n; ++i) out[i] = (int)p[i]; }
  else if (dt == "<u4") { auto p = a.as<uint32_t>(); for (int64_t i = 0; i < n; ++i) out[i] = (int)p[i]; }
  else throw std::runtime_error("nd_integers: unsupported dtype " + dt);
  return out;
}

writable::strings to_strings(const std::vector<std::string>& v) {
  writable::strings out(v.size());
  for (size_t i = 0; i < v.size(); ++i) out[i] = v[i];
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

}  // namespace

[[cpp11::register]]
list lstar_cpp_read(std::string path) {
  lstar::Dataset ds = lstar::read(path);

  writable::list axes(ds.axes.size());
  writable::strings axnames(ds.axes.size());
  for (size_t k = 0; k < ds.axes.size(); ++k) {
    const auto& a = ds.axes[k];
    axes[k] = writable::list({"labels"_nm = to_strings(a.labels),
                              "origin"_nm = a.origin, "role"_nm = a.role});
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
    } else {
      fl = writable::list({"role"_nm = f.role, "span"_nm = to_strings(f.span),
                           "encoding"_nm = f.encoding, "state"_nm = f.state, "subtype"_nm = f.subtype,
                           "dense"_nm = nd_doubles(f.dense), "shape"_nm = to_ints(f.dense.shape)});
    }
    fields[k] = fl;
    fnames[k] = f.name;
  }
  fields.attr("names") = fnames;

  return writable::list({"kind"_nm = ds.kind, "spec_version"_nm = ds.spec_version,
                         "profiles"_nm = to_strings(ds.profiles),
                         "dropped"_nm = to_strings(ds.dropped),
                         "axes"_nm = axes, "fields"_nm = fields});
}

// Write: R list -> libstar Dataset -> Zarr store. The R layer disassembles Matrix objects
// into (data, indices, indptr) / dense vectors before calling this.
[[cpp11::register]]
void lstar_cpp_write(list ds, std::string path) {
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
    if (fl.encoding == "csc" || fl.encoding == "csr") {
      integers shp = f["shape"];
      for (R_xlen_t j = 0; j < shp.size(); ++j) fl.shape.push_back((int64_t)shp[j]);
      doubles dat = f["data"], ind = f["indices"], ptr = f["indptr"];
      fl.data    = nd_from_doubles(dat, {(int64_t)dat.size()}, "<f8");
      fl.indices = nd_from_doubles(ind, {(int64_t)ind.size()}, "<i4");
      fl.indptr  = nd_from_doubles(ptr, {(int64_t)ptr.size()}, "<i4");
    } else if (fl.encoding == "utf8") {
      fl.strings = as_cpp<std::vector<std::string>>(f["strings"]);
    } else {
      integers shp = f["shape"];
      std::vector<int64_t> sh;
      for (R_xlen_t j = 0; j < shp.size(); ++j) sh.push_back((int64_t)shp[j]);
      doubles dn = f["dense"];
      fl.dense = nd_from_doubles(dn, sh, "<f8");
    }
    out.fields.push_back(std::move(fl));
  }
  lstar::write(out, path);
}
