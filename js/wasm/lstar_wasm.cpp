// lstar._wasm — Emscripten/embind binding over libstar's translation primitives.
//
// The browser/Node binding of the same header-only C++ core that backs the R (cpp11) and Python
// (pybind11) packages: the compute kernels run in WebAssembly so a TypeScript viewer gets the exact
// same numbers without a server. Only the compute primitives are exposed here (no IO / filesystem /
// zlib) — Zarr I/O is the TypeScript layer's job (zarrita.js); this module receives assembled typed
// arrays and returns results.
//
// Inputs are JS typed arrays (copied in); outputs are fresh JS typed arrays. Single-threaded for now
// (no OpenMP/pthreads — that needs cross-origin isolation; a later step).
#include <cstdint>
#include <string>
#include <vector>

#include <emscripten/bind.h>
#include <emscripten/val.h>

#include "lstar/lstar.hpp"

using namespace emscripten;

// ---- marshalling helpers ----
static std::vector<int64_t> to_i64(const val& a) {
    std::vector<double> d = convertJSArrayToNumberVector<double>(a);
    std::vector<int64_t> out(d.size());
    for (size_t i = 0; i < d.size(); ++i) out[i] = static_cast<int64_t>(d[i]);
    return out;
}
// copy a C++ vector into a JS-owned typed array (the .slice() copies before the C++ vector dies)
static val to_f64(const std::vector<double>& v) {
    return val(typed_memory_view(v.size(), v.data())).call<val>("slice");
}
static val to_i32(const std::vector<int64_t>& v) {
    std::vector<int32_t> t(v.begin(), v.end());
    return val(typed_memory_view(t.size(), t.data())).call<val>("slice");
}
static std::vector<int> to_int(const val& a) {
    std::vector<double> d = convertJSArrayToNumberVector<double>(a);
    std::vector<int> out(d.size());
    for (size_t i = 0; i < d.size(); ++i) out[i] = static_cast<int>(d[i]);
    return out;
}

// Zero-aware per-column mean/variance of a CSC measure (optionally over log1p values).
// data: Float32Array|Float64Array (nnz); indptr: Int32Array|Int32 (ncols+1). -> {mean,var,nnz}.
static val colMeanVar(val data_js, val indptr_js, int nrows, int n_threads, bool lognorm) {
    std::vector<double> data = convertJSArrayToNumberVector<double>(data_js);
    std::vector<int64_t> indptr = to_i64(indptr_js);
    int64_t ncols = static_cast<int64_t>(indptr.size()) - 1;
    auto s = lstar::csc_col_mean_var(data.data(), indptr.data(), ncols, nrows, n_threads, lognorm);
    std::vector<int64_t> nnz(s.nnz.begin(), s.nnz.end());
    val out = val::object();
    out.set("mean", to_f64(s.mean));
    out.set("var", to_f64(s.var));
    out.set("nnz", to_i32(nnz));
    return out;
}

// CSC -> CSR storage transpose (orientation flip, e.g. gene-major <-> cell-major). -> {data,indices,indptr}.
static val cscToCsr(val data_js, val indices_js, val indptr_js, int nrows, int ncols) {
    std::vector<double> data = convertJSArrayToNumberVector<double>(data_js);
    std::vector<int64_t> indices = to_i64(indices_js);
    std::vector<int64_t> indptr = to_i64(indptr_js);
    auto r = lstar::csc_to_csr(data.data(), indices.data(), indptr.data(), nrows, ncols);
    val out = val::object();
    out.set("data", to_f64(r.data));
    out.set("indices", to_i32(r.indices));
    out.set("indptr", to_i32(r.indptr));
    return out;
}

// Per-group sufficient stats over a CSC measure (cells x genes). group: Int32Array (length nrows),
// cell -> group in [0,ngroups) or <0 to skip. -> {sum,sumsq,n_expr} flat (ngroups x ncols), ngenes.
static val colSumByGroup(val data_js, val indptr_js, val indices_js, int nrows, int ncols, val group_js, int ngroups, bool lognorm) {
    std::vector<double> data = convertJSArrayToNumberVector<double>(data_js);
    std::vector<int64_t> indptr = to_i64(indptr_js), indices = to_i64(indices_js);
    std::vector<int> grp = to_int(group_js);
    auto s = lstar::csc_col_sum_by_group(data.data(), indptr.data(), indices.data(), nrows, ncols, grp.data(), ngroups, lognorm, 1);
    val out = val::object();
    out.set("sum", to_f64(s.sum)); out.set("sumsq", to_f64(s.sumsq)); out.set("n_expr", to_f64(s.n_expr));
    out.set("ngroups", ngroups); out.set("ngenes", (double)s.ngenes);
    return out;
}

// Subsample DE ranker over a CSR submatrix (sampled cells x genes). membership: Int32Array
// (length nrows), 0=A, 1=B, <0=skip. -> {meanA,meanB,lfc, nA,nB}. Caller ranks by |lfc|.
static val subsampleDeRank(val data_js, val indptr_js, val indices_js, int nrows, int ngenes, val membership_js, bool lognorm) {
    std::vector<double> data = convertJSArrayToNumberVector<double>(data_js);
    std::vector<int64_t> indptr = to_i64(indptr_js), indices = to_i64(indices_js);
    std::vector<int> mem = to_int(membership_js);
    auto r = lstar::subsample_de_rank(data.data(), indptr.data(), indices.data(), nrows, ngenes, mem.data(), lognorm);
    val out = val::object();
    out.set("meanA", to_f64(r.meanA)); out.set("meanB", to_f64(r.meanB)); out.set("lfc", to_f64(r.lfc));
    out.set("nA", (double)r.nA); out.set("nB", (double)r.nB);
    return out;
}

static std::string version() { return "lstar-wasm 0.0.2"; }

EMSCRIPTEN_BINDINGS(lstar_wasm) {
    function("colMeanVar", &colMeanVar);
    function("cscToCsr", &cscToCsr);
    function("colSumByGroup", &colSumByGroup);
    function("subsampleDeRank", &subsampleDeRank);
    function("version", &version);
}
