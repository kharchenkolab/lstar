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

// gzip-compress raw bytes (RFC1952) via the core's deflate_stream -- the write-side codec for chunked
// stores. Produces exactly what the Python/C++ writers emit (.zarray compressor {"id":"gzip"}), so the
// bytes are decodable by zarrita's GzipCodec, Python numcodecs.GZip, and the C++ reader unchanged.
// (Built with -sUSE_ZLIB=1 -DLSTAR_HAVE_ZLIB; deflate/inflate are otherwise excluded from the kernels.)
static val gzipCompress(val bytes_js, int level) {
    std::vector<uint8_t> src = convertJSArrayToNumberVector<uint8_t>(bytes_js);
    std::vector<uint8_t> out = lstar::deflate_stream(src.data(), src.size(), level, /*gzip=*/true);
    return val(typed_memory_view(out.size(), out.data())).call<val>("slice");
}

// viewer@0.1: 1-vs-rest markers from per-(group,gene) stats. S, NE are Float64Array flat group-major
// (g*ngenes + gene); nper = group sizes; ncells = total. -> {lfc, padj} flat GENE-major (gene*ngroups
// + g) -- the spec's ng x K orientation.
static val markersOneVsRest(val S_js, val NE_js, val nper_js, int ngroups, int ngenes, double ncells) {
    std::vector<double> S = convertJSArrayToNumberVector<double>(S_js);
    std::vector<double> NE = convertJSArrayToNumberVector<double>(NE_js);
    std::vector<int64_t> nper = to_i64(nper_js);
    auto m = lstar::markers_one_vs_rest(S.data(), NE.data(), nper.data(), ngroups, (int64_t)ngenes, (int64_t)ncells);
    val out = val::object();
    out.set("lfc", to_f64(m.lfc)); out.set("padj", to_f64(m.padj));
    out.set("ngenes", (double)ngenes); out.set("ngroups", (double)ngroups);
    return out;
}

// viewer@0.1: per-gene overdispersion score (pagoda2 lowess + F-test). mean/var Float64Array, nobs
// per-gene (expressing cells). -> Float64Array od. The same kernel the prep uses, so live == prepped.
static val overdispersion(val mean_js, val var_js, val nobs_js) {
    std::vector<double> mean = convertJSArrayToNumberVector<double>(mean_js);
    std::vector<double> var = convertJSArrayToNumberVector<double>(var_js);
    std::vector<int64_t> nobs = to_i64(nobs_js);
    auto od = lstar::overdispersion(mean.data(), var.data(), nobs.data(), (int64_t)mean.size());
    return to_f64(od);
}

static std::string version() { return "lstar-wasm 0.0.4"; }

EMSCRIPTEN_BINDINGS(lstar_wasm) {
    function("colMeanVar", &colMeanVar);
    function("cscToCsr", &cscToCsr);
    function("colSumByGroup", &colSumByGroup);
    function("subsampleDeRank", &subsampleDeRank);
    function("markersOneVsRest", &markersOneVsRest);
    function("overdispersion", &overdispersion);
    function("gzipCompress", &gzipCompress);
    function("version", &version);
}
