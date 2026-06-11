// Cross-impl test for the upgraded reader: read a *chunked, gzip-compressed* store written by
// Python, check the heavy fields decode correctly, exercise the csc<->csr transpose primitive,
// and write the store back (now with a consolidated .zmetadata) for Python to re-open.
//
//   usage: test_chunked <in.lstar.zarr> <out.lstar.zarr> <expected_nnz> <expected_sum>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iostream>

#include "lstar/lstar.hpp"

static int fail(const std::string& m) { std::cerr << "FAIL: " << m << "\n"; return 1; }

int main(int argc, char** argv) {
    std::string in = argc > 1 ? argv[1] : "/tmp/cx_chunked.lstar.zarr";
    std::string out = argc > 2 ? argv[2] : "/tmp/cx_chunked_cpp.lstar.zarr";
    int64_t exp_nnz = argc > 3 ? std::atoll(argv[3]) : -1;
    double exp_sum = argc > 4 ? std::atof(argv[4]) : 0.0;

    lstar::Dataset ds = lstar::read(in);
    auto* counts = ds.field("counts");
    if (!counts || counts->encoding != "csc") return fail("counts field/enc");
    int64_t nnz = counts->data.nelem();
    // index arrays normalize to int64 (small for indptr; indices only needed by the transpose);
    // the value array stays in its stored dtype (float32 here) -- no widening copy.
    std::vector<int64_t> idx = lstar::as_i64(counts->indices);
    std::vector<int64_t> ip = lstar::as_i64(counts->indptr);
    auto* emb = ds.field("pca");
    if (!emb || emb->dense.shape.size() != 2) return fail("pca dense");
    int many = 1;
#ifdef _OPENMP
    many = std::min(16, omp_get_max_threads());   // capture before a serial run pins threads to 1
#endif

    // Everything below runs on the native value pointer T* (float or double) -- memory-lean, with
    // double accumulation inside the kernels.
    auto run = [&](const auto* d) -> int {
        double s = 0;
        for (int64_t i = 0; i < nnz; ++i) s += static_cast<double>(d[i]);
        std::cout << "  [c++] read chunked+gzip: counts csc " << counts->shape[0] << "x"
                  << counts->shape[1] << " nnz=" << nnz << " sum=" << s << "\n";
        if (exp_nnz >= 0 && nnz != exp_nnz) return fail("nnz mismatch");
        if (exp_nnz >= 0 && std::abs(s - exp_sum) > 1e-5 * std::abs(exp_sum)) return fail("sum mismatch");

        auto csr = lstar::csc_to_csr(d, idx.data(), ip.data(), counts->shape[0], counts->shape[1]);
        double s2 = 0;
        for (auto x : csr.data) s2 += static_cast<double>(x);
        if ((int64_t)csr.data.size() != nnz) return fail("transpose nnz");
        if (std::abs(s2 - s) > 1e-5 * std::abs(s)) return fail("transpose sum");
        if (csr.indptr.front() != 0 || csr.indptr.back() != nnz) return fail("transpose indptr ends");
        for (size_t i = 1; i < csr.indptr.size(); ++i)
            if (csr.indptr[i] < csr.indptr[i - 1]) return fail("transpose indptr monotone");
        std::cout << "  [c++] csc_to_csr (lean, keeps value dtype): nnz=" << csr.data.size()
                  << " preserved, indptr[" << csr.indptr.size() << "] -> " << csr.indptr.back() << "\n";

        auto timed = [&](int nt) {
            auto t0 = std::chrono::high_resolution_clock::now();
            auto st = lstar::csc_col_mean_var(d, ip.data(), counts->shape[1], counts->shape[0], nt, true);
            auto t1 = std::chrono::high_resolution_clock::now();
            return std::make_pair(std::chrono::duration<double, std::milli>(t1 - t0).count(), st);
        };
        auto one = timed(1);
        auto all = timed(many);
        bool same = one.second.mean == all.second.mean && one.second.var == all.second.var;
        std::cout << "  [c++] csc_col_mean_var (log1p) over " << counts->shape[1] << " genes: 1 thread "
                  << one.first << " ms, " << many << " threads " << all.first << " ms ("
                  << (one.first / all.first) << "x); results identical=" << (same ? "yes" : "no") << "\n";
        return same ? 0 : fail("threaded reduction result differs from serial");
    };

    int rc = (counts->data.dtype == "<f4") ? run(counts->data.as<float>())
                                           : run(counts->data.as<double>());
    if (rc != 0) return rc;

    lstar::write(ds, out);
    std::cout << "  [c++] wrote " << out << " (+ consolidated .zmetadata)\n";
    return 0;
}
