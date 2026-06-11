// Benchmark: multithreaded per-column mean/variance over a synthetic CSC matrix.
// Demonstrates the OpenMP threading model of libstar's translation primitives.
//
// Two kernels are timed at 1 thread vs all threads:
//   plain   - raw zero-aware moments; trivial arithmetic per nonzero -> memory-bandwidth-bound,
//             barely scales with threads (the bottleneck is the streaming load of `data`).
//   lognorm - the same moments but over log1p(value), the pagoda2 lazy-normalized-view pattern.
//             The transcendental per nonzero makes it compute-bound, so the column-parallel
//             schedule actually scales toward the core count.
// The point is to show *which* primitives benefit from threading and why, not to claim every
// pass does.
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <vector>

#include "lstar/lstar.hpp"

int main() {
  const int64_t nrows = 100000, ncols = 8000;
  const int64_t per_col = 2000;  // ~16M nonzeros
  std::vector<int64_t> indptr(ncols + 1, 0);
  for (int64_t j = 0; j < ncols; ++j) indptr[j + 1] = indptr[j] + per_col;
  const int64_t nnz = indptr[ncols];
  std::vector<double> data(nnz);
  uint64_t st = 88172645463325252ULL;  // xorshift, deterministic
  auto rnd = [&]() { st ^= st << 13; st ^= st >> 7; st ^= st << 17;
                     return (double)((st >> 11) * (1.0 / 9007199254740992.0)); };
  for (int64_t i = 0; i < nnz; ++i) data[i] = rnd() * 10.0;  // raw-count-like magnitudes

  auto run = [&](int nt, bool lognorm) {
    auto t0 = std::chrono::high_resolution_clock::now();
    auto s = lstar::csc_col_mean_var(data.data(), indptr.data(), ncols, nrows, nt, lognorm);
    auto t1 = std::chrono::high_resolution_clock::now();
    double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
    return std::make_pair(ms, s.var[0]);
  };

  int max_threads = 1;
#ifdef _OPENMP
  max_threads = omp_get_max_threads();
#endif
  // Cap at 16 so the figure is meaningful on big machines (per-column work saturates well
  // before the full core count, and a 56-way run is dominated by scheduling, not compute).
  int many = max_threads < 16 ? max_threads : 16;
  printf("OpenMP max threads = %d (benchmarking 1 vs %d); matrix %lldx%lld (%lld nnz)\n",
         max_threads, many, (long long)nrows, (long long)ncols, (long long)nnz);

  bool ok = true;
  for (int variant = 0; variant < 2; ++variant) {
    bool lognorm = (variant == 1);
    auto one = run(1, lognorm);
    auto all = run(many, lognorm);  // explicit thread count; n_threads=0 would inherit prior state
    printf("  %-8s: 1 thread = %7.1f ms, %2d threads = %7.1f ms (%.2fx); var[0]=%.6g\n",
           lognorm ? "lognorm" : "plain", one.first, many, all.first,
           one.first / all.first, one.second);
    ok = ok && (one.second == all.second);  // identical result regardless of thread count
  }
  return ok ? 0 : 1;
}
