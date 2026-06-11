// lstar._accel — pybind11 binding over libstar's translation primitives.
//
// Exposes the compute kernels (not the IO) so the Python package can run reductions at C++/OpenMP
// speed when this extension is present, falling back to the pure-Python implementation when it is
// not. Inputs are taken zero-copy in their stored dtype (float32 measures stay float32); the GIL
// is released around the kernel so OpenMP scales and other threads run.
#include <cstdint>

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include "lstar/lstar.hpp"

namespace py = pybind11;

// Per-column zero-aware mean/variance of a CSC block/matrix. `data` is float32 or float64 (kept in
// place); `indptr` is int32 or int64 (normalized to int64). Returns (mean f64, var f64, nnz i64).
static py::tuple col_mean_var(py::array data, py::array indptr, int64_t nrows,
                              int n_threads, bool lognorm) {
    auto ip = py::array_t<int64_t, py::array::c_style | py::array::forcecast>::ensure(indptr);
    if (!ip || ip.ndim() != 1 || ip.size() < 1)
        throw std::runtime_error("indptr must be a 1-D integer array of length ncols+1");
    const int64_t ncols = ip.size() - 1;
    const int64_t* ipp = ip.data();

    lstar::ColStats s;
    const auto dt = data.dtype();
    if (dt.kind() == 'f' && dt.itemsize() == 4) {
        auto d = py::array_t<float, py::array::c_style>::ensure(data);   // zero-copy if contiguous
        const float* dp = d.data();
        py::gil_scoped_release rel;
        s = lstar::csc_col_mean_var(dp, ipp, ncols, nrows, n_threads, lognorm);
    } else {
        auto d = py::array_t<double, py::array::c_style | py::array::forcecast>::ensure(data);
        const double* dp = d.data();
        py::gil_scoped_release rel;
        s = lstar::csc_col_mean_var(dp, ipp, ncols, nrows, n_threads, lognorm);
    }
    return py::make_tuple(py::array_t<double>(s.mean.size(), s.mean.data()),
                          py::array_t<double>(s.var.size(), s.var.data()),
                          py::array_t<int64_t>(s.nnz.size(), s.nnz.data()));
}

// CSC -> CSR storage transpose (orientation flip), value dtype preserved.
static py::tuple csc_to_csr(py::array data, py::array indices, py::array indptr,
                            int64_t nrows, int64_t ncols) {
    auto idx = py::array_t<int64_t, py::array::c_style | py::array::forcecast>::ensure(indices);
    auto ip = py::array_t<int64_t, py::array::c_style | py::array::forcecast>::ensure(indptr);
    const int64_t* idxp = idx.data();
    const int64_t* ipp = ip.data();
    const auto dt = data.dtype();
    if (dt.kind() == 'f' && dt.itemsize() == 4) {
        auto d = py::array_t<float, py::array::c_style>::ensure(data);
        const float* dp = d.data();
        lstar::CsxArrays<float> r;
        { py::gil_scoped_release rel; r = lstar::csc_to_csr(dp, idxp, ipp, nrows, ncols); }
        return py::make_tuple(py::array_t<float>(r.data.size(), r.data.data()),
                              py::array_t<int64_t>(r.indices.size(), r.indices.data()),
                              py::array_t<int64_t>(r.indptr.size(), r.indptr.data()));
    }
    auto d = py::array_t<double, py::array::c_style | py::array::forcecast>::ensure(data);
    const double* dp = d.data();
    lstar::CsxArrays<double> r;
    { py::gil_scoped_release rel; r = lstar::csc_to_csr(dp, idxp, ipp, nrows, ncols); }
    return py::make_tuple(py::array_t<double>(r.data.size(), r.data.data()),
                          py::array_t<int64_t>(r.indices.size(), r.indices.data()),
                          py::array_t<int64_t>(r.indptr.size(), r.indptr.data()));
}

static int max_threads() {
#ifdef _OPENMP
    return omp_get_max_threads();
#else
    return 1;
#endif
}

PYBIND11_MODULE(_accel, m) {
    m.doc() = "libstar compute kernels (OpenMP) for the lstar Python package";
    m.def("col_mean_var", &col_mean_var, py::arg("data"), py::arg("indptr"), py::arg("nrows"),
          py::arg("n_threads") = 0, py::arg("lognorm") = false);
    m.def("csc_to_csr", &csc_to_csr, py::arg("data"), py::arg("indices"), py::arg("indptr"),
          py::arg("nrows"), py::arg("ncols"));
    m.def("max_threads", &max_threads);
#ifdef _OPENMP
    m.attr("openmp") = true;
#else
    m.attr("openmp") = false;
#endif
}
