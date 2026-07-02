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

// Per-(group, gene) sufficient stats over a CSC measure: sum, sumsq, nnz of each group's cells.
// `group_of_cell` (length nrows) maps a cell/row to a group id in [0, ngroups). Returns three
// (ngroups, ngenes) f64 arrays — the viewer profile's cluster stats, computed on the shared core.
static py::tuple col_sum_by_group(py::array data, py::array indptr, py::array indices,
                                  int64_t nrows, int64_t ncols, py::array group_of_cell,
                                  int ngroups, bool lognorm, int n_threads) {
    auto ip = py::array_t<int64_t, py::array::c_style | py::array::forcecast>::ensure(indptr);
    auto idx = py::array_t<int64_t, py::array::c_style | py::array::forcecast>::ensure(indices);
    auto grp = py::array_t<int32_t, py::array::c_style | py::array::forcecast>::ensure(group_of_cell);
    if (!ip || ip.size() != ncols + 1) throw std::runtime_error("indptr must have length ncols+1");
    if (!grp || grp.size() != nrows) throw std::runtime_error("group_of_cell must have length nrows");
    const int64_t* ipp = ip.data();
    const int64_t* idxp = idx.data();
    const int* gp = grp.data();

    lstar::GroupStats s;
    const auto dt = data.dtype();
    if (dt.kind() == 'f' && dt.itemsize() == 4) {
        auto d = py::array_t<float, py::array::c_style>::ensure(data);
        const float* dp = d.data();
        py::gil_scoped_release rel;
        s = lstar::csc_col_sum_by_group(dp, ipp, idxp, nrows, ncols, gp, ngroups, lognorm, n_threads);
    } else {
        auto d = py::array_t<double, py::array::c_style | py::array::forcecast>::ensure(data);
        const double* dp = d.data();
        py::gil_scoped_release rel;
        s = lstar::csc_col_sum_by_group(dp, ipp, idxp, nrows, ncols, gp, ngroups, lognorm, n_threads);
    }
    const std::vector<py::ssize_t> shape{(py::ssize_t)s.ngroups, (py::ssize_t)s.ngenes};
    return py::make_tuple(py::array_t<double>(shape, s.sum.data()),
                          py::array_t<double>(shape, s.sumsq.data()),
                          py::array_t<double>(shape, s.n_expr.data()));
}

// Two-group mean(log1p) difference over a CSR *cell-major* measure (e.g. the viewer de_panel):
// indptr is over rows/cells (length nrows+1), indices are gene ids. membership (length nrows) is
// 0=A, 1=B, <0=excluded. Returns (meanA, meanB, lfc) f64 over genes plus the per-side cell counts.
static py::tuple subsample_de_rank(py::array data, py::array indptr, py::array indices,
                                   int64_t nrows, int64_t ngenes, py::array membership, bool lognorm) {
    auto ip = py::array_t<int64_t, py::array::c_style | py::array::forcecast>::ensure(indptr);
    auto idx = py::array_t<int64_t, py::array::c_style | py::array::forcecast>::ensure(indices);
    auto mem = py::array_t<int32_t, py::array::c_style | py::array::forcecast>::ensure(membership);
    if (!ip || ip.size() != nrows + 1) throw std::runtime_error("indptr must have length nrows+1 (CSR cell-major)");
    if (!mem || mem.size() != nrows) throw std::runtime_error("membership must have length nrows");
    const int64_t* ipp = ip.data();
    const int64_t* idxp = idx.data();
    const int* mp = mem.data();

    lstar::DERank r;
    const auto dt = data.dtype();
    if (dt.kind() == 'f' && dt.itemsize() == 4) {
        auto d = py::array_t<float, py::array::c_style>::ensure(data);
        const float* dp = d.data();
        py::gil_scoped_release rel;
        r = lstar::subsample_de_rank(dp, ipp, idxp, nrows, ngenes, mp, lognorm);
    } else {
        auto d = py::array_t<double, py::array::c_style | py::array::forcecast>::ensure(data);
        const double* dp = d.data();
        py::gil_scoped_release rel;
        r = lstar::subsample_de_rank(dp, ipp, idxp, nrows, ngenes, mp, lognorm);
    }
    return py::make_tuple(py::array_t<double>(r.meanA.size(), r.meanA.data()),
                          py::array_t<double>(r.meanB.size(), r.meanB.data()),
                          py::array_t<double>(r.lfc.size(), r.lfc.data()), r.nA, r.nB);
}

// viewer@0.1: 1-vs-rest markers from per-(group,gene) stats. S, NE are (ngroups, ngenes) f64
// (group-major, over log1p); nper = group sizes (int64, ngroups); ncells = total. Returns
// (lfc, padj) each (ngenes, ngroups) f64 — the spec's gene-major ng×K orientation.
static py::tuple markers_one_vs_rest(py::array S_, py::array NE_, py::array nper_, int64_t ncells) {
    auto S = py::array_t<double, py::array::c_style | py::array::forcecast>::ensure(S_);
    auto NE = py::array_t<double, py::array::c_style | py::array::forcecast>::ensure(NE_);
    auto np = py::array_t<int64_t, py::array::c_style | py::array::forcecast>::ensure(nper_);
    if (!S || S.ndim() != 2) throw std::runtime_error("S must be a 2-D (ngroups, ngenes) array");
    const int ngroups = (int)S.shape(0); const int64_t ngenes = S.shape(1);
    if (NE.ndim() != 2 || NE.shape(0) != ngroups || NE.shape(1) != ngenes)
        throw std::runtime_error("NE must match S shape");
    if (np.size() != ngroups) throw std::runtime_error("nper must have length ngroups");
    lstar::Markers mk;
    { py::gil_scoped_release rel; mk = lstar::markers_one_vs_rest(S.data(), NE.data(), np.data(), ngroups, ngenes, ncells); }
    const std::vector<py::ssize_t> shape{(py::ssize_t)ngenes, (py::ssize_t)ngroups};
    return py::make_tuple(py::array_t<double>(shape, mk.lfc.data()),
                          py::array_t<double>(shape, mk.padj.data()));
}

// viewer@0.1: per-gene overdispersion score (pagoda2 lowess + F-test). mean/var f64 and nobs i64,
// each length ngenes (e.g. from col_mean_var over log1p). Returns od f64 (ngenes).
static py::array_t<double> overdispersion(py::array mean_, py::array var_, py::array nobs_) {
    auto mean = py::array_t<double, py::array::c_style | py::array::forcecast>::ensure(mean_);
    auto var  = py::array_t<double, py::array::c_style | py::array::forcecast>::ensure(var_);
    auto nobs = py::array_t<int64_t, py::array::c_style | py::array::forcecast>::ensure(nobs_);
    const int64_t ng = mean.size();
    if (var.size() != ng || nobs.size() != ng) throw std::runtime_error("mean/var/nobs length mismatch");
    std::vector<double> od;
    { py::gil_scoped_release rel; od = lstar::overdispersion(mean.data(), var.data(), nobs.data(), ng); }
    return py::array_t<double>(od.size(), od.data());
}

// viewer@0.1: canonical cell order (the single source all bindings share). primary_code int32 (ncells);
// emb float64 (ncells, >=2; first 2 cols used) or None -> cluster-only. Returns pos_of int64 (ncells).
static py::array_t<int64_t> viewer_cell_order(py::array primary_code, py::object emb, int64_t grid) {
    auto pc = py::array_t<int32_t, py::array::c_style | py::array::forcecast>::ensure(primary_code);
    if (!pc || pc.ndim() != 1) throw std::runtime_error("primary_code must be a 1-D int array");
    const int64_t ncells = pc.size();
    std::vector<int64_t> pos;
    if (emb.is_none()) {
        py::gil_scoped_release rel;
        pos = lstar::viewer_cell_order(pc.data(), nullptr, ncells, grid);
    } else {
        auto e = py::array_t<double, py::array::c_style | py::array::forcecast>::ensure(emb);
        if (!e || e.ndim() != 2 || e.shape(1) < 2) throw std::runtime_error("emb must be (ncells, >=2)");
        if (e.shape(0) != ncells) throw std::runtime_error("emb rows must equal ncells");
        const int64_t ecols = e.shape(1);
        const double* ep = e.data();
        std::vector<double> emb2((size_t)ncells * 2);          // pack first 2 cols row-major (matches Python emb[:, :2])
        for (int64_t i = 0; i < ncells; i++) { emb2[(size_t)(2 * i)] = ep[i * ecols]; emb2[(size_t)(2 * i + 1)] = ep[i * ecols + 1]; }
        py::gil_scoped_release rel;
        pos = lstar::viewer_cell_order(pc.data(), emb2.data(), ncells, grid);
    }
    return py::array_t<int64_t>(pos.size(), pos.data());
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
    m.def("col_sum_by_group", &col_sum_by_group, py::arg("data"), py::arg("indptr"),
          py::arg("indices"), py::arg("nrows"), py::arg("ncols"), py::arg("group_of_cell"),
          py::arg("ngroups"), py::arg("lognorm") = false, py::arg("n_threads") = 0);
    m.def("subsample_de_rank", &subsample_de_rank, py::arg("data"), py::arg("indptr"),
          py::arg("indices"), py::arg("nrows"), py::arg("ngenes"), py::arg("membership"),
          py::arg("lognorm") = true);
    m.def("markers_one_vs_rest", &markers_one_vs_rest, py::arg("S"), py::arg("NE"),
          py::arg("nper"), py::arg("ncells"));
    m.def("overdispersion", &overdispersion, py::arg("mean"), py::arg("var"), py::arg("nobs"));
    m.def("max_threads", &max_threads);
#ifdef _OPENMP
    m.attr("openmp") = true;
#else
    m.attr("openmp") = false;
#endif
}
