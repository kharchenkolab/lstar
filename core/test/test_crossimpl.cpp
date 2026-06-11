// Cross-implementation test: read an L* store written by Python, validate it, write it
// back out. A Python script then reads the C++-written store and asserts equality.
//
//   usage: test_crossimpl <in.lstar.zarr> <out.lstar.zarr>
#include <cmath>
#include <iostream>

#include "lstar/lstar.hpp"

static int fail(const std::string& msg) {
    std::cerr << "FAIL: " << msg << "\n";
    return 1;
}

int main(int argc, char** argv) {
    std::string in = argc > 1 ? argv[1] : "/tmp/cx.lstar.zarr";
    std::string out = argc > 2 ? argv[2] : "/tmp/cx_cpp.lstar.zarr";

    lstar::Dataset ds = lstar::read(in);
    std::cout << "read " << in << ": kind=" << ds.kind << " axes=" << ds.axes.size()
              << " fields=" << ds.fields.size() << "\n";

    auto* cells = ds.axis("cells");
    if (!cells || cells->labels.size() != 100) return fail("cells axis");
    if (cells->labels.front() != "cell0" || cells->labels.back() != "cell99")
        return fail("cell labels: " + cells->labels.front() + ".." + cells->labels.back());

    auto* pcaAx = ds.axis("pca");
    if (!pcaAx || pcaAx->origin != "derived" || pcaAx->labels.size() != 10) return fail("pca axis");

    auto* counts = ds.field("counts");
    if (!counts || counts->encoding != "csc" || counts->state != "raw") return fail("counts field");
    const double* d = counts->data.as<double>();
    double s = 0;
    for (int64_t i = 0; i < counts->data.nelem(); ++i) s += d[i];
    std::cout << "  counts: csc " << counts->shape[0] << "x" << counts->shape[1]
              << " nnz=" << counts->data.nelem() << " sum=" << s << "\n";

    auto* pca = ds.field("pca");
    if (!pca || pca->role != "embedding" || pca->dense.dtype != "<f4") return fail("pca field");
    if (pca->dense.shape != std::vector<int64_t>{100, 10}) return fail("pca shape");

    auto* leiden = ds.field("leiden");
    if (!leiden || leiden->encoding != "utf8" || leiden->strings.size() != 100)
        return fail("leiden field");
    std::cout << "  leiden: utf8 label, [0]=" << leiden->strings[0] << "\n";

    lstar::write(ds, out);
    std::cout << "wrote " << out << " (C++ round-trip)\n";
    return 0;
}
