// STORED single-file .lstar.zarr.zip — C++ read/write round-trip + cross-surface converter.
//
//   usage: test_zip <in> <out>
// Reads <in> (a directory store OR a .lstar.zarr.zip) and writes <out> (either form). The shell
// driver (conformance/zip.sh) chains this with the Python surface: Python writes a .zip, C++ reads
// it and writes a .zip, Python reads that back and asserts field-for-field equality. Also self-checks
// a pure C++ round-trip (read <in>, write <out>, re-read <out>) and, for a .zip <out>, that the
// result is a single file whose read-back has the same axes/fields.
#include <iostream>

#include "lstar/lstar.hpp"

namespace fs = std::filesystem;

static int fail(const std::string& m) { std::cerr << "FAIL: " << m << "\n"; return 1; }

int main(int argc, char** argv) {
    if (argc < 3) return fail("usage: test_zip <in> <out>");
    std::string in = argv[1], out = argv[2];

    lstar::Dataset ds = lstar::read(in);
    std::cout << "read " << in << ": kind=" << ds.kind << " axes=" << ds.axes.size()
              << " fields=" << ds.fields.size() << "\n";

    lstar::write(ds, out);
    bool out_is_zip = fs::path(out).extension() == ".zip";
    if (out_is_zip && !fs::is_regular_file(out)) return fail(".zip output must be a single file");
    std::cout << "wrote " << out << (out_is_zip ? " (single-file zip)\n" : " (directory)\n");

    // read the freshly-written store back and check the shape survived the round-trip
    lstar::Dataset rt = lstar::read(out);
    if (rt.axes.size() != ds.axes.size()) return fail("axis count changed on round-trip");
    if (rt.fields.size() != ds.fields.size()) return fail("field count changed on round-trip");
    for (auto& f : ds.fields) {
        auto* g = rt.field(f.name);
        if (!g) return fail("field missing on round-trip: " + f.name);
        if (g->encoding != f.encoding) return fail("encoding changed on round-trip: " + f.name);
    }
    std::cout << "round-trip OK: " << rt.axes.size() << " axes, " << rt.fields.size() << " fields\n";
    return 0;
}
