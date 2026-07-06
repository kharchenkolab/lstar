// v3-format conformance: read a store (any format), and either (a) re-emit it as v3 and confirm the
// read-back is value-identical (write mode), or (b) compare it to another store (compare mode). The
// comparator checks every axis/field/aux value across all encodings, so a format that drops or reshapes
// anything fails loudly. Read auto-probes v2 (.zarray/.zgroup) vs v3 (zarr.json); write mode emits v3.
#include "lstar/lstar.hpp"
#include <iostream>
#include <string>

using namespace lstar;
static int fails = 0;
static void ck(const std::string& n, bool ok) {
    if (!ok) { std::cout << "  FAIL  " << n << "\n"; fails++; }
}
static bool eqnd(const NdArray& a, const NdArray& b) {
    return a.dtype == b.dtype && a.shape == b.shape && a.bytes == b.bytes;
}
static void cmp(const Dataset& A, const Dataset& B) {
    ck("kind", A.kind == B.kind);
    ck("spec_version", A.spec_version == B.spec_version);
    ck("profiles", A.profiles == B.profiles);
    ck("dropped", A.dropped == B.dropped);
    ck("axis count", A.axes.size() == B.axes.size());
    for (size_t i = 0; i < std::min(A.axes.size(), B.axes.size()); ++i) {
        auto& x = A.axes[i]; auto& y = B.axes[i];
        ck("axis[" + x.name + "]", x.name == y.name && x.origin == y.origin && x.role == y.role
            && x.induced_by == y.induced_by && x.labels == y.labels && x.provenance == y.provenance);
    }
    ck("field count", A.fields.size() == B.fields.size());
    for (size_t i = 0; i < std::min(A.fields.size(), B.fields.size()); ++i) {
        auto& x = A.fields[i]; auto& y = B.fields[i];
        std::string t = "field[" + x.name + "]";
        ck(t + " meta", x.name == y.name && x.role == y.role && x.encoding == y.encoding
            && x.state == y.state && x.subtype == y.subtype && x.coverage == y.coverage
            && x.span == y.span && x.index_axis == y.index_axis);
        ck(t + " flags", x.has_directed == y.has_directed && x.directed == y.directed
            && x.has_weighted == y.has_weighted && x.weighted == y.weighted
            && x.has_ordered == y.has_ordered && x.ordered == y.ordered);
        ck(t + " shape", x.shape == y.shape);
        if (x.encoding == "csr" || x.encoding == "csc")
            ck(t + " csc arrays", eqnd(x.data, y.data) && eqnd(x.indices, y.indices) && eqnd(x.indptr, y.indptr));
        else if (x.encoding == "utf8") ck(t + " strings", x.strings == y.strings);
        else if (x.encoding == "categorical") ck(t + " codes+cats", eqnd(x.codes, y.codes) && x.categories == y.categories);
        else ck(t + " dense", eqnd(x.dense, y.dense));
        ck(t + " mask", x.has_mask == y.has_mask && (!x.has_mask || eqnd(x.mask, y.mask)));
        ck(t + " index", x.has_index == y.has_index && (!x.has_index || eqnd(x.index, y.index)));
        ck(t + " prov/uncert", x.provenance == y.provenance && x.uncertainty == y.uncertainty);
    }
    ck("aux count", A.aux.size() == B.aux.size());
    for (size_t i = 0; i < std::min(A.aux.size(), B.aux.size()); ++i) {
        auto& x = A.aux[i]; auto& y = B.aux[i];
        ck("aux[" + x.ns + "]", x.ns == y.ns && x.attrs == y.attrs && x.leaves.size() == y.leaves.size());
        for (size_t j = 0; j < std::min(x.leaves.size(), y.leaves.size()); ++j) {
            auto& p = x.leaves[j]; auto& q = y.leaves[j];
            ck("aux leaf[" + p.id + "]", p.id == q.id && p.kind == q.kind && p.strings == q.strings && eqnd(p.dense, q.dense));
        }
    }
}
int main(int argc, char** argv) {
    if (argc < 4) { std::cerr << "usage: test_v3 write|shard|compare <a> <b>\n"; return 2; }
    std::string mode = argv[1], a = argv[2], b = argv[3];
    json gz; gz["id"] = "gzip"; gz["level"] = 5;
    if (mode == "write") {                               // read a (any fmt) -> write b as v3 -> compare
        Dataset A = read(a);
        write(A, b, 0, gz, zarr::ZarrFormat::v3);
        Dataset B = read(b);
        cmp(A, B);
    } else if (mode == "shard") {                        // read a -> write b as chunked+SHARDED v3 -> compare
        Dataset A = read(a);
        write(A, b, /*chunk_elems*/1000, gz, zarr::ZarrFormat::v3, /*shard_elems*/4000);
        Dataset B = read(b);
        cmp(A, B);
    } else if (mode == "compare") {                      // read both -> compare (either format)
        cmp(read(a), read(b));
    } else { std::cerr << "unknown mode: " << mode << "\n"; return 2; }
    std::cout << (fails ? "  v3 " + mode + ": " + std::to_string(fails) + " FAILURES\n"
                        : "  v3 " + mode + ": value-identical\n");
    return fails ? 1 : 0;
}
