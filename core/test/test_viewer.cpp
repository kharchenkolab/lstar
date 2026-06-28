// Unit test for the viewer@0.1 recipe kernels (markers_one_vs_rest, overdispersion + its incomplete-
// beta F-tail, hilbert_xy2d, cell_order_pos). Numeric cross-language agreement is checked separately
// by conformance/viewer.sh; this pins the C++ math against hand-rolled references.
#include <cmath>
#include <iostream>
#include <vector>

#include "lstar/lstar.hpp"

static int failures = 0;
static void check(const char* name, bool ok) {
    std::cout << "  " << (ok ? "OK  " : "FAIL") << "  " << name << "\n";
    if (!ok) failures++;
}

int main() {
    using namespace lstar;

    // --- markers_one_vs_rest: orientation (ng x K, gene-major) + the 1-vs-rest formula ---------
    {
        const int K = 2; const int64_t ng = 3, ncells = 10;
        std::vector<double> S  = {2.0, 0.0, 4.0,   3.0, 6.0, 1.0};   // group-major g*ng+gene
        std::vector<double> NE = {3.0, 0.0, 4.0,   5.0, 6.0, 1.0};
        std::vector<int64_t> nper = {4, 6};
        Markers m = markers_one_vs_rest(S.data(), NE.data(), nper.data(), K, ng, ncells);
        bool ok = (m.ngenes == ng && m.ngroups == K
                   && (int64_t)m.lfc.size() == ng * K && (int64_t)m.padj.size() == ng * K);
        // reference: lfc[gene,g] = S[g,gene]/nper[g] - (grand[gene]-S[g,gene])/(ncells-nper[g])
        for (int64_t j = 0; j < ng && ok; ++j) {
            double grand = S[0 * ng + j] + S[1 * ng + j];
            for (int g = 0; g < K; ++g) {
                double sg = S[(size_t)g * ng + j];
                double exp_lfc = sg / nper[g] - (grand - sg) / (ncells - nper[g]);
                double exp_padj = std::exp(-std::fabs(exp_lfc * std::sqrt(NE[(size_t)g * ng + j] + 1.0)));
                exp_padj = exp_padj < 1e-12 ? 1e-12 : (exp_padj > 1.0 ? 1.0 : exp_padj);
                ok = ok && std::fabs(m.lfc[(size_t)j * K + g] - exp_lfc) < 1e-12
                        && std::fabs(m.padj[(size_t)j * K + g] - exp_padj) < 1e-12;
            }
        }
        check("markers_one_vs_rest: gene-major orientation + formula", ok);
    }

    // --- reg_incbeta_log: F-tail at known points (P(F>1; d,d) = 0.5) ----------------------------
    {
        // upper-tail p = I_{1/(1+f)}(d/2, d/2); f=1, d=10 -> I_0.5(5,5) = 0.5 -> -log = ln 2
        double lp = reg_incbeta_log(5.0, 5.0, 0.5);
        check("reg_incbeta_log(5,5,0.5) == log(0.5)", std::fabs(lp - std::log(0.5)) < 1e-9);
        // overdispersed (f>1 -> x<0.5) gives a smaller tail -> more negative log -> larger -log
        double lp_hi = reg_incbeta_log(5.0, 5.0, 1.0 / (1.0 + 4.0));   // f=4
        check("reg_incbeta_log monotone in f (more OD -> larger score)", -lp_hi > -lp);
    }

    // --- overdispersion: a gene off the mean-variance trend scores highest ----------------------
    {
        const int64_t ng = 41;
        std::vector<double> mean(ng), var(ng); std::vector<int64_t> nobs(ng, 60);
        for (int64_t j = 0; j < ng; ++j) {                 // 40 genes on the line log(var)=1.2*log(mean)+0.1
            double m = 0.2 + 0.1 * (double)j;
            mean[j] = m; var[j] = std::exp(1.2 * std::log(m) + 0.1);
        }
        var[20] *= 6.0;                                    // one strongly overdispersed gene
        std::vector<double> od = overdispersion(mean.data(), var.data(), nobs.data(), ng);
        bool argmax_is_20 = true; for (int64_t j = 0; j < ng; ++j) if (od[j] > od[20]) argmax_is_20 = false;
        check("overdispersion: off-trend gene scores highest", argmax_is_20 && od[20] > 0.0);
        // a gene that is on the trend has a near-baseline score (residual ~0 -> f~1 -> ~ln2)
        check("overdispersion: on-trend gene near baseline", od[5] < od[20] && od[5] >= 0.0);
    }

    // --- hilbert_xy2d: the canonical 4x4 order (matches reorder.mjs / viewer.py test) -----------
    {
        const int64_t expected[16] = {0, 1, 14, 15, 3, 2, 13, 12, 4, 7, 8, 11, 5, 6, 9, 10};
        bool ok = true;
        for (int64_t y = 0; y < 4; ++y) for (int64_t x = 0; x < 4; ++x)
            ok = ok && hilbert_xy2d(4, x, y) == expected[y * 4 + x];
        check("hilbert_xy2d: canonical 4x4 curve", ok);
    }

    // --- cell_order_pos: cluster-contiguous valid permutation -----------------------------------
    {
        std::vector<int> code = {1, 0, 1, 0, 2, 1, 0};     // 7 cells, 3 clusters
        std::vector<int64_t> pos = cell_order_pos(code.data(), nullptr, (int64_t)code.size());
        // valid permutation
        std::vector<int> seen(7, 0); bool perm = true;
        for (int64_t p : pos) { if (p < 0 || p >= 7 || seen[(size_t)p]) perm = false; else seen[(size_t)p] = 1; }
        // cluster-major: cells of cluster 0 occupy the lowest physical rows, then 1, then 2
        bool contiguous = true; int64_t prev = -1;
        for (int c = 0; c <= 2; ++c) for (size_t i = 0; i < code.size(); ++i) if (code[i] == c) {
            if (pos[i] < prev) contiguous = false; prev = pos[i];
        }
        check("cell_order_pos: valid cluster-contiguous permutation", perm && contiguous);
    }

    std::cout << (failures == 0 ? "\ntest_viewer OK\n" : "\ntest_viewer FAIL\n");
    return failures == 0 ? 0 : 1;
}
