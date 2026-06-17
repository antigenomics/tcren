// Fitting-alignment discrimination for MHC pseudosequence (MPS) matching.
//
// A NetMHCpan pseudosequence is 34 groove residues scattered along the MHC chain. To find
// which allele's pseudosequence a chain carries we thread each candidate 34-mer fully through
// the chain, leaving chain residues free to be skipped (a "fitting" / semi-global alignment).
// The best-scoring candidate identifies the allele; the alignment of that candidate marks the
// (scattered) residues. This file is the hot path: scoring ~4k candidates per chain.
//
// Scoring matches Bio.Align's fitting configuration exactly (BLOSUM62, query/placed gaps
// open -11 / extend -1, target/chain gaps free), so the Python fallback in tcren.mhc.pseudo
// gives identical results.

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <algorithm>
#include <array>
#include <string>
#include <utility>
#include <vector>

namespace py = pybind11;

namespace {

// BLOSUM62, alphabet "ARNDCQEGHILKMFPSTWYVBZX*" (Bio.Align.substitution_matrices order).
constexpr int NA = 24;
constexpr char ALPHABET[NA + 1] = "ARNDCQEGHILKMFPSTWYVBZX*";
constexpr int BLOSUM62[NA][NA] = {
    {4,-1,-2,-2,0,-1,-1,0,-2,-1,-1,-1,-1,-2,-1,1,0,-3,-2,0,-2,-1,0,-4},
    {-1,5,0,-2,-3,1,0,-2,0,-3,-2,2,-1,-3,-2,-1,-1,-3,-2,-3,-1,0,-1,-4},
    {-2,0,6,1,-3,0,0,0,1,-3,-3,0,-2,-3,-2,1,0,-4,-2,-3,3,0,-1,-4},
    {-2,-2,1,6,-3,0,2,-1,-1,-3,-4,-1,-3,-3,-1,0,-1,-4,-3,-3,4,1,-1,-4},
    {0,-3,-3,-3,9,-3,-4,-3,-3,-1,-1,-3,-1,-2,-3,-1,-1,-2,-2,-1,-3,-3,-2,-4},
    {-1,1,0,0,-3,5,2,-2,0,-3,-2,1,0,-3,-1,0,-1,-2,-1,-2,0,3,-1,-4},
    {-1,0,0,2,-4,2,5,-2,0,-3,-3,1,-2,-3,-1,0,-1,-3,-2,-2,1,4,-1,-4},
    {0,-2,0,-1,-3,-2,-2,6,-2,-4,-4,-2,-3,-3,-2,0,-2,-2,-3,-3,-1,-2,-1,-4},
    {-2,0,1,-1,-3,0,0,-2,8,-3,-3,-1,-2,-1,-2,-1,-2,-2,2,-3,0,0,-1,-4},
    {-1,-3,-3,-3,-1,-3,-3,-4,-3,4,2,-3,1,0,-3,-2,-1,-3,-1,3,-3,-3,-1,-4},
    {-1,-2,-3,-4,-1,-2,-3,-4,-3,2,4,-2,2,0,-3,-2,-1,-2,-1,1,-4,-3,-1,-4},
    {-1,2,0,-1,-3,1,1,-2,-1,-3,-2,5,-1,-3,-1,0,-1,-3,-2,-2,0,1,-1,-4},
    {-1,-1,-2,-3,-1,0,-2,-3,-2,1,2,-1,5,0,-2,-1,-1,-1,-1,1,-3,-1,-1,-4},
    {-2,-3,-3,-3,-2,-3,-3,-3,-1,0,0,-3,0,6,-4,-2,-2,1,3,-1,-3,-3,-1,-4},
    {-1,-2,-2,-1,-3,-1,-1,-2,-2,-3,-3,-1,-2,-4,7,-1,-1,-4,-3,-2,-2,-1,-2,-4},
    {1,-1,1,0,-1,0,0,0,-1,-2,-2,0,-1,-2,-1,4,1,-3,-2,-2,0,0,0,-4},
    {0,-1,0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1,1,5,-2,-2,0,-1,-1,0,-4},
    {-3,-3,-4,-4,-2,-2,-3,-2,-2,-3,-2,-3,-1,1,-4,-3,-2,11,2,-3,-4,-3,-2,-4},
    {-2,-2,-2,-3,-2,-1,-2,-3,2,-1,-1,-2,-1,3,-3,-2,-2,2,7,-1,-3,-2,-1,-4},
    {0,-3,-3,-3,-1,-2,-2,-3,-3,3,1,-2,1,-1,-2,-2,0,-3,-1,4,-3,-2,-1,-4},
    {-2,-1,3,4,-3,0,1,-1,0,-3,-4,0,-3,-3,-2,0,-1,-4,-3,-3,4,1,-1,-4},
    {-1,0,0,1,-3,3,4,-2,0,-3,-3,1,-1,-3,-1,0,-1,-3,-2,-2,1,4,-1,-4},
    {0,-1,-1,-1,-2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-2,0,0,-2,-1,-1,-1,-1,-1,-4},
    {-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,1},
};

constexpr double GAP_OPEN = -11.0;     // gap in the placed (pseudo) sequence
constexpr double GAP_EXTEND = -1.0;
constexpr double NEG = -1e9;

std::array<int, 256> build_index() {
    std::array<int, 256> idx{};
    idx.fill(22);  // 'X' — unknown residues
    for (int i = 0; i < NA; ++i) idx[static_cast<unsigned char>(ALPHABET[i])] = i;
    return idx;
}
const std::array<int, 256> CODE = build_index();

inline int sub(char a, char b) {
    return BLOSUM62[CODE[static_cast<unsigned char>(a)]][CODE[static_cast<unsigned char>(b)]];
}

// Fitting score: `placed` is consumed except for free END gaps (matching Bio's
// query_end_gap_score = 0); `freeq` residues may be skipped at no cost anywhere. Internal
// gaps in `placed` are affine-penalised.
double fitting_score(const std::string& placed, const std::string& freeq) {
    const int m = static_cast<int>(placed.size());
    const int n = static_cast<int>(freeq.size());
    if (m == 0) return 0.0;
    std::vector<double> prevM(n + 1, 0.0), prevG(n + 1, NEG);  // row 0: free leading skips
    std::vector<double> curM(n + 1), curG(n + 1);
    double answer = prevM[n];  // i = 0: whole `placed` is a free trailing end gap → 0
    for (int i = 1; i <= m; ++i) {
        curM[0] = 0.0;     // free leading end gap in `placed`
        curG[0] = NEG;
        const char pa = placed[i - 1];
        for (int j = 1; j <= n; ++j) {
            const double s = sub(pa, freeq[j - 1]);
            curG[j] = std::max(prevM[j] + GAP_OPEN, prevG[j] + GAP_EXTEND);
            curM[j] = std::max(std::max(prevM[j - 1] + s, prevG[j - 1] + s),
                               std::max(curM[j - 1], curG[j - 1]));
        }
        answer = std::max(answer, curM[n]);  // free trailing end gap: may stop after row i
        std::swap(prevM, curM);
        std::swap(prevG, curG);
    }
    return answer;
}

// (best index, best score) over candidates; ties resolved to the first (max). GIL released.
std::pair<int, double> best_hit(const std::string& freeq,
                                const std::vector<std::string>& placed_candidates) {
    int best = -1;
    double best_score = NEG;
    {
        py::gil_scoped_release release;
        for (int k = 0; k < static_cast<int>(placed_candidates.size()); ++k) {
            const double s = fitting_score(placed_candidates[k], freeq);
            if (s > best_score) { best_score = s; best = k; }
        }
    }
    return {best, best_score};
}

// Matched (placed_pos, free_pos) column pairs of the best fitting alignment. Full-matrix
// traceback (m*n) — only called once for the chosen candidate, so memory is small.
std::vector<std::pair<int, int>> align(const std::string& placed, const std::string& freeq) {
    const int m = static_cast<int>(placed.size());
    const int n = static_cast<int>(freeq.size());
    std::vector<std::pair<int, int>> out;
    if (m == 0 || n == 0) return out;

    // M[i][j], G[i][j] with backpointers (0=match/diag, 1=from-G-diag, 2=free-skip, 3=gap).
    std::vector<std::vector<double>> M(m + 1, std::vector<double>(n + 1, NEG));
    std::vector<std::vector<double>> G(m + 1, std::vector<double>(n + 1, NEG));
    std::vector<std::vector<char>> bM(m + 1, std::vector<char>(n + 1, 0));
    for (int j = 0; j <= n; ++j) M[0][j] = 0.0;          // free leading target skips
    for (int i = 1; i <= m; ++i) {
        M[i][0] = 0.0;                                   // free leading end gap in `placed`
        const char pa = placed[i - 1];
        for (int j = 1; j <= n; ++j) {
            G[i][j] = std::max(M[i - 1][j] + GAP_OPEN, G[i - 1][j] + GAP_EXTEND);
            const double s = sub(pa, freeq[j - 1]);
            const double diagM = M[i - 1][j - 1] + s;
            const double diagG = G[i - 1][j - 1] + s;
            const double skipM = M[i][j - 1];
            const double skipG = G[i][j - 1];
            double best = diagM; char bp = 0;
            if (diagG > best) { best = diagG; bp = 1; }
            if (skipM > best) { best = skipM; bp = 2; }
            if (skipG > best) { best = skipG; bp = 3; }
            M[i][j] = best; bM[i][j] = bp;
        }
    }
    // Free trailing end gap in `placed`: start from the row with the best M[i][n].
    int i = 0;
    for (int k = 0; k <= m; ++k) if (M[k][n] >= M[i][n]) i = k;
    int j = n;
    bool inGap = false;
    while (i > 0 && j > 0) {
        if (inGap) {                                    // placed residue i is an internal gap
            if (G[i][j] == M[i - 1][j] + GAP_OPEN) inGap = false;
            --i;
            continue;
        }
        const char bp = bM[i][j];
        if (bp == 2) { --j; }                           // free skip of a free-seq residue
        else if (bp == 3) { inGap = true; }             // entering a gap run
        else {                                          // match column (placed i-1, free j-1)
            out.emplace_back(i - 1, j - 1);
            --i; --j;
        }
    }
    std::reverse(out.begin(), out.end());
    return out;
}

}  // namespace

PYBIND11_MODULE(_align, m) {
    m.doc() = "Fitting-alignment discrimination for MHC pseudosequence matching (BLOSUM62).";
    m.def("fitting_score", &fitting_score, py::arg("placed"), py::arg("free"),
          "Fitting-alignment score: `placed` consumed in full, `free` skips are cost-free.");
    m.def("best_hit", &best_hit, py::arg("free"), py::arg("candidates"),
          "(index, score) of the best-fitting candidate over `free` (GIL released).");
    m.def("align", &align, py::arg("placed"), py::arg("free"),
          "Matched (placed_pos, free_pos) column pairs of the best fitting alignment.");
    m.attr("__version__") = "0.1.0";
}
