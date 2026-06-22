// Potential-guided peptide refinement (knowledge-based, NOT physics MD).
//
// A rigid-body Metropolis Monte-Carlo local refinement of a peptide against fixed partner atoms
// (TCR + MHC), scored by a residue-level statistical potential (the TCRen/MJ contact energy) plus a
// soft heavy-atom clash penalty. This is the "easy in C" refinement for tcren — it reuses the dense
// potential matrix and needs no force field. Real ref2015 relaxation = Rosetta (subprocess), not this.
//
// Energy(peptide pose) =
//     Σ over peptide-residue × partner-residue pairs in contact (min heavy-atom dist ≤ cutoff)
//         potential[aa_pep][aa_partner]                       (favourable contacts are negative)
//   + clash_w · Σ over heavy-atom pairs with d < clash_d0  (clash_d0 − d)²
//
// The Python wrapper (tcren.refine.refine_peptide) pre-filters partner atoms to the interface
// shell, so the brute-force O(Np·Nq) energy is small; MC runs thousands of steps in well under a second.

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <algorithm>
#include <cmath>
#include <random>
#include <tuple>
#include <vector>

namespace py = pybind11;

namespace {

struct Vec3 { double x, y, z; };

inline double dist2(const Vec3& a, const Vec3& b) {
    const double dx = a.x - b.x, dy = a.y - b.y, dz = a.z - b.z;
    return dx * dx + dy * dy + dz * dz;
}

// Energy of a peptide pose against the fixed partner atoms. A harmonic restraint to the starting
// coordinates (`pep0`) keeps the refinement LOCAL — without it, minimising Σ(potential over
// contacts) by rigid moves trivially ejects the peptide from its pocket.
double energy(const std::vector<Vec3>& pep, const std::vector<Vec3>& pep0,
              const std::vector<int>& pep_atom_res, const std::vector<int>& pep_res_aa,
              const std::vector<Vec3>& par, const std::vector<int>& par_atom_res,
              const std::vector<int>& par_res_aa,
              const std::vector<double>& pot, int n_aa,
              double cutoff, double clash_d0, double clash_w, double restraint_w,
              int n_pep_res, int n_par_res) {
    const double cut2 = cutoff * cutoff, d0_2 = clash_d0 * clash_d0;
    std::vector<char> in_contact(static_cast<size_t>(n_pep_res) * n_par_res, 0);
    double clash = 0.0, restraint = 0.0;
    for (size_t a = 0; a < pep.size(); ++a) {
        const int ri = pep_atom_res[a];
        restraint += restraint_w * dist2(pep[a], pep0[a]);
        for (size_t b = 0; b < par.size(); ++b) {
            const double d2 = dist2(pep[a], par[b]);
            if (d2 < d0_2) { const double d = std::sqrt(d2); clash += clash_w * (clash_d0 - d) * (clash_d0 - d); }
            if (d2 <= cut2) in_contact[static_cast<size_t>(ri) * n_par_res + par_atom_res[b]] = 1;
        }
    }
    double e = clash + restraint;
    for (int ri = 0; ri < n_pep_res; ++ri)
        for (int rj = 0; rj < n_par_res; ++rj)
            if (in_contact[static_cast<size_t>(ri) * n_par_res + rj])
                e += pot[static_cast<size_t>(pep_res_aa[ri]) * n_aa + par_res_aa[rj]];
    return e;
}

// Apply a rigid transform (rotation about the peptide centroid, then translation) in place.
void rigid_move(std::vector<Vec3>& pep, const Vec3& centroid,
                double ax, double ay, double az, double angle, const Vec3& t) {
    // Rodrigues rotation matrix for axis (ax,ay,az) (assumed ~unit) and `angle`.
    const double c = std::cos(angle), s = std::sin(angle), C = 1.0 - c;
    const double R00 = c + ax * ax * C,      R01 = ax * ay * C - az * s, R02 = ax * az * C + ay * s;
    const double R10 = ay * ax * C + az * s, R11 = c + ay * ay * C,      R12 = ay * az * C - ax * s;
    const double R20 = az * ax * C - ay * s, R21 = az * ay * C + ax * s, R22 = c + az * az * C;
    for (auto& p : pep) {
        const double x = p.x - centroid.x, y = p.y - centroid.y, z = p.z - centroid.z;
        p.x = centroid.x + R00 * x + R01 * y + R02 * z + t.x;
        p.y = centroid.y + R10 * x + R11 * y + R12 * z + t.y;
        p.z = centroid.z + R20 * x + R21 * y + R22 * z + t.z;
    }
}

py::tuple refine(py::array_t<double> pep_xyz, py::array_t<int> pep_atom_res,
                 py::array_t<int> pep_res_aa,
                 py::array_t<double> par_xyz, py::array_t<int> par_atom_res,
                 py::array_t<int> par_res_aa,
                 py::array_t<double> potential,
                 double cutoff, double clash_d0, double clash_w, double restraint_w,
                 int n_steps, double trans_sigma, double rot_sigma,
                 double temp0, double temp1, unsigned int seed) {
    auto load_xyz = [](py::array_t<double> a) {
        auto r = a.unchecked<2>();
        std::vector<Vec3> v(r.shape(0));
        for (py::ssize_t i = 0; i < r.shape(0); ++i) v[i] = {r(i, 0), r(i, 1), r(i, 2)};
        return v;
    };
    auto load_int = [](py::array_t<int> a) {
        auto r = a.unchecked<1>();
        std::vector<int> v(r.shape(0));
        for (py::ssize_t i = 0; i < r.shape(0); ++i) v[i] = r(i);
        return v;
    };
    std::vector<Vec3> pep = load_xyz(pep_xyz), par = load_xyz(par_xyz);
    std::vector<int> pep_ar = load_int(pep_atom_res), pep_ra = load_int(pep_res_aa);
    std::vector<int> par_ar = load_int(par_atom_res), par_ra = load_int(par_res_aa);
    const int n_aa = static_cast<int>(potential.shape(0));
    std::vector<double> pot(potential.size());
    std::copy(potential.data(), potential.data() + potential.size(), pot.begin());
    // `*_res_aa` is indexed by residue (its length is the residue count); `*_atom_res` maps each
    // atom to a residue index in [0, n_res).
    const int n_pep_res = static_cast<int>(pep_ra.size());
    const int n_par_res = static_cast<int>(par_ra.size());

    const std::vector<Vec3> pep0 = pep;  // restraint reference (the input pose)
    auto E = [&](const std::vector<Vec3>& p) {
        return energy(p, pep0, pep_ar, pep_ra, par, par_ar, par_ra, pot, n_aa,
                      cutoff, clash_d0, clash_w, restraint_w, n_pep_res, n_par_res);
    };

    std::mt19937 rng(seed);
    std::normal_distribution<double> gauss(0.0, 1.0);
    std::uniform_real_distribution<double> uni(0.0, 1.0);

    double cur_e = E(pep), best_e = cur_e;
    std::vector<Vec3> best = pep;
    int n_accept = 0;
    {
        py::gil_scoped_release release;
        for (int step = 0; step < n_steps; ++step) {
            const double T = temp0 + (temp1 - temp0) * (n_steps > 1 ? double(step) / (n_steps - 1) : 1.0);
            Vec3 centroid{0, 0, 0};
            for (const auto& p : pep) { centroid.x += p.x; centroid.y += p.y; centroid.z += p.z; }
            const double inv = pep.empty() ? 0.0 : 1.0 / pep.size();
            centroid.x *= inv; centroid.y *= inv; centroid.z *= inv;

            std::vector<Vec3> trial = pep;
            double ax = gauss(rng), ay = gauss(rng), az = gauss(rng);
            const double norm = std::sqrt(ax * ax + ay * ay + az * az) + 1e-12;
            ax /= norm; ay /= norm; az /= norm;
            const double angle = gauss(rng) * rot_sigma;
            const Vec3 t{gauss(rng) * trans_sigma, gauss(rng) * trans_sigma, gauss(rng) * trans_sigma};
            rigid_move(trial, centroid, ax, ay, az, angle, t);

            const double trial_e = E(trial);
            const double dE = trial_e - cur_e;
            if (dE <= 0.0 || uni(rng) < std::exp(-dE / std::max(T, 1e-9))) {
                pep = std::move(trial);
                cur_e = trial_e;
                ++n_accept;
                if (cur_e < best_e) { best_e = cur_e; best = pep; }
            }
        }
    }

    py::array_t<double> out({static_cast<py::ssize_t>(best.size()), py::ssize_t(3)});
    auto w = out.mutable_unchecked<2>();
    for (size_t i = 0; i < best.size(); ++i) { w(i, 0) = best[i].x; w(i, 1) = best[i].y; w(i, 2) = best[i].z; }
    return py::make_tuple(out, best_e, n_accept);
}

}  // namespace

PYBIND11_MODULE(_refine, m) {
    m.doc() = "Potential-guided rigid-body peptide refinement (statistical potential + soft clash).";
    m.def("refine", &refine,
          py::arg("pep_xyz"), py::arg("pep_atom_res"), py::arg("pep_res_aa"),
          py::arg("par_xyz"), py::arg("par_atom_res"), py::arg("par_res_aa"),
          py::arg("potential"), py::arg("cutoff") = 5.0, py::arg("clash_d0") = 3.0,
          py::arg("clash_w") = 1.0, py::arg("restraint_w") = 1.0, py::arg("n_steps") = 2000,
          py::arg("trans_sigma") = 0.2, py::arg("rot_sigma") = 0.05, py::arg("temp0") = 1.0,
          py::arg("temp1") = 0.05, py::arg("seed") = 0,
          "Rigid-body Metropolis MC of the peptide; returns (best_xyz, best_energy, n_accept).");
    m.attr("__version__") = "0.1.0";
}
