// Potential-guided peptide refinement (knowledge-based, NOT physics MD).
//
// A rigid-body Metropolis Monte-Carlo local refinement of a peptide against fixed partner atoms
// (TCR + MHC), scored by the DOPE atom-level distance-dependent statistical potential (Shen & Sali,
// Protein Science 2006) plus a harmonic restraint to the input pose. DOPE is used here *only* for
// refinement, independently of the TCRen/MJ potentials tcren uses for epitope scoring, so the pose
// is not optimised against the same quantity it is later scored with. DOPE's short-range bins are
// strongly repulsive, so it supplies its own clash term -- there is no separate clash penalty.
//
// Energy(pose) =  Σ_{p∈peptide, q∈partner, d_pq < d_max}  φ_DOPE(c_p, c_q, d_pq)
//              +  w_r Σ_p ‖x_p − x_p^0‖²
// where c_p, c_q are MODELLER mean-force atom classes and φ_DOPE is the (linearly interpolated)
// tabulated potential. The Python wrapper pre-filters partner atoms to the interface shell, so the
// brute-force O(Np·Nq) energy is small; MC runs thousands of steps in well under a second.

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

// DOPE energy of a peptide pose against the fixed partner atoms, plus a harmonic restraint to the
// starting coordinates `pep0` (which keeps the refinement local -- DOPE alone could still drift the
// peptide to another favourable pocket).
double energy(const std::vector<Vec3>& pep, const std::vector<Vec3>& pep0,
              const std::vector<int>& pep_cls, const std::vector<Vec3>& par,
              const std::vector<int>& par_cls, const float* table,
              int n_cls, int n_bins, double x_start, double dx, double d_max, double restraint_w) {
    const double d_max2 = d_max * d_max;
    double e = 0.0;
    for (size_t a = 0; a < pep.size(); ++a) {
        e += restraint_w * dist2(pep[a], pep0[a]);
        const int cp = pep_cls[a];
        if (cp < 0) continue;
        const float* row = table + static_cast<size_t>(cp) * n_cls * n_bins;
        for (size_t b = 0; b < par.size(); ++b) {
            const int cq = par_cls[b];
            if (cq < 0) continue;
            const double d2 = dist2(pep[a], par[b]);
            if (d2 >= d_max2) continue;
            const float* knots = row + static_cast<size_t>(cq) * n_bins;
            const double t = (std::sqrt(d2) - x_start) / dx;
            if (t <= 0.0) {
                e += knots[0];                                   // short-range repulsive cap
            } else {
                const int k = static_cast<int>(t);
                if (k >= n_bins - 1) {
                    e += knots[n_bins - 1];
                } else {
                    const double f = t - k;
                    e += knots[k] * (1.0 - f) + knots[k + 1] * f;  // linear interpolation
                }
            }
        }
    }
    return e;
}

// Apply a rigid transform (rotation about the peptide centroid, then translation) in place.
void rigid_move(std::vector<Vec3>& pep, const Vec3& centroid,
                double ax, double ay, double az, double angle, const Vec3& t) {
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

py::tuple refine(py::array_t<double> pep_xyz, py::array_t<int> pep_class,
                 py::array_t<double> par_xyz, py::array_t<int> par_class,
                 py::array_t<float> dope_table, int n_cls, int n_bins,
                 double x_start, double dx, double restraint_w,
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
    std::vector<int> pep_cls = load_int(pep_class), par_cls = load_int(par_class);
    std::vector<float> table(dope_table.data(), dope_table.data() + dope_table.size());
    const double d_max = x_start + (n_bins - 1) * dx;  // beyond the last knot: no interaction

    const std::vector<Vec3> pep0 = pep;  // restraint reference (the input pose)
    auto E = [&](const std::vector<Vec3>& p) {
        return energy(p, pep0, pep_cls, par, par_cls, table.data(), n_cls, n_bins,
                      x_start, dx, d_max, restraint_w);
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
    m.doc() = "Potential-guided rigid-body peptide refinement (DOPE atom-level statistical potential).";
    m.def("refine", &refine,
          py::arg("pep_xyz"), py::arg("pep_class"), py::arg("par_xyz"), py::arg("par_class"),
          py::arg("dope_table"), py::arg("n_cls"), py::arg("n_bins"),
          py::arg("x_start"), py::arg("dx"), py::arg("restraint_w") = 0.5,
          py::arg("n_steps") = 2000, py::arg("trans_sigma") = 0.2, py::arg("rot_sigma") = 0.05,
          py::arg("temp0") = 1.0, py::arg("temp1") = 0.05, py::arg("seed") = 0,
          "Rigid-body Metropolis MC of the peptide; returns (best_xyz, best_energy, n_accept).");
    m.attr("__version__") = "0.2.0";
}
