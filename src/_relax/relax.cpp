// Interface energy + (later) sampling for native ΔΔG — the FlexPepDock/InterfaceAnalyzer replacement.
//
// repack: the rotamer packer — every side chain placed in the chi conformer DOPE likes best, with the
// Boltzmann weights over conformers so a caller can average instead of choose. This is what the native
// refiner was missing: `dope` returned 44 heavy atoms where the crystal peptide has 77, so it could not
// be compared with OpenMM/FlexPepDock on anything side-chain-sensitive (see refine/CPP_REWRITE.md).
//
// interface_energy: the DOPE atom-level statistical-potential interaction energy across an interface
// (peptide <-> partner heavy-atom pairs only — the peptide-internal and partner-internal terms are not
// summed, so this is already E_bound − E_separated for the cross terms, i.e. the interaction energy).
// A 20-line lift of energy() from src/_refine/refine.cpp with the harmonic restraint dropped. Stdlib
// only, pybind11, C++17 — the house pattern. The rotamer repack + flexible-backbone relax that turn
// this into a full ΔΔG land here later (relax_interface); this is the energy core they all call.

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <cmath>
#include <vector>

#ifndef M_PI  // MSVC does not define M_PI without _USE_MATH_DEFINES
#define M_PI 3.14159265358979323846
#endif

namespace py = pybind11;

namespace {

struct Vec3 {
    double x, y, z;
};
inline double dist2(const Vec3& a, const Vec3& b) {
    const double dx = a.x - b.x, dy = a.y - b.y, dz = a.z - b.z;
    return dx * dx + dy * dy + dz * dz;
}

// Sum the tabulated (linearly interpolated) DOPE potential over peptide<->partner heavy-atom pairs
// within range. `pep_class`/`par_class` are MODELLER mean-force atom classes (-1 = skip). Identical
// binning/interpolation to refine.cpp so the two agree.
double interface_energy(py::array_t<double> pep_xyz, py::array_t<int> pep_class,
                        py::array_t<double> par_xyz, py::array_t<int> par_class,
                        py::array_t<float> dope_table, int n_cls, int n_bins,
                        double x_start, double dx) {
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
    std::vector<int> pcl = load_int(pep_class), qcl = load_int(par_class);
    std::vector<float> table(dope_table.data(), dope_table.data() + dope_table.size());
    const double d_max = x_start + (n_bins - 1) * dx;
    const double d_max2 = d_max * d_max;

    double e = 0.0;
    py::gil_scoped_release release;
    for (size_t a = 0; a < pep.size(); ++a) {
        const int cp = pcl[a];
        if (cp < 0) continue;
        const float* row = table.data() + static_cast<size_t>(cp) * n_cls * n_bins;
        for (size_t b = 0; b < par.size(); ++b) {
            const int cq = qcl[b];
            if (cq < 0) continue;
            const double d2 = dist2(pep[a], par[b]);
            if (d2 >= d_max2) continue;
            const float* knots = row + static_cast<size_t>(cq) * n_bins;
            const double t = (std::sqrt(d2) - x_start) / dx;
            if (t <= 0.0) {
                e += knots[0];
            } else {
                const int k = static_cast<int>(t);
                if (k >= n_bins - 1) {
                    e += knots[n_bins - 1];
                } else {
                    const double f = t - k;
                    e += knots[k] * (1.0 - f) + knots[k + 1] * f;
                }
            }
        }
    }
    return e;
}

// Rotamer repack: place every side chain in the chi conformer DOPE likes best, and report the
// Boltzmann weights over the conformers so a caller can average rather than pick.
//
// Rotating every atom past C-beta about the CA-CB axis IS a chi1 change (deeper torsions ride along
// unchanged), and the same holds at each successive depth, so the geometry here is exact rather than
// interpolated. Which bonds rotate is decided in Python (tcren.rotamers.chi_axes) and handed over as
// a flat CSR-style description, because that needs residue chemistry and this needs speed.
//
// Torsions are enumerated as a product over a residue's chi angles: `n_steps` conformers per angle,
// `n_steps ^ n_chi` per residue. The first conformer of every residue is always the input one, so a
// repack can never return a pose worse than the one it was given.
//
// Mean field, deliberately: each residue is weighted against every other residue held at its INPUT
// conformation. Two side chains that would have to move together are not coupled, which is the
// approximation a dead-end-elimination packer removes and this one does not.
// The layout the caller hands over, all flat int arrays:
//   res_lo[r], res_hi[r]                -> this residue's atoms in `xyz`
//   res_chi_ptr[r],  res_chi_ptr[r+1]   -> its torsions in chi_a / chi_b / chi_mov_ptr
//   chi_a[i], chi_b[i]                  -> the two atoms defining torsion i's axis
//   moving[chi_mov_ptr[i] : ptr[i+1]]   -> the atoms torsion i rotates
//   env_atom[env_ptr[r] : ptr[r+1]]     -> the atoms residue r is scored against
//
// `xyz` holds EVERY heavy atom of the structure, not only the residues being repacked: a side chain
// has to be scored against the partner chains it packs against, and indexing the environment into a
// peptide-only array repacked it in vacuum.
inline void rotate_about(std::vector<Vec3>& xyz, const Vec3& a, const Vec3& b,
                         const int* mov, int n_mov, double angle) {
    const double ux = b.x - a.x, uy = b.y - a.y, uz = b.z - a.z;
    const double norm = std::sqrt(ux * ux + uy * uy + uz * uz);
    if (norm < 1e-9) return;
    const double kx = ux / norm, ky = uy / norm, kz = uz / norm;
    const double c = std::cos(angle), s = std::sin(angle), c1 = 1.0 - c;
    for (int m = 0; m < n_mov; ++m) {
        Vec3& p = xyz[mov[m]];
        const double vx = p.x - b.x, vy = p.y - b.y, vz = p.z - b.z;
        const double cx = ky * vz - kz * vy, cy = kz * vx - kx * vz, cz = kx * vy - ky * vx;
        const double dot = kx * vx + ky * vy + kz * vz;
        p.x = b.x + vx * c + cx * s + kx * dot * c1;
        p.y = b.y + vy * c + cy * s + ky * dot * c1;
        p.z = b.z + vz * c + cz * s + kz * dot * c1;
    }
}

// DOPE energy of one residue's atoms against a fixed environment, ignoring its own atoms.
inline double residue_energy(const std::vector<Vec3>& xyz, const std::vector<int>& cls,
                             int lo, int hi, const std::vector<Vec3>& env,
                             const std::vector<int>& env_cls, const std::vector<float>& table,
                             int n_cls, int n_bins, double x_start, double dx, double d_max2) {
    double e = 0.0;
    for (int a = lo; a < hi; ++a) {
        const int cp = cls[a];
        if (cp < 0) continue;
        const float* row = table.data() + static_cast<size_t>(cp) * n_cls * n_bins;
        for (size_t b = 0; b < env.size(); ++b) {
            const int cq = env_cls[b];
            if (cq < 0) continue;
            const double d2 = dist2(xyz[a], env[b]);
            if (d2 >= d_max2) continue;
            const float* knots = row + static_cast<size_t>(cq) * n_bins;
            const double t = (std::sqrt(d2) - x_start) / dx;
            if (t <= 0.0) {
                e += knots[0];
            } else {
                const int k = static_cast<int>(t);
                e += (k >= n_bins - 1) ? knots[n_bins - 1]
                                       : knots[k] * (1.0 - (t - k)) + knots[k + 1] * (t - k);
            }
        }
    }
    return e;
}

py::dict repack(py::array_t<double> xyz_in, py::array_t<int> atom_class,
                py::array_t<int> res_lo, py::array_t<int> res_hi, py::array_t<int> res_chi_ptr,
                py::array_t<int> chi_a, py::array_t<int> chi_b, py::array_t<int> chi_mov_ptr,
                py::array_t<int> moving, py::array_t<int> env_atom,
                py::array_t<int> env_ptr, py::array_t<float> dope_table, int n_cls, int n_bins,
                double x_start, double dx, int n_steps, double temperature) {
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

    std::vector<Vec3> xyz = load_xyz(xyz_in);
    std::vector<int> cls = load_int(atom_class);
    std::vector<int> rlo = load_int(res_lo), rhi = load_int(res_hi), rcp = load_int(res_chi_ptr);
    std::vector<int> ca = load_int(chi_a), cb = load_int(chi_b), cmp_ = load_int(chi_mov_ptr);
    std::vector<int> mov = load_int(moving);
    std::vector<int> env = load_int(env_atom), envp = load_int(env_ptr);
    std::vector<float> table(dope_table.data(), dope_table.data() + dope_table.size());

    const int n_res = static_cast<int>(rlo.size());
    const double d_max = x_start + (n_bins - 1) * dx;
    const double d_max2 = d_max * d_max;
    const double two_pi = 2.0 * M_PI;

    std::vector<double> best_e(n_res, 0.0);
    std::vector<int> n_conf(n_res, 1);
    std::vector<double> weights;          // ragged: n_conf[r] entries per residue, concatenated
    std::vector<int> weight_ptr(n_res + 1, 0);
    std::vector<Vec3> out = xyz;          // repacked coordinates (best conformer per residue)

    {
        py::gil_scoped_release release;
        std::vector<Vec3> work, env_xyz;
        std::vector<int> env_cls;
        std::vector<double> energies;
        for (int r = 0; r < n_res; ++r) {
            const int a_lo = rlo[r], a_hi = rhi[r];
            const int c_lo = rcp[r], c_hi = rcp[r + 1];
            const int n_chi = c_hi - c_lo;

            // The environment for this residue: its neighbours, at their input conformation.
            env_xyz.clear();
            env_cls.clear();
            for (int k = envp[r]; k < envp[r + 1]; ++k) {
                env_xyz.push_back(xyz[env[k]]);
                env_cls.push_back(cls[env[k]]);
            }

            int total = 1;
            for (int c = 0; c < n_chi; ++c) total *= n_steps;
            n_conf[r] = total;
            energies.assign(total, 0.0);

            int best = 0;
            double best_energy = 0.0;
            std::vector<Vec3> best_xyz(xyz.begin() + a_lo, xyz.begin() + a_hi);
            for (int t = 0; t < total; ++t) {
                work.assign(xyz.begin(), xyz.end());
                // Decode t as a mixed-radix index over the residue's chi angles; step 0 is the
                // input conformation, so t == 0 reproduces the input exactly.
                int rem = t;
                for (int c = 0; c < n_chi; ++c) {
                    const int step = rem % n_steps;
                    rem /= n_steps;
                    if (step == 0) continue;
                    const int ci = c_lo + c;
                    rotate_about(work, work[ca[ci]], work[cb[ci]], mov.data() + cmp_[ci],
                                 cmp_[ci + 1] - cmp_[ci], two_pi * step / n_steps);
                }
                const double e = residue_energy(work, cls, a_lo, a_hi, env_xyz, env_cls, table,
                                                n_cls, n_bins, x_start, dx, d_max2);
                energies[t] = e;
                if (t == 0 || e < best_energy) {
                    best_energy = e;
                    best = t;
                    best_xyz.assign(work.begin() + a_lo, work.begin() + a_hi);
                }
            }
            for (int a = a_lo; a < a_hi; ++a) out[a] = best_xyz[a - a_lo];
            best_e[r] = best_energy;

            // Boltzmann weights, shifted by the minimum so exp() cannot overflow.
            double z = 0.0;
            const double beta = 1.0 / (temperature > 1e-9 ? temperature : 1e-9);
            weight_ptr[r + 1] = weight_ptr[r];
            for (int t = 0; t < total; ++t) {
                const double w = std::exp(-(energies[t] - best_energy) * beta);
                weights.push_back(w);
                ++weight_ptr[r + 1];
                z += w;
            }
            for (int k = weight_ptr[r]; k < weight_ptr[r + 1]; ++k) weights[k] /= z;
        }
    }

    py::array_t<double> out_xyz(std::vector<py::ssize_t>{static_cast<py::ssize_t>(out.size()), 3});
    auto w = out_xyz.mutable_unchecked<2>();
    for (size_t i = 0; i < out.size(); ++i) {
        w(i, 0) = out[i].x;
        w(i, 1) = out[i].y;
        w(i, 2) = out[i].z;
    }
    py::dict res;
    res["xyz"] = out_xyz;
    res["energy"] = py::array_t<double>(best_e.size(), best_e.data());
    res["n_conformers"] = py::array_t<int>(n_conf.size(), n_conf.data());
    res["weights"] = py::array_t<double>(weights.size(), weights.data());
    res["weight_ptr"] = py::array_t<int>(weight_ptr.size(), weight_ptr.data());
    return res;
}

}  // namespace

PYBIND11_MODULE(_relax, m) {
    m.doc() = "Native interface energy + sampling for ΔΔG (DOPE interaction energy; repack/relax later).";
    m.def("interface_energy", &interface_energy,
          py::arg("pep_xyz"), py::arg("pep_class"), py::arg("par_xyz"), py::arg("par_class"),
          py::arg("dope_table"), py::arg("n_cls"), py::arg("n_bins"), py::arg("x_start"), py::arg("dx"),
          "DOPE interaction energy over peptide<->partner heavy-atom pairs (the interface ΔG core).");
    m.def("repack", &repack, py::arg("xyz"), py::arg("atom_class"), py::arg("res_lo"),
          py::arg("res_hi"), py::arg("res_chi_ptr"), py::arg("chi_a"), py::arg("chi_b"),
          py::arg("chi_mov_ptr"),
          py::arg("moving"), py::arg("env_atom"), py::arg("env_ptr"), py::arg("dope_table"),
          py::arg("n_cls"), py::arg("n_bins"), py::arg("x_start"), py::arg("dx"),
          py::arg("n_steps") = 3, py::arg("temperature") = 1.0,
          "Rotamer repack: best chi conformer per residue under DOPE, plus Boltzmann weights.");
    m.attr("__version__") = "0.2.0";
}
