// Interface energy + (later) sampling for native ΔΔG — the FlexPepDock/InterfaceAnalyzer replacement.
//
// relax_interface: flexible-backbone Metropolis MC of the peptide — how stable is the conformation
// it was handed? The readout is RMSF/drift under thermal sampling, not a better pose, and the
// intra-peptide term is switchable so its contribution to that stability can be measured directly.
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

#include <algorithm>
#include <atomic>
#include <cmath>
#include <limits>
#include <random>
#include <thread>
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

// Validate public native inputs while the GIL is held, before unchecked indexing.
void check_vector(const py::array_t<int>& a, py::ssize_t size = -1) {
    if (a.ndim() != 1 || a.size() > std::numeric_limits<int>::max() ||
        (size >= 0 && a.size() != size))
        throw py::value_error("invalid integer vector shape or length");
}
void check_range(const py::array_t<int>& a, py::ssize_t lo, py::ssize_t hi) {
    check_vector(a);
    auto v = a.unchecked<1>();
    for (py::ssize_t i = 0; i < a.size(); ++i)
        if (v(i) < lo || v(i) >= hi) throw py::value_error("atom index or class out of range");
}
void check_xyz(const py::array_t<double>& a, const py::array_t<int>& cls, int n_cls) {
    if (a.ndim() != 2 || a.shape(1) != 3 || a.shape(0) > std::numeric_limits<int>::max())
        throw py::value_error("coordinates must have shape (N, 3)");
    check_vector(cls, a.shape(0));
    check_range(cls, -1, n_cls);
    auto v = a.unchecked<2>();
    for (py::ssize_t i = 0; i < a.shape(0); ++i)
        for (int j = 0; j < 3; ++j)
            if (!std::isfinite(v(i, j))) throw py::value_error("coordinates must be finite");
}
void check_table(const py::array_t<float>& a, int nc, int nb, double x0, double dx) {
    if (nc <= 0 || nb <= 0 || !std::isfinite(x0) || x0 < 0 ||
        !std::isfinite(dx) || dx <= 0 || !std::isfinite(x0 + (nb - 1) * dx) ||
        static_cast<size_t>(nc) > std::numeric_limits<size_t>::max() / nc / nb ||
        !(a.flags() & py::array::c_style) ||
        (a.ndim() != 1 && a.ndim() != 3) ||
        static_cast<size_t>(a.size()) != static_cast<size_t>(nc) * nc * nb)
        throw py::value_error("invalid DOPE table dimensions, layout or distance grid");
    if (a.ndim() == 3 && (a.shape(0) != nc || a.shape(1) != nc || a.shape(2) != nb))
        throw py::value_error("invalid DOPE table shape");
    for (py::ssize_t i = 0; i < a.size(); ++i)
        if (!std::isfinite(a.data()[i])) throw py::value_error("DOPE table must be finite");
}
void check_csr(const py::array_t<int>& a, py::ssize_t rows, py::ssize_t size) {
    check_vector(a, rows + 1);
    auto v = a.unchecked<1>();
    if (v(0) != 0 || v(rows) != size) throw py::value_error("invalid CSR endpoints");
    for (py::ssize_t i = 1; i <= rows; ++i)
        if (v(i) < v(i - 1)) throw py::value_error("CSR pointers must be monotone");
}

// Sum the tabulated (linearly interpolated) DOPE potential over peptide<->partner heavy-atom pairs
// within range. `pep_class`/`par_class` are MODELLER mean-force atom classes (-1 = skip). Identical
// binning/interpolation to refine.cpp so the two agree.
double interface_energy(py::array_t<double> pep_xyz, py::array_t<int> pep_class,
                        py::array_t<double> par_xyz, py::array_t<int> par_class,
                        py::array_t<float> dope_table, int n_cls, int n_bins,
                        double x_start, double dx) {
    check_table(dope_table, n_cls, n_bins, x_start, dx);
    check_xyz(pep_xyz, pep_class, n_cls);
    check_xyz(par_xyz, par_class, n_cls);
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
// `n_steps ^ n_chi` per residue. Including the input conformer guarantees each chosen residue's
// energy cannot worsen against its fixed input environment. Combining those choices does not
// guarantee a non-increasing total energy or absence of new clashes between moved residues.
//
// Mean field, deliberately: each residue is weighted against every other residue held at its INPUT
// conformation. Two side chains that would have to move together are not coupled, which is the
// approximation a dead-end-elimination packer removes and this one does not.
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
py::dict repack(py::array_t<double> xyz_in, py::array_t<int> atom_class,
                py::array_t<int> res_lo, py::array_t<int> res_hi, py::array_t<int> res_chi_ptr,
                py::array_t<int> chi_a, py::array_t<int> chi_b, py::array_t<int> chi_mov_ptr,
                py::array_t<int> moving, py::array_t<int> env_atom,
                py::array_t<int> env_ptr, py::array_t<float> dope_table, int n_cls, int n_bins,
                double x_start, double dx, int n_steps, double temperature) {
    check_table(dope_table, n_cls, n_bins, x_start, dx);
    check_xyz(xyz_in, atom_class, n_cls);
    if (n_steps <= 0 || !std::isfinite(temperature) || temperature <= 0)
        throw py::value_error("repack steps and temperature must be positive");
    const auto n = xyz_in.shape(0);
    check_vector(res_lo); check_vector(res_hi, res_lo.size());
    check_range(res_lo, 0, n + 1); check_range(res_hi, 0, n + 1);
    check_range(chi_a, 0, n); check_vector(chi_b, chi_a.size()); check_range(chi_b, 0, n);
    check_range(moving, 0, n); check_range(env_atom, 0, n);
    check_csr(res_chi_ptr, res_lo.size(), chi_a.size());
    check_csr(chi_mov_ptr, chi_a.size(), moving.size());
    check_csr(env_ptr, res_lo.size(), env_atom.size());
    auto lo = res_lo.unchecked<1>(), hi = res_hi.unchecked<1>();
    auto ptr = res_chi_ptr.unchecked<1>();
    int count = 0;
    for (py::ssize_t r = 0; r < res_lo.size(); ++r) {
        if (lo(r) > hi(r)) throw py::value_error("residue range is reversed");
        int total = 1;
        for (int c = ptr(r); c < ptr(r + 1); ++c) {
            if (total > std::numeric_limits<int>::max() / n_steps)
                throw py::value_error("rotamer grid exceeds integer capacity");
            total *= n_steps;
        }
        if (count > std::numeric_limits<int>::max() - total)
            throw py::value_error("rotamer weights exceed integer capacity");
        count += total;
    }
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

// Flexible-backbone Metropolis MC of the peptide — how stable is the conformation it was given?
//
// The question this exists to answer is not "what is the best pose" but "does the pose hold". A
// contact potential scores whichever conformation it is handed; it cannot tell a peptide that its
// own side chains hold in the TCR-facing conformation from one that merely happens to be drawn that
// way. Sampling backbone torsions at temperature and measuring how far the peptide wanders does.
//
// Moves are torsional: perturb one backbone phi or psi and rotate everything downstream of it. On a
// peptide both of whose termini are free that is exact and needs no loop closure — the same
// rotate_about primitive the packer uses, applied to N-CA and CA-C instead of CA-CB.
//
// Energy(pose) = DOPE(peptide <-> partner)                     the interface
//              + intra_w * DOPE(peptide <-> peptide, |i-j| >= min_sep)   the peptide with itself
//              + anchor_w * sum ||x - x0||^2  over the anchor atoms      it stays in the groove
//
// `intra_w` is the switch the hypothesis needs. Sewell (2026-08): intra-peptide interactions
// stabilise the productive bulge, and a poor binder can make many contacts yet fail to hold that
// conformation — so an additive contact model should fail exactly where the intra term matters.
// Running the same MC with intra_w = 1 and intra_w = 0 measures that directly.
py::dict relax_interface(py::array_t<double> xyz_in, py::array_t<int> atom_class,
                         int pep_lo, int pep_hi, py::array_t<int> pep_res,
                         py::array_t<int> tor_a, py::array_t<int> tor_b,
                         py::array_t<int> tor_mov_ptr, py::array_t<int> tor_mov,
                         py::array_t<int> par_atom, py::array_t<int> anchor_atom,
                         py::array_t<float> dope_table, int n_cls, int n_bins,
                         double x_start, double dx, int n_steps, double temperature,
                         double sigma_deg, double intra_w, double anchor_w, int min_sep,
                         int burn_in, int seed) {
    check_table(dope_table, n_cls, n_bins, x_start, dx);
    check_xyz(xyz_in, atom_class, n_cls);
    const auto n = xyz_in.shape(0);
    if (pep_lo < 0 || pep_hi <= pep_lo || pep_hi > n || n_steps <= 0 || burn_in < 0 ||
        min_sep < 0 || !std::isfinite(temperature) || temperature <= 0 ||
        !std::isfinite(sigma_deg) || sigma_deg <= 0 ||
        !std::isfinite(intra_w) || !std::isfinite(anchor_w))
        throw py::value_error("invalid peptide range or relaxation parameters");
    check_vector(pep_res, pep_hi - pep_lo);
    check_range(pep_res, 0, std::numeric_limits<int>::max());
    check_range(tor_a, 0, n); check_vector(tor_b, tor_a.size()); check_range(tor_b, 0, n);
    check_range(tor_mov, 0, n); check_range(par_atom, 0, n);
    check_range(anchor_atom, pep_lo, pep_hi);
    check_csr(tor_mov_ptr, tor_a.size(), tor_mov.size());
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
    std::vector<int> cls = load_int(atom_class), pres = load_int(pep_res);
    std::vector<int> ta = load_int(tor_a), tb = load_int(tor_b), tmp = load_int(tor_mov_ptr);
    std::vector<int> tmov = load_int(tor_mov), par = load_int(par_atom);
    std::vector<int> anc = load_int(anchor_atom);
    std::vector<float> table(dope_table.data(), dope_table.data() + dope_table.size());

    const int n_pep = pep_hi - pep_lo;
    const int n_tor = static_cast<int>(ta.size());
    const double d_max = x_start + (n_bins - 1) * dx;
    const double d_max2 = d_max * d_max;
    const std::vector<Vec3> start(xyz.begin() + pep_lo, xyz.begin() + pep_hi);

    // Tabulated pair term, used for both the interface and the intra-peptide sum.
    auto pair_e = [&](const Vec3& p, int cp, const Vec3& q, int cq) -> double {
        if (cp < 0 || cq < 0) return 0.0;
        const double d2 = dist2(p, q);
        if (d2 >= d_max2) return 0.0;
        const float* knots = table.data() + (static_cast<size_t>(cp) * n_cls + cq) * n_bins;
        const double t = (std::sqrt(d2) - x_start) / dx;
        if (t <= 0.0) return knots[0];
        const int k = static_cast<int>(t);
        if (k >= n_bins - 1) return knots[n_bins - 1];
        return knots[k] * (1.0 - (t - k)) + knots[k + 1] * (t - k);
    };

    auto total_energy = [&](const std::vector<Vec3>& c) -> double {
        double e = 0.0;
        for (int a = pep_lo; a < pep_hi; ++a) {
            for (size_t b = 0; b < par.size(); ++b) e += pair_e(c[a], cls[a], c[par[b]], cls[par[b]]);
        }
        if (intra_w != 0.0) {
            // Sequence-local pairs are in contact by covalent geometry, not by folding, so they
            // carry no information about whether the conformation is held together.
            for (int a = pep_lo; a < pep_hi; ++a) {
                for (int b = a + 1; b < pep_hi; ++b) {
                    if (std::abs(pres[a - pep_lo] - pres[b - pep_lo]) < min_sep) continue;
                    e += intra_w * pair_e(c[a], cls[a], c[b], cls[b]);
                }
            }
        }
        if (anchor_w != 0.0) {
            for (size_t k = 0; k < anc.size(); ++k) {
                const int a = anc[k];
                e += anchor_w * dist2(c[a], start[a - pep_lo]);
            }
        }
        return e;
    };

    std::mt19937 rng(static_cast<unsigned>(seed));
    std::uniform_real_distribution<double> uni(0.0, 1.0);
    std::normal_distribution<double> gauss(0.0, sigma_deg * M_PI / 180.0);
    std::uniform_int_distribution<int> pick(0, std::max(n_tor - 1, 0));

    const double e_start = total_energy(xyz);
    double e_cur = e_start, e_best = e_cur;
    std::vector<Vec3> best(xyz.begin() + pep_lo, xyz.begin() + pep_hi);
    // Running mean and sum-of-squares per peptide atom, over the post-burn-in samples: RMSF is the
    // spread of the ensemble, drift the distance of its mean from where it started.
    std::vector<Vec3> mean(n_pep, {0.0, 0.0, 0.0});
    std::vector<double> sq(n_pep, 0.0);
    long n_sample = 0, n_accept = 0;

    if (n_tor > 0) {
        py::gil_scoped_release release;
        std::vector<Vec3> work;
        for (int step = 0; step < n_steps; ++step) {
            work = xyz;
            const int t = pick(rng);
            rotate_about(work, work[ta[t]], work[tb[t]], tmov.data() + tmp[t],
                         tmp[t + 1] - tmp[t], gauss(rng));
            const double e_new = total_energy(work);
            if (e_new <= e_cur || uni(rng) < std::exp(-(e_new - e_cur) / temperature)) {
                xyz.swap(work);
                e_cur = e_new;
                ++n_accept;
                if (e_cur < e_best) {
                    e_best = e_cur;
                    best.assign(xyz.begin() + pep_lo, xyz.begin() + pep_hi);
                }
            }
            if (step >= burn_in) {
                ++n_sample;
                for (int i = 0; i < n_pep; ++i) {
                    const Vec3& p = xyz[pep_lo + i];
                    mean[i].x += p.x;
                    mean[i].y += p.y;
                    mean[i].z += p.z;
                    sq[i] += p.x * p.x + p.y * p.y + p.z * p.z;
                }
            }
        }
    }

    double rmsf = 0.0, drift = 0.0;
    if (n_sample > 0) {
        for (int i = 0; i < n_pep; ++i) {
            const double mx = mean[i].x / n_sample, my = mean[i].y / n_sample,
                         mz = mean[i].z / n_sample;
            // Var(x)+Var(y)+Var(z) = <r.r> - <r>.<r>; clamped because catastrophic cancellation can
            // push a genuinely rigid atom microscopically below zero.
            rmsf += std::max(0.0, sq[i] / n_sample - (mx * mx + my * my + mz * mz));
            const double dxm = mx - start[i].x, dym = my - start[i].y, dzm = mz - start[i].z;
            drift += dxm * dxm + dym * dym + dzm * dzm;
        }
        rmsf = std::sqrt(rmsf / n_pep);
        drift = std::sqrt(drift / n_pep);
    }

    py::array_t<double> out_xyz(std::vector<py::ssize_t>{n_pep, 3});
    auto w = out_xyz.mutable_unchecked<2>();
    for (int i = 0; i < n_pep; ++i) {
        w(i, 0) = best[i].x;
        w(i, 1) = best[i].y;
        w(i, 2) = best[i].z;
    }
    py::dict res;
    res["xyz"] = out_xyz;
    res["energy"] = e_best;
    res["energy_start"] = e_start;
    res["energy_final"] = e_cur;
    res["rmsf"] = rmsf;
    res["drift"] = drift;
    res["accept_rate"] = n_tor > 0 ? static_cast<double>(n_accept) / n_steps : 0.0;
    res["n_samples"] = static_cast<long>(n_sample);
    return res;
}


// ---------------------------------------------------------------------------------------------
// substitution_scan: the amino-acid background taken in 3D rather than by swapping a letter.
//
// `smoothed_reference` sums phi_i(a) = sum_j sigma_ij e(a, y_j) over a FROZEN contact map, so a
// contact an arginine made is still charged to the alanine that replaces it. This builds every
// candidate residue at every rotamer of a chi grid, and reports three outcomes per (position,
// candidate, partner) instead of one: no contact, contact with a probability, or -- when every
// rotamer of that candidate collides with the environment -- an inadmissible substitution that must
// leave the background average rather than sit at its least-bad rotamer.
//
// One call per structure. The loop is positions x 20 candidates x 3^n_chi rotamers, each needing a
// DOPE energy against the environment and a heavy-atom distance against every partner, so a Python
// pool over residues is the failure mode to avoid: the GIL is released and the positions are threaded
// here instead.
// ---------------------------------------------------------------------------------------------

// Natural extension reference frame: place an atom from three placed ones plus (bond, angle, torsion).
// The exact inverse of tcren.refine.sidechain.dihedral -- getting the sign backwards leaves every
// bond and angle right and every atom 180 degrees out, which reads as a plausible RMSD, not an error.
inline Vec3 nerf(const Vec3& a, const Vec3& b, const Vec3& c,
                 double bond, double angle, double torsion) {
    double bx = c.x - b.x, by = c.y - b.y, bz = c.z - b.z;
    const double bn = std::sqrt(bx * bx + by * by + bz * bz);
    bx /= bn; by /= bn; bz /= bn;
    const double ax = b.x - a.x, ay = b.y - a.y, az = b.z - a.z;
    double nx = ay * bz - az * by, ny = az * bx - ax * bz, nz = ax * by - ay * bx;
    const double nn = std::sqrt(nx * nx + ny * ny + nz * nz);
    nx /= nn; ny /= nn; nz /= nn;
    const double mx = ny * bz - nz * by, my = nz * bx - nx * bz, mz = nx * by - ny * bx;
    const double d0 = -bond * std::cos(angle);
    const double d1 = bond * std::sin(angle) * std::cos(torsion);
    const double d2 = bond * std::sin(angle) * std::sin(torsion);
    return {c.x + d0 * bx + d1 * mx + d2 * nx,
            c.y + d0 * by + d1 * my + d2 * ny,
            c.z + d0 * bz + d1 * mz + d2 * nz};
}

struct ScanIn {
    const double* bb;          // (P, 4, 3) N, CA, C, CB per scan position
    const int* t_a;            // (20, A) slot of the torsion's first reference atom
    const int* t_b;
    const int* t_c;
    const int* t_chi;          // (20, A) 0 = fixed torsion, else the 1-based chi index
    const double* t_bond;      // (20, A)
    const double* t_ang;
    const double* t_off;
    const int* t_cls;          // (20, A) DOPE atom class of the built atom
    const double* t_rad;       // (20, A) van der Waals radius, for the clash test
    const int* n_atom;         // (20)
    const int* n_chi;          // (20)
    const double* chi_vals;    // (G) the chi grid, shared by every chi
    int n_grid;
    const double* env_xyz;     // (E, 3)
    const int* env_cls;        // (E)
    const double* env_rad;     // (E)
    const int* env_res;        // (E) scan position this atom belongs to, or -1
    const int* env_ptr;        // (P+1) CSR into env_idx: the atoms near each scan position
    const int* env_idx;        // (nnz)
    int n_env;
    const double* par_xyz;     // (M, 3)
    const int* par_ptr;        // (n_par + 1) CSR: par_xyz rows belonging to each partner residue
    int n_par_atom;
    const int* pair_i;         // (K) scan position
    const int* pair_j;         // (K) partner residue
    int n_pair;
    const float* table;        // DOPE
    int n_cls, n_bins;
    double x_start, dx;
    double cutoff, clash_tol, temperature;
    int max_atom;
};

// One (position, candidate) cell: build every rotamer, weight it, and accumulate the contact
// probability against each partner residue this position is paired with.
void scan_cell(const ScanIn& in, int pos, int aa, double* p_out, unsigned char* adm,
               int* n_clash_out, int* n_rot_out,
               const std::vector<std::vector<int>>& pairs_of_pos) {
    const int na = in.n_atom[aa];
    const int nc = in.n_chi[aa];
    const int A = in.max_atom;
    int n_rot = 1;
    for (int k = 0; k < nc; ++k) n_rot *= in.n_grid;
    n_rot_out[pos * 20 + aa] = n_rot;
    if (na == 0) {                      // alanine and glycine: nothing past CB, nothing to place
        adm[pos * 20 + aa] = 1;
        return;
    }

    const Vec3 N = {in.bb[(pos * 4 + 0) * 3 + 0], in.bb[(pos * 4 + 0) * 3 + 1], in.bb[(pos * 4 + 0) * 3 + 2]};
    const Vec3 CA = {in.bb[(pos * 4 + 1) * 3 + 0], in.bb[(pos * 4 + 1) * 3 + 1], in.bb[(pos * 4 + 1) * 3 + 2]};
    const Vec3 C = {in.bb[(pos * 4 + 2) * 3 + 0], in.bb[(pos * 4 + 2) * 3 + 1], in.bb[(pos * 4 + 2) * 3 + 2]};
    const Vec3 CB = {in.bb[(pos * 4 + 3) * 3 + 0], in.bb[(pos * 4 + 3) * 3 + 1], in.bb[(pos * 4 + 3) * 3 + 2]};

    const double d_max = in.x_start + (in.n_bins - 1) * in.dx;
    const double d_max2 = d_max * d_max;
    const double cut2 = in.cutoff * in.cutoff;

    std::vector<Vec3> slot(4 + na);
    std::vector<double> chi(4, 0.0);
    std::vector<double> energy(n_rot, 0.0);
    std::vector<unsigned char> clashed(n_rot, 0);
    // hit[r * n_pair_here + q] -- did rotamer r put this candidate within cutoff of partner q
    const std::vector<int>& mypairs = pairs_of_pos[pos];
    const int np = static_cast<int>(mypairs.size());
    std::vector<unsigned char> hit(static_cast<size_t>(n_rot) * (np > 0 ? np : 1), 0);

    for (int r = 0; r < n_rot; ++r) {
        int q = r;
        for (int k = 0; k < nc; ++k) { chi[k] = in.chi_vals[q % in.n_grid]; q /= in.n_grid; }
        slot[0] = N; slot[1] = CA; slot[2] = C; slot[3] = CB;
        for (int t = 0; t < na; ++t) {
            const int idx = aa * A + t;
            const double tor = in.t_off[idx] + (in.t_chi[idx] ? chi[in.t_chi[idx] - 1] : 0.0);
            slot[4 + t] = nerf(slot[in.t_a[idx]], slot[in.t_b[idx]], slot[in.t_c[idx]],
                               in.t_bond[idx], in.t_ang[idx], tor);
        }
        // Environment pass: DOPE energy and the clash test in one sweep. The residue's own atoms
        // are skipped -- the built geometry is ideal by construction, so an intra-residue overlap
        // is a property of the table and not of the pose.
        double e = 0.0;
        bool bad = false;
        for (int t = 0; t < na && !bad; ++t) {
            const int idx = aa * A + t;
            const int cp = in.t_cls[idx];
            const double rad = in.t_rad[idx];
            const Vec3& x = slot[4 + t];
            const float* row = (cp >= 0) ? in.table + static_cast<size_t>(cp) * in.n_cls * in.n_bins
                                         : nullptr;
            for (int e_i = in.env_ptr[pos]; e_i < in.env_ptr[pos + 1]; ++e_i) {
                const int b = in.env_idx[e_i];
                if (in.env_res[b] == pos) continue;
                const Vec3 y = {in.env_xyz[b * 3], in.env_xyz[b * 3 + 1], in.env_xyz[b * 3 + 2]};
                const double d2 = dist2(x, y);
                if (d2 >= d_max2) continue;
                const double d = std::sqrt(d2);
                // A clashed rotamer is discarded whole -- its energy is never weighted and its
                // contacts are never counted -- so stop the moment one overlap is found rather
                // than finishing the sweep to compute a number nothing reads.
                if (rad + in.env_rad[b] - d > in.clash_tol) { bad = true; break; }
                if (row && in.env_cls[b] >= 0) {
                    int bin = static_cast<int>((d - in.x_start) / in.dx);
                    if (bin < 0) bin = 0;
                    if (bin >= in.n_bins) continue;
                    e += row[static_cast<size_t>(in.env_cls[b]) * in.n_bins + bin];
                }
            }
        }
        energy[r] = e;
        clashed[r] = bad ? 1 : 0;
        if (bad || np == 0) continue;
        // Partner pass: closest heavy-atom distance from this rotamer to each paired partner
        // residue. CB counts -- it is part of the side chain the candidate presents.
        for (int qi = 0; qi < np; ++qi) {
            const int jres = in.pair_j[mypairs[qi]];
            bool touch = false;
            for (int m = in.par_ptr[jres]; m < in.par_ptr[jres + 1] && !touch; ++m) {
                const Vec3 y = {in.par_xyz[m * 3], in.par_xyz[m * 3 + 1], in.par_xyz[m * 3 + 2]};
                if (dist2(CB, y) < cut2) { touch = true; break; }
                for (int t = 0; t < na; ++t)
                    if (dist2(slot[4 + t], y) < cut2) { touch = true; break; }
            }
            hit[static_cast<size_t>(r) * np + qi] = touch ? 1 : 0;
        }
    }

    int n_bad = 0;
    for (int r = 0; r < n_rot; ++r) n_bad += clashed[r];
    n_clash_out[pos * 20 + aa] = n_bad;
    const bool admissible = n_bad < n_rot;
    adm[pos * 20 + aa] = admissible ? 1 : 0;
    if (!admissible || np == 0) return;

    // Boltzmann weights over the SURVIVING rotamers only. A candidate that fits nowhere has already
    // left; one that fits somewhere is averaged over the places it fits, not over the places it does not.
    double emin = 0.0;
    bool first = true;
    for (int r = 0; r < n_rot; ++r) {
        if (clashed[r]) continue;
        if (first || energy[r] < emin) { emin = energy[r]; first = false; }
    }
    double zsum = 0.0;
    std::vector<double> w(n_rot, 0.0);
    for (int r = 0; r < n_rot; ++r) {
        if (clashed[r]) continue;
        w[r] = std::exp(-(energy[r] - emin) / std::max(in.temperature, 1e-6));
        zsum += w[r];
    }
    if (zsum <= 0.0) return;
    for (int qi = 0; qi < np; ++qi) {
        double p = 0.0;
        for (int r = 0; r < n_rot; ++r)
            if (!clashed[r] && hit[static_cast<size_t>(r) * np + qi]) p += w[r];
        p_out[static_cast<size_t>(mypairs[qi]) * 20 + aa] = p / zsum;  // row = pair index
    }
}

// pybind entry point. Everything is marshalled flat by tcren.energetics.background, which owns the
// packing and the neighbour lists; this owns the arithmetic and the threads.
py::dict substitution_scan(
        py::array_t<double> bb,
        py::array_t<int> t_a, py::array_t<int> t_b, py::array_t<int> t_c, py::array_t<int> t_chi,
        py::array_t<double> t_bond, py::array_t<double> t_ang, py::array_t<double> t_off,
        py::array_t<int> t_cls, py::array_t<double> t_rad,
        py::array_t<int> n_atom, py::array_t<int> n_chi, py::array_t<double> chi_vals,
        py::array_t<double> env_xyz, py::array_t<int> env_cls, py::array_t<double> env_rad,
        py::array_t<int> env_res, py::array_t<int> env_ptr, py::array_t<int> env_idx,
        py::array_t<double> par_xyz, py::array_t<int> par_ptr,
        py::array_t<int> pair_i, py::array_t<int> pair_j,
        py::array_t<float> dope_table, int n_cls, int n_bins, double x_start, double dx,
        double cutoff, double clash_tol, double temperature, int n_threads) {
    // Validate the Python trust boundary before any raw-pointer indexing or GIL release.
    auto shape = [](const py::array& a, std::initializer_list<py::ssize_t> dims) {
        if (!(a.flags() & py::array::c_style) || a.ndim() != static_cast<int>(dims.size()))
            throw py::value_error("scan arrays must be C-contiguous and have the expected rank");
        int k = 0;
        for (auto d : dims) {
            if ((d >= 0 && a.shape(k) != d) || a.shape(k) > std::numeric_limits<int>::max() / 20)
                throw py::value_error("scan array shape is invalid or too large");
            ++k;
        }
    };
    shape(bb, {-1, 4, 3}); shape(t_a, {20, -1});
    const auto width = t_a.shape(1);
    for (const auto& a : {t_b, t_c, t_chi, t_cls}) shape(a, {20, width});
    for (const auto& a : {t_bond, t_ang, t_off, t_rad}) shape(a, {20, width});
    shape(n_atom, {20}); shape(n_chi, {20}); shape(chi_vals, {-1});
    shape(env_xyz, {-1, 3}); shape(par_xyz, {-1, 3});
    const auto ne = env_xyz.shape(0);
    shape(env_cls, {ne}); shape(env_rad, {ne}); shape(env_res, {ne});
    shape(env_ptr, {bb.shape(0) + 1}); shape(env_idx, {-1});
    shape(par_ptr, {-1}); shape(pair_i, {-1}); shape(pair_j, {pair_i.shape(0)});
    if (n_cls <= 0 || n_bins <= 0 || n_threads < 0 || chi_vals.size() == 0 ||
        !std::isfinite(x_start) || x_start < 0 || !std::isfinite(dx) || dx <= 0 ||
        !std::isfinite(cutoff) || cutoff <= 0 || !std::isfinite(clash_tol) || clash_tol < 0 ||
        !std::isfinite(temperature) || temperature <= 0)
        throw py::value_error("invalid scan dimensions, grid, thread count, or scalar parameters");
    // The Python packer supplies a flattened table; accept its equivalent tensor as well.
    if (!(dope_table.flags() & py::array::c_style) ||
        (dope_table.ndim() != 1 && dope_table.ndim() != 3) ||
        static_cast<size_t>(n_cls) > std::numeric_limits<size_t>::max() / n_cls / n_bins ||
        static_cast<size_t>(dope_table.size()) != static_cast<size_t>(n_cls) * n_cls * n_bins)
        throw py::value_error("invalid DOPE table size or layout");
    if (dope_table.ndim() == 3) shape(dope_table, {n_cls, n_cls, n_bins});
    auto finite = [](const auto& a) {
        for (py::ssize_t i = 0; i < a.size(); ++i)
            if (!std::isfinite(a.data()[i])) throw py::value_error("scan arrays must be finite");
    };
    for (const auto& a : {bb, t_bond, t_ang, t_off, t_rad, chi_vals, env_xyz, env_rad, par_xyz}) finite(a);
    finite(dope_table);
    auto range = [](const auto& a, int lo, py::ssize_t hi) {
        for (py::ssize_t i = 0; i < a.size(); ++i)
            if (a.data()[i] < lo || a.data()[i] >= hi)
                throw py::value_error("scan index or topology value is out of range");
    };
    range(n_atom, 0, width + 1); range(n_chi, 0, 5);
    range(env_cls, -1, n_cls); range(env_res, -1, bb.shape(0));
    range(env_idx, 0, ne); range(pair_i, 0, bb.shape(0));
    auto csr = [](const auto& a, py::ssize_t size) {
        if (a.size() == 0 || a.data()[0] != 0 || a.data()[a.size() - 1] != size)
            throw py::value_error("invalid scan CSR endpoints");
        for (py::ssize_t i = 1; i < a.size(); ++i)
            if (a.data()[i] < a.data()[i - 1]) throw py::value_error("scan CSR must be monotone");
    };
    csr(env_ptr, env_idx.size()); csr(par_ptr, par_xyz.shape(0));
    range(pair_j, 0, par_ptr.size() - 1);
    for (int aa = 0; aa < 20; ++aa) {
        int nr = 1;
        for (int k = 0; k < n_chi.data()[aa]; ++k) {
            if (nr > std::numeric_limits<int>::max() / chi_vals.size())
                throw py::value_error("rotamer grid exceeds integer capacity");
            nr *= static_cast<int>(chi_vals.size());
        }
        for (int t = 0; t < n_atom.data()[aa]; ++t) {
            const auto i = aa * width + t;
            for (const auto& a : {t_a, t_b, t_c})
                if (a.data()[i] < 0 || a.data()[i] >= 4 + t)
                    throw py::value_error("topology must reference an already placed atom");
            if (t_chi.data()[i] < 0 || t_chi.data()[i] > n_chi.data()[aa] ||
                t_cls.data()[i] < -1 || t_cls.data()[i] >= n_cls ||
                t_bond.data()[i] <= 0 || t_rad.data()[i] < 0)
                throw py::value_error("invalid active topology parameters");
        }
    }
    const int P = static_cast<int>(bb.shape(0));
    const int K = static_cast<int>(pair_i.shape(0));
    const int A = static_cast<int>(t_a.shape(1));

    ScanIn in{};
    in.bb = bb.data(); in.t_a = t_a.data(); in.t_b = t_b.data(); in.t_c = t_c.data();
    in.t_chi = t_chi.data(); in.t_bond = t_bond.data(); in.t_ang = t_ang.data();
    in.t_off = t_off.data(); in.t_cls = t_cls.data(); in.t_rad = t_rad.data();
    in.n_atom = n_atom.data(); in.n_chi = n_chi.data();
    in.chi_vals = chi_vals.data(); in.n_grid = static_cast<int>(chi_vals.shape(0));
    in.env_xyz = env_xyz.data(); in.env_cls = env_cls.data(); in.env_rad = env_rad.data();
    in.env_res = env_res.data(); in.env_ptr = env_ptr.data(); in.env_idx = env_idx.data();
    in.n_env = static_cast<int>(env_cls.shape(0));
    in.par_xyz = par_xyz.data(); in.par_ptr = par_ptr.data();
    in.n_par_atom = static_cast<int>(par_xyz.shape(0));
    in.pair_i = pair_i.data(); in.pair_j = pair_j.data(); in.n_pair = K;
    in.table = dope_table.data(); in.n_cls = n_cls; in.n_bins = n_bins;
    in.x_start = x_start; in.dx = dx;
    in.cutoff = cutoff; in.clash_tol = clash_tol; in.temperature = temperature;
    in.max_atom = A;

    // Which pairs each position owns, so a cell writes only its own output rows.
    std::vector<std::vector<int>> pairs_of_pos(P);
    for (int k = 0; k < K; ++k) pairs_of_pos[in.pair_i[k]].push_back(k);

    auto p_out = py::array_t<double>({static_cast<py::ssize_t>(K), static_cast<py::ssize_t>(20)});
    auto adm = py::array_t<unsigned char>({static_cast<py::ssize_t>(P), static_cast<py::ssize_t>(20)});
    auto n_clash = py::array_t<int>({static_cast<py::ssize_t>(P), static_cast<py::ssize_t>(20)});
    auto n_rot = py::array_t<int>({static_cast<py::ssize_t>(P), static_cast<py::ssize_t>(20)});
    double* pp = p_out.mutable_data();
    unsigned char* ap = adm.mutable_data();
    int* cp = n_clash.mutable_data();
    int* rp = n_rot.mutable_data();
    std::fill(pp, pp + static_cast<size_t>(K) * 20, 0.0);
    std::fill(ap, ap + static_cast<size_t>(P) * 20, 0);
    std::fill(cp, cp + static_cast<size_t>(P) * 20, 0);
    std::fill(rp, rp + static_cast<size_t>(P) * 20, 0);

    {
        py::gil_scoped_release release;
        int nt = n_threads > 0 ? n_threads : static_cast<int>(std::thread::hardware_concurrency());
        if (nt < 1) nt = 1;
        // Thread over the FLATTENED (position, candidate) cells, not over positions. On the
        // peptide side there are only nine positions, so a per-position split caps the pool at
        // nine however many cores are free; the cells are 20x more numerous and equally
        // independent, since each writes its own output slots and reads nothing another writes.
        // Atomic cursor rather than a static chunk because the cells are not equal work -- a
        // tryptophan is nine atoms and an alanine is none.
        const int n_cell = P * 20;
        std::atomic<int> next{0};
        auto worker = [&]() {
            for (;;) {
                const int c = next.fetch_add(1, std::memory_order_relaxed);
                if (c >= n_cell) return;
                scan_cell(in, c / 20, c % 20, pp, ap, cp, rp, pairs_of_pos);
            }
        };
        if (nt > n_cell) nt = n_cell > 0 ? n_cell : 1;
        if (nt == 1 || n_cell == 0) {
            worker();
        } else {
            std::vector<std::thread> pool;
            pool.reserve(nt);
            for (int t = 0; t < nt; ++t) pool.emplace_back(worker);
            for (auto& th : pool) th.join();
        }
    }

    py::dict out;
    out["p"] = p_out;
    out["admissible"] = adm;
    out["n_clash"] = n_clash;
    out["n_rotamer"] = n_rot;
    return out;
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
    m.def("relax_interface", &relax_interface, py::arg("xyz"), py::arg("atom_class"),
          py::arg("pep_lo"), py::arg("pep_hi"), py::arg("pep_res"), py::arg("tor_a"),
          py::arg("tor_b"), py::arg("tor_mov_ptr"), py::arg("tor_mov"), py::arg("par_atom"),
          py::arg("anchor_atom"), py::arg("dope_table"), py::arg("n_cls"), py::arg("n_bins"),
          py::arg("x_start"), py::arg("dx"), py::arg("n_steps") = 4000,
          py::arg("temperature") = 4.0, py::arg("sigma_deg") = 6.0, py::arg("intra_w") = 1.0,
          py::arg("anchor_w") = 1.0, py::arg("min_sep") = 3, py::arg("burn_in") = 500,
          py::arg("seed") = 0,
          "Flexible-backbone Metropolis MC of the peptide; returns the best pose plus RMSF/drift.");
    m.def("substitution_scan", &substitution_scan,
          py::arg("bb"), py::arg("t_a"), py::arg("t_b"), py::arg("t_c"), py::arg("t_chi"),
          py::arg("t_bond"), py::arg("t_ang"), py::arg("t_off"), py::arg("t_cls"), py::arg("t_rad"),
          py::arg("n_atom"), py::arg("n_chi"), py::arg("chi_vals"),
          py::arg("env_xyz"), py::arg("env_cls"), py::arg("env_rad"), py::arg("env_res"),
          py::arg("env_ptr"), py::arg("env_idx"), py::arg("par_xyz"), py::arg("par_ptr"),
          py::arg("pair_i"), py::arg("pair_j"), py::arg("dope_table"), py::arg("n_cls"),
          py::arg("n_bins"), py::arg("x_start"), py::arg("dx"), py::arg("cutoff"),
          py::arg("clash_tol"), py::arg("temperature"), py::arg("n_threads") = 0,
          "Every amino acid built at every rotamer: contact probability, admissibility and the "
          "clashing-rotamer count, one call per structure.");
    m.attr("__version__") = "0.4.0";
}
