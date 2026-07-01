// Native interface-geometry descriptors for TCR binder/non-binder ranking.
//
// The AF-orthogonal structural signal that beats AlphaFold/TCRmodel2 confidence: interface burial
// (buried SASA), shape/packing, size, dual-chain balance and H-bonds. AlphaFold cannot manufacture a
// large, well-packed, buried interface for a non-cognate TCR, so these discriminate binders even though
// the pairwise contact *energy* on forced AF poses does not. All compute is here in C++ (stdlib only,
// pybind11) — Python is thin glue (chain annotation + the final small logistic). Matches the manuscript
// feature definitions (scripts/tcrvdb_features/{geom_sasa,geom_cov,geom_fast}.py) so the trained model
// reproduces.
//
// Kernels: shrake_rupley (per-atom SASA, grid neighbour list), interface_hbonds (polar-pair count),
// contact_descriptors (interface size + dual-chain balance). Lawrence-Colman shape complementarity is
// added later (the one hard descriptor).

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <algorithm>
#include <cmath>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace py = pybind11;

namespace {

struct Vec3 {
    double x, y, z;
};
inline double d2(const Vec3& a, const Vec3& b) {
    const double dx = a.x - b.x, dy = a.y - b.y, dz = a.z - b.z;
    return dx * dx + dy * dy + dz * dz;
}

std::vector<Vec3> load_xyz(py::array_t<double> a) {
    auto r = a.unchecked<2>();
    std::vector<Vec3> v(r.shape(0));
    for (py::ssize_t i = 0; i < r.shape(0); ++i) v[i] = {r(i, 0), r(i, 1), r(i, 2)};
    return v;
}
std::vector<int> load_int(py::array_t<int> a) {
    auto r = a.unchecked<1>();
    std::vector<int> v(r.shape(0));
    for (py::ssize_t i = 0; i < r.shape(0); ++i) v[i] = r(i);
    return v;
}

// Golden-spiral (Fibonacci) unit sphere of n points — the sampling ShrakeRupley uses.
std::vector<Vec3> unit_sphere(int n) {
    std::vector<Vec3> pts(n);
    const double off = 2.0 / n, inc = M_PI * (3.0 - std::sqrt(5.0));
    for (int k = 0; k < n; ++k) {
        const double y = k * off - 1.0 + off / 2.0;
        const double r = std::sqrt(std::max(0.0, 1.0 - y * y));
        const double phi = k * inc;
        pts[k] = {std::cos(phi) * r, y, std::sin(phi) * r};
    }
    return pts;
}

// Per-atom SASA via Shrake-Rupley with a uniform-grid neighbour list (fast on full complexes).
py::array_t<double> shrake_rupley(py::array_t<double> xyz, py::array_t<double> radii,
                                  double probe, int n_points) {
    std::vector<Vec3> pos = load_xyz(xyz);
    auto rr = radii.unchecked<1>();
    const int n = static_cast<int>(pos.size());
    std::vector<double> rad(n);
    double rmax = 0.0;
    for (int i = 0; i < n; ++i) { rad[i] = rr(i) + probe; rmax = std::max(rmax, rad[i]); }

    py::array_t<double> out(n);
    auto w = out.mutable_unchecked<1>();
    if (n == 0) return out;

    // Grid: cell size = 2*rmax so any contacting neighbour lies in the 27-cell stencil.
    const double cell = std::max(2.0 * rmax, 1e-3);
    auto key = [cell](double v) { return static_cast<long>(std::floor(v / cell)); };
    std::unordered_map<long long, std::vector<int>> grid;
    auto hash = [](long a, long b, long c) {
        return (static_cast<long long>(a) * 73856093LL) ^ (static_cast<long long>(b) * 19349663LL) ^
               (static_cast<long long>(c) * 83492791LL);
    };
    for (int i = 0; i < n; ++i)
        grid[hash(key(pos[i].x), key(pos[i].y), key(pos[i].z))].push_back(i);

    const std::vector<Vec3> sphere = unit_sphere(n_points);
    {
        py::gil_scoped_release release;
        for (int i = 0; i < n; ++i) {
            // Collect neighbours from the 27-cell stencil.
            std::vector<int> nb;
            const long cx = key(pos[i].x), cy = key(pos[i].y), cz = key(pos[i].z);
            for (long dx = -1; dx <= 1; ++dx)
                for (long dy = -1; dy <= 1; ++dy)
                    for (long dz = -1; dz <= 1; ++dz) {
                        auto it = grid.find(hash(cx + dx, cy + dy, cz + dz));
                        if (it == grid.end()) continue;
                        for (int j : it->second) {
                            if (j == i) continue;
                            const double rr2 = (rad[i] + rad[j]) * (rad[i] + rad[j]);
                            if (d2(pos[i], pos[j]) < rr2) nb.push_back(j);
                        }
                    }
            int acc = 0;
            for (const Vec3& s : sphere) {
                const Vec3 p{pos[i].x + rad[i] * s.x, pos[i].y + rad[i] * s.y, pos[i].z + rad[i] * s.z};
                bool buried = false;
                for (int j : nb) {
                    if (d2(p, pos[j]) < rad[j] * rad[j]) { buried = true; break; }
                }
                if (!buried) ++acc;
            }
            w(i) = (4.0 * M_PI * rad[i] * rad[i]) * (static_cast<double>(acc) / n_points);
        }
    }
    return out;
}

// Interface H-bond proxy: polar (N/O) heavy-atom pairs between the two sides within dist_cutoff.
// Matches the manuscript's crude n_hbond (no angle term) that trained the binder model.
int interface_hbonds(py::array_t<double> donor_xyz, py::array_t<double> acceptor_xyz,
                     double dist_cutoff) {
    std::vector<Vec3> a = load_xyz(donor_xyz), b = load_xyz(acceptor_xyz);
    const double c2 = dist_cutoff * dist_cutoff;
    int count = 0;
    py::gil_scoped_release release;
    for (const Vec3& pa : a)
        for (const Vec3& pb : b)
            if (d2(pa, pb) < c2) ++count;
    return count;
}

// Interface size (# distinct TCR residues contacting peptide+MHC) and dual-chain balance.
// pm_cov_ntcr: distinct TCR residues (TRA/TRB, disambiguated by chain flag) with a heavy atom within
//   contact_cut of any peptide OR MHC atom.
// chain_balance: na = # peptide atoms whose nearest TRA atom < bal_cut, nb likewise for TRB;
//   min(na,nb)/max(na+nb,1)  (== geom_cov.py).
py::dict contact_descriptors(py::array_t<double> tcra_xyz, py::array_t<int> tcra_res,
                             py::array_t<double> tcrb_xyz, py::array_t<int> tcrb_res,
                             py::array_t<double> pep_xyz, py::array_t<double> mhc_xyz,
                             double contact_cut, double bal_cut) {
    std::vector<Vec3> ta = load_xyz(tcra_xyz), tb = load_xyz(tcrb_xyz);
    std::vector<Vec3> pep = load_xyz(pep_xyz), mhc = load_xyz(mhc_xyz);
    std::vector<int> ra = load_int(tcra_res), rb = load_int(tcrb_res);
    const double cc2 = contact_cut * contact_cut, bc2 = bal_cut * bal_cut;

    // pm_cov_ntcr: distinct engaged TCR residues (a<res>, b<res> namespaced by chain).
    std::unordered_set<long long> engaged;
    auto scan_side = [&](const std::vector<Vec3>& t, const std::vector<int>& res, long long tag) {
        for (size_t i = 0; i < t.size(); ++i) {
            bool hit = false;
            for (const Vec3& q : pep) { if (d2(t[i], q) < cc2) { hit = true; break; } }
            if (!hit) for (const Vec3& q : mhc) { if (d2(t[i], q) < cc2) { hit = true; break; } }
            if (hit) engaged.insert(tag * 100000LL + res[i]);
        }
    };
    {
        py::gil_scoped_release release;
        scan_side(ta, ra, 1);
        scan_side(tb, rb, 2);
    }
    const int pm_cov_ntcr = static_cast<int>(engaged.size());

    // chain_balance from per-peptide-atom nearest-TRA / nearest-TRB distances.
    int na = 0, nb = 0;
    {
        py::gil_scoped_release release;
        for (const Vec3& p : pep) {
            double mina = 1e30, minb = 1e30;
            for (const Vec3& q : ta) mina = std::min(mina, d2(p, q));
            for (const Vec3& q : tb) minb = std::min(minb, d2(p, q));
            if (mina < bc2) ++na;
            if (minb < bc2) ++nb;
        }
    }
    const double balance = static_cast<double>(std::min(na, nb)) / std::max(na + nb, 1);

    py::dict out;
    out["pm_cov_ntcr"] = pm_cov_ntcr;
    out["chain_balance"] = balance;
    out["n_pep_near_tra"] = na;
    out["n_pep_near_trb"] = nb;
    return out;
}

}  // namespace

PYBIND11_MODULE(_geom, m) {
    m.doc() = "Native interface-geometry descriptors (SASA, H-bonds, size/balance) for binder ranking.";
    m.def("shrake_rupley", &shrake_rupley, py::arg("xyz"), py::arg("radii"),
          py::arg("probe") = 1.4, py::arg("n_points") = 100,
          "Per-atom solvent-accessible surface area (Shrake-Rupley, grid neighbour list).");
    m.def("interface_hbonds", &interface_hbonds, py::arg("donor_xyz"), py::arg("acceptor_xyz"),
          py::arg("dist_cutoff") = 3.5, "Count polar-polar heavy-atom pairs across an interface.");
    m.def("contact_descriptors", &contact_descriptors, py::arg("tcra_xyz"), py::arg("tcra_res"),
          py::arg("tcrb_xyz"), py::arg("tcrb_res"), py::arg("pep_xyz"), py::arg("mhc_xyz"),
          py::arg("contact_cut") = 5.0, py::arg("bal_cut") = 4.5,
          "Interface size (# engaged TCR residues) + dual-chain balance.");
    m.attr("__version__") = "0.1.0";
}
