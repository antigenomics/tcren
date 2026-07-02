// Interface energy + (later) sampling for native ΔΔG — the FlexPepDock/InterfaceAnalyzer replacement.
//
// interface_energy: the DOPE atom-level statistical-potential interaction energy across an interface
// (peptide <-> partner heavy-atom pairs only — the peptide-internal and partner-internal terms are not
// summed, so this is already E_bound − E_separated for the cross terms, i.e. the interaction energy).
// A 20-line lift of energy() from src/_refine/refine.cpp with the harmonic restraint dropped. Stdlib
// only, pybind11, C++17 — the house pattern. The rotamer repack + flexible-backbone relax that turn
// this into a full ΔΔG land here later (functions repack / relax_interface); this is the energy core
// they all call.

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <cmath>
#include <vector>

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

}  // namespace

PYBIND11_MODULE(_relax, m) {
    m.doc() = "Native interface energy + sampling for ΔΔG (DOPE interaction energy; repack/relax later).";
    m.def("interface_energy", &interface_energy,
          py::arg("pep_xyz"), py::arg("pep_class"), py::arg("par_xyz"), py::arg("par_class"),
          py::arg("dope_table"), py::arg("n_cls"), py::arg("n_bins"), py::arg("x_start"), py::arg("dx"),
          "DOPE interaction energy over peptide<->partner heavy-atom pairs (the interface ΔG core).");
    m.attr("__version__") = "0.1.0";
}
