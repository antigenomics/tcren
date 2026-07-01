// Cyclic Coordinate Descent (CCD) loop closure for peptide backbone modelling.
//
// An open-source, license-free geometric kernel for driving a kinematic backbone chain so that a set
// of "anchor" atoms reach target positions (e.g. predicted MHC-groove pocket centroids), while the
// rest of the chain follows as a linkage. This is the analytic core of the open replacement for
// Rosetta FlexPepDock's loop refinement / MODELLER's loopmodel anchor restraints -- it uses ONLY the
// C++ standard library (no Eigen/Boost), matching the existing tcren._refine / tcren._align kernels.
//
// Method: CCD (Canutescu & Dunbrack, Protein Science 2003). For each rotatable backbone bond, taken
// in order from the chain base toward the moving tip, we rotate the downstream atoms about the bond
// axis by the angle that MINIMISES the weighted sum of squared distances between the anchor atoms and
// their targets. That optimal angle has a closed form,
//
//     theta* = atan2( Sum_k w_k * ( s_k . (axis x r_k) ) ,  Sum_k w_k * ( s_k . r_k ) )
//
// where r_k is the anchor's position relative to the bond axis (perpendicular component) and s_k is
// the target's. Because each per-bond step is locally optimal, the anchor RMSD decreases monotonically
// and the sweep converges. KIC (kinematic / robotics closure, Coutsias 2004) is the higher-accuracy
// alternative for the 3-pivot exact-closure case and is exposed as a stub for a future companion.

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace py = pybind11;

namespace {

struct Vec3 {
    double x, y, z;
};

inline Vec3 sub(const Vec3& a, const Vec3& b) { return {a.x - b.x, a.y - b.y, a.z - b.z}; }
inline Vec3 add(const Vec3& a, const Vec3& b) { return {a.x + b.x, a.y + b.y, a.z + b.z}; }
inline Vec3 scale(const Vec3& a, double s) { return {a.x * s, a.y * s, a.z * s}; }
inline double dot(const Vec3& a, const Vec3& b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
inline Vec3 cross(const Vec3& a, const Vec3& b) {
    return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}
inline double norm(const Vec3& a) { return std::sqrt(dot(a, a)); }

// Rotate point p about the axis (unit `u`) through `o` by `angle` (Rodrigues' rotation).
inline Vec3 rotate_about(const Vec3& p, const Vec3& o, const Vec3& u, double angle) {
    const Vec3 v = sub(p, o);
    const double c = std::cos(angle), s = std::sin(angle);
    const Vec3 term1 = scale(v, c);
    const Vec3 term2 = scale(cross(u, v), s);
    const Vec3 term3 = scale(u, dot(u, v) * (1.0 - c));
    return add(o, add(add(term1, term2), term3));
}

// Weighted anchor RMSD to targets.
double anchor_rmsd(const std::vector<Vec3>& coords, const std::vector<int>& moving,
                   const std::vector<Vec3>& targets, const std::vector<double>& w) {
    double num = 0.0, den = 0.0;
    for (size_t k = 0; k < moving.size(); ++k) {
        const Vec3 d = sub(coords[moving[k]], targets[k]);
        num += w[k] * dot(d, d);
        den += w[k];
    }
    return den > 0.0 ? std::sqrt(num / den) : 0.0;
}

// CCD closure. `bonds` is (M,2): rotatable bond m has axis from atom bonds[m][0] to bonds[m][1];
// rotating it moves every atom with index > bonds[m][1] (downstream of the bond). `moving`/`targets`
// are the anchor atom indices and their target xyz; `w` per-anchor weights. Returns the closed coords,
// the final anchor RMSD, and the iterations used.
py::tuple ccd_close(py::array_t<double> coords_in, py::array_t<int> bonds_in,
                    py::array_t<int> moving_in, py::array_t<double> targets_in,
                    py::array_t<double> weights_in, int max_iter, double tol) {
    auto cr = coords_in.unchecked<2>();
    std::vector<Vec3> coords(cr.shape(0));
    for (py::ssize_t i = 0; i < cr.shape(0); ++i) coords[i] = {cr(i, 0), cr(i, 1), cr(i, 2)};

    auto br = bonds_in.unchecked<2>();
    std::vector<std::pair<int, int>> bonds(br.shape(0));
    for (py::ssize_t i = 0; i < br.shape(0); ++i) bonds[i] = {br(i, 0), br(i, 1)};

    auto mr = moving_in.unchecked<1>();
    std::vector<int> moving(mr.shape(0));
    for (py::ssize_t i = 0; i < mr.shape(0); ++i) moving[i] = mr(i);

    auto tr = targets_in.unchecked<2>();
    std::vector<Vec3> targets(tr.shape(0));
    for (py::ssize_t i = 0; i < tr.shape(0); ++i) targets[i] = {tr(i, 0), tr(i, 1), tr(i, 2)};

    auto wr = weights_in.unchecked<1>();
    std::vector<double> w(wr.shape(0));
    for (py::ssize_t i = 0; i < wr.shape(0); ++i) w[i] = wr(i);

    const int n = static_cast<int>(coords.size());

    // Validate every index/length BEFORE releasing the GIL: the inner loops index std::vector with
    // operator[] (no bounds check), so an out-of-range index from Python would be undefined behaviour.
    if (targets.size() < moving.size() || w.size() < moving.size())
        throw std::invalid_argument("targets and weights must each have >= len(moving) rows");
    for (const auto& bond : bonds) {
        if (bond.first < 0 || bond.first >= n || bond.second < 0 || bond.second >= n)
            throw std::invalid_argument("bond index out of range [0, n_atoms)");
    }
    for (int mi : moving) {
        if (mi < 0 || mi >= n)
            throw std::invalid_argument("moving (anchor) index out of range [0, n_atoms)");
    }

    double rmsd = anchor_rmsd(coords, moving, targets, w);
    int it = 0;
    {
        py::gil_scoped_release release;
        for (; it < max_iter && rmsd > tol; ++it) {
            for (const auto& bond : bonds) {
                const Vec3 o = coords[bond.first];
                Vec3 axis = sub(coords[bond.second], o);
                const double an = norm(axis);
                if (an < 1e-9) continue;
                axis = scale(axis, 1.0 / an);

                // Closed-form optimal angle over the anchors downstream of this bond.
                double num = 0.0, den = 0.0;  // num = Sum w (s . (axis x r)), den = Sum w (s . r)
                for (size_t k = 0; k < moving.size(); ++k) {
                    if (moving[k] <= bond.second) continue;  // upstream anchor: unaffected by this bond
                    const Vec3 p = coords[moving[k]];
                    // r = current anchor perpendicular to axis (relative to o)
                    const Vec3 vp = sub(p, o);
                    const Vec3 r = sub(vp, scale(axis, dot(axis, vp)));
                    // s = target perpendicular to axis (relative to o)
                    const Vec3 vt = sub(targets[k], o);
                    const Vec3 s = sub(vt, scale(axis, dot(axis, vt)));
                    num += w[k] * dot(s, cross(axis, r));
                    den += w[k] * dot(s, r);
                }
                if (num == 0.0 && den == 0.0) continue;
                const double theta = std::atan2(num, den);
                if (std::abs(theta) < 1e-12) continue;
                for (int j = bond.second + 1; j < n; ++j)
                    coords[j] = rotate_about(coords[j], o, axis, theta);
            }
            rmsd = anchor_rmsd(coords, moving, targets, w);
        }
    }

    py::array_t<double> out({static_cast<py::ssize_t>(coords.size()), py::ssize_t(3)});
    auto ow = out.mutable_unchecked<2>();
    for (size_t i = 0; i < coords.size(); ++i) { ow(i, 0) = coords[i].x; ow(i, 1) = coords[i].y; ow(i, 2) = coords[i].z; }
    return py::make_tuple(out, rmsd, it);
}

}  // namespace

PYBIND11_MODULE(_fold, m) {
    m.doc() = "Open-source CCD loop-closure kernel for anchor-restrained peptide backbone modelling.";
    m.def("ccd_close", &ccd_close,
          py::arg("coords"), py::arg("bonds"), py::arg("moving"), py::arg("targets"),
          py::arg("weights"), py::arg("max_iter") = 1000, py::arg("tol") = 0.08,
          "Cyclic Coordinate Descent: rotate backbone bonds to drive anchor atoms onto targets.\n"
          "Returns (closed_coords (N,3), final_anchor_rmsd, iterations).");
    m.attr("__version__") = "0.1.0";
}
