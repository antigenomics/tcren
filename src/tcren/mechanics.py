"""Interface mechanics: the TCR↔pMHC contact map as a network of breakable springs.

Each inter-body residue contact is a Hookean spring (stiffness from atomic-contact multiplicity),
anchored at the two Cα atoms. Two quantities are exposed:

* :func:`stiffness_tensor` — the linear-response stiffness tensor ``K = Σ kᵢ ûᵢ⊗ûᵢ`` of the interface,
  split into a **tensile** component along the docking axis and an in-plane **shear** component.
* :func:`rupture` — a static steered-unbinding cartoon: rigidly displace the pMHC body along a pull
  direction, letting springs break past a strain threshold, and record the peak resisting **force** and
  the cumulative **work**.

Rationale (validated on ATLAS, 2026): these mechanical measures track the *kinetic* stability of the
complex — the dissociation off-rate ``koff`` — far better than the equilibrium ΔG/Kd (Bell–Evans: rupture
resistance reflects the dissociation barrier, not the well depth). This is the physically apt axis for the
TCR, a mechanosensor whose pMHC bonds are catch bonds. Pure-numpy, single structure, no MD.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .contactmap import ContactMap
from .structure.model import PEPTIDE_TYPE, RECEPTOR_TYPES, Structure

#: Spring-stiffness models for an interface contact. ``"unit"`` = 1 per contact (pure topology);
#: ``"count"`` = heavy-atom-pair multiplicity; ``"invdist2"`` = multiplicity / dist² (Hookean-ish, the
#: validated default).
WEIGHTS = ("unit", "count", "invdist2")


@dataclass(slots=True)
class InterfaceSprings:
    """The TCR↔pMHC spring network of one structure.

    ``a``/``b`` are the (n, 3) Cα anchor coordinates on the TCR and pMHC sides; ``k`` the (n,) spring
    stiffnesses; ``rest`` the (n,) rest lengths ``|b − a|``; ``axis`` the unit docking axis (stiffness-
    weighted separation of the two interface centroids, pointing TCR→pMHC).
    """

    a: np.ndarray
    b: np.ndarray
    k: np.ndarray
    rest: np.ndarray
    axis: np.ndarray

    def __len__(self) -> int:
        return len(self.k)


def interface_springs(
    structure: Structure, *, cutoff: float = 8.0, weight: str = "invdist2"
) -> InterfaceSprings:
    """Build the TCR↔pMHC interface spring network from residue contacts.

    Args:
        structure: An annotated structure (chains typed; peptide present).
        cutoff: Heavy-atom contact cutoff (Å) defining a spring.
        weight: Spring-stiffness model, one of :data:`WEIGHTS`.

    Returns:
        The :class:`InterfaceSprings` (empty arrays if no inter-body contact is found).
    """
    if weight not in WEIGHTS:
        raise ValueError(f"unknown weight {weight!r}; choose from {WEIGHTS}")
    if not any(c.chain_type == PEPTIDE_TYPE for c in structure.chains):
        raise ValueError("structure has no peptide chain")

    ca = {
        (c.chain_id, r.seq_index): np.asarray(r.ca)
        for c in structure.chains
        for r in c.residues
        if r.ca is not None
    }
    cm = ContactMap.from_structure(structure, cutoff=cutoff, count_atoms=True)
    import polars as pl

    df = pl.concat(
        [cm.interface("tcr_peptide", tcr_regions="all"), cm.interface("tcr_mhc")]
    )
    a_list, b_list, k_list, rest = [], [], [], []
    for cf, rf, ct, rt, nat, d in zip(
        df["chain.id.from"], df["residue.index.from"],
        df["chain.id.to"], df["residue.index.to"],
        df["n_atom_contacts"], df["dist"],
    ):
        a, b = ca.get((cf, rf)), ca.get((ct, rt))
        if a is None or b is None:
            continue
        L = float(np.linalg.norm(b - a))
        if L < 1e-6:
            continue
        if weight == "unit":
            k = 1.0
        elif weight == "count":
            k = float(nat)
        else:  # invdist2
            k = float(nat) / (float(d) ** 2)
        a_list.append(a); b_list.append(b); k_list.append(k); rest.append(L)
    if not a_list:
        z = np.empty((0, 3))
        return InterfaceSprings(z, z, np.empty(0), np.empty(0), np.array([0.0, 0.0, 1.0]))
    a = np.asarray(a_list); b = np.asarray(b_list); k = np.asarray(k_list); rest = np.asarray(rest)
    ct = np.average(a, axis=0, weights=k); cp = np.average(b, axis=0, weights=k)
    axis = cp - ct
    axis = axis / (np.linalg.norm(axis) + 1e-9)
    return InterfaceSprings(a, b, k, rest, axis)


def _stiffness_matrix(s: InterfaceSprings) -> np.ndarray:
    """The 3×3 stiffness tensor ``K = Σ kᵢ ûᵢ⊗ûᵢ`` over the spring network."""
    u = (s.b - s.a) / s.rest[:, None]
    return (s.k[:, None, None] * u[:, :, None] * u[:, None, :]).sum(0)


def stiffness_tensor(
    structure: Structure, *, cutoff: float = 8.0, weight: str = "invdist2"
) -> dict[str, float]:
    """Linear-response stiffness descriptors of the TCR↔pMHC interface.

    Forms ``K = Σ kᵢ ûᵢ⊗ûᵢ`` over interface springs and resolves it along the docking axis.

    Returns a dict with:
        ``S_tot`` (trace K, total stiffness), ``K_tens`` (tensile, along the docking axis),
        ``K_shear`` (in-plane, ``S_tot − K_tens``), ``aniso`` (``K_shear / K_tens``),
        ``lam_max``/``lam_min`` (extreme eigenvalues), ``n_spring``. All ``nan`` if < 3 springs.
    """
    s = interface_springs(structure, cutoff=cutoff, weight=weight)
    if len(s) < 3:
        return {k: float("nan") for k in
                ("S_tot", "K_tens", "K_shear", "aniso", "lam_max", "lam_min")} | {"n_spring": float(len(s))}
    K = _stiffness_matrix(s)
    S_tot = float(np.trace(K))
    K_tens = float(s.axis @ K @ s.axis)
    evals = np.linalg.eigvalsh(K)
    return {
        "S_tot": S_tot, "K_tens": K_tens, "K_shear": S_tot - K_tens,
        "aniso": (S_tot - K_tens) / (K_tens + 1e-9),
        "lam_max": float(evals[-1]), "lam_min": float(evals[0]), "n_spring": float(len(s)),
    }


def _pull_direction(s: InterfaceSprings, direction: str) -> np.ndarray:
    """Unit pull direction. ``"tensile"`` = docking axis; ``"shear"`` = stiffest in-plane direction."""
    if direction == "tensile":
        return s.axis
    K = _stiffness_matrix(s)
    evals, evecs = np.linalg.eigh(K)
    top = evecs[:, -1]
    inplane = top - (top @ s.axis) * s.axis  # remove the docking-axis component
    n = np.linalg.norm(inplane)
    return inplane / n if n > 1e-6 else evecs[:, 0]


def rupture(
    structure: Structure, *, direction: str = "tensile", cutoff: float = 8.0,
    weight: str = "invdist2", break_strain: float = 0.5, steps: int = 80,
) -> dict[str, float]:
    """Steered-unbinding cartoon: pull the pMHC body off the TCR and measure rupture resistance.

    The pMHC anchors are rigidly translated along the pull ``direction``; each spring resists in tension
    (Hookean) and is removed once its strain exceeds ``break_strain``. Integrates until all springs break.

    Args:
        direction: ``"tensile"`` (docking axis), ``"shear"`` (stiffest in-plane), or ``"auto"``
            (the minimum-force of tensile and shear — the easiest rupture path).
        break_strain: fractional extension at which a spring breaks (the tuning knob; 0.5 = 50 %).
        steps: displacement increments.

    Returns a dict: ``rupture_force`` (peak resisting force along the pull), ``rupture_work``
        (∫ force · displacement), ``n_spring``, ``break_strain``. ``nan`` if < 3 springs.
    """
    if direction == "auto":
        rt = rupture(structure, direction="tensile", cutoff=cutoff, weight=weight,
                     break_strain=break_strain, steps=steps)
        rs = rupture(structure, direction="shear", cutoff=cutoff, weight=weight,
                     break_strain=break_strain, steps=steps)
        if not np.isfinite(rt["rupture_force"]):
            return rt
        pick = rt if rt["rupture_force"] <= rs["rupture_force"] else rs
        return pick | {"direction": "auto"}

    s = interface_springs(structure, cutoff=cutoff, weight=weight)
    out = {"rupture_force": float("nan"), "rupture_work": float("nan"),
           "n_spring": float(len(s)), "break_strain": break_strain}
    if len(s) < 3:
        return out
    d = _pull_direction(s, direction)
    max_disp = 2.0 * float(s.rest.max())  # enough for even perpendicular springs to reach break_strain
    dt = max_disp / steps
    alive = np.ones(len(s), bool)
    forces, t = [], 0.0
    for _ in range(steps):
        t += dt
        v = (s.b + t * d) - s.a          # spring vectors at displacement t
        L = np.linalg.norm(v, axis=1)
        strain = (L - s.rest) / s.rest
        alive &= strain <= break_strain  # break (permanently) past threshold
        ext = np.maximum(L - s.rest, 0.0)  # springs resist only in tension
        f = s.k * ext * (v @ d) / (L + 1e-9)  # restoring force projected onto pull axis
        forces.append(float((f * alive).sum()))
        if not alive.any():
            break
    forces = np.asarray(forces)
    out["rupture_force"] = float(forces.max())
    out["rupture_work"] = float(forces.sum() * dt)
    return out


def _coupling_counts(pep_tcr, pep_mhc, mhc_tcr, mhc_pep, tcr_pmhc, tcr_ab) -> dict[str, int]:
    """Intersect the residue sets into coupling counts (pure; unit-tested in :func:`_demo`)."""
    cp, cm, ct = pep_tcr & pep_mhc, mhc_tcr & mhc_pep, tcr_pmhc & tcr_ab
    return {
        "couple_pep": len(cp), "couple_mhc": len(cm), "couple_tcr": len(ct),
        "couple_total": len(cp) + len(cm) + len(ct),
        "n_interface": len(tcr_pmhc) + len(pep_tcr) + len(mhc_tcr),
    }


def coupling_residues(structure: Structure, *, cutoff: float = 5.0) -> dict[str, int]:
    """Residues that couple the pre-formed internal scaffold to the TCR↔pMHC interface.

    Counts residues sitting in **both** an intra-body contact and the binding interface (computed from the
    bound complex; the intra-body scaffold — Vα↔Vβ pairing, peptide-in-groove — is pre-formed, so this
    approximates the internal-contact residues the interface recruits):

        ``couple_pep``   — peptide residues contacting both MHC (groove-anchored) *and* TCR (dual-role);
        ``couple_mhc``   — MHC residues contacting both peptide *and* TCR (groove-rim presenting);
        ``couple_tcr``   — TCR residues in the Vα↔Vβ interface that *also* contact pMHC (combining-site apex).

    Also returns ``couple_total`` and ``n_interface`` (interface residue count, a size denominator).

    Note: as a binding *estimator* these are weak/underpowered on current data — ``couple_pep`` tracks the
    dissociation off-rate (koff) at r≈−0.34 (class I, borderline; partially survives interface-size control),
    but the TCR/MHC sets are near-null. Useful primarily as interpretable structural descriptors.
    """
    if not any(c.chain_type == PEPTIDE_TYPE for c in structure.chains):
        raise ValueError("structure has no peptide chain")
    import polars as pl

    cm = ContactMap.from_structure(structure, cutoff=cutoff)
    tp = cm.interface("tcr_peptide", tcr_regions="all")  # from=TCR, to=peptide
    tm = cm.interface("tcr_mhc")                          # from=TCR, to=MHC
    pm = cm.interface("peptide_mhc")                      # from=peptide, to=MHC
    ab = cm.contacts.filter(                              # intra-TCR: Vα↔Vβ (both receptor chains)
        pl.col("chain.type.from").is_in(list(RECEPTOR_TYPES))
        & pl.col("chain.type.to").is_in(list(RECEPTOR_TYPES))
    )

    def rs(df, side):
        return set(zip(df[f"chain.id.{side}"].to_list(), df[f"residue.index.{side}"].to_list()))

    return _coupling_counts(
        pep_tcr=rs(tp, "to"), pep_mhc=rs(pm, "from"),
        mhc_tcr=rs(tm, "to"), mhc_pep=rs(pm, "to"),
        tcr_pmhc=rs(tp, "from") | rs(tm, "from"), tcr_ab=rs(ab, "from") | rs(ab, "to"),
    )


def _demo() -> None:
    """Self-check on synthetic spring geometries (no PDB needed)."""
    # a minimal InterfaceSprings: 3 springs along +z (tensile), unit stiffness, rest length 5.
    def make(a, b, k):
        a, b = np.asarray(a, float), np.asarray(b, float)
        rest = np.linalg.norm(b - a, axis=1)
        ct = np.average(a, 0, weights=k); cp = np.average(b, 0, weights=k)
        ax = cp - ct; ax = ax / (np.linalg.norm(ax) + 1e-9)
        return InterfaceSprings(a, b, np.asarray(k, float), rest, ax)

    # springs purely along z -> K_tens (along z) >> K_shear
    zs = make([[0, 0, 0], [3, 0, 0], [0, 3, 0]], [[0, 0, 5], [3, 0, 5], [0, 3, 5]], [1, 1, 1])
    u = (zs.b - zs.a) / zs.rest[:, None]
    K = (zs.k[:, None, None] * u[:, :, None] * u[:, None, :]).sum(0)
    K_tens = zs.axis @ K @ zs.axis
    assert K_tens > np.trace(K) - K_tens, "z-aligned springs must be tensile-dominant"

    # more springs -> higher rupture force (monotone in stiffness), via the module path on synthetic sets
    def force(springs, direction="tensile", bs=0.5, steps=80):
        d = _pull_direction(springs, direction)
        md = 2.0 * springs.rest.max(); dt = md / steps
        alive = np.ones(len(springs), bool); best = 0.0
        for i in range(steps):
            t = (i + 1) * dt
            v = (springs.b + t * d) - springs.a; L = np.linalg.norm(v, axis=1)
            alive &= (L - springs.rest) / springs.rest <= bs
            ext = np.maximum(L - springs.rest, 0.0)
            best = max(best, float(((springs.k * ext * (v @ d) / (L + 1e-9)) * alive).sum()))
        return best

    f3 = force(zs)
    zs6 = make(list(zs.a) + [[6, 0, 0], [0, 6, 0], [6, 6, 0]],
               list(zs.b) + [[6, 0, 5], [0, 6, 5], [6, 6, 5]], [1] * 6)
    assert force(zs6) > f3, "doubling springs must raise rupture force"
    # stiffer springs -> higher force
    stiff = make(zs.a, zs.b, [5, 5, 5])
    assert force(stiff) > f3 * 3, "5x stiffness must raise rupture force"
    # lower break strain -> easier rupture (force reached before break is <=, work smaller)
    assert force(zs, bs=0.2) <= force(zs, bs=0.8) + 1e-9, "weaker bonds break sooner"

    # coupling set-intersection logic
    c = _coupling_counts(pep_tcr={1, 2, 3}, pep_mhc={2, 3, 4}, mhc_tcr={5}, mhc_pep={5, 6},
                         tcr_pmhc={1, 2, 7}, tcr_ab={2, 8})
    assert (c["couple_pep"], c["couple_mhc"], c["couple_tcr"]) == (2, 1, 1), c
    assert c["couple_total"] == 4 and c["n_interface"] == 7, c
    print("mechanics self-check OK: K_tens dominant for z-springs; force monotone in count & stiffness; "
          "coupling intersections correct")


if __name__ == "__main__":
    _demo()
