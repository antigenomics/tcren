"""Contact stability / fragility at the TCR:peptide interface.

The 5 Å contact cutoff is a hard threshold on a continuous distance: a contact at 4.9 Å has almost
no positional slack (a ~1 Å shift, or tightening the cutoff, kills it), while one at 3.5 Å is robust.
:func:`contact_stability` reads that slack directly off the contact map. For each TCR:peptide
*contact* --- a receptor-residue / peptide-residue pair whose closest heavy-atom pair is within
``cutoff`` (the same contact definition as :func:`tcren.contacts.geometry.all_atom_contacts`) --- the
*margin* ``m = cutoff - dmin`` is the contact's positional tolerance, and a rigid isotropic shift of
size ``delta`` loses it with probability ``clip((delta - m) / (2*delta), 0, 1)``. Aggregated over the
interface these give ``mean_margin`` (how deep the contacts sit), ``frac_robust`` (fraction with
``m >= delta``) and ``exp_lost`` (expected contacts lost under a ``delta`` shift) --- a coordinate-only
read of interface positional confidence, the physical analogue of an interface PAE.

The per-residue-pair minimum-distance scan is a native ``_geom`` kernel; the numpy implementation
behind it (:func:`_contact_stability_numpy`) is the reference and a fallback where the extension is
unavailable.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..structure.model import PEPTIDE_TYPE, RECEPTOR_TYPES, Structure

try:
    from .. import _geom  # native contact-stability kernel (built by scikit-build-core)
except ImportError:  # pragma: no cover - pure-Python fallback if the extension is unavailable
    _geom = None


@dataclass(frozen=True, slots=True)
class StabilityReport:
    """TCR:peptide contact-stability summary.

    Attributes:
        n_contacts: Receptor-residue / peptide-residue contacts within ``cutoff``.
        mean_margin: Mean ``cutoff - dmin`` over contacts (Å); larger = contacts sit deeper.
        frac_marg_lt1: Fraction of contacts with margin below ``delta`` (fragile).
        frac_robust: Fraction of contacts with margin at least ``delta`` (robust to a ``delta`` shift).
        exp_lost: Expected number of contacts lost under a ``delta`` Å isotropic shift.
    """

    n_contacts: int
    mean_margin: float
    frac_marg_lt1: float
    frac_robust: float
    exp_lost: float


def _interface_atoms(structure: Structure, chain_types):
    """``(xyz (N,3) float64, res_id (N,) int32)`` for the heavy atoms of the given chain types.

    ``res_id`` is a globally-unique residue key ``chain_position * 100000 + seq_index`` --- ``seq_index``
    alone is per-chain, so it would collide between the peptide and a TCR chain (and between TRA/TRB).
    """
    xyz, res = [], []
    for ci, chain in enumerate(structure.chains):
        if chain.chain_type not in chain_types:
            continue
        for r in chain.residues:
            uid = ci * 100000 + r.seq_index
            for a in r.atoms:
                if a.element.strip().upper() == "H":
                    continue
                xyz.append(a.coord)
                res.append(uid)
    return np.asarray(xyz, float).reshape(-1, 3), np.asarray(res, np.int32)


def contact_stability(
    structure: Structure, *, cutoff: float = 5.0, delta: float = 1.0
) -> StabilityReport:
    """Contact stability / fragility of the TCR:peptide interface.

    Args:
        structure: A chain-typed complex with a peptide chain and at least one receptor chain.
        cutoff: Contact distance cutoff (Å); the closest heavy-atom pair of a residue pair must be
            within this to count as a contact (matches ``all_atom_contacts``).
        delta: Positional shift (Å) used for the fragility metrics (``frac_marg_lt1``, ``exp_lost``).

    Returns:
        A :class:`StabilityReport`.

    Raises:
        ValueError: If the structure has no peptide chain or no receptor chain.
    """
    pep_xyz, pep_res = _interface_atoms(structure, (PEPTIDE_TYPE,))
    tcr_xyz, tcr_res = _interface_atoms(structure, RECEPTOR_TYPES)
    if len(pep_xyz) == 0:
        raise ValueError(f"no peptide heavy atoms (chain_type == {PEPTIDE_TYPE!r}) in {structure.pdb_id!r}")
    if len(tcr_xyz) == 0:
        raise ValueError(f"no receptor heavy atoms (chain_type in {RECEPTOR_TYPES}) in {structure.pdb_id!r}")

    if _geom is not None:
        d = _geom.contact_stability(pep_xyz, pep_res, tcr_xyz, tcr_res, cutoff, delta)
    else:
        d = _contact_stability_numpy(pep_xyz, pep_res, tcr_xyz, tcr_res, cutoff, delta)
    return StabilityReport(
        n_contacts=int(d["n5"]),
        mean_margin=float(d["mean_margin"]),
        frac_marg_lt1=float(d["frac_marg_lt1"]),
        frac_robust=float(d["frac_robust"]),
        exp_lost=float(d["exp_lost"]),
    )


def _contact_stability_numpy(pep_xyz, pep_res, tcr_xyz, tcr_res, cutoff, delta):
    """Pure-numpy reference for the native ``contact_stability`` kernel."""
    d = np.sqrt(((tcr_xyz[:, None, :] - pep_xyz[None, :, :]) ** 2).sum(-1))  # (T, P)
    ti, pj = np.where(d <= cutoff)
    mind: dict[tuple[int, int], float] = {}
    for t, p in zip(ti, pj):
        key = (int(tcr_res[t]), int(pep_res[p]))
        dd = float(d[t, p])
        if key not in mind or dd < mind[key]:
            mind[key] = dd
    n5 = len(mind)
    if n5 == 0:
        return {"n5": 0, "mean_margin": 0.0, "frac_marg_lt1": 0.0, "frac_robust": 0.0, "exp_lost": 0.0}
    margins = cutoff - np.asarray(list(mind.values()), float)
    ploss = np.clip((delta - margins) / (2.0 * delta), 0.0, 1.0)
    return {
        "n5": n5,
        "mean_margin": float(margins.mean()),
        "frac_marg_lt1": float((margins < delta).mean()),
        "frac_robust": float((margins >= delta).mean()),
        "exp_lost": float(ploss.sum()),
    }


def _selfcheck() -> None:
    """Tiny synthetic assertion so ``python -m tcren.stability`` fails loudly if the logic breaks."""
    from ..structure.model import Atom, Chain, Residue

    def atom(name, el, xyz):
        return Atom(name, el, np.asarray(xyz, float))

    def res(i, resname, aa, atoms):
        return Residue(i, i + 1, "", aa, resname, tuple(atoms))

    # One receptor residue at the origin; two peptide residues. Contact A at 3.0 Å (margin 2.0,
    # robust), contact B at 4.5 Å (margin 0.5 < delta, fragile); a 3rd peptide residue at 6 Å is
    # beyond the cutoff and must not count.
    pep = Chain("C", [res(0, "PHE", "F", [atom("CA", "C", [3.0, 0, 0])]),
                      res(1, "GLY", "G", [atom("CA", "C", [4.5, 0, 0])]),
                      res(2, "ALA", "A", [atom("CA", "C", [6.0, 0, 0])])], chain_type=PEPTIDE_TYPE)
    trb = Chain("B", [res(0, "TYR", "Y", [atom("CA", "C", [0, 0, 0])])], chain_type="TRB")
    rep = contact_stability(Structure("synth", [pep, trb]))
    assert rep.n_contacts == 2, rep
    assert abs(rep.mean_margin - (2.0 + 0.5) / 2) < 1e-9, rep.mean_margin
    assert abs(rep.frac_robust - 0.5) < 1e-9, rep.frac_robust           # A robust, B fragile
    assert abs(rep.frac_marg_lt1 - 0.5) < 1e-9, rep.frac_marg_lt1
    # exp_lost = clip((1-2)/2,0,1) + clip((1-0.5)/2,0,1) = 0 + 0.25
    assert abs(rep.exp_lost - 0.25) < 1e-9, rep.exp_lost
    # numpy reference must agree with whatever backend ran above.
    pep_xyz, pep_res = _interface_atoms(Structure("synth", [pep, trb]), (PEPTIDE_TYPE,))
    tcr_xyz, tcr_res = _interface_atoms(Structure("synth", [pep, trb]), RECEPTOR_TYPES)
    ref = _contact_stability_numpy(pep_xyz, pep_res, tcr_xyz, tcr_res, 5.0, 1.0)
    assert ref["n5"] == 2 and abs(ref["exp_lost"] - 0.25) < 1e-9, ref
    print("tcren.stability self-check OK")


if __name__ == "__main__":
    _selfcheck()
