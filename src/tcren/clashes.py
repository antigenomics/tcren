"""Steric-clash detection at the peptide interface.

:func:`interface_clashes` counts heavy-atom van der Waals overlaps between the peptide chain and its
TCR/MHC partners --- the signature of a forced or wrong-register pose, e.g. an AlphaFold / TCRmodel
peptide swap that seats the peptide non-physically. A *clash* is a non-bonded heavy-atom pair whose
separation is shorter than the sum of their Bondi vdW radii by more than a tolerance (Molprobity uses
``0.4`` Å for a "bad" clash; ``0.6`` Å here marks a *severe* one).

This is a structure-quality check: a generated complex with a heavy clash burden is geometrically
non-physical, so its contact energy is read off a distorted interface (see
:mod:`tcren.refine.register` for the register-specific diagnostic and correction). numpy-only --- no
scipy, no compiled kernel.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .structure.model import PEPTIDE_TYPE, Structure

#: Bondi van der Waals radii (Å) by element symbol; ``_DEFAULT_RADIUS`` covers anything unlisted.
BONDI_RADII: dict[str, float] = {
    "H": 1.20, "C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80, "P": 1.80,
    "F": 1.47, "CL": 1.75, "BR": 1.85, "I": 1.98, "SE": 1.90,
}
_DEFAULT_RADIUS = 1.70
_MAX_RADIUS = max(BONDI_RADII.values())


@dataclass(frozen=True, slots=True)
class ClashPair:
    """A single peptide↔partner clashing residue pair."""

    peptide_residue: str  # e.g. "PHE5"
    partner_residue: str  # e.g. "TYR84"
    partner_chain_type: str  # TRA / TRB / MHCa / …
    overlap: float  # Å the two atoms interpenetrate (vdw_i + vdw_j − distance)


@dataclass(frozen=True, slots=True)
class ClashReport:
    """Interface steric-clash summary for a chain-typed complex.

    Attributes:
        n_clashes: Heavy-atom pairs overlapping by more than ``tolerance``.
        n_severe: Subset overlapping by more than ``severe``.
        max_overlap: Largest single overlap (Å); ``0.0`` if clash-free.
        clash_score: Sum of all overlaps > ``tolerance`` (Å) — a total-burden measure.
        by_partner: Clash count per partner ``chain_type`` (TRA/TRB/MHCa/…).
        worst: Up to ``top`` worst clashing residue pairs, largest overlap first.
        n_peptide_atoms: Heavy atoms on the peptide chain (denominator context).
    """

    n_clashes: int
    n_severe: int
    max_overlap: float
    clash_score: float
    by_partner: dict[str, int]
    worst: tuple[ClashPair, ...]
    n_peptide_atoms: int

    @property
    def clashing(self) -> bool:
        """True if any pair overlaps by more than the tolerance."""
        return self.n_clashes > 0


def _radius(element: str) -> float:
    return BONDI_RADII.get(element.strip().upper(), _DEFAULT_RADIUS)


def _heavy_atoms(residues, chain_label: str):
    """Yield ``(coord, radius, residue_label, chain_label)`` for non-hydrogen atoms."""
    for res in residues:
        label = f"{res.resname}{res.pdb_index}{res.insertion_code}".strip()
        for a in res.atoms:
            if a.element.strip().upper() == "H":
                continue
            yield a.coord, _radius(a.element), label, chain_label


def interface_clashes(
    structure: Structure,
    *,
    tolerance: float = 0.4,
    severe: float = 0.6,
    top: int = 8,
) -> ClashReport:
    """Detect steric clashes between the peptide chain and its partners.

    Two heavy atoms on different chains clash when their separation is shorter than the sum of their
    Bondi vdW radii by more than ``tolerance`` Å. Only peptide↔partner pairs are examined (partners =
    every non-peptide chain); intra-chain and peptide-internal pairs are ignored, so no bonded pair is
    ever counted.

    Args:
        structure: A chain-typed complex with a peptide chain (``chain_type == 'PEPTIDE'``).
        tolerance: vdW overlap (Å) above which a pair counts as a clash (Molprobity ``0.4``).
        severe: overlap (Å) above which a clash is also counted as *severe*.
        top: How many worst residue pairs to return.

    Returns:
        A :class:`ClashReport`.

    Raises:
        ValueError: If the structure has no peptide chain.
    """
    pep = next((c for c in structure.chains if c.chain_type == PEPTIDE_TYPE), None)
    if pep is None:
        raise ValueError(f"no peptide chain (chain_type == {PEPTIDE_TYPE!r}) in {structure.pdb_id!r}")

    pep_xyz, pep_rad, pep_lab = [], [], []
    for coord, rad, label, _ in _heavy_atoms(pep.residues, pep.chain_id):
        pep_xyz.append(coord); pep_rad.append(rad); pep_lab.append(label)
    if not pep_xyz:
        raise ValueError(f"peptide chain of {structure.pdb_id!r} has no heavy atoms")
    P = np.asarray(pep_xyz, float)
    Pr = np.asarray(pep_rad, float)

    n = n_sev = 0
    max_ov = 0.0
    score = 0.0
    by_partner: dict[str, int] = {}
    pairs: list[ClashPair] = []
    # A clash needs distance < r_i + r_j − tol ≤ 2·max_radius − tol; prefilter partner atoms to that shell.
    shell = 2.0 * _MAX_RADIUS - tolerance

    for chain in structure.chains:
        if chain is pep:
            continue
        ctype = chain.chain_type or chain.chain_id
        qx, qr, ql = [], [], []
        for coord, rad, label, _ in _heavy_atoms(chain.residues, chain.chain_id):
            qx.append(coord); qr.append(rad); ql.append(label)
        if not qx:
            continue
        Q = np.asarray(qx, float)
        Qr = np.asarray(qr, float)
        # pairwise distances, peptide (rows) × partner (cols)
        d = np.sqrt(((P[:, None, :] - Q[None, :, :]) ** 2).sum(-1))
        overlap = (Pr[:, None] + Qr[None, :]) - d
        hit = overlap > tolerance
        if not hit.any():
            continue
        cnt = int(hit.sum())
        n += cnt
        by_partner[ctype] = by_partner.get(ctype, 0) + cnt
        n_sev += int((overlap > severe).sum())
        score += float(overlap[hit].sum())
        max_ov = max(max_ov, float(overlap[hit].max()))
        for i, j in zip(*np.where(hit)):
            pairs.append(ClashPair(pep_lab[i], ql[j], ctype, float(overlap[i, j])))

    pairs.sort(key=lambda p: p.overlap, reverse=True)
    return ClashReport(
        n_clashes=n, n_severe=n_sev, max_overlap=max_ov, clash_score=score,
        by_partner=by_partner, worst=tuple(pairs[:top]), n_peptide_atoms=len(pep_xyz),
    )


def has_clashes(structure: Structure, *, tolerance: float = 0.4) -> bool:
    """Convenience predicate: does the peptide interface have any clash beyond ``tolerance``?"""
    return interface_clashes(structure, tolerance=tolerance).clashing


def _selfcheck() -> None:
    """Tiny synthetic assertion so ``python -m tcren.clashes`` fails loudly if the logic breaks."""
    from .structure.model import Atom, Chain, Residue, Structure

    def atom(name, el, xyz):
        return Atom(name, el, np.asarray(xyz, float))

    def res(i, resname, aa, atoms):
        return Residue(i, i + 1, "", aa, resname, tuple(atoms))

    # One peptide C atom at origin; one MHC C atom 2.0 Å away → overlap 1.70+1.70−2.0 = 1.40 (severe).
    pep = Chain("C", [res(0, "PHE", "F", [atom("CZ", "C", [0, 0, 0])])], chain_type=PEPTIDE_TYPE)
    mhc = Chain("D", [res(0, "TYR", "Y", [atom("OH", "O", [2.0, 0, 0])])], chain_type="MHCa")
    s = Structure("synth", [pep, mhc])
    rep = interface_clashes(s)
    assert rep.n_clashes == 1, rep
    assert rep.n_severe == 1, rep
    assert abs(rep.max_overlap - (1.70 + 1.52 - 2.0)) < 1e-9, rep.max_overlap
    assert rep.by_partner == {"MHCa": 1}, rep.by_partner
    # Move the partner far away → clash-free.
    mhc2 = Chain("D", [res(0, "TYR", "Y", [atom("OH", "O", [50.0, 0, 0])])], chain_type="MHCa")
    assert not interface_clashes(Structure("synth", [pep, mhc2])).clashing
    print("tcren.clashes self-check OK")


if __name__ == "__main__":
    _selfcheck()
