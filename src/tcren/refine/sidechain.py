"""Build a side chain from a backbone and a chi vector.

:func:`tcren.refine.substitute.substitute_residues` re-types a residue in 3D but truncates it at
Cbeta, because the atoms past Cbeta are identity-specific and it has no way to place them. That is
exact for alanine -- alanine's heavy atoms *are* N, CA, C, O and CB -- and understates the reach of
every other type, which is why the poly-alanine reference is the one the shipped scores use. This
module supplies the missing half: given the backbone and a chi vector, place every heavy atom of any
of the nineteen non-glycine types.

**Chemistry is declared, geometry is measured.** :data:`TOPOLOGY` names, per residue type, the atoms
past Cbeta in placement order, the three already-placed atoms whose frame each one hangs off, and
which chi angle drives its torsion. Nothing numeric lives there. The bond lengths, bond angles and
the fixed torsion offsets are fitted from the Canonical2026 crystals by
``scripts/fit_sidechain_geometry.py`` and shipped as ``data/sidechain_geometry.json`` -- the same
arrangement as :data:`tcren.geometry.BOND_CA` and :data:`tcren.refine.substitute.virtual_cb`'s
offset, and for the same reason: a table measured on the set the rest of the library is derived from
needs no external source and cannot drift away from it.

The torsion of an atom is ``chi[k] + offset`` when :data:`TOPOLOGY` names a chi index ``k``, and the
bare ``offset`` when it does not. Both cases are one measured number, so a branch atom (LEU:CD2 off
CD1, ASP:OD2 off OD1) and a ring atom (PHE:CE1, TRP:CH2) are placed by the same rule -- the branch
carries an offset near 120 or 180 degrees that rotates with its chi, the ring one that does not
rotate at all.

Placement is the natural extension reference frame (NeRF): three placed atoms plus a bond, an angle
and a torsion give the fourth, which is the standard forward step of an internal-coordinate chain
and the exact inverse of the dihedral :func:`chi_angles` reads back.
"""
from __future__ import annotations

import json
import math
from functools import lru_cache

import numpy as np

from ..structure.model import Atom, Residue

#: Element of every side-chain atom name, keyed by its first character. Heavy atoms only; the
#: structure model carries no hydrogens and the contact definitions are heavy-atom throughout.
_ELEMENT = {"C": "C", "N": "N", "O": "O", "S": "S"}

#: ``residue -> ((atom, a, b, c, chi), ...)``, in placement order.
#:
#: ``atom`` is placed at a measured bond length from ``c``, a measured angle ``(b, c, atom)`` and a
#: torsion ``(a, b, c, atom)``. ``chi`` is the 1-based chi index driving that torsion, or ``None``
#: when the torsion is fixed. Every one of ``a``, ``b``, ``c`` is either a backbone atom or an atom
#: placed earlier in the same tuple, so one forward pass suffices.
#:
#: Proline is present and is a special case: its CD closes a ring onto the backbone N, which the
#: open chain here does not enforce. A proline built from a chi vector is therefore approximate,
#: and :func:`build_sidechain` returns that approximation without a ring-closure constraint.
TOPOLOGY: dict[str, tuple[tuple[str, str, str, str, int | None], ...]] = {
    "A": (),
    "G": (),
    "S": (("OG",  "N", "CA", "CB", 1),),
    "C": (("SG",  "N", "CA", "CB", 1),),
    "T": (("OG1", "N", "CA", "CB", 1),
          ("CG2", "N", "CA", "CB", 1)),
    "V": (("CG1", "N", "CA", "CB", 1),
          ("CG2", "N", "CA", "CB", 1)),
    "I": (("CG1", "N", "CA", "CB", 1),
          ("CG2", "N", "CA", "CB", 1),
          ("CD1", "CA", "CB", "CG1", 2)),
    "L": (("CG",  "N", "CA", "CB", 1),
          ("CD1", "CA", "CB", "CG", 2),
          ("CD2", "CA", "CB", "CG", 2)),
    "P": (("CG",  "N", "CA", "CB", 1),
          ("CD",  "CA", "CB", "CG", 2)),
    "M": (("CG",  "N", "CA", "CB", 1),
          ("SD",  "CA", "CB", "CG", 2),
          ("CE",  "CB", "CG", "SD", 3)),
    "F": (("CG",  "N", "CA", "CB", 1),
          ("CD1", "CA", "CB", "CG", 2),
          ("CD2", "CA", "CB", "CG", 2),
          ("CE1", "CB", "CG", "CD1", None),
          ("CE2", "CB", "CG", "CD2", None),
          ("CZ",  "CG", "CD1", "CE1", None)),
    "Y": (("CG",  "N", "CA", "CB", 1),
          ("CD1", "CA", "CB", "CG", 2),
          ("CD2", "CA", "CB", "CG", 2),
          ("CE1", "CB", "CG", "CD1", None),
          ("CE2", "CB", "CG", "CD2", None),
          ("CZ",  "CG", "CD1", "CE1", None),
          ("OH",  "CD1", "CE1", "CZ", None)),
    "W": (("CG",  "N", "CA", "CB", 1),
          ("CD1", "CA", "CB", "CG", 2),
          ("CD2", "CA", "CB", "CG", 2),
          ("NE1", "CB", "CG", "CD1", None),
          ("CE2", "CB", "CG", "CD2", None),
          ("CE3", "CB", "CG", "CD2", None),
          ("CZ2", "CG", "CD2", "CE2", None),
          ("CZ3", "CG", "CD2", "CE3", None),
          ("CH2", "CD2", "CE2", "CZ2", None)),
    "H": (("CG",  "N", "CA", "CB", 1),
          ("ND1", "CA", "CB", "CG", 2),
          ("CD2", "CA", "CB", "CG", 2),
          ("CE1", "CB", "CG", "ND1", None),
          ("NE2", "CB", "CG", "CD2", None)),
    "D": (("CG",  "N", "CA", "CB", 1),
          ("OD1", "CA", "CB", "CG", 2),
          ("OD2", "CA", "CB", "CG", 2)),
    "N": (("CG",  "N", "CA", "CB", 1),
          ("OD1", "CA", "CB", "CG", 2),
          ("ND2", "CA", "CB", "CG", 2)),
    "E": (("CG",  "N", "CA", "CB", 1),
          ("CD",  "CA", "CB", "CG", 2),
          ("OE1", "CB", "CG", "CD", 3),
          ("OE2", "CB", "CG", "CD", 3)),
    "Q": (("CG",  "N", "CA", "CB", 1),
          ("CD",  "CA", "CB", "CG", 2),
          ("OE1", "CB", "CG", "CD", 3),
          ("NE2", "CB", "CG", "CD", 3)),
    "K": (("CG",  "N", "CA", "CB", 1),
          ("CD",  "CA", "CB", "CG", 2),
          ("CE",  "CB", "CG", "CD", 3),
          ("NZ",  "CG", "CD", "CE", 4)),
    "R": (("CG",  "N", "CA", "CB", 1),
          ("CD",  "CA", "CB", "CG", 2),
          ("NE",  "CB", "CG", "CD", 3),
          ("CZ",  "CG", "CD", "NE", 4),
          ("NH1", "CD", "NE", "CZ", None),
          ("NH2", "CD", "NE", "CZ", None)),
}

#: How many chi angles each type actually has, from :data:`TOPOLOGY`.
N_CHI: dict[str, int] = {
    aa: max((t[4] for t in rows if t[4] is not None), default=0) for aa, rows in TOPOLOGY.items()
}


def nerf(a: np.ndarray, b: np.ndarray, c: np.ndarray,
         bond: float, angle: float, torsion: float) -> np.ndarray:
    """The fourth atom, from three placed ones plus an internal-coordinate triple.

    Natural extension reference frame. ``angle`` is the ``(b, c, new)`` bond angle and ``torsion``
    the ``(a, b, c, new)`` dihedral, both in radians; ``bond`` is the ``c``-to-new distance in
    angstrom. Exactly inverts the dihedral convention :func:`chi_angles` reads with, which is what
    lets the fitted table be measured by read-back and used by construction.
    """
    bc = c - b
    bc /= np.linalg.norm(bc)
    n = np.cross(b - a, bc)
    n /= np.linalg.norm(n)
    m = np.cross(n, bc)
    d = np.array([-bond * math.cos(angle),
                  bond * math.sin(angle) * math.cos(torsion),
                  bond * math.sin(angle) * math.sin(torsion)])
    return c + d[0] * bc + d[1] * m + d[2] * n


def dihedral(p0: np.ndarray, p1: np.ndarray, p2: np.ndarray, p3: np.ndarray) -> float:
    """The ``(p0, p1, p2, p3)`` dihedral in radians, in the convention :func:`nerf` builds with.

    Exactly inverts :func:`nerf`: feeding this back with the measured bond and angle reproduces
    ``p3`` to 5.6e-14 over 200 random quadruples. The sign of ``b0`` is the whole convention and
    getting it backwards costs a silent 180 degrees, which reads as a plausible 4 angstrom RMSD
    rather than as an error -- :func:`_selfcheck` exists for that reason.
    """
    b0, b1, b2 = p0 - p1, p2 - p1, p3 - p2
    b1 = b1 / np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1
    return math.atan2(np.dot(np.cross(b1, v), w), np.dot(v, w))


def bond_angle(p0: np.ndarray, p1: np.ndarray, p2: np.ndarray) -> float:
    """The ``(p0, p1, p2)`` bond angle in radians."""
    u, v = p0 - p1, p2 - p1
    c = float(np.dot(u, v) / (np.linalg.norm(u) * np.linalg.norm(v)))
    return math.acos(min(1.0, max(-1.0, c)))


@lru_cache(maxsize=1)
def geometry() -> dict[str, dict[str, dict[str, float]]]:
    """The fitted internal-coordinate table, ``aa -> atom -> {bond, angle, offset, ...}``.

    Measured from the 374 Canonical2026 crystals over 176,475 residues by
    ``scripts/fit_sidechain_geometry.py``; angles and offsets are radians. Cached, because it is
    read once per placed atom in the hot path.
    """
    from importlib import resources

    text = resources.files("tcren.data").joinpath("sidechain_geometry.json").read_text()
    return json.loads(text)["geometry"]


def build_sidechain(n: np.ndarray, ca: np.ndarray, c: np.ndarray, cb: np.ndarray,
                    aa: str, chis) -> dict[str, np.ndarray]:
    """Every heavy atom past Cbeta of residue type ``aa`` at chi vector ``chis``.

    Args:
        n, ca, c, cb: backbone coordinates, and the Cbeta the side chain grows from. Build Cbeta
            with :func:`tcren.refine.substitute.virtual_cb` when the residue does not carry one.
        aa: one-letter code. ``"A"`` and ``"G"`` return ``{}``, correctly -- alanine's heavy atoms
            are exactly the backbone plus Cbeta and glycine has neither.
        chis: chi angles in radians, at least :data:`N_CHI` of them. Extra entries are ignored, so
            one grid can be handed to every type.

    Returns:
        ``atom name -> coordinate``. Empty for A and G. Proline is approximate:
        the open-chain construction does not constrain CD to close onto backbone N.

    Raises:
        KeyError: for a non-standard one-letter code.
        ValueError: if ``chis`` is shorter than the type needs.
    """
    rows = TOPOLOGY[aa]
    if not rows:
        return {}
    need = N_CHI[aa]
    chis = list(chis)
    if len(chis) < need:
        raise ValueError(f"{aa} needs {need} chi angles, got {len(chis)}")
    g = geometry()[aa]
    xyz: dict[str, np.ndarray] = {"N": n, "CA": ca, "C": c, "CB": cb}
    for atom, a, b, cc, chi in rows:
        p = g[atom]
        torsion = p["offset"] + (chis[chi - 1] if chi is not None else 0.0)
        xyz[atom] = nerf(xyz[a], xyz[b], xyz[cc], p["bond"], p["angle"], torsion)
    return {k: v for k, v in xyz.items() if k not in ("N", "CA", "C", "CB")}


def chi_angles(residue: Residue) -> list[float]:
    """The chi angles a residue as deposited actually adopts, in radians.

    Read through the same frames :func:`build_sidechain` builds with, so
    ``build_sidechain(..., chi_angles(r))`` reconstructs ``r`` up to the fitted bond lengths and
    angles. Returns ``[]`` for a type with no chi, and stops at the first chi whose atoms are
    missing -- a truncated side chain yields the chi angles it can still support.
    """
    xyz = {a.name: np.asarray(a.coord, float) for a in residue.atoms}
    out: list[float] = []
    for k in range(1, N_CHI.get(residue.aa, 0) + 1):
        hit = next(((atom, a, b, c) for atom, a, b, c, chi in TOPOLOGY[residue.aa] if chi == k),
                   None)
        if hit is None or not set(hit) <= xyz.keys():
            break
        atom, a, b, c = hit
        out.append(dihedral(xyz[a], xyz[b], xyz[c], xyz[atom]))
    return out


def rebuild_residue(residue: Residue, aa: str, chis) -> Residue:
    """``residue`` re-typed to ``aa`` at chi vector ``chis``, with a full side chain.

    The backbone is untouched; Cbeta is kept when present and built by ideal geometry when not, the
    same rule :func:`tcren.refine.substitute._rewrite_residue` applies. Unlike that function this
    one does **not** stop at Cbeta, so the result can be scored at atom level.
    """
    from .substitute import _KEEP, _ONE_TO_THREE, virtual_cb

    resname = _ONE_TO_THREE.get(aa)
    if resname is None:
        raise ValueError(f"non-standard amino acid {aa!r}")
    bb = {a.name: np.asarray(a.coord, float) for a in residue.atoms}
    if not {"N", "CA", "C"} <= bb.keys():
        raise ValueError(f"residue {residue.seq_index} has no complete backbone")
    keep = _KEEP - ({"CB"} if aa == "G" else set())
    atoms = tuple(a for a in residue.atoms if a.name in keep)
    if aa == "G":
        return Residue(residue.seq_index, residue.pdb_index, residue.insertion_code,
                       aa, resname, atoms)
    cb = bb.get("CB")
    if cb is None:
        cb = virtual_cb(bb["N"], bb["CA"], bb["C"])
        atoms += (Atom(name="CB", element="C", coord=cb),)
    built = build_sidechain(bb["N"], bb["CA"], bb["C"], cb, aa, chis)
    atoms += tuple(Atom(name=k, element=_ELEMENT[k[0]], coord=v) for k, v in built.items())
    return Residue(residue.seq_index, residue.pdb_index, residue.insertion_code,
                   aa, resname, atoms)


def _selfcheck() -> None:
    """Assert the two directions invert each other, on random points and on a real geometry.

    Run as ``python -m tcren.refine.sidechain``. The first assertion is the one that matters: a
    sign error in :func:`dihedral` leaves every bond length and bond angle correct and every atom
    180 degrees around its parent, which reads back as a plausible-looking 4 angstrom RMSD instead
    of an exception. Measured on 120 Canonical2026 crystals the round trip is 0.0784 angstrom mean
    over 54,954 residues; anything near 4 means the convention has drifted again.
    """
    rng = np.random.default_rng(0)
    for _ in range(200):
        a, b, c, d = rng.normal(size=(4, 3)) * 3.0
        back = nerf(a, b, c, float(np.linalg.norm(d - c)), bond_angle(b, c, d),
                    dihedral(a, b, c, d))
        assert np.linalg.norm(d - back) < 1e-9, "nerf and dihedral disagree"

    # A tryptophan built at a chi vector reads that vector back.
    n, ca, c = np.array([0.0, 1.4, 0.0]), np.zeros(3), np.array([1.5, 0.0, 0.0])
    from .substitute import virtual_cb

    cb = virtual_cb(n, ca, c)
    want = [-1.1, 0.7]
    built = build_sidechain(n, ca, c, cb, "W", want)
    assert set(built) == {t[0] for t in TOPOLOGY["W"]}, "W built the wrong atom set"
    xyz = {"N": n, "CA": ca, "C": c, "CB": cb, **built}
    got = []
    for k in range(1, N_CHI["W"] + 1):
        at, a, b, cc = next((r[:4] for r in TOPOLOGY["W"] if r[4] == k))
        got.append(dihedral(xyz[a], xyz[b], xyz[cc], xyz[at]))
    assert np.allclose(got, want, atol=1e-9), f"chi round trip failed: {got} != {want}"

    assert build_sidechain(n, ca, c, cb, "A", []) == {}, "alanine has nothing past CB"
    assert build_sidechain(n, ca, c, cb, "G", []) == {}, "glycine has nothing past CA"
    print("sidechain: ok")


if __name__ == "__main__":
    _selfcheck()
