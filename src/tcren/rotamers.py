"""Rotamer-averaged contacts — a contact map that does not depend on one side-chain guess.

A modelled side chain is a guess, and the contact map reads it as fact. Rotate a Tyr by one χ1 step
and a contact appears or vanishes; the energy moves with it, though nothing about the two residues
changed. On forced poses (AlphaFold peptide swaps, threaded mutants) that is a large part of why a
pairwise contact energy stops discriminating.

Rather than pick a better single rotamer, this samples the χ angles of each interface side chain,
weights each rotamer by its Boltzmann factor under DOPE, and returns a **contact probability** per
residue pair instead of a 0/1 indicator. The score becomes ``sum_ij p_ij * e(a_i, b_j)`` — the same
sum with a softer indicator — via ``score_peptides(..., weights=...)``.

**What is exact and what is not.**

* Rotating every atom beyond Cβ about the Cα–Cβ axis *is* a χ1 change: deeper torsions ride along
  unchanged. The same holds at each subsequent depth. So the geometry is exact, not interpolated.
* Rotamers are enumerated on a uniform grid anchored on the **native** χ (so the input pose is
  always in the set), not drawn from a backbone-dependent library. A Dunbrack-style library would
  give better priors; DOPE supplies the energy here instead.
* Residues are weighted **independently**, each against the rest of the structure held at its input
  conformation (mean field). Two side chains that would have to move together are not coupled.
* ``max_chi`` defaults to 2. χ1 and χ2 carry the reorientation; χ3/χ4 on Arg, Lys, Met and Glu are
  reachable by raising it, at 3× the rotamers per extra angle.

DOPE is used for the weights and nothing else, deliberately: the potential tcren *scores* with must
not also be the one that decided which contacts exist.

**Measured (six crystals, χ1 of every peptide side chain rotated 120° — a deliberately wrong
guess).** The hard 5 Å contact set keeps a Jaccard of only **0.66** against its unperturbed self;
the rotamer-averaged map keeps **0.95**. The energy is the sharper reading: mean ``|ΔΦ|`` falls from
**0.524 to 0.054**, a factor of ten, against interface energies whose own magnitude is 0.4–2.2. On
1ao7 the hard map's error under one wrong rotamer (−0.64) is *larger than the energy itself*
(−0.50), and on 2ckb it flips a +0.43 to +1.66. That is the failure mode this module exists to
remove.

Because most alternative rotamers clash, DOPE separates them by ~200 units and the weights are
sharp; ``temperature`` around 1 is therefore close to a repack, and raising it hedges further. The
default is not tuned for softness — 1.0 measured best on ``|ΔΦ|`` and near-best on Jaccard.
"""

from __future__ import annotations

import numpy as np
import polars as pl

from .structure.model import PEPTIDE_TYPE

#: Rotatable side-chain torsions per residue (one-letter). Ala/Gly have none; Pro's ring is closed.
N_CHI = {
    "S": 1, "C": 1, "T": 1, "V": 1,
    "I": 2, "L": 2, "D": 2, "N": 2, "F": 2, "Y": 2, "W": 2, "H": 2,
    "M": 3, "E": 3, "Q": 3,
    "K": 4, "R": 4,
    "A": 0, "G": 0, "P": 0,
}

#: Greek-letter depth encoded in a PDB atom name's second character (CB -> 0, CG1 -> 1, ...).
_GREEK = {"B": 0, "G": 1, "D": 2, "E": 3, "Z": 4, "H": 5}
_BACKBONE = {"N", "CA", "C", "O", "OXT"}

DEFAULT_STEP = 120.0        #: degrees between sampled rotamers (the g-/t/g+ spacing)
DEFAULT_TEMPERATURE = 1.0   #: Boltzmann temperature in DOPE energy units
#: How far a side chain can sweep about Cβ (Arg, the longest). Used only to reject residue pairs
#: that no rotation could bring into contact.
_MAX_SIDECHAIN_REACH = 7.0


def _depth(name: str) -> int:
    """Greek depth of a side-chain atom name, or ``-1`` for backbone/hydrogen."""
    if name in _BACKBONE or name.startswith("H"):
        return -1
    return _GREEK.get(name[1:2], -1)


def chi_axes(residue) -> list[tuple[int, int, np.ndarray]]:
    """Rotatable torsions of one residue as ``(axis_start, axis_end, moving_atom_indices)``.

    χ_n rotates about the bond from the depth ``n-2`` atom to the depth ``n-1`` atom, moving every
    atom deeper than that (Cα–Cβ for χ1 moves everything from Cγ out). Indices are into
    ``residue.atoms``. Branch points resolve to the alphabetically first name at the depth, which is
    the convention that makes Ile's χ2 run along Cβ–Cγ1 rather than Cβ–Cγ2.
    """
    names = [a.name for a in residue.atoms]
    depth = [_depth(n) for n in names]
    by_depth: dict[int, list[int]] = {}
    for i, d in enumerate(depth):
        if d >= 0:
            by_depth.setdefault(d, []).append(i)
    ca = names.index("CA") if "CA" in names else None
    if ca is None:
        return []

    out = []
    for chi in range(1, N_CHI.get(residue.aa, 0) + 1):
        end = by_depth.get(chi - 1)
        moving = [i for i, d in enumerate(depth) if d >= chi]
        if not end or not moving:
            break
        end_i = min(end, key=lambda i: names[i])
        start_i = ca if chi == 1 else min(by_depth[chi - 2], key=lambda i: names[i])
        out.append((start_i, end_i, np.asarray(moving, dtype=np.int64)))
    return out


def _rotate(coords: np.ndarray, start: np.ndarray, end: np.ndarray,
            moving: np.ndarray, degrees: float) -> np.ndarray:
    """Rotate ``coords[moving]`` about the ``start``->``end`` axis by ``degrees`` (Rodrigues)."""
    axis = end - start
    norm = np.linalg.norm(axis)
    if norm < 1e-9:
        return coords
    k = axis / norm
    theta = np.radians(degrees)
    v = coords[moving] - end
    rotated = (v * np.cos(theta) + np.cross(k, v) * np.sin(theta)
               + np.outer(v @ k, k) * (1 - np.cos(theta)))
    out = coords.copy()
    out[moving] = rotated + end
    return out


def residue_rotamers(residue, max_chi: int = 2, step: float = DEFAULT_STEP) -> np.ndarray:
    """All sampled conformers of one residue as ``(n_rotamers, n_atoms, 3)``.

    The first entry is always the input conformation, so a caller that keeps only the best rotamer
    can never do worse than the pose it was given.
    """
    base = np.asarray([a.coord for a in residue.atoms], dtype=np.float64)
    axes = chi_axes(residue)[:max_chi]
    if not axes:
        return base[None, :, :]

    n_steps = max(int(round(360.0 / step)), 1)
    offsets = [i * step for i in range(1, n_steps)]      # 0 is the conformer we already hold
    conformers = [base]
    for start_i, end_i, moving in axes:
        conformers = conformers + [
            _rotate(conf, conf[start_i], conf[end_i], moving, off)
            for conf in conformers for off in offsets
        ]
    return np.asarray(conformers)


def _dope_weights(coords_by_rotamer: np.ndarray, atom_class: np.ndarray,
                  partner_xyz: np.ndarray, partner_class: np.ndarray,
                  temperature: float) -> np.ndarray:
    """Boltzmann weights over one residue's rotamers, from DOPE against fixed partners."""
    from . import _relax
    from .refine import _dope

    table, _amap, x_start, dx, nbins = _dope()
    n_cls = table.shape[0]
    energies = np.array([
        _relax.interface_energy(np.ascontiguousarray(xyz), atom_class, partner_xyz, partner_class,
                                table, n_cls, nbins, x_start, dx)
        for xyz in coords_by_rotamer
    ], dtype=np.float64)
    e = energies - energies.min()
    w = np.exp(-e / max(temperature, 1e-6))
    return w / w.sum()


def contact_probabilities(structure, interface: str = "tcr_peptide", *, cutoff: float = 5.0,
                          max_chi: int = 2, step: float = DEFAULT_STEP,
                          temperature: float = DEFAULT_TEMPERATURE,
                          shell: float = 12.0) -> pl.DataFrame:
    """Rotamer-averaged contact probability for every residue pair near an interface.

    Each interface side chain is rotated through its χ grid, weighted by ``exp(-E_DOPE / T)`` against
    the rest of the structure, and a pair's probability is the chance that *some* pair of their
    rotamers is within ``cutoff`` — ``p_ij = sum_r sum_s w_r w_s [d(r, s) <= cutoff]``.

    Pairs the input pose does not contact can still acquire probability, and pairs it does contact
    can fall below 1. Both are the point: a contact map built from one modelled rotamer asserts
    certainty it does not have.

    Args:
        structure: chain-typed, annotated structure.
        interface: ``"tcr_peptide"``, ``"tcr_mhc"`` or ``"peptide_mhc"``.
        cutoff: contact distance (Å), matching the hard contact map.
        max_chi: how many χ angles to sample per residue (see the module docstring).
        step: degrees between sampled rotamers.
        temperature: Boltzmann temperature in DOPE units. Larger = flatter weights; as it
            approaches 0 the result converges on the single best rotamer per residue.
        shell: only partner atoms within this distance enter a residue's DOPE weighting.

    Returns:
        ``chain.id.from``, ``residue.index.from``, ``chain.id.to``, ``residue.index.to``,
        ``residue.aa.from``, ``residue.aa.to``, ``dist`` (the input pose's closest heavy-atom
        distance) and ``p`` (the rotamer-averaged contact probability), one row per pair with
        ``p > 0``. ``from`` is the interface's first side (the TCR for ``"tcr_peptide"``), as
        :meth:`~tcren.contactmap.ContactMap.interface` orients it, so a directed potential indexes
        the right way round.

    Raises:
        ValueError: for an unknown ``interface``, or if the structure is not chain-typed.
    """
    from scipy.spatial import cKDTree

    from .contactmap import MHC_TYPES, RECEPTOR_TYPES
    from .refine import _dope

    sides = {"tcr_peptide": (RECEPTOR_TYPES, (PEPTIDE_TYPE,)),
             "tcr_mhc": (RECEPTOR_TYPES, MHC_TYPES),
             "peptide_mhc": ((PEPTIDE_TYPE,), MHC_TYPES)}
    if interface not in sides:
        raise ValueError(f"unknown interface {interface!r}")
    from_types, to_types = sides[interface]

    _table, atom_class_map, _x0, _dx, _nb = _dope()

    def classes(res):
        return np.asarray([atom_class_map.get(f"{res.resname}:{a.name}", -1) for a in res.atoms],
                          dtype=np.int32)

    # `all_xyz` spans every chain-typed residue (the DOPE environment); `on_interface` picks out the
    # ones whose rotamers we sample, keyed by the same global index so the neighbour search can map
    # an atom straight back to a residue.
    all_xyz, all_cls, on_interface = [], [], {}
    for chain in structure.chains:
        if chain.chain_type is None:
            continue
        for res in chain.residues:
            all_xyz.append(np.asarray([a.coord for a in res.atoms], dtype=np.float64))
            all_cls.append(classes(res))
            if chain.chain_type in from_types or chain.chain_type in to_types:
                on_interface[len(all_xyz) - 1] = (chain, res)
    if not on_interface:
        raise ValueError("no chain-typed residues on this interface; run classify_chains first")

    flat_xyz = np.vstack(all_xyz)
    flat_cls = np.concatenate(all_cls)
    owner = np.concatenate([np.full(len(x), i) for i, x in enumerate(all_xyz)])
    tree = cKDTree(flat_xyz)

    def neighbours(idx: int, radius: float) -> np.ndarray:
        """Indices of atoms within ``radius`` of any atom of residue ``idx``, excluding its own."""
        hits = tree.query_ball_point(all_xyz[idx], radius)
        if len(hits) == 0:
            return np.zeros(0, dtype=int)
        # np.concatenate over the per-atom hit arrays, not a generator over every hit: flattening
        # ~1.6M indices one at a time cost more than the whole DOPE evaluation.
        near = np.unique(np.concatenate([np.asarray(h, dtype=int) for h in hits if h]
                                        or [np.zeros(0, dtype=int)]))
        return near[owner[near] != idx]

    # Per-residue rotamer sets and their Boltzmann weights, each against everything else held fixed.
    packed = {}
    for idx, (_chain, res) in on_interface.items():
        confs = residue_rotamers(res, max_chi=max_chi, step=step)
        if len(confs) == 1:
            packed[idx] = (confs, np.ones(1))
            continue
        near = neighbours(idx, shell)
        w = _dope_weights(confs, all_cls[idx], np.ascontiguousarray(flat_xyz[near]),
                          np.ascontiguousarray(flat_cls[near]), temperature)
        packed[idx] = (confs, w)

    # Only pairs that could touch after rotation are worth the O(rotamers x atoms) distance work.
    # A side chain sweeps at most ~7 A about Cbeta, so two residues further apart than the cutoff
    # plus twice that can never contact whatever their rotamers do.
    reach = cutoff + 2 * _MAX_SIDECHAIN_REACH

    rows = []
    for idx_a, (chain_a, res_a) in on_interface.items():
        if chain_a.chain_type not in from_types:
            continue
        ca, wa = packed[idx_a]
        for idx_b in sorted({int(owner[j]) for j in neighbours(idx_a, reach)} & set(on_interface)):
            chain_b, res_b = on_interface[idx_b]
            if chain_b.chain_type not in to_types or chain_a.chain_id == chain_b.chain_id:
                continue
            cb, wb = packed[idx_b]
            # (n_rot_a, n_rot_b) closest-atom distance over every rotamer pair.
            d = np.sqrt(((ca[:, None, :, None, :] - cb[None, :, None, :, :]) ** 2).sum(-1)
                        ).min(axis=(2, 3))
            p = float(wa @ (d <= cutoff).astype(np.float64) @ wb)
            if p <= 0.0:
                continue
            # `from` is the interface's first side, exactly as ContactMap.interface orients it.
            # NOT the (chain.id, residue.index) canonical order all_atom_contacts uses: TCRen is a
            # **directed** TCR->peptide matrix, so orienting the pair by chain letter would index it
            # backwards whenever the peptide chain happens to sort first.
            rows.append({"chain.id.from": chain_a.chain_id, "residue.index.from": res_a.seq_index,
                         "chain.id.to": chain_b.chain_id, "residue.index.to": res_b.seq_index,
                         "residue.aa.from": res_a.aa, "residue.aa.to": res_b.aa,
                         "dist": float(d[0, 0]), "p": p})

    schema = {"chain.id.from": pl.Utf8, "residue.index.from": pl.Int64, "chain.id.to": pl.Utf8,
              "residue.index.to": pl.Int64, "residue.aa.from": pl.Utf8, "residue.aa.to": pl.Utf8,
              "dist": pl.Float64, "p": pl.Float64}
    return pl.DataFrame(rows, schema=schema).sort("p", descending=True)


def soft_energy(structure, potential, interface: str = "tcr_peptide", **kwargs) -> float:
    """Interface energy summed over rotamer-averaged contact probabilities.

    ``sum_ij p_ij * e(a_i, b_j)`` — the same sum as
    :func:`tcren.pipeline._interface_energy`, with the hard 0/1 contact indicator replaced by
    :func:`contact_probabilities`.

    Args:
        structure: chain-typed, annotated structure.
        potential: the pairwise potential (TCRen for TCR:peptide).
        interface: which interface.
        **kwargs: passed to :func:`contact_probabilities`.

    Returns:
        The summed energy.
    """
    from .pipeline import _interface_energy

    p = contact_probabilities(structure, interface, **kwargs)
    if p.is_empty():
        return 0.0
    return _interface_energy(p, potential, weights=p["p"].to_numpy())
