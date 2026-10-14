"""The amino-acid background taken in 3D: no contact, contact, or collision.

:func:`tcren.energetics.mutation.smoothed_reference` is a *virtual* scan. It holds the contact map
fixed and re-indexes the potential, so ``phi_i(a) = sum_j sigma_ij e(a, y_j)`` charges every
candidate for whatever contacts the deposited side chain happened to make -- a contact an arginine
reaches across is still counted after that arginine is notionally an alanine. Over the 374 reference
crystals only 54 % of 5 angstrom TCR:peptide residue pairs have both side chains in range at all,
so that is not a small correction.

This module takes the same background structurally. Every candidate residue is *built* at every
rotamer of a chi grid (:mod:`tcren.refine.sidechain`), weighted by DOPE against the environment, and
each (position, candidate, partner) resolves to one of three outcomes rather than one number:

no contact
    No surviving rotamer of that candidate comes within ``cutoff`` of the partner. ``p = 0``, and
    the pair contributes nothing.
contact
    Some surviving rotamer does. ``p`` is the Boltzmann-weighted fraction that do, and the pair
    contributes ``p * e(a, y_j)``.
collision
    *Every* rotamer of that candidate overlaps the environment beyond ``clash_tol``. The candidate
    is **inadmissible** and leaves the background sum entirely, its prior mass renormalised over the
    survivors. It must not enter at its least-bad rotamer, which is what a plain
    ``w / w.sum()`` over all rotamers does -- that reports an impossible substitution as the
    least-bad one and lets it set the reference free energy.

The reference free energy is therefore taken over the admissible set ``A_i``:

.. math:: F_i = -\\beta^{-1}\\,\\ln \\sum_{a \\in A_i} \\tilde p_a\\, e^{-\\beta \\varphi_i(a)}

with :math:`\\tilde p_a` the background renormalised over :math:`A_i`. ``|A_i|`` is itself a
read-out: a position where seventeen of twenty residues collide is a different statement about the
interface from one where all twenty fit.

The work is positions x 20 candidates x ``n_grid ** n_chi`` rotamers, each needing a DOPE energy
against the environment and a heavy-atom distance to every partner, so it runs in one batched call
into :func:`tcren._relax.substitution_scan` per structure -- never one call per residue, and never a
Python pool around it. :func:`scan_numpy` is the reference implementation the kernel is checked
against, in the same relation ``clashes._clash_pairs_numpy`` stands in to the native clash kernel.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from functools import lru_cache

import numpy as np

from ..refine.sidechain import N_CHI, TOPOLOGY, geometry, nerf

#: The twenty one-letter codes, in the order every ``(..., 20)`` axis below is indexed by.
AA: tuple[str, ...] = tuple("ACDEFGHIKLMNPQRSTVWY")

#: Degrees between sampled chi values, from -60. Chosen with :data:`DEFAULT_CLASH_TOL` against one
#: invariant -- a residue the crystal deposited must be admissible at its own position -- by
#: ``scripts/calibrate_background_scan.py`` over 400 deposited peptide residues in 40 Native2026
#: crystals. At 60 degrees and 0.8 angstrom that reads **0.9725**, mean |A| 17.26 of 20, at 23 ms
#: a structure. Halving the step again does **not** help (0.9700 at 30 degrees, 15x the cost), so
#: the residual 2.75 % is strained crystal geometry rather than a coarse grid, and tightening
#: further only buys time. 120 degrees is the g-/t/g+ spacing and reads 0.9300.
DEFAULT_STEP = 60.0

#: Overlap, angstrom, beyond which two heavy atoms are called a collision. Looser than
#: :func:`tcren.clashes.has_clashes`\' 0.4, because that one judges a pose whose side chains are
#: where the crystal put them and this one judges a side chain built from ideal internal
#: coordinates on a grid. See :data:`DEFAULT_STEP` for the measurement that picked it.
DEFAULT_CLASH_TOL = 0.8

#: Boltzmann temperature in DOPE units, as :mod:`tcren.energetics.rotamers` uses.
DEFAULT_TEMPERATURE = 1.0

#: Van der Waals radii, angstrom, keyed by element -- the same table :func:`tcren.clashes._radius`
#: reads, so a collision here and a collision there are the same event.
_RADIUS = {"C": 1.7, "N": 1.55, "O": 1.52, "S": 1.8}

_SLOT = {"N": 0, "CA": 1, "C": 2, "CB": 3}


@dataclass(frozen=True)
class ScanResult:
    """What one structure's scan returns.

    Attributes:
        p: ``(n_pairs, 20)`` contact probability of each candidate at ``pair_i`` against the partner
            residue ``pair_j``. Zero for a no-contact and for an inadmissible candidate alike --
            read ``admissible`` to tell those apart.
        admissible: ``(n_positions, 20)`` bool. False when every rotamer collided.
        n_clash: ``(n_positions, 20)`` how many rotamers collided.
        n_rotamer: ``(n_positions, 20)`` how many were enumerated. ``n_clash == n_rotamer`` is the
            collision case.
        pair_i, pair_j: the pair index arrays the ``p`` rows are keyed by.
    """

    p: np.ndarray
    admissible: np.ndarray
    n_clash: np.ndarray
    n_rotamer: np.ndarray
    pair_i: np.ndarray
    pair_j: np.ndarray

    def admissible_fraction(self) -> np.ndarray:
        """``|A_i| / 20`` per position -- how much of the background that position can hold."""
        return self.admissible.mean(axis=1)


@lru_cache(maxsize=1)
def _tables(step: float = DEFAULT_STEP):
    """The topology and fitted geometry packed dense over the twenty types, for the kernel.

    ``(t_a, t_b, t_c, t_chi, t_bond, t_ang, t_off, t_cls, t_rad, n_atom, n_chi, chi_vals)``, every
    ``(20, A)`` with ``A`` the widest side chain (nine, tryptophan). Cached: it is the same for
    every structure and rebuilding it per call would dominate a small scan.
    """
    from ..refine import _dope
    from ..refine.substitute import _ONE_TO_THREE

    _table, atom_class, _x0, _dx, _nb = _dope()
    g = geometry()
    A = max(len(TOPOLOGY[a]) for a in AA)
    t_a = np.zeros((20, A), np.int32); t_b = np.zeros((20, A), np.int32)
    t_c = np.zeros((20, A), np.int32); t_chi = np.zeros((20, A), np.int32)
    t_bond = np.zeros((20, A)); t_ang = np.zeros((20, A)); t_off = np.zeros((20, A))
    t_cls = np.full((20, A), -1, np.int32); t_rad = np.zeros((20, A))
    n_atom = np.zeros(20, np.int32); n_chi = np.zeros(20, np.int32)
    for i, aa in enumerate(AA):
        rows = TOPOLOGY[aa]
        n_atom[i] = len(rows); n_chi[i] = N_CHI[aa]
        order = {name: 4 + k for k, (name, *_rest) in enumerate(rows)}
        for k, (name, a, b, c, chi) in enumerate(rows):
            p = g[aa][name]
            t_a[i, k] = _SLOT.get(a, order.get(a, -1))
            t_b[i, k] = _SLOT.get(b, order.get(b, -1))
            t_c[i, k] = _SLOT.get(c, order.get(c, -1))
            t_chi[i, k] = 0 if chi is None else chi
            t_bond[i, k] = p["bond"]; t_ang[i, k] = p["angle"]; t_off[i, k] = p["offset"]
            t_cls[i, k] = atom_class.get(f"{_ONE_TO_THREE[aa]}:{name}", -1)
            t_rad[i, k] = _RADIUS[name[0]]
    # Staggered, not eclipsed: at the 120 degree default this is g-, g+, t (-60, 60, 180), the
    # three wells a chi actually occupies. Starting the grid at 0 would sample the eclipsed
    # conformation every type avoids and miss all three.
    n = max(1, int(round(360.0 / step)))
    chi_vals = np.radians(-60.0 + np.arange(n) * step)
    return (t_a, t_b, t_c, t_chi, t_bond, t_ang, t_off, t_cls, t_rad, n_atom, n_chi, chi_vals)


#: How far a side chain can sweep about Cbeta. Two residues further apart than ``cutoff`` plus twice
#: this can never contact whatever their rotamers do -- the bound
#: :func:`tcren.energetics.rotamers.contact_probabilities` prunes with.
MAX_REACH = 7.0


def _interface_sides(interface: str):
    from ..contactmap import MHC_TYPES, RECEPTOR_TYPES
    from ..structure.model import PEPTIDE_TYPE

    sides = {"tcr_peptide": (RECEPTOR_TYPES, (PEPTIDE_TYPE,)),
             "tcr_mhc": (RECEPTOR_TYPES, MHC_TYPES),
             "peptide_mhc": ((PEPTIDE_TYPE,), MHC_TYPES)}
    if interface not in sides:
        raise ValueError(f"unknown interface {interface!r}, expected one of {sorted(sides)}")
    return sides[interface]


def _pack(structure, interface: str, side: str, cutoff: float, shell: float):
    """Flatten one structure into the arrays the kernel reads.

    Returns ``(inputs, positions, partners)`` where ``positions`` and ``partners`` are the
    ``(chain, residue)`` lists the two output axes are keyed by. The varying side is ``side``:
    ``"from"`` scans the interface's first side (the receptor for ``tcr_peptide``) and ``"to"`` the
    second, which is the difference between referencing the receptor and referencing the peptide.
    """
    from scipy.spatial import cKDTree

    from ..refine import _dope

    from_types, to_types = _interface_sides(interface)
    if side not in ("from", "to"):
        raise ValueError(f"side must be 'from' or 'to', got {side!r}")
    scan_types = from_types if side == "from" else to_types
    par_types = to_types if side == "from" else from_types

    _table, atom_class, _x0, _dx, _nb = _dope()

    env_xyz, env_cls, env_rad, env_res = [], [], [], []
    positions, partners = [], []
    par_atoms: list[np.ndarray] = []

    for chain in structure.chains:
        if chain.chain_type is None:
            continue
        for res in chain.residues:
            names = {a.name for a in res.atoms}
            scan_here = chain.chain_type in scan_types and {"N", "CA", "C"} <= names
            k = -1
            if scan_here and res.aa in TOPOLOGY:
                k = len(positions)
                positions.append((chain, res))
            for a in res.atoms:
                env_xyz.append(np.asarray(a.coord, float))
                env_cls.append(atom_class.get(f"{res.resname}:{a.name}", -1))
                env_rad.append(_RADIUS.get(a.element.strip().upper()[:1], 1.7))
                env_res.append(k)
            if chain.chain_type in par_types:
                partners.append((chain, res))
                par_atoms.append(np.asarray([a.coord for a in res.atoms], float))
    if not positions or not partners:
        raise ValueError(f"interface {interface!r} side {side!r} has no residues to scan")

    env_xyz = np.ascontiguousarray(env_xyz)
    env_cls = np.asarray(env_cls, np.int32)
    env_rad = np.asarray(env_rad, float)
    env_res = np.asarray(env_res, np.int32)

    # Backbone frame per scan position; Cbeta is built when the residue has none, the same rule
    # `substitute._rewrite_residue` applies, so glycine is scanned like anything else.
    from ..refine.substitute import virtual_cb

    bb = np.empty((len(positions), 4, 3))
    for i, (_c, res) in enumerate(positions):
        xyz = {a.name: np.asarray(a.coord, float) for a in res.atoms}
        cb = xyz.get("CB")
        if cb is None:
            cb = virtual_cb(xyz["N"], xyz["CA"], xyz["C"])
        bb[i] = [xyz["N"], xyz["CA"], xyz["C"], cb]

    par_ptr = np.zeros(len(partners) + 1, np.int32)
    for j, a in enumerate(par_atoms):
        par_ptr[j + 1] = par_ptr[j] + len(a)
    par_xyz = np.ascontiguousarray(np.vstack(par_atoms))

    # Pairs worth reporting: a candidate can only reach a partner within cutoff + 2 * MAX_REACH of
    # its own Cbeta. Everything beyond that is a no-contact for every rotamer of every candidate.
    # Both neighbour searches are ONE tree query over every Cbeta at once, not one per position:
    # cKDTree takes an (n, 3) array and does the descent in C. The per-position Python loop this
    # replaces was the whole of `_pack`'s cost.
    reach = cutoff + 2.0 * MAX_REACH
    par_owner = np.concatenate([np.full(len(a), j) for j, a in enumerate(par_atoms)])
    cb_all = np.ascontiguousarray(bb[:, 3])
    par_hits = cKDTree(par_xyz).query_ball_point(cb_all, reach)
    pair_i, pair_j = [], []
    for i, hits in enumerate(par_hits):
        if not len(hits):
            continue
        for j in np.unique(par_owner[np.asarray(hits, int)]):
            pair_i.append(i); pair_j.append(int(j))

    # Environment neighbour list, CSR. Off Cbeta at the DOPE shell plus the reach, so a rotamer
    # that swings out still sees everything it could touch.
    env_hits = cKDTree(env_xyz).query_ball_point(cb_all, shell + MAX_REACH)
    idx = [np.asarray(h, np.int32) for h in env_hits]
    ptr = np.zeros(len(positions) + 1, np.int32)
    ptr[1:] = np.cumsum([len(h) for h in idx])
    env_idx = np.concatenate(idx).astype(np.int32) if idx else np.zeros(0, np.int32)

    return (dict(bb=np.ascontiguousarray(bb.reshape(-1, 3)).reshape(len(positions), 4, 3),
                 env_xyz=env_xyz, env_cls=env_cls, env_rad=env_rad, env_res=env_res,
                 env_ptr=ptr, env_idx=env_idx, par_xyz=par_xyz, par_ptr=par_ptr,
                 pair_i=np.asarray(pair_i, np.int32), pair_j=np.asarray(pair_j, np.int32)),
            positions, partners)


def scan_structure(structure, interface: str = "tcr_peptide", side: str = "to", *,
                   cutoff: float = 5.0, step: float = DEFAULT_STEP,
                   clash_tol: float = DEFAULT_CLASH_TOL,
                   temperature: float = DEFAULT_TEMPERATURE, shell: float = 12.0,
                   n_threads: int = 0, engine: str = "native"):
    """Scan every amino acid at every scannable position of one interface, in 3D.

    One batched call into the native kernel per structure. Never call this per residue and never
    wrap it in a process pool -- the kernel releases the GIL and threads over positions itself.

    Args:
        structure: chain-typed, and annotated if the interface names the MHC.
        interface: ``"tcr_peptide"``, ``"tcr_mhc"`` or ``"peptide_mhc"``.
        side: which side varies. ``"to"`` (default) scans the second side -- the peptide for
            ``tcr_peptide``, which is the peptide-ranking direction; ``"from"`` scans the first.
        cutoff: contact distance, angstrom, matching the hard contact map.
        step: degrees between sampled chi values.
        clash_tol: overlap beyond which two heavy atoms collide.
        temperature: Boltzmann temperature in DOPE units for the rotamer average.
        shell: DOPE environment radius about Cbeta.
        n_threads: kernel threads; 0 (default) uses every core.
        engine: ``"native"`` or ``"numpy"`` -- the reference implementation, for the agreement test.

    Returns:
        :class:`ScanResult`, plus ``positions`` and ``partners`` as ``(chain, residue)`` lists.
    """
    packed, positions, partners = _pack(structure, interface, side, cutoff, shell)
    t = _tables(step)
    from ..refine import _dope

    table, _amap, x_start, dx, nbins = _dope()
    kw = dict(t_a=t[0], t_b=t[1], t_c=t[2], t_chi=t[3], t_bond=t[4], t_ang=t[5], t_off=t[6],
              t_cls=t[7], t_rad=t[8], n_atom=t[9], n_chi=t[10], chi_vals=t[11],
              dope_table=table, n_cls=table.shape[0], n_bins=nbins, x_start=x_start, dx=dx,
              cutoff=cutoff, clash_tol=clash_tol, temperature=temperature)
    if engine == "numpy":
        out = scan_numpy(**packed, **kw)
    elif engine == "native":
        from .. import _relax

        out = _relax.substitution_scan(**packed, **kw, n_threads=n_threads)
    else:
        raise ValueError(f"engine must be 'native' or 'numpy', got {engine!r}")
    return (ScanResult(p=np.asarray(out["p"]), admissible=np.asarray(out["admissible"], bool),
                       n_clash=np.asarray(out["n_clash"]), n_rotamer=np.asarray(out["n_rotamer"]),
                       pair_i=packed["pair_i"], pair_j=packed["pair_j"]),
            positions, partners)


def scan_numpy(*, bb, t_a, t_b, t_c, t_chi, t_bond, t_ang, t_off, t_cls, t_rad, n_atom, n_chi,
               chi_vals, env_xyz, env_cls, env_rad, env_res, env_ptr, env_idx, par_xyz, par_ptr,
               pair_i, pair_j, dope_table, n_cls, n_bins, x_start, dx, cutoff, clash_tol,
               temperature):
    """Pure-numpy reference for :func:`tcren._relax.substitution_scan`.

    Same signature, same arithmetic, no threads. Kept so the native kernel has something to be
    checked against elementwise -- the relation ``clashes._clash_pairs_numpy`` stands in to
    ``_geom.interface_clashes``. Far too slow for a benchmark; it exists to be right, not fast.
    """
    P = bb.shape[0]
    K = len(pair_i)
    d_max = x_start + (n_bins - 1) * dx
    p_out = np.zeros((K, 20))
    adm = np.zeros((P, 20), np.uint8)
    n_clash = np.zeros((P, 20), np.int32)
    n_rot_out = np.zeros((P, 20), np.int32)
    pairs_of_pos: dict[int, list[int]] = {}
    for k in range(K):
        pairs_of_pos.setdefault(int(pair_i[k]), []).append(k)

    for pos in range(P):
        near = env_idx[env_ptr[pos]:env_ptr[pos + 1]]
        near = near[env_res[near] != pos]
        e_xyz, e_cls, e_rad = env_xyz[near], env_cls[near], env_rad[near]
        mypairs = pairs_of_pos.get(pos, [])
        for aa in range(20):
            na, nc = int(n_atom[aa]), int(n_chi[aa])
            n_rot = int(chi_vals.size) ** nc
            n_rot_out[pos, aa] = n_rot
            if na == 0:
                adm[pos, aa] = 1
                continue
            energy = np.zeros(n_rot)
            clashed = np.zeros(n_rot, bool)
            hit = np.zeros((n_rot, len(mypairs)), bool)
            for r in range(n_rot):
                q, chi = r, [0.0] * 4
                for k in range(nc):
                    chi[k] = float(chi_vals[q % chi_vals.size]); q //= chi_vals.size
                slot = [bb[pos, 0], bb[pos, 1], bb[pos, 2], bb[pos, 3]]
                for t in range(na):
                    tor = t_off[aa, t] + (chi[int(t_chi[aa, t]) - 1] if t_chi[aa, t] else 0.0)
                    slot.append(nerf(slot[int(t_a[aa, t])], slot[int(t_b[aa, t])],
                                     slot[int(t_c[aa, t])], t_bond[aa, t], t_ang[aa, t], tor))
                built = np.asarray(slot[4:])
                d = np.linalg.norm(built[:, None, :] - e_xyz[None, :, :], axis=-1)
                clashed[r] = bool(((t_rad[aa, :na, None] + e_rad[None, :] - d) > clash_tol).any())
                ok = (d < d_max) & (e_cls[None, :] >= 0) & (t_cls[aa, :na, None] >= 0)
                if ok.any():
                    ai, bi = np.nonzero(ok)
                    bins = np.clip(((d[ai, bi] - x_start) / dx).astype(int), 0, n_bins - 1)
                    # float64 accumulation, matching the kernel: the table is float32 and
                    # summing it in float32 costs 5e-7 against the C++ double accumulator.
                    energy[r] = float(
                        dope_table[t_cls[aa, ai], e_cls[bi], bins].astype(np.float64).sum())
                for qi, k in enumerate(mypairs):
                    j = int(pair_j[k])
                    y = par_xyz[par_ptr[j]:par_ptr[j + 1]]
                    side = np.vstack([bb[pos, 3][None, :], built])
                    hit[r, qi] = bool(
                        (np.linalg.norm(side[:, None, :] - y[None, :, :], axis=-1) < cutoff).any())
            n_clash[pos, aa] = int(clashed.sum())
            live = ~clashed
            adm[pos, aa] = 1 if live.any() else 0
            if not live.any() or not mypairs:
                continue
            w = np.zeros(n_rot)
            w[live] = np.exp(-(energy[live] - energy[live].min()) / max(temperature, 1e-6))
            z = w.sum()
            if z <= 0:
                continue
            for qi, k in enumerate(mypairs):
                p_out[k, aa] = float(w[live & hit[:, qi]].sum() / z)
    return {"p": p_out, "admissible": adm, "n_clash": n_clash, "n_rotamer": n_rot_out}


def reference_free_energy(result: ScanResult, positions, partners, potential, *,
                          side: str = "to", background: np.ndarray | None = None,
                          beta: float = 1.0, scale: float = 1.0):
    """The per-position free energy over the **admissible** amino-acid background, and its spread.

    .. math::
        \\varphi_i(a) = \\text{scale}\\sum_j p_{ij}(a)\\, e(a, y_j), \\qquad
        F_i = -\\beta^{-1}\\ln\\!\\sum_{a \\in A_i} \\tilde p_a\\, e^{-\\beta\\varphi_i(a)}

    with :math:`A_i` the candidates no rotamer of which collides and :math:`\\tilde p` the
    background renormalised over it. A candidate that fits nowhere is **absent** from the sum, not
    present at its least-bad rotamer: a substitution the interface cannot accommodate is not a very
    unfavourable one, it is not a substitution.

    Args:
        result: what :func:`scan_structure` returned.
        positions, partners: its two ``(chain, residue)`` lists.
        potential: the residue-pair potential. TCRen2 is **directional** -- ``value(from, to)`` --
            so which slot the candidate takes follows ``side``, and getting it backwards silently
            transposes a matrix that is not symmetric. A pair the potential leaves undefined (it
            carries no receptor-side cysteine) contributes nothing, the rule
            :func:`tcren.energetics.scoring._interface_energy` applies with ``nansum``.
        side: which side was scanned, matching the :func:`scan_structure` call.
        background: 20-vector of amino-acid frequencies in :data:`AA` order, or ``None`` (default)
            for equimolar.
        beta: inverse temperature in the potential's units.
        scale: the interface's own coefficient, ``1 / sd_Native2026(Phi_I)``, so several interfaces
            can be summed commensurably. See :func:`tcren.energetics.scoring._phi_scale`.

    Returns:
        ``{"phi": (n_positions, 20), "F": (n_positions,), "dPhi": (n_positions,),
        "varPhi": (n_positions,), "n_admissible": (n_positions,), "native": (n_positions,)}``.
        ``phi`` is ``nan`` for an inadmissible candidate. ``dPhi`` and ``native`` are ``nan`` at a
        position whose own residue is inadmissible -- which should not happen on a deposited
        structure and is the invariant ``scripts/calibrate_background_scan.py`` picks the grid by.
    """
    P = len(positions)
    phi = np.zeros((P, 20))
    if side not in ("from", "to"):
        raise ValueError(f"side must be 'from' or 'to', got {side!r}")

    def _e(cand: str, part: str) -> float:
        a, b = (cand, part) if side == "from" else (part, cand)
        try:
            return float(potential.value(a, b))
        except KeyError:
            return 0.0

    cache: dict[str, np.ndarray] = {}
    for k, (i, j) in enumerate(zip(result.pair_i, result.pair_j)):
        y = partners[int(j)][1].aa
        e = cache.get(y)
        if e is None:
            e = cache[y] = np.nan_to_num(np.array([_e(a, y) for a in AA], float), nan=0.0)
        phi[int(i)] += result.p[k] * e
    phi *= scale
    phi[~result.admissible] = np.nan

    p0 = np.full(20, 1.0 / 20.0) if background is None else (
        np.asarray(background, float) / np.sum(background))
    F = np.full(P, np.nan)
    var = np.full(P, np.nan)
    native = np.full(P, np.nan)
    n_adm = result.admissible.sum(axis=1)
    for i in range(P):
        live = result.admissible[i]
        if not live.any():
            continue
        w = p0[live] / p0[live].sum()
        e = phi[i, live]
        shift = e.min()
        u = w * np.exp(-beta * (e - shift))
        z = u.sum()
        if z <= 0:
            continue
        F[i] = shift - math.log(z) / beta
        q = u / z
        m1 = float((q * e).sum())
        var[i] = float((q * e * e).sum() - m1 * m1)
        aa = positions[i][1].aa
        if aa in AA and live[AA.index(aa)]:
            native[i] = phi[i, AA.index(aa)]
    return {"phi": phi, "F": F, "dPhi": native - F, "varPhi": var,
            "n_admissible": n_adm, "native": native}
