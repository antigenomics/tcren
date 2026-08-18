"""Peptide conformational stability — does the pose hold, or was it merely drawn that way?

A contact potential scores whichever conformation it is handed. It cannot tell a peptide whose own
side chains hold it in the TCR-facing conformation from one that happens to have been modelled
there, because both present the same contact list. That blind spot is a specific, testable
hypothesis rather than a general complaint:

    Sewell (2026-08, ``manuscripts/2026-tcren/suggestions/sewell.txt``): Supplementary Fig. 4 of
    Dolton et al. (*J Clin Invest* 2024) shows a substantial **intra-peptide** interaction between
    P3 and P6 that stabilises the central bulge the TCR reads. "Poor binders could perhaps still
    make many contacts but fail to stabilise the productive peptide conformation" — which would
    explain why an additive contact model describes some systems well and others badly.

:func:`peptide_stability` measures it. Backbone φ/ψ are sampled by Metropolis Monte Carlo at
temperature against the DOPE potential, and the readout is not a better pose but how far the peptide
wanders: ``rmsf`` (the spread of the sampled ensemble) and ``drift`` (how far its mean moves from the
conformation it was given). A peptide held by its own interactions stays; one that is not, does not.

**The intra-peptide term is a switch, not a fixture.** ``intra_weight=0`` reruns the identical
sampling with the peptide's contacts with itself removed, so ``delta_rmsf`` between the two is the
stabilisation those interactions actually provide. That difference is the quantity Sewell's
hypothesis is about, and :func:`stability_table` computes it for a set of structures.

Moves are torsional and exact: perturb one φ or ψ and rotate every atom downstream of that bond. A
peptide has two free termini, so no loop closure is needed — the same rotation primitive
:mod:`tcren.rotamers` uses for χ, applied to N–Cα and Cα–C.

**What this is not.** It is not molecular dynamics: there is no solvent, no force field, no time, and
DOPE is a knowledge-based potential rather than an energy surface with physical units. It samples a
knowledge-based conformational basin, and ``rmsf`` is comparable *between* structures run with the
same settings, not against an MD RMSF in Å.
"""

from __future__ import annotations

from dataclasses import dataclass, asdict

import numpy as np
import polars as pl

from .structure.model import PEPTIDE_TYPE

DEFAULT_STEPS = 4000
DEFAULT_TEMPERATURE = 4.0     #: MC temperature in DOPE units
DEFAULT_SIGMA = 6.0           #: st. dev. of a single φ/ψ perturbation, degrees
DEFAULT_ANCHOR_W = 1.0        #: harmonic weight holding the anchor Cα in their pockets


@dataclass(slots=True)
class Stability:
    """How well one peptide holds the conformation it was given."""

    structure_id: str
    peptide: str
    rmsf: float            #: ensemble spread, Å — larger = floppier
    drift: float           #: distance of the ensemble mean from the input pose, Å
    energy: float          #: best DOPE energy seen
    energy_start: float    #: DOPE energy of the pose as given
    energy_gap: float      #: ``energy_start - energy``: how much better MC could do than the input
    accept_rate: float
    n_samples: int
    intra_weight: float

    def to_dict(self) -> dict:
        return asdict(self)


def _peptide_chain(structure):
    pep = next((c for c in structure.chains if c.chain_type == PEPTIDE_TYPE), None)
    if pep is None:
        raise ValueError(f"no peptide chain in {structure.pdb_id!r}; is the structure chain-typed?")
    return pep


def backbone_torsions(peptide_atoms):
    """φ/ψ torsions of a peptide as ``(axis_start, axis_end, moving_atom_indices)``.

    ``peptide_atoms`` is ``[(residue_index, atom_name), ...]`` in the flat order the coordinates are
    packed in. Indices returned are into that same array.

    A torsional rotation splits the chain at the bond and turns one side, so:

    * **φ** rotates about N–Cα, moving everything of that residue *except* N and Cα — the side chain
      included, since Cβ hangs off Cα and is not on the axis — plus every later residue;
    * **ψ** rotates about Cα–C, moving only the carbonyl O and every later residue. The side chain
      stays: it is on the Cα side of the bond.
    """
    res_of = np.array([r for r, _ in peptide_atoms])
    name_of = [n for _, n in peptide_atoms]
    out = []
    for res in sorted(set(res_of.tolist())):
        here = np.flatnonzero(res_of == res)
        later = np.flatnonzero(res_of > res)
        idx = {name_of[i]: int(i) for i in here}
        if not {"N", "CA", "C"} <= set(idx):
            continue                                   # incomplete backbone: no torsion to define
        phi_mov = np.concatenate([[i for i in here if name_of[i] not in ("N", "CA")], later])
        if len(phi_mov):
            out.append((idx["N"], idx["CA"], phi_mov.astype(np.int64)))
        psi_mov = np.concatenate([[i for i in here if name_of[i] == "O"], later])
        if len(psi_mov):
            out.append((idx["CA"], idx["C"], psi_mov.astype(np.int64)))
    return out


def _pack(structure, shell, anchors):
    """Flatten a structure for ``_relax.relax_interface``."""
    from scipy.spatial import cKDTree

    from .refine import _dope

    _table, atom_class, _x0, _dx, _nb = _dope()
    pep = _peptide_chain(structure)

    xyz, cls, pep_atoms, is_pep = [], [], [], []
    for chain in structure.chains:
        for res in chain.residues:
            for a in res.atoms:
                if a.element == "H":
                    continue
                xyz.append(a.coord)
                cls.append(atom_class.get(f"{res.resname}:{a.name}", -1))
                is_pep.append(chain is pep)
                if chain is pep:
                    pep_atoms.append((res.seq_index, a.name))
    xyz = np.asarray(xyz, dtype=np.float64)
    is_pep = np.asarray(is_pep)
    pep_idx = np.flatnonzero(is_pep)
    if pep_idx.size == 0:
        raise ValueError("peptide chain has no heavy atoms")
    # The kernel takes the peptide as one contiguous range; it is, because chains are packed whole.
    lo, hi = int(pep_idx[0]), int(pep_idx[-1]) + 1
    if hi - lo != pep_idx.size:
        raise ValueError("peptide atoms are not contiguous in the packed array")

    tors = backbone_torsions(pep_atoms)
    tor_a = np.array([lo + t[0] for t in tors], dtype=np.int32)
    tor_b = np.array([lo + t[1] for t in tors], dtype=np.int32)
    tor_mov, tor_mov_ptr = [], [0]
    for _s, _e, mov in tors:
        tor_mov.extend((lo + mov).tolist())
        tor_mov_ptr.append(len(tor_mov))

    # Partner atoms: only those near the peptide can contribute within DOPE's range, and the energy
    # is O(n_pep x n_par) per MC step, so this shell is what makes thousands of steps affordable.
    tree = cKDTree(xyz[~is_pep])
    other = np.flatnonzero(~is_pep)
    hits = tree.query_ball_point(xyz[lo:hi], shell)
    near = np.unique(np.concatenate([np.asarray(h, dtype=int) for h in hits if h]
                                    or [np.zeros(0, dtype=int)]))

    anchor_atom = [lo + i for i, (r, n) in enumerate(pep_atoms)
                   if n == "CA" and r in anchors]
    return {
        "xyz": xyz,
        "atom_class": np.asarray(cls, dtype=np.int32),
        "pep_lo": lo, "pep_hi": hi,
        "pep_res": np.asarray([r for r, _ in pep_atoms], dtype=np.int32),
        "tor_a": tor_a, "tor_b": tor_b,
        "tor_mov_ptr": np.asarray(tor_mov_ptr, dtype=np.int32),
        "tor_mov": np.asarray(tor_mov, dtype=np.int32),
        "par_atom": np.asarray(other[near], dtype=np.int32),
        "anchor_atom": np.asarray(anchor_atom, dtype=np.int32),
    }


def peptide_stability(structure, *, intra_weight: float = 1.0, n_steps: int = DEFAULT_STEPS,
                      temperature: float = DEFAULT_TEMPERATURE, sigma_deg: float = DEFAULT_SIGMA,
                      anchor_weight: float = DEFAULT_ANCHOR_W, min_sep: int = 3,
                      shell: float = 12.0, burn_in: int | None = None,
                      seed: int = 0) -> Stability:
    """Sample the peptide backbone and report how far it wanders from the pose it was given.

    Args:
        structure: chain-typed, MHC-annotated structure.
        intra_weight: weight on the peptide's DOPE contacts **with itself**. ``1.0`` includes them,
            ``0.0`` removes them — the comparison Sewell's hypothesis turns on.
        n_steps: Metropolis steps.
        temperature: MC temperature in DOPE units. Higher = more exploration; the value matters only
            in that every structure being compared must use the same one.
        sigma_deg: st. dev. of one φ/ψ perturbation, in degrees.
        anchor_weight: harmonic weight pinning the anchor Cα, so the peptide samples conformations
            rather than falling out of the groove. Anchors come from
            :func:`tcren.refine.predict_anchors`.
        min_sep: sequence separation below which intra-peptide pairs are ignored — neighbours are in
            contact by covalent geometry, not by folding.
        shell: partner atoms within this distance of the peptide enter the energy.
        burn_in: steps discarded before sampling; defaults to ``n_steps // 8``.
        seed: RNG seed. The result is deterministic given one.

    Returns:
        A :class:`Stability`.

    Raises:
        ValueError: if the structure has no chain-typed peptide, or its atoms are not contiguous.
    """
    from . import _relax
    from .refine import _dope
    from .refine.anchors import native_peptide, predict_anchors

    table, _amap, x_start, dx, nbins = _dope()
    pep = _peptide_chain(structure)
    sequence = native_peptide(structure)
    anchors = {pep.residues[i].seq_index for i in predict_anchors(sequence, structure).anchors
               if 0 <= i < len(pep.residues)}

    spec = _pack(structure, shell, anchors)
    out = _relax.relax_interface(
        dope_table=table, n_cls=table.shape[0], n_bins=nbins, x_start=x_start, dx=dx,
        n_steps=n_steps, temperature=temperature, sigma_deg=sigma_deg, intra_w=intra_weight,
        anchor_w=anchor_weight, min_sep=min_sep,
        burn_in=n_steps // 8 if burn_in is None else burn_in, seed=seed, **spec)

    return Stability(
        structure_id=getattr(structure, "pdb_id", "") or "", peptide=sequence,
        rmsf=float(out["rmsf"]), drift=float(out["drift"]), energy=float(out["energy"]),
        energy_start=float(out["energy_start"]),
        energy_gap=float(out["energy_start"]) - float(out["energy"]),
        accept_rate=float(out["accept_rate"]), n_samples=int(out["n_samples"]),
        intra_weight=float(intra_weight),
    )


def stability_table(structures, *, intra_weights=(1.0, 0.0), **kwargs) -> pl.DataFrame:
    """Run :func:`peptide_stability` at each ``intra_weight`` and report the difference.

    The difference is the point. ``delta_rmsf = rmsf(intra=0) - rmsf(intra=1)`` is how much the
    peptide's own contacts steady it: positive means removing them lets the backbone wander further,
    i.e. those interactions are holding the conformation together.

    Args:
        structures: an iterable of chain-typed, MHC-annotated structures.
        intra_weights: the weights to run. The default pair gives ``delta_rmsf``/``delta_drift``.
        **kwargs: passed to :func:`peptide_stability`.

    Returns:
        One row per structure: ``structure.id``, ``peptide``, ``rmsf``/``drift``/``energy`` at each
        weight (suffixed ``_intra1``/``_intra0`` for the default pair), and the deltas.
    """
    rows = []
    for s in structures:
        row = {"structure.id": getattr(s, "pdb_id", "") or "", "peptide": ""}
        by_w = {}
        for w in intra_weights:
            st = peptide_stability(s, intra_weight=w, **kwargs)
            by_w[w] = st
            row["peptide"] = st.peptide
            tag = "intra1" if w == 1.0 else ("intra0" if w == 0.0 else f"intra{w:g}")
            row[f"rmsf_{tag}"] = st.rmsf
            row[f"drift_{tag}"] = st.drift
            row[f"energy_{tag}"] = st.energy
            row[f"gap_{tag}"] = st.energy_gap
            row[f"accept_{tag}"] = st.accept_rate
        if 1.0 in by_w and 0.0 in by_w:
            row["delta_rmsf"] = by_w[0.0].rmsf - by_w[1.0].rmsf
            row["delta_drift"] = by_w[0.0].drift - by_w[1.0].drift
            row["delta_gap"] = by_w[0.0].energy_gap - by_w[1.0].energy_gap
        rows.append(row)
    return pl.DataFrame(rows)
