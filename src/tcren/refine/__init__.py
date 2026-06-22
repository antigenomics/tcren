"""Peptide substitution + potential-guided refinement.

:func:`substitute_peptide` threads a new sequence onto the peptide backbone; :func:`refine_peptide`
runs a knowledge-based rigid-body Monte-Carlo refinement of the peptide pose (statistical potential
+ soft clash) via the compiled ``tcren._refine`` kernel. This is a lightweight, in-character refine
for a knowledge-based method — NOT physics relaxation (that is Rosetta, as a subprocess).
"""

from __future__ import annotations

import numpy as np

from ..structure.model import PEPTIDE_TYPE, Atom, Chain, Residue, Structure
from .substitute import substitute_peptide

__all__ = ["substitute_peptide", "refine_peptide"]


def _padded_potential(potential):
    """Dense potential matrix + aa→index map, padded with a neutral (zero) row/col for unknowns."""
    matrix, index = potential.as_matrix()
    n = matrix.shape[0]
    padded = np.zeros((n + 1, n + 1), dtype=np.float64)
    padded[:n, :n] = np.nan_to_num(np.asarray(matrix, dtype=np.float64), nan=0.0)
    return padded, index, n  # n = index of the neutral "unknown" row/col


def _rebuild_peptide(structure: Structure, pep: Chain, coords: np.ndarray) -> Structure:
    """Return a copy of ``structure`` with the peptide atoms moved to ``coords`` (same order)."""
    k = 0
    new_res = []
    for res in pep.residues:
        atoms = tuple(Atom(a.name, a.element, np.asarray(coords[k + j], dtype=np.float64))
                      for j, a in enumerate(res.atoms))
        k += len(res.atoms)
        new_res.append(Residue(res.seq_index, res.pdb_index, res.insertion_code,
                               res.aa, res.resname, atoms))
    new_pep = Chain(pep.chain_id, new_res, chain_type=pep.chain_type,
                    chain_supertype=pep.chain_supertype)
    chains = [new_pep if c is pep else c for c in structure.chains]
    return Structure(structure.pdb_id, chains, complex_species=structure.complex_species,
                     cell_type=structure.cell_type)


def refine_peptide(structure: Structure, potential=None, *, shell: float = 12.0,
                   cutoff: float = 5.0, clash_d0: float = 3.0, clash_w: float = 1.0,
                   restraint_w: float = 1.0, n_steps: int = 2000, trans_sigma: float = 0.2,
                   rot_sigma: float = 0.05, temp0: float = 1.0, temp1: float = 0.05, seed: int = 0):
    """Rigid-body refine the peptide pose against its TCR+MHC partners; ``(structure, energy)``.

    Scores with ``potential`` (default the MJ general contact potential — appropriate for the
    peptide's groove/TCR packing) over inter-chain residue contacts within ``cutoff`` Å, plus a soft
    heavy-atom clash penalty and a **harmonic restraint to the input pose** (``restraint_w``) that
    keeps the search local (without it a rigid contact-energy minimisation just ejects the peptide).
    Only partner residues within ``shell`` Å of the peptide are considered. The structure must be
    chain-typed (peptide = chain of ``chain_type == 'PEPTIDE'``). Returns a copy with the refined
    peptide coordinates and the final energy. Requires the compiled ``_refine`` ext.
    """
    from .. import _refine  # built by scikit-build-core; refinement requires it (no Python fallback)
    from ..potential import mj

    pot = potential if potential is not None else mj()
    padded, index, unk = _padded_potential(pot)
    aa_idx = lambda aa: index.get(aa, unk)  # noqa: E731

    pep = next((c for c in structure.chains if c.chain_type == PEPTIDE_TYPE), None)
    if pep is None:
        raise ValueError(f"no peptide chain in structure {structure.pdb_id!r}")

    pep_xyz, pep_atom_res, pep_res_aa = [], [], []
    for ri, res in enumerate(pep.residues):
        pep_res_aa.append(aa_idx(res.aa))
        for a in res.atoms:
            pep_xyz.append(a.coord)
            pep_atom_res.append(ri)
    pep_xyz = np.asarray(pep_xyz, dtype=np.float64)
    if not len(pep_xyz):
        raise ValueError("peptide chain has no atoms")

    par_xyz, par_atom_res, par_res_aa = [], [], []
    rj = 0
    for chain in structure.chains:
        if chain is pep:  # partners = every non-peptide chain (TCR + MHC + β2m); shell filters distant ones
            continue
        for res in chain.residues:
            coords = np.array([a.coord for a in res.atoms], dtype=np.float64)
            if not len(coords):
                continue
            if np.sqrt(((coords[:, None, :] - pep_xyz[None, :, :]) ** 2).sum(-1).min()) > shell:
                continue
            par_res_aa.append(aa_idx(res.aa))
            for a in res.atoms:
                par_xyz.append(a.coord)
                par_atom_res.append(rj)
            rj += 1
    par_xyz = np.asarray(par_xyz, dtype=np.float64).reshape(-1, 3)

    best, energy, _n_accept = _refine.refine(
        pep_xyz, np.asarray(pep_atom_res, np.int32), np.asarray(pep_res_aa, np.int32),
        par_xyz, np.asarray(par_atom_res, np.int32), np.asarray(par_res_aa, np.int32),
        padded, cutoff, clash_d0, clash_w, restraint_w, n_steps, trans_sigma, rot_sigma,
        temp0, temp1, seed,
    )
    return _rebuild_peptide(structure, pep, np.asarray(best)), float(energy)
