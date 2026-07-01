"""Peptide substitution + potential-guided refinement.

:func:`substitute_peptide` threads a new sequence onto the peptide backbone; :func:`refine_peptide`
runs a knowledge-based rigid-body Monte-Carlo refinement of the peptide pose via the compiled
``tcren._refine`` kernel. The refinement energy is the **DOPE** atom-level distance-dependent
statistical potential (Shen & Sali, *Protein Science* 2006), used here *independently* of the
TCRen/MJ potentials tcren scores epitopes with — so the pose is not optimised against the same
quantity it is later scored with. This is a lightweight, knowledge-based refine, NOT physics
relaxation (that is Rosetta FlexPepDock, as a subprocess).
"""

from __future__ import annotations

from functools import lru_cache
from importlib import resources

import numpy as np

from ..structure.model import PEPTIDE_TYPE, Atom, Chain, Residue, Structure
from .anchors import Decomposition, native_peptide, predict_anchors
from .rmsd import PeptideRMSD, peptide_rmsd
from .substitute import substitute_peptide

__all__ = [
    "substitute_peptide", "refine_peptide",
    "predict_anchors", "native_peptide", "Decomposition",
    "peptide_rmsd", "PeptideRMSD",
    "model_peptide", "available_engines",
]


def __getattr__(name):
    # Lazy re-export of the engine-orchestration layer (keeps `import tcren.refine` cheap and
    # avoids any import cycle through the engine registry).
    if name == "model_peptide":
        from .model import model_peptide
        return model_peptide
    if name == "available_engines":
        from .engines import available_engines
        return available_engines
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


@lru_cache(maxsize=1)
def _dope():
    """Load the bundled DOPE potential: ``(table (Nc,Nc,Nbins) float32, atom_class, x_start, dx)``."""
    path = resources.files("tcren.data").joinpath("dope_potential.npz")
    with resources.as_file(path) as p:
        d = np.load(p)
        table = np.ascontiguousarray(d["table"], dtype=np.float32)
        atom_class = {str(k): int(v) for k, v in zip(d["keys"].tolist(), d["vals"].tolist())}
        return table, atom_class, float(d["x_start"]), float(d["dx"]), int(d["nbins"])


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


def refine_peptide(structure: Structure, *, shell: float = 12.0, restraint_w: float = 0.5,
                   n_steps: int = 2000, trans_sigma: float = 0.2, rot_sigma: float = 0.05,
                   temp0: float = 1.0, temp1: float = 0.05, seed: int = 0):
    """Rigid-body refine the peptide pose against its TCR+MHC partners; ``(structure, energy)``.

    The energy is the DOPE atom-level distance-dependent statistical potential summed over all
    peptide$\\leftrightarrow$partner heavy-atom pairs within DOPE's range (its short-range bins are
    repulsive, so it provides its own clash term), plus a harmonic restraint to the input pose
    (``restraint_w``) that keeps the search local. Only partner atoms within ``shell`` Å of the
    peptide are considered. The structure must be chain-typed (peptide = chain of
    ``chain_type == 'PEPTIDE'``). Requires the compiled ``_refine`` ext + the bundled DOPE table.
    """
    from .. import _refine  # built by scikit-build-core; refinement requires it (no Python fallback)

    table, atom_class, x_start, dx, nbins = _dope()
    n_cls = table.shape[0]

    def cls(resname: str, atom: str) -> int:
        return atom_class.get(f"{resname}:{atom}", -1)

    pep = next((c for c in structure.chains if c.chain_type == PEPTIDE_TYPE), None)
    if pep is None:
        raise ValueError(f"no peptide chain in structure {structure.pdb_id!r}")

    # Peptide: ALL atoms (moved rigidly and written back); class -1 (e.g. H / non-standard) is
    # skipped in the energy but still transformed so the chain stays intact.
    pep_xyz, pep_class = [], []
    for res in pep.residues:
        for a in res.atoms:
            pep_xyz.append(a.coord)
            pep_class.append(cls(res.resname, a.name))
    pep_xyz = np.asarray(pep_xyz, dtype=np.float64)
    if not len(pep_xyz):
        raise ValueError("peptide chain has no atoms")

    # Partner: only DOPE-typed atoms, within `shell` of the peptide.
    par_xyz, par_class = [], []
    for chain in structure.chains:
        if chain is pep:
            continue
        for res in chain.residues:
            for a in res.atoms:
                c = cls(res.resname, a.name)
                if c >= 0:
                    par_xyz.append(a.coord)
                    par_class.append(c)
    par_xyz = np.asarray(par_xyz, dtype=np.float64).reshape(-1, 3)
    par_class = np.asarray(par_class, dtype=np.int32)
    if len(par_xyz):
        keep = np.sqrt(((par_xyz[:, None, :] - pep_xyz[None, :, :]) ** 2).sum(-1).min(1)) <= shell
        par_xyz, par_class = par_xyz[keep], par_class[keep]

    best, energy, _n_accept = _refine.refine(
        pep_xyz, np.asarray(pep_class, np.int32), par_xyz, par_class,
        table.reshape(-1), n_cls, nbins, x_start, dx, restraint_w,
        n_steps, trans_sigma, rot_sigma, temp0, temp1, seed,
    )
    return _rebuild_peptide(structure, pep, np.asarray(best)), float(energy)
