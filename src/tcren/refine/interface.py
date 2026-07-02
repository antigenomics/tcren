"""Interface interaction energy — the ΔΔG core (native, DOPE, no PyRosetta).

``interface_energy(structure)`` sums the DOPE atom-level statistical potential over peptide↔partner
heavy-atom pairs — the interaction energy across the TCR/MHC:peptide interface. Because only the cross
(peptide↔partner) terms are summed, it is already ``E_bound − E_separated`` for the interaction, i.e.
the interface ΔG contribution of the peptide. This is the quantity a sampling ΔΔG differences between a
mutant and the native (each first repacked/relaxed); see the ``_relax`` kernel roadmap in CPP_REWRITE.md.

Reuses the bundled DOPE table + atom-class map from :mod:`tcren.refine` (``_dope``), and the compiled
``tcren._relax`` kernel.
"""

from __future__ import annotations

import numpy as np

from ..structure.model import PEPTIDE_TYPE, Structure


def interface_energy(structure: Structure, *, shell: float = 12.0) -> float:
    """DOPE interaction energy across the peptide↔partner interface (lower = more favourable).

    The structure must be chain-typed (peptide = ``chain_type == 'PEPTIDE'``). Only partner heavy atoms
    within ``shell`` Å of the peptide are considered (beyond DOPE's range they contribute nothing).
    Requires the compiled ``tcren._relax`` ext + the bundled DOPE table.
    """
    from . import _dope  # bundled DOPE table + atom-class map (tcren.refine)
    from .. import _relax

    table, atom_class, x_start, dx, nbins = _dope()
    n_cls = table.shape[0]

    def cls(resname: str, atom: str) -> int:
        return atom_class.get(f"{resname}:{atom}", -1)

    pep = next((c for c in structure.chains if c.chain_type == PEPTIDE_TYPE), None)
    if pep is None:
        raise ValueError(f"no peptide chain in {structure.pdb_id!r}")

    pep_xyz, pep_class = [], []
    for res in pep.residues:
        for a in res.atoms:
            pep_xyz.append(a.coord)
            pep_class.append(cls(res.resname, a.name))
    pep_xyz = np.asarray(pep_xyz, dtype=np.float64).reshape(-1, 3)
    pep_class = np.asarray(pep_class, dtype=np.int32)
    if len(pep_xyz) == 0:
        raise ValueError("peptide chain has no atoms")

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

    return float(_relax.interface_energy(
        pep_xyz, pep_class, par_xyz, par_class, table.reshape(-1), n_cls, nbins, x_start, dx))
