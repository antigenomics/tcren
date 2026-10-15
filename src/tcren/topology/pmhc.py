"""Receptor-independent pMHC topography and peptide solvent exposure.

The supplied coordinates are unchanged; removing a receptor does not undo induced
fit. SASA is Shrake–Rupley (100 sphere points, 1.4 Å probe) in the groove frame.
"""
from __future__ import annotations

import numpy as np

from .surface import KYTE_DOOLITTLE, _groove_frame, surface_map, surface_stats


def pmhc_features(structure) -> dict:
    """Return pMHC-only surface scalars and peptide/side-chain SASA in Å².

    Receptor coordinates enter neither the surface normal nor the occluding atoms.
    ``sasa_hydropathy`` weights Kyte–Doolittle values by exposed side-chain area;
    glycine has zero side-chain area. No class-II core/register is guessed.
    All MHC chains retain their physical identities. Nonstandard peptides fail.
    """
    from Bio.PDB.Atom import Atom
    from Bio.PDB.Chain import Chain
    from Bio.PDB.Model import Model
    from Bio.PDB.Residue import Residue
    from Bio.PDB.SASA import ShrakeRupley
    from scipy.spatial import cKDTree

    from ..structure.model import Structure

    peptide = structure.by_type("PEPTIDE")
    mhc = structure.by_type("MHCa", "MHCb")
    if len(peptide) != 1 or not mhc or set(peptide[0].sequence()) - set(KYTE_DOOLITTLE):
        raise ValueError("pMHC features require one standard-amino-acid peptide and annotated MHC")
    classes = {c.chain_supertype for c in mhc if c.chain_supertype in {"MHCI", "MHCII"}}
    if len(classes) != 1:
        raise ValueError("pMHC requires one unambiguous annotated MHC class")
    if classes == {"MHCII"} and {c.chain_type for c in mhc} != {"MHCa", "MHCb"}:
        raise ValueError("class II requires both MHC groove chains")
    selected = Structure(structure.pdb_id, structure.by_type("PEPTIDE", "MHCa", "MHCb", "B2M"))
    origin, basis = _groove_frame(selected)
    values = surface_stats(surface_map(selected, scale="kd"))
    values.update(peptide=peptide[0].sequence(), peptide_length=len(peptide[0].residues),
                  mhc_class="II" if classes == {"MHCII"} else "I")
    sr = ShrakeRupley(n_points=100, probe_radius=1.4)
    # Atoms farther than the sum of expanded radii cannot occlude peptide SASA.
    # Filter before invoking the existing native-neighbour SASA implementation.
    peptide_xyz = np.array([a.coord for r in peptide[0].residues for a in r.atoms
                            if a.element not in {"H", "D"}])
    if not len(peptide_xyz):
        raise ValueError("peptide has no atoms")
    tree = cKDTree(peptide_xyz)
    atoms = [a for c in selected.chains for r in c.residues for a in r.atoms
             if a.element not in {"H", "D"}]
    radius = max(sr.radii_dict[a.element] for a in atoms) + sr.probe_radius
    model = Model(0)
    neighbours = tree.query(np.array([a.coord for a in atoms]))[0] <= 2 * radius
    nearby = {id(a) for a, keep in zip(atoms, neighbours) if keep}
    for chain in selected.chains:
        biochain = Chain(chain.chain_id)
        model.add(biochain)
        for residue in chain.residues:
            biores = Residue((" ", residue.seq_index, " "), residue.resname, " ")
            biochain.add(biores)
            for atom in residue.atoms:
                if atom.element in {"H", "D"}:
                    continue
                if chain.chain_type != "PEPTIDE" and id(atom) not in nearby:
                    continue
                if atom.name in biores:  # Match the writer's single-conformer atom convention.
                    continue
                biores.add(Atom(atom.name, (atom.coord - origin) @ basis.T, 0., 1., " ",
                                atom.name, len(biores), element=atom.element))
    sr.compute(model, level="A")
    backbone = {"N", "CA", "C", "O", "OXT", "OT1", "OT2", "OC1", "OC2"}
    total, side, weighted = 0., 0., 0.
    for residue in peptide[0].residues:
        atoms = model[peptide[0].chain_id][(" ", residue.seq_index, " ")]
        area = sum(a.sasa for a in atoms if a.name not in backbone)
        total += sum(a.sasa for a in atoms)
        side += area
        weighted += area * KYTE_DOOLITTLE[residue.aa]
    values.update(peptide_sasa_A2=float(total), sidechain_sasa_A2=float(side),
                  sasa_hydropathy=float(weighted / side) if side else float("nan"))
    return values
