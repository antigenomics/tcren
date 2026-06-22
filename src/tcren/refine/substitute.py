"""Backbone-preserving peptide substitution.

``score_peptides`` scores a candidate peptide *virtually* (it re-indexes the potential matrix over
the native contact map — no atoms move). When you want to actually re-dock / refine a candidate you
first need its coordinates: :func:`substitute_peptide` threads an equal-length sequence onto the
existing peptide backbone, keeping N/Cα/C/O(+Cβ) and dropping the old side-chain atoms beyond Cβ
(a refiner / rotamer repack rebuilds them). Pure data-model manipulation; returns a new structure.
"""

from __future__ import annotations

from ..structure.model import PEPTIDE_TYPE, Chain, Residue, Structure

# Atoms retained on substitution: backbone + Cβ (the rest of the side chain is identity-specific).
_KEEP = {"N", "CA", "C", "O", "CB"}

_ONE_TO_THREE = {
    "A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP", "C": "CYS", "Q": "GLN", "E": "GLU",
    "G": "GLY", "H": "HIS", "I": "ILE", "L": "LEU", "K": "LYS", "M": "MET", "F": "PHE",
    "P": "PRO", "S": "SER", "T": "THR", "W": "TRP", "Y": "TYR", "V": "VAL",
}


def substitute_peptide(structure: Structure, new_peptide: str,
                       chain_type: str = PEPTIDE_TYPE) -> Structure:
    """Return a copy of ``structure`` with the peptide chain threaded to ``new_peptide``.

    The peptide backbone (and Cβ) is preserved; side-chain atoms beyond Cβ are dropped (and Cβ
    too for any position mutated to glycine). ``new_peptide`` must equal the peptide length and use
    the 20 standard one-letter amino acids.

    Raises:
        ValueError: if there is no peptide chain, the length differs, or a code is non-standard.
    """
    new_peptide = new_peptide.upper()
    pep = next((c for c in structure.chains if c.chain_type == chain_type), None)
    if pep is None:
        raise ValueError(f"no {chain_type} chain in structure {structure.pdb_id!r}")
    if len(new_peptide) != len(pep.residues):
        raise ValueError(
            f"length mismatch: peptide has {len(pep.residues)} residues, got {len(new_peptide)}"
        )

    new_residues: list[Residue] = []
    for res, aa in zip(pep.residues, new_peptide):
        resname = _ONE_TO_THREE.get(aa)
        if resname is None:
            raise ValueError(f"non-standard amino acid {aa!r} in peptide")
        keep = _KEEP - ({"CB"} if aa == "G" else set())
        atoms = tuple(a for a in res.atoms if a.name in keep)
        new_residues.append(Residue(res.seq_index, res.pdb_index, res.insertion_code,
                                    aa, resname, atoms))

    new_pep = Chain(chain_id=pep.chain_id, residues=new_residues,
                    chain_type=pep.chain_type, chain_supertype=pep.chain_supertype)
    chains = [new_pep if c is pep else c for c in structure.chains]
    return Structure(pdb_id=structure.pdb_id, chains=chains,
                     complex_species=structure.complex_species, cell_type=structure.cell_type)
