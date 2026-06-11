"""Parse PDB / mmCIF files into the :mod:`tcren.structure.model` data model."""

from __future__ import annotations

from pathlib import Path

import numpy as np
from Bio.Data.PDBData import protein_letters_3to1_extended
from Bio.PDB import MMCIFParser, PDBParser

from .model import Atom, Chain, Residue, Structure

# Extended three→one map covers modified residues (MSE→M, SEC→U, …).
_THREE_TO_ONE = dict(protein_letters_3to1_extended)
_WATER = {"HOH", "WAT", "DOD"}


def _one_letter(resname: str) -> str | None:
    """Map a three-letter residue name to one letter, or ``None`` if not an amino acid."""
    return _THREE_TO_ONE.get(resname.strip().upper())


def _select_atoms(residue, keep_hydrogens: bool) -> tuple[Atom, ...]:
    """Collect atoms, keeping *every* alternate conformer.

    The legacy mir contact definition takes the minimum inter-atomic distance over all
    alternate locations, so each altloc position is retained as a separate atom.
    """
    atoms: list[Atom] = []
    for atom in residue.get_atoms():
        children = atom.disordered_get_list() if atom.is_disordered() else [atom]
        for child in children:
            element = (child.element or child.get_name()[0]).strip().upper()
            if not keep_hydrogens and element == "H":
                continue
            atoms.append(
                Atom(
                    name=child.get_name().strip(),
                    element=element,
                    coord=np.asarray(child.get_coord(), dtype=np.float64),
                )
            )
    return tuple(atoms)


def parse_structure(
    path: str | Path,
    pdb_id: str | None = None,
    model: int = 0,
    keep_hydrogens: bool = True,
) -> Structure:
    """Parse a structure file into a :class:`Structure`.

    Residues are taken in author order; only amino-acid residues (standard or modified,
    via the extended three→one table) are kept — waters, ions and ligands are dropped.
    Each kept residue receives a 0-based sequential ``seq_index`` per chain, matching the
    legacy ``mir`` ``residue.index``.

    Args:
        path: Path to a ``.pdb``/``.ent`` or ``.cif``/``.mmcif`` file.
        pdb_id: Structure identifier; defaults to the file stem.
        model: Model index to read (default 0 — the first model).
        keep_hydrogens: Keep hydrogen atoms (default ``True`` — the legacy mir contact
            definition counts hydrogens when a structure provides them).

    Returns:
        The parsed :class:`Structure`.
    """
    path = Path(path)
    pdb_id = pdb_id or path.stem
    suffix = path.suffix.lower()
    if suffix in (".cif", ".mmcif"):
        parser = MMCIFParser(QUIET=True)
    else:
        parser = PDBParser(QUIET=True)

    bio = parser.get_structure(pdb_id, str(path))
    bio_model = list(bio)[model]

    chains: list[Chain] = []
    for bio_chain in bio_model:
        residues: list[Residue] = []
        seq_index = 0
        for bio_res in bio_chain:
            hetflag, resseq, icode = bio_res.id
            resname = bio_res.get_resname().strip().upper()
            if resname in _WATER:
                continue
            # The legacy mir indexes only ATOM records (blank het flag); it skips every
            # HETATM — ligands, ions, and even modified residues such as CIR
            # (citrulline) or MSE that sit inside a polymer chain. Unknown ATOM
            # residues (e.g. the AMN chain cap) are kept and labelled 'X'.
            if hetflag.strip():
                continue
            aa = _one_letter(resname)
            if aa is None:
                aa = "X"
            atoms = _select_atoms(bio_res, keep_hydrogens)
            if not atoms:
                continue
            residues.append(
                Residue(
                    seq_index=seq_index,
                    pdb_index=int(resseq),
                    insertion_code=icode.strip(),
                    aa=aa if len(aa) == 1 else "X",
                    resname=resname,
                    atoms=atoms,
                )
            )
            seq_index += 1
        if residues:
            chains.append(Chain(chain_id=bio_chain.id, residues=residues))

    return Structure(pdb_id=pdb_id, chains=chains)
