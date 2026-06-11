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


def _trim_constant_regions(structure: Structure, min_score: float) -> None:
    """Drop each chain's C-terminal TCR constant domain in place (V-domain preserved).

    The constant region is C-terminal, so trimming removes trailing residues and leaves
    the variable-domain ``seq_index`` values unchanged (contacts/scoring unaffected). A
    no-op for chains without a constant domain (e.g. variable-only or non-TCR chains).
    """
    from ..annotation.cgene import constant_span

    for chain in structure.chains:
        span = constant_span(chain.sequence(), min_score=min_score)
        if span is None:
            continue
        start, _end = span
        if 0 < start < len(chain.residues):
            chain.residues = chain.residues[:start]


def import_structure(
    path: str | Path,
    pdb_id: str | None = None,
    model: int = 0,
    keep_hydrogens: bool = True,
    trim_c_gene: bool = True,
    keep_c_gene: bool = False,
    min_constant_score: float = 80.0,
) -> Structure:
    """Parse a structure and prepare it for interface analysis.

    Wraps :func:`parse_structure`, records the αβ/γδ cell type from the TCR constant
    region, and — by default — trims that constant region so downstream analysis works on
    the variable domains and the interface.

    Args:
        path, pdb_id, model, keep_hydrogens: as in :func:`parse_structure`.
        trim_c_gene: Trim the TCR constant domain (default ``True``).
        keep_c_gene: Retain the constant domain even if ``trim_c_gene`` is set. **Use this
            for molecular-dynamics / FlexPepDock and any workflow that needs the full
            chain** — those depend on the presence of the C-gene.
        min_constant_score: Minimum constant-region alignment score to trim on.

    Returns:
        The imported :class:`Structure` with ``cell_type`` set.
    """
    # TODO: molecular dynamics, FlexPepDock, and full-chain workflows depend on the
    # presence of the C-gene — pass keep_c_gene=True there.
    from ..annotation.cgene import cell_type as _cell_type

    structure = parse_structure(path, pdb_id=pdb_id, model=model, keep_hydrogens=keep_hydrogens)
    structure.cell_type = _cell_type(structure, min_score=min_constant_score)
    if trim_c_gene and not keep_c_gene:
        _trim_constant_regions(structure, min_score=min_constant_score)
    return structure
