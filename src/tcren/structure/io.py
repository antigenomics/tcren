"""Parse PDB / mmCIF files into the :mod:`tcren.structure.model` data model.

Accepts plain ``.pdb``/``.ent``/``.cif``/``.mmcif`` files, their gzip-compressed forms
(``.pdb.gz``/``.cif.gz`` …), and — for batches — directories or ``.tar``/``.tar.gz`` archives
of any of those (see :func:`iter_structures`). Structure identifiers are resolved from the file
name by :func:`structure_id_from_path`.
"""

from __future__ import annotations

import gzip
import tarfile
import tempfile
from collections.abc import Callable, Iterator
from pathlib import Path

import numpy as np
from Bio.Data.PDBData import protein_letters_3to1_extended
from Bio.PDB import MMCIFParser, PDBParser

from .model import Atom, Chain, Residue, Structure

# Extended three→one map covers modified residues (MSE→M, SEC→U, …).
_THREE_TO_ONE = dict(protein_letters_3to1_extended)
_WATER = {"HOH", "WAT", "DOD"}

# Recognised structure extensions (mmCIF first). A file also matches with a trailing ``.gz``.
_CIF_SUFFIXES = (".cif", ".mmcif")
_PDB_SUFFIXES = (".pdb", ".ent")
STRUCTURE_SUFFIXES = _CIF_SUFFIXES + _PDB_SUFFIXES
_TAR_SUFFIXES = (".tar", ".tar.gz", ".tgz")


def _strip_gz(name: str) -> tuple[str, bool]:
    """``(name_without_trailing_.gz, was_gzipped)``."""
    return (name[:-3], True) if name.lower().endswith(".gz") else (name, False)


def is_structure_file(name: str | Path) -> bool:
    """True if ``name`` is a (optionally gzipped) PDB/mmCIF structure file."""
    inner, _ = _strip_gz(Path(name).name)
    return inner.lower().endswith(STRUCTURE_SUFFIXES)


def structure_id_from_path(path: str | Path) -> str:
    """Resolve a structure identifier from a file name.

    Strips a trailing ``.gz`` and the structure extension, then takes the part before the
    first ``_`` (so ``4x5w_renumbered.cif`` and ``1ao7.pdb.gz`` and
    ``6uk4_TCRpMHCmodels.pdb`` all resolve to their PDB id).
    """
    inner, _ = _strip_gz(Path(path).name)
    stem = inner.rsplit(".", 1)[0] if "." in inner else inner
    return stem.split("_")[0]


def _structure_format(name: str) -> str:
    """``"cif"`` or ``"pdb"`` for a (possibly gzipped) structure file name."""
    inner, _ = _strip_gz(Path(name).name)
    return "cif" if inner.lower().endswith(_CIF_SUFFIXES) else "pdb"


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
    inner, gzipped = _strip_gz(path.name)
    pdb_id = pdb_id or inner.rsplit(".", 1)[0]
    parser = MMCIFParser(QUIET=True) if _structure_format(path.name) == "cif" else PDBParser(QUIET=True)

    if gzipped:
        with gzip.open(path, "rt") as handle:
            bio = parser.get_structure(pdb_id, handle)
    else:
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


def structure_paths(src: str | Path) -> list[Path]:
    """List structure files for ``src`` (a single file or a directory), sorted.

    Recognises plain and gzipped PDB/mmCIF (``.pdb``, ``.cif.gz``, …). For archives or
    streaming, use :func:`iter_structures`.
    """
    src = Path(src)
    if src.is_dir():
        return sorted(p for p in src.iterdir() if is_structure_file(p))
    return [src]


def iter_structures(
    src: str | Path,
    importer: Callable[..., Structure] = import_structure,
    on_error: str = "raise",
    **kwargs,
) -> Iterator[tuple[str, Structure]]:
    """Yield ``(pdb_id, Structure)`` for a file, directory, or ``.tar``/``.tar.gz`` archive.

    Handles plain and gzipped PDB/mmCIF (``.pdb``/``.cif``/``.pdb.gz``/``.cif.gz`` …); a
    directory is scanned for those; a tar archive is streamed member-by-member (each member
    materialised to a temp file so the path-based ``importer`` works unchanged). The
    identifier is resolved per file by :func:`structure_id_from_path`.

    Args:
        src: structure file, directory, or tar archive.
        importer: per-file parser — :func:`import_structure` (default, trims the C-gene) or
            :func:`parse_structure` (parity-pure). Extra ``kwargs`` are forwarded to it.
        on_error: ``"raise"`` (default) or ``"skip"`` to ignore files that fail to parse.
    """
    src = Path(src)
    name = src.name.lower()

    def _safe(path: Path, pdb_id: str) -> Structure | None:
        try:
            return importer(path, pdb_id=pdb_id, **kwargs)
        except Exception:
            if on_error == "raise":
                raise
            return None

    if src.is_file() and name.endswith(_TAR_SUFFIXES):
        with tarfile.open(src) as tar:
            for member in tar.getmembers():
                if not (member.isfile() and is_structure_file(member.name)):
                    continue
                inner, _ = _strip_gz(Path(member.name).name)
                ext = "." + inner.rsplit(".", 1)[-1] + (".gz" if member.name.lower().endswith(".gz") else "")
                fh = tar.extractfile(member)
                if fh is None:
                    continue
                with tempfile.NamedTemporaryFile(suffix=ext) as tmp:
                    tmp.write(fh.read())
                    tmp.flush()
                    s = _safe(Path(tmp.name), structure_id_from_path(member.name))
                if s is not None:
                    yield structure_id_from_path(member.name), s
        return

    for path in structure_paths(src):
        s = _safe(path, structure_id_from_path(path))
        if s is not None:
            yield structure_id_from_path(path), s


def _atom_name_field(name: str) -> str:
    """PDB columns 13-16 for an atom name (the standard left/right justification rule)."""
    return f"{name:<4}" if len(name) >= 4 else f" {name:<3}"


def pdb_lines(structure: Structure, transform=None, keep_hydrogens: bool = True) -> list[str]:
    """ATOM/TER/END record lines for ``structure`` (optionally coordinate-transformed).

    One conformer per atom name per residue (drops duplicate altlocs). ``transform`` is an
    optional ``coord -> coord`` callable (e.g. for an oriented frame); identity if ``None``.
    Author residue numbers + insertion codes are preserved.
    """
    lines: list[str] = []
    serial = 1
    for chain in structure.chains:
        chain_id = (chain.chain_id or " ")[0]
        last = None
        for res in chain.residues:
            seen: set[str] = set()
            icode = (res.insertion_code or " ")[:1] or " "
            for atom in res.atoms:
                element = (atom.element or atom.name[:1]).strip().upper()
                if (not keep_hydrogens and element == "H") or atom.name in seen:
                    continue
                seen.add(atom.name)
                x, y, z = transform(atom.coord) if transform else atom.coord
                lines.append(
                    f"ATOM  {serial % 100000:>5} {_atom_name_field(atom.name)} "
                    f"{res.resname:>3} {chain_id}{res.pdb_index:>4}{icode}   "
                    f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {element:>2}"
                )
                serial += 1
                last = res
        if last is not None:
            lines.append(
                f"TER   {serial % 100000:>5}      {last.resname:>3} "
                f"{chain_id}{last.pdb_index:>4}{(last.insertion_code or ' ')[:1] or ' '}"
            )
            serial += 1
    lines.append("END")
    return lines


def write_pdb(structure: Structure, path: str | Path, transform=None,
              keep_hydrogens: bool = True) -> Path:
    """Write ``structure`` to a PDB file; return the path.

    A ``.gz`` suffix (``foo.pdb.gz``) transparently gzip-compresses the output.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    text = "\n".join(pdb_lines(structure, transform, keep_hydrogens)) + "\n"
    if path.name.lower().endswith(".gz"):
        with gzip.open(path, "wt") as fh:
            fh.write(text)
    else:
        path.write_text(text)
    return path


def write_structure(structure: Structure, path: str | Path, **kwargs) -> Path:
    """Format-dispatch writer (PDB / gzipped PDB only for now; mmCIF deferred)."""
    path = Path(path)
    inner, _ = _strip_gz(path.name)
    if not inner.lower().endswith(_PDB_SUFFIXES):
        raise ValueError(f"only PDB output is supported, got {path.name!r}")
    return write_pdb(structure, path, **kwargs)
