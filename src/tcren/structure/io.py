"""Parse PDB / mmCIF files into the :mod:`tcren.structure.model` data model.

Accepts plain ``.pdb``/``.ent``/``.cif``/``.mmcif`` files, their gzip-compressed forms
(``.pdb.gz``/``.cif.gz`` …), and — for batches — directories or ``.tar``/``.tar.gz`` archives
of any of those (see :func:`iter_structures`). Structure identifiers are resolved from the file
name by :func:`structure_id_from_path`.
"""

from __future__ import annotations

import glob
import gzip
import tarfile
import tempfile
from collections.abc import Callable, Iterable, Iterator
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
#: a plain-text list of structure paths, one per line (see structure_paths)
_MANIFEST_SUFFIXES = (".txt", ".list", ".lst")


def _strip_gz(name: str) -> tuple[str, bool]:
    """``(name_without_trailing_.gz, was_gzipped)``."""
    return (name[:-3], True) if name.lower().endswith(".gz") else (name, False)


def is_structure_file(name: str | Path) -> bool:
    """True if ``name`` is a (optionally gzipped) PDB/mmCIF structure file.

    Rejects macOS AppleDouble sidecars (``._4x5w.pdb``), which carry the extension of the
    file they shadow but hold a binary resource fork -- tarring a structure set on HFS+
    writes one beside every member, and feeding one to a parser raises a decode error
    several frames from the cause.
    """
    inner, _ = _strip_gz(Path(name).name)
    if inner.startswith("._"):
        return False
    return inner.lower().endswith(STRUCTURE_SUFFIXES)


def structure_stem(path: str | Path) -> str:
    """The file name with a trailing ``.gz`` and the structure extension removed."""
    inner, _ = _strip_gz(Path(path).name)
    return inner.rsplit(".", 1)[0] if "." in inner else inner


def structure_id_from_path(path: str | Path) -> str:
    """Resolve a structure identifier from a file name.

    Strips a trailing ``.gz`` and the structure extension, then takes the part before the
    first ``_`` (so ``4x5w_renumbered.cif`` and ``1ao7.pdb.gz`` and
    ``6uk4_TCRpMHCmodels.pdb`` all resolve to their PDB id).

    .. warning::
       Lossy by design, and *not* unique for cohorts whose file names encode metadata after an
       underscore -- ``VDJdb_Model_603_min.pdb`` and ``VDJdb_Model_604_min.pdb`` both give
       ``VDJdb``. Prefer :func:`iter_structures`, which detects that case for a whole set and
       falls back to the full stem rather than silently collapsing rows.
    """
    return structure_stem(path).split("_")[0]


def resolve_structure_ids(paths) -> dict[str, str]:
    """Map each path to an id: the PDB-id prefix when that is unique over the set, else the stem.

    The prefix rule is what makes ``4x5w_renumbered.cif`` come back as ``4x5w``, and it is right
    for RCSB-derived files. It is wrong -- silently, and in a way that destroys rows downstream --
    for any set whose names carry metadata after the first underscore. Deciding per SET rather than
    per file keeps the convenience where it is unambiguous and refuses it where it is not.
    """
    paths = list(paths)
    short = [structure_id_from_path(p) for p in paths]
    if len(set(short)) == len(paths):
        return {str(p): s for p, s in zip(paths, short)}
    return {str(p): structure_stem(p) for p in paths}


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


# =================================================================================================
# fast PDB path
# =================================================================================================
# Biopython's PDBParser is 86% of the wall clock of any dataset-scale pass through tcren: profiled
# over Native2026 it spends 19.0 s of a 30.3 s run in ``_parse_coordinates`` alone, building 2.1 M
# Atom and 2.4 M Entity objects we immediately throw away, and a further 5.3 s is spent in
# ``_select_atoms`` unwrapping them again. A PDB ATOM record is fixed-column, so the whole file can
# be sliced as one uint8 array and the residue boundaries found with a single ``np.flatnonzero``.
#
# This path is exact or it does not run. It bails to Biopython -- returning ``None`` -- on anything
# it is not certain of: a blank element column (Biopython infers the element with rules worth not
# duplicating), a short or ragged ATOM line, or a coordinate field that will not parse. Equality
# with the Biopython result on every ``tests/assets/pdb`` file is asserted by the unit tests.
_ATOM = b"ATOM  "


def _pdb_bytes(path: Path, gzipped: bool) -> bytes:
    opener = gzip.open if gzipped else open
    with opener(path, "rb") as fh:
        return fh.read()


def _parse_pdb_fast(raw: bytes, pdb_id: str, keep_hydrogens: bool) -> Structure | None:
    """Vectorised ATOM-record parse, or ``None`` if this file needs the reference parser."""
    end = raw.find(b"\nENDMDL")                      # first model only, as ``model=0`` asks
    if end != -1:
        raw = raw[:end]
    lines = [ln for ln in raw.split(b"\n") if ln[:6] == _ATOM]
    n = len(lines)
    if not n:
        return None
    lens = list(map(len, lines))                     # one pass; three separate min/max cost more
    width, narrowest = max(lens), min(lens)
    if narrowest < 78:                               # no element column to read
        return None
    # Padding costs an ljust per line, so skip it on the usual file where every record is the
    # same length -- which is what the format specifies and what every writer here emits.
    flat = b"".join(lines) if narrowest == width else b"".join(ln.ljust(width) for ln in lines)
    buf = np.frombuffer(flat, dtype=np.uint8).reshape(n, width)

    def col(a, b, dtype):                            # one fixed column range, one value per line
        return np.frombuffer(buf[:, a:b].tobytes(), dtype=dtype)

    resnames, chain_ids = col(17, 20, "S3"), col(21, 22, "S1")
    resseqs, icodes = col(22, 26, "S4"), col(26, 27, "S1")
    elements = np.char.upper(np.char.strip(col(76, 78, "S2")))
    if (elements == b"").any():                      # Biopython would infer these; do not guess
        return None
    try:
        # float32 first, deliberately. Biopython reads coordinates into a float32 array before
        # anything widens them, so a PDB's three decimals come back as 59.42599869, not 59.426.
        # Parsing the text straight to float64 is the more faithful reading of the file, but it
        # would move every downstream number in its last digits against everything computed
        # before. The gain is ~2e-5 A on a coordinate; matching bit for bit is worth more.
        xyz = np.stack([col(a, a + 8, "S8").astype(np.float32) for a in (30, 38, 46)],
                       axis=1).astype(np.float64)
        resnum = resseqs.astype(np.int64)
    except ValueError:
        return None

    # Atom names and elements repeat: a whole file holds ~40 distinct names over ~5,000 records,
    # so decoding the uniques and indexing turns 1.4 M bytes.decode() calls into a few dozen.
    def decoded(arr, strip=False):
        uniq, inv = np.unique(arr, return_inverse=True)
        table = [(b.decode().strip() if strip else b.decode()) for b in uniq]
        return table, inv

    name_tab, name_ix = decoded(col(12, 16, "S4"), strip=True)
    el_tab, el_ix = decoded(elements)
    keep = np.ones(n, bool) if keep_hydrogens else (elements != b"H")

    # A residue break is any change in (chain, resseq, icode, resname) -- the same test
    # StructureBuilder applies. Waters and any non-amino resname are dropped afterwards.
    same = ((chain_ids[1:] == chain_ids[:-1]) & (resseqs[1:] == resseqs[:-1])
            & (icodes[1:] == icodes[:-1]) & (resnames[1:] == resnames[:-1]))
    starts = np.flatnonzero(~same) + 1
    bounds = list(zip(np.r_[0, starts].tolist(), np.r_[starts, n].tolist()))

    name_ix, el_ix, keep_l = name_ix.tolist(), el_ix.tolist(), keep.tolist()
    by_chain: dict[str, list[Residue]] = {}
    for lo, hi in bounds:
        resname = resnames[lo].decode().strip().upper()
        if resname in _WATER:
            continue
        atoms = tuple(Atom(name=name_tab[name_ix[i]], element=el_tab[el_ix[i]], coord=xyz[i])
                      for i in range(lo, hi) if keep_l[i])
        if not atoms:
            continue
        residues = by_chain.setdefault(chain_ids[lo].decode(), [])
        aa = _one_letter(resname) or "X"
        residues.append(Residue(seq_index=len(residues), pdb_index=int(resnum[lo]),
                                insertion_code=icodes[lo].decode().strip(),
                                aa=aa if len(aa) == 1 else "X", resname=resname, atoms=atoms))
    return Structure(pdb_id=pdb_id,
                     chains=[Chain(chain_id=c, residues=r) for c, r in by_chain.items() if r])


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
    is_cif = _structure_format(path.name) == "cif"
    if not is_cif and model == 0:
        fast = _parse_pdb_fast(_pdb_bytes(path, gzipped), pdb_id, keep_hydrogens)
        if fast is not None:
            return fast
    parser = MMCIFParser(QUIET=True) if is_cif else PDBParser(QUIET=True)

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


def mean_bfactor(path: str | Path, chain: str | None = None) -> float:
    """Mean B-factor over a structure file, or over one chain of it.

    **In a model written by AlphaFold or TCRmodel2 the B-factor column IS the per-residue pLDDT**,
    so this is how a generated structure's own confidence is read back off disk. In a crystal it is
    the crystallographic B-factor and means something entirely different; the caller has to know
    which kind of file it is holding.

    This is *supplied* data — the generator's read-out, not a quantity tcren computes — which is why
    it is a file reader here rather than a descriptor in
    :data:`tcren.recognition.DESCRIPTORS`. :func:`parse_structure` deliberately drops B-factors,
    since they are not part of the geometry the rest of the package reasons about; this exists so
    that reading them does not require a second PDB parser.

    Args:
        path: a ``.pdb`` / ``.pdb.gz`` file. mmCIF is not supported.
        chain: restrict to one author chain id, or ``None`` for every atom in the file.

    Returns:
        The mean, or ``nan`` if the file holds no atom line matching ``chain``.
    """
    import gzip

    path = Path(path)
    opener = gzip.open if path.suffix == ".gz" else open
    total, n = 0.0, 0
    with opener(path, "rt") as fh:
        for line in fh:
            if not line.startswith(("ATOM  ", "HETATM")):
                continue
            if chain is not None and line[21] != chain:
                continue
            try:
                total += float(line[60:66])
                n += 1
            except ValueError:            # a malformed or absent B-factor column
                continue
    return total / n if n else float("nan")


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
    """List structure files for ``src``, sorted.

    ``src`` may be a single structure file, a directory (scanned for structure files), a glob
    pattern (``models/*.pdb.gz``), or a **manifest**: a ``.txt``/``.list``/``.lst`` file holding
    one path per line, ``#`` comments and blank lines ignored, relative paths resolved against
    the manifest's own directory. Recognises plain and gzipped PDB/mmCIF (``.pdb``, ``.cif.gz``,
    …). For archives or streaming, use :func:`iter_structures`.
    """
    text = str(src)
    if any(ch in text for ch in "*?[") and not Path(text).exists():
        return sorted(Path(p) for p in glob.glob(text, recursive=True) if is_structure_file(p))
    src = Path(src)
    if src.is_dir():
        return sorted(p for p in src.iterdir() if is_structure_file(p))
    if src.suffix.lower() in _MANIFEST_SUFFIXES:
        out: list[Path] = []
        for raw in src.read_text().splitlines():
            line = raw.split("#", 1)[0].strip()
            if line:
                out.append(Path(line) if Path(line).is_absolute() else src.parent / line)
        return out
    return [src]


def resolve_sources(sources: str | Path | Iterable[str | Path]) -> list[str]:
    """Split a CLI ``-s`` spec (or several) into individual sources for :func:`iter_structures`.

    Accepts one spec or an iterable of them, and splits each on commas, so
    ``-s a.pdb.gz,b.pdb.gz``, a repeated ``-s``, and a shell glob all mean the same thing.
    Directories, tar archives, globs and manifests are passed through untouched — they are
    expanded downstream, where the tar streaming lives.
    """
    if isinstance(sources, (str, Path)):
        sources = [sources]
    return [part for spec in sources for part in str(spec).split(",") if part.strip()]


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
            members = [m for m in tar.getmembers() if m.isfile() and is_structure_file(m.name)]
            ids = resolve_structure_ids(m.name for m in members)
            for member in members:
                inner, _ = _strip_gz(Path(member.name).name)
                ext = "." + inner.rsplit(".", 1)[-1] + (".gz" if member.name.lower().endswith(".gz") else "")
                fh = tar.extractfile(member)
                if fh is None:
                    continue
                pdb_id = ids[member.name]
                with tempfile.NamedTemporaryFile(suffix=ext) as tmp:
                    tmp.write(fh.read())
                    tmp.flush()
                    s = _safe(Path(tmp.name), pdb_id)
                if s is not None:
                    yield pdb_id, s
        return

    paths = list(structure_paths(src))
    ids = resolve_structure_ids(paths)
    for path in paths:
        s = _safe(path, ids[str(path)])
        if s is not None:
            yield ids[str(path)], s


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


_CIF_COLUMNS = (
    "group_PDB", "id", "type_symbol", "label_atom_id", "label_alt_id", "label_comp_id",
    "label_asym_id", "label_seq_id", "pdbx_PDB_ins_code", "Cartn_x", "Cartn_y", "Cartn_z",
    "occupancy", "B_iso_or_equiv", "auth_seq_id", "auth_asym_id", "pdbx_PDB_model_num",
)


def cif_lines(structure: Structure, transform=None, keep_hydrogens: bool = True) -> list[str]:
    """Minimal mmCIF ``atom_site`` loop for ``structure`` (optionally transformed).

    Same atom selection as :func:`pdb_lines` (one conformer per atom name per residue). Only
    the ``_atom_site`` category is written — enough to round-trip coordinates + chain/residue
    identity through the Biopython MMCIF parser, which is all tcren consumes.
    """
    lines = [f"data_{structure.pdb_id or 'structure'}", "#", "loop_"]
    lines += [f"_atom_site.{c}" for c in _CIF_COLUMNS]
    serial = 1
    for chain in structure.chains:
        chain_id = (chain.chain_id or "A")[0]
        for res in chain.residues:
            seen: set[str] = set()
            for atom in res.atoms:
                element = (atom.element or atom.name[:1]).strip().upper()
                if (not keep_hydrogens and element == "H") or atom.name in seen:
                    continue
                seen.add(atom.name)
                x, y, z = transform(atom.coord) if transform else atom.coord
                icode = (res.insertion_code or "?")[:1] or "?"
                lines.append(
                    f"ATOM {serial} {element} {atom.name} . {res.resname} {chain_id} "
                    f"{res.pdb_index} {icode} {x:.3f} {y:.3f} {z:.3f} 1.00 0.00 "
                    f"{res.pdb_index} {chain_id} 1"
                )
                serial += 1
    lines.append("#")
    return lines


def write_pdb(structure: Structure, path: str | Path, transform=None,
              keep_hydrogens: bool = True) -> Path:
    """Write ``structure`` to a PDB file; return the path.

    A ``.gz`` suffix (``foo.pdb.gz``) transparently gzip-compresses the output.
    """
    return _write_text(path, "\n".join(pdb_lines(structure, transform, keep_hydrogens)) + "\n")


def _write_text(path: str | Path, text: str) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.name.lower().endswith(".gz"):
        with gzip.open(path, "wt") as fh:
            fh.write(text)
    else:
        path.write_text(text)
    return path


def structure_output_path(directory: str | Path, pdb_id: str, mmcif: bool = False,
                          compress: bool = False) -> Path:
    """Build an output path ``<directory>/<pdb_id>.<ext>`` from format flags.

    ``.pdb`` by default, ``.cif`` if ``mmcif``, with a trailing ``.gz`` if ``compress``.
    """
    ext = ".cif" if mmcif else ".pdb"
    return Path(directory) / f"{pdb_id}{ext}{'.gz' if compress else ''}"


def write_structure(structure: Structure, path: str | Path, transform=None,
                    keep_hydrogens: bool = True) -> Path:
    """Format-dispatch writer: PDB or mmCIF, optionally gzipped (by the path suffix)."""
    inner, _ = _strip_gz(Path(path).name)
    if inner.lower().endswith(_CIF_SUFFIXES):
        return _write_text(path, "\n".join(cif_lines(structure, transform, keep_hydrogens)) + "\n")
    return write_pdb(structure, path, transform=transform, keep_hydrogens=keep_hydrogens)
