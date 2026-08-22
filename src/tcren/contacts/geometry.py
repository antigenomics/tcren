"""Atom-level contact and Cα-distance computation.

Ports the legacy ``mir`` ``compute-pdb-contacts`` / ``compute-pdb-geom`` steps using a
``scipy.spatial.cKDTree`` for the all-atom neighbour search.
"""

from __future__ import annotations

import numpy as np
import polars as pl
from scipy.spatial import cKDTree

from ..structure.model import Structure

#: Main-chain heavy atoms. Everything else (CB onward) belongs to the side chain, which is what a
#: residue's one-letter identity stands for -- and what a residue-level potential prices.
BACKBONE_ATOMS = frozenset({"N", "CA", "C", "O", "OXT"})


def _atom_arrays(structure: Structure):
    """Flatten the structure into per-atom coordinate and metadata arrays.

    Hydrogens are dropped. Every contact definition in the package is a *heavy-atom* one — the
    thresholds, the potentials derived from them, and the chemical typing in
    :mod:`tcren.contact_types` all assume it — but ``parse_structure``/``import_structure`` keep
    hydrogens by default, so a model that carries them (AlphaFold relaxed, OpenMM, PDBFixer) used
    to report an H as the closest atom of a pair and defeat the typing entirely.
    :mod:`tcren.clashes` and :mod:`tcren.stability` already filtered here; this brings the contact
    layer in line.
    """
    coords: list[np.ndarray] = []
    chain_ids: list[str] = []
    res_idx: list[int] = []
    res_aa: list[str] = []
    atom_names: list[str] = []
    for chain in structure.chains:
        for res in chain.residues:
            for atom in res.atoms:
                if atom.element == "H":
                    continue
                coords.append(atom.coord)
                chain_ids.append(chain.chain_id)
                res_idx.append(res.seq_index)
                res_aa.append(res.aa)
                atom_names.append(atom.name)
    if not coords:
        return (np.empty((0, 3)), np.array([]), np.array([]), np.array([]), np.array([]))
    return (
        np.asarray(coords, dtype=np.float64),
        np.asarray(chain_ids, dtype=object),
        np.asarray(res_idx, dtype=np.int64),
        np.asarray(res_aa, dtype=object),
        np.asarray(atom_names, dtype=object),
    )


def all_atom_contacts(
    structure: Structure,
    cutoff: float = 5.0,
    count_atoms: bool = False,
    scope: str = "inter",
    atom_pairs: bool = False,
    sidechain: bool = False,
) -> pl.DataFrame:
    """Closest atom contact for each residue pair within ``cutoff`` Å.

    For every pair of residues that have at least one heavy-atom pair within ``cutoff``
    (inclusive, matching the legacy ``dist <= 5``), the row with the minimum atom–atom
    distance is kept.

    Args:
        structure: The (parsed) structure.
        cutoff: Contact distance threshold (Å, inclusive).
        count_atoms: When ``True`` add an extra ``n_atom_contacts`` column = the count
            of heavy-atom pairs within ``cutoff`` for that residue pair (always ``>= 1``
            for any kept row, and ``>=`` the single closest-atom row this function keeps).
            Default ``False`` keeps the schema and every value byte-identical to the
            legacy output.
        scope: Which residue pairs to keep. ``"inter"`` (default) keeps only pairs on
            different chains, which is what every interface score in the package is
            built on and what the legacy output contained. ``"intra"`` keeps only pairs
            within one chain, and ``"all"`` keeps both. Intra-chain pairs include
            sequence neighbours, which are in contact by covalent geometry rather than
            by folding — filter on sequence separation before interpreting them
            (:func:`peptide_internal_contacts` does this).
        sidechain: When ``True`` add ``sc.from``, ``sc.to`` and ``n_sc_pairs``: whether each side
            puts a **side-chain** heavy atom within ``cutoff`` of the other residue, and how many
            of the residue pair's atom pairs are side-chain to side-chain. A 5 Å residue pair whose
            only atoms in range are the two backbones is not an interaction between those two
            residue *identities*, and a residue-level potential that prices it is charging for
            chemistry that is not happening. The flags are computed over **all** atom pairs, so
            they survive the collapse to the closest one.
        atom_pairs: When ``True`` return **every** heavy-atom pair within ``cutoff`` instead of
            collapsing each residue pair to its closest one. Chemical typing needs this: a salt
            bridge whose nearest atom pair happens to be Arg CG–Asp OD1 is invisible in the
            collapsed table, because the charged NH1–OD2 pair it also makes was discarded.
            Changes the row count, not the schema.

    Returns:
        Columns: ``chain.id.from``, ``residue.index.from``, ``chain.id.to``,
        ``residue.index.to``, ``residue.aa.from``, ``residue.aa.to``,
        ``atom.from``, ``atom.to``, ``dist`` (plus ``n_atom_contacts`` when
        ``count_atoms`` is set). Each unordered residue pair appears once,
        in ``(chain.id, residue.index)`` lexicographic order.
    """
    coords, chain_ids, res_idx, res_aa, atom_names = _atom_arrays(structure)
    schema = {
        "chain.id.from": pl.Utf8,
        "residue.index.from": pl.Int64,
        "chain.id.to": pl.Utf8,
        "residue.index.to": pl.Int64,
        "residue.aa.from": pl.Utf8,
        "residue.aa.to": pl.Utf8,
        "atom.from": pl.Utf8,
        "atom.to": pl.Utf8,
        "dist": pl.Float64,
    }
    if count_atoms:
        schema["n_atom_contacts"] = pl.UInt32
    if sidechain:
        schema["sc.from"] = pl.Boolean
        schema["sc.to"] = pl.Boolean
        schema["n_sc_pairs"] = pl.UInt32
    if len(coords) == 0:
        return pl.DataFrame(schema=schema)

    tree = cKDTree(coords)
    pairs = tree.query_pairs(r=cutoff, output_type="ndarray")
    if len(pairs) == 0:
        return pl.DataFrame(schema=schema)

    i, j = pairs[:, 0], pairs[:, 1]
    if scope not in ("inter", "intra", "all"):
        raise ValueError(f"scope must be 'inter', 'intra' or 'all', got {scope!r}")
    if scope != "all":
        same_chain = chain_ids[i] == chain_ids[j]
        keep = same_chain if scope == "intra" else ~same_chain
        i, j = i[keep], j[keep]
    if scope != "inter":
        # Two atoms of one residue are not a contact between residues.
        distinct = (chain_ids[i] != chain_ids[j]) | (res_idx[i] != res_idx[j])
        i, j = i[distinct], j[distinct]
    if len(i) == 0:
        return pl.DataFrame(schema=schema)

    dist = np.linalg.norm(coords[i] - coords[j], axis=1)

    # Orient each pair canonically by (chain.id, residue.index) so the from/to labelling
    # is deterministic and symmetric duplicates collapse.
    ci, cj = chain_ids[i], chain_ids[j]
    ri, rj = res_idx[i], res_idx[j]
    swap = (ci > cj) | ((ci == cj) & (ri > rj))

    def pick(a, b):
        return np.where(swap, b, a)

    df = pl.DataFrame(
        {
            "chain.id.from": pick(ci, cj),
            "residue.index.from": pick(ri, rj),
            "chain.id.to": pick(cj, ci),
            "residue.index.to": pick(rj, ri),
            "residue.aa.from": pick(res_aa[i], res_aa[j]),
            "residue.aa.to": pick(res_aa[j], res_aa[i]),
            "atom.from": pick(atom_names[i], atom_names[j]),
            "atom.to": pick(atom_names[j], atom_names[i]),
            "dist": dist,
        }
    )
    group_keys = [
        "chain.id.from", "residue.index.from", "chain.id.to", "residue.index.to",
    ]
    if count_atoms:
        # One row per atom-atom pair, so the group size is the heavy-atom-pair count.
        df = df.with_columns(pl.len().over(group_keys).alias("n_atom_contacts"))
    if sidechain:
        # Aggregated over the group before any collapse, so the flags describe the whole residue
        # pair rather than whichever atom pair happens to be closest.
        bb = list(BACKBONE_ATOMS)
        scf, sct = ~pl.col("atom.from").is_in(bb), ~pl.col("atom.to").is_in(bb)
        df = df.with_columns(
            scf.any().over(group_keys).alias("sc.from"),
            sct.any().over(group_keys).alias("sc.to"),
            (scf & sct).sum().over(group_keys).cast(pl.UInt32).alias("n_sc_pairs"),
        )
    if not atom_pairs:
        # Keep the closest atom pair per residue pair.
        df = df.sort("dist").group_by(group_keys, maintain_order=True).first()
    return df.select(list(schema.keys()))


def peptide_internal_contacts(
    structure: Structure,
    cutoff: float = 5.0,
    min_seq_sep: int = 3,
    count_atoms: bool = True,
) -> pl.DataFrame:
    """Contacts *within* the peptide chain — the term every interface score omits.

    Every energy in this package sums over inter-chain contacts only, so a peptide that
    is held in a particular conformation by its own side chains scores the same as one
    that is not. This returns those omitted pairs, so that assumption can be tested
    rather than inherited.

    Sequence neighbours are in contact because they are bonded, not because the peptide
    folded that way, so pairs closer than ``min_seq_sep`` in sequence are dropped; the
    default of 3 keeps ``i``/``i+3`` and beyond, which is the shortest separation that
    can carry a side-chain-to-side-chain interaction across a turn. The 5.0 Å default is
    the same contact definition the rest of the package uses, so an internal contact and
    an interface contact mean the same thing.

    Args:
        structure: The (parsed, annotated) structure; the chain typed ``PEPTIDE`` is used.
        cutoff: Contact distance threshold (Å, inclusive).
        min_seq_sep: Minimum ``|i - j|`` in residue index for a pair to be kept.
        count_atoms: Carry the ``n_atom_contacts`` heavy-atom-pair count (default ``True``
            here, since the count is the quantity of interest for an internal contact).

    Returns:
        The :func:`all_atom_contacts` schema, restricted to the peptide chain.
        Empty (with schema) when the structure has no ``PEPTIDE`` chain.
    """
    peptide = next((c for c in structure.chains if c.chain_type == "PEPTIDE"), None)
    contacts = all_atom_contacts(
        structure, cutoff=cutoff, count_atoms=count_atoms, scope="intra"
    )
    if peptide is None:
        return contacts.clear()
    return contacts.filter(
        (pl.col("chain.id.from") == peptide.chain_id)
        & (pl.col("chain.id.to") == peptide.chain_id)
        & ((pl.col("residue.index.to") - pl.col("residue.index.from")).abs() >= min_seq_sep)
    )


def _representative_arrays(structure: Structure, kind: str):
    """Per-residue representative-atom arrays (one atom per residue).

    ``kind`` is ``"ca"`` (Cα) or ``"cb"`` (Cβ, falling back to Cα for glycine / missing Cβ).
    Residues with no representative atom are skipped.
    """
    coords, chain_ids, res_idx, res_aa = [], [], [], []
    for chain in structure.chains:
        for res in chain.residues:
            rep = res.cb_or_ca if kind == "cb" else res.ca
            if rep is None:
                continue
            coords.append(rep)
            chain_ids.append(chain.chain_id)
            res_idx.append(res.seq_index)
            res_aa.append(res.aa)
    if not coords:
        return np.empty((0, 3)), np.array([]), np.array([]), np.array([])
    return (
        np.asarray(coords, dtype=np.float64),
        np.asarray(chain_ids, dtype=object),
        np.asarray(res_idx, dtype=np.int64),
        np.asarray(res_aa, dtype=object),
    )


def representative_atom_contacts(
    structure: Structure, kind: str = "ca", cutoff: float = 12.0
) -> pl.DataFrame:
    """Inter-chain residue contacts by a single representative atom per residue.

    ``kind="ca"`` uses Cα (default cutoff 12 Å); ``kind="cb"`` uses Cβ with a glycine/
    missing-Cβ fallback to Cα (default cutoff 8 Å). Mirrors :func:`all_atom_contacts`'
    residue-pair schema (``atom.from``/``atom.to`` carry the representative atom kind).
    """
    coords, chain_ids, res_idx, res_aa = _representative_arrays(structure, kind)
    rep = "CB" if kind == "cb" else "CA"
    schema = {
        "chain.id.from": pl.Utf8, "residue.index.from": pl.Int64,
        "chain.id.to": pl.Utf8, "residue.index.to": pl.Int64,
        "residue.aa.from": pl.Utf8, "residue.aa.to": pl.Utf8,
        "atom.from": pl.Utf8, "atom.to": pl.Utf8, "dist": pl.Float64,
    }
    if len(coords) == 0:
        return pl.DataFrame(schema=schema)
    tree = cKDTree(coords)
    pairs = tree.query_pairs(r=cutoff, output_type="ndarray")
    if len(pairs) == 0:
        return pl.DataFrame(schema=schema)
    i, j = pairs[:, 0], pairs[:, 1]
    inter = chain_ids[i] != chain_ids[j]
    i, j = i[inter], j[inter]
    if len(i) == 0:
        return pl.DataFrame(schema=schema)
    dist = np.linalg.norm(coords[i] - coords[j], axis=1)
    ci, cj = chain_ids[i], chain_ids[j]
    ri, rj = res_idx[i], res_idx[j]
    swap = (ci > cj) | ((ci == cj) & (ri > rj))

    def pick(a, b):
        return np.where(swap, b, a)

    return pl.DataFrame({
        "chain.id.from": pick(ci, cj), "residue.index.from": pick(ri, rj),
        "chain.id.to": pick(cj, ci), "residue.index.to": pick(rj, ri),
        "residue.aa.from": pick(res_aa[i], res_aa[j]),
        "residue.aa.to": pick(res_aa[j], res_aa[i]),
        "atom.from": rep, "atom.to": rep, "dist": dist,
    }).select(list(schema.keys()))


def ca_distance_matrix(structure: Structure) -> tuple[np.ndarray, list[tuple[str, int]]]:
    """Pairwise Cα–Cα distance matrix over all residues with a Cα atom.

    Returns:
        ``(matrix, keys)`` where ``matrix[a, b]`` is the Cα distance and ``keys[a]`` is
        the ``(chain_id, seq_index)`` of row/column ``a``.
    """
    coords: list[np.ndarray] = []
    keys: list[tuple[str, int]] = []
    for chain in structure.chains:
        for res in chain.residues:
            ca = res.ca
            if ca is not None:
                coords.append(ca)
                keys.append((chain.chain_id, res.seq_index))
    if not coords:
        return np.empty((0, 0)), []
    arr = np.asarray(coords)
    diff = arr[:, None, :] - arr[None, :, :]
    return np.linalg.norm(diff, axis=2), keys
