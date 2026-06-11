"""Canonical polars tables for complementarity maps: residue markup + contacts.

Both tables are tidy and keyed for joins on ``(structure_id, structure_chain, aa_index)``
where ``aa_index`` is the 0-based position in the chain's protein sequence (the legacy
mir ``residue.index``); ``residue_index`` is the author/PDB numbering (a display label).
The contacts table wraps the parity-preserving :func:`all_atom_contacts`, so it is
compatible with the original TCRen contact calculation.
"""

from __future__ import annotations

import polars as pl

from ..contacts.geometry import all_atom_contacts, ca_distance_matrix
from ..structure.model import Structure
from .frame import ProjectionResult

# chain_type -> complex_chain label
_COMPLEX_CHAIN = {
    "TRA": "tra", "TRB": "trb", "TRD": "trd", "TRG": "trg",
    "IGH": "trb", "IGK": "tra", "IGL": "tra",  # antibody mimics mapped to the TCR slots
    "PEPTIDE": "peptide", "MHCa": "mhca", "MHCb": "mhcb", "B2M": "b2m",
}
_COMPLEX_CHAINS = set(_COMPLEX_CHAIN)

# Backbone atom names (for the backbone vs sidechain flag).
_BACKBONE = {"N", "CA", "C", "O"}

# Charged side-chain atoms for the salt-bridge test.
_ACIDIC = {("D", "OD1"), ("D", "OD2"), ("E", "OE1"), ("E", "OE2")}
_BASIC = {("K", "NZ"), ("R", "NH1"), ("R", "NH2"), ("R", "NE"),
          ("H", "ND1"), ("H", "NE2")}
_AROMATIC_AA = {"F", "W", "Y", "H"}


def _complex_region(region_type: str | None) -> str | None:
    """Normalise a region type to the complex_region vocabulary."""
    if region_type is None:
        return None
    if region_type.startswith("HELIX_"):
        return "mhc_" + region_type.lower()  # HELIX_A1 -> mhc_helix_a1
    return region_type.lower()  # CDR1->cdr1, FR1->fr1, GROOVE_FLOOR->groove_floor, PEPTIDE->peptide


def _region_of_chain(chain) -> dict[int, str]:
    return {r.seq_index: reg.region_type for reg in chain.regions for r in reg.residues}


def residue_markup_table(
    structure: Structure, projection: ProjectionResult | None = None
) -> pl.DataFrame:
    """Per-residue markup with chain/region annotation and groove-frame coordinates.

    Args:
        structure: A chain-typed, annotated structure.
        projection: Optional :class:`ProjectionResult`; its groove-frame coordinates are
            joined on ``(chain_id, seq_index)`` (``x/y/z`` and in-plane ``u/v``/``height``).

    Returns:
        Columns: ``structure_id, structure_chain, complex_chain, complex_region,
        residue_index, aa_index, aa_len, aa, x, y, z, u, v, height``.
    """
    coords = {}
    if projection is not None:
        for key, xyz in zip(projection.keys, projection.coords3d):
            coords[key] = xyz

    rows = []
    for chain in structure.chains:
        if chain.chain_type not in _COMPLEX_CHAINS:
            continue
        region_of = _region_of_chain(chain)
        aa_len = len(chain.residues)
        for res in chain.residues:
            xyz = coords.get((chain.chain_id, res.seq_index))
            rows.append(
                {
                    "structure_id": structure.pdb_id,
                    "structure_chain": chain.chain_id,
                    "complex_chain": _COMPLEX_CHAIN[chain.chain_type],
                    "complex_region": _complex_region(region_of.get(res.seq_index)),
                    "residue_index": res.pdb_index,
                    "aa_index": res.seq_index,
                    "aa_len": aa_len,
                    "aa": res.aa,
                    "x": float(xyz[0]) if xyz is not None else None,
                    "y": float(xyz[1]) if xyz is not None else None,
                    "z": float(xyz[2]) if xyz is not None else None,
                    "u": float(xyz[0]) if xyz is not None else None,
                    "v": float(xyz[1]) if xyz is not None else None,
                    "height": float(xyz[2]) if xyz is not None else None,
                }
            )
    schema = {
        "structure_id": pl.Utf8, "structure_chain": pl.Utf8, "complex_chain": pl.Utf8,
        "complex_region": pl.Utf8, "residue_index": pl.Int64, "aa_index": pl.Int64,
        "aa_len": pl.Int64, "aa": pl.Utf8, "x": pl.Float64, "y": pl.Float64,
        "z": pl.Float64, "u": pl.Float64, "v": pl.Float64, "height": pl.Float64,
    }
    return pl.DataFrame(rows, schema=schema)


def ca_contacts_table(structure: Structure, threshold: float = 8.0) -> pl.DataFrame:
    """Inter-chain Cα–Cα "chain contacts" within ``threshold`` Å (bold edges in the map).

    Columns: ``structure_id, structure_chain_1, structure_chain_2, aa_index_1, aa_index_2,
    ca_dist``. Complements the all-atom :func:`contacts_table` (which gives the dashed
    inter-residue edges).
    """
    matrix, keys = ca_distance_matrix(structure)
    rows = []
    n = len(keys)
    for i in range(n):
        c1, r1 = keys[i]
        for j in range(i + 1, n):
            c2, r2 = keys[j]
            if c1 == c2:
                continue
            d = matrix[i, j]
            if d <= threshold:
                rows.append(
                    {
                        "structure_id": structure.pdb_id,
                        "structure_chain_1": c1, "structure_chain_2": c2,
                        "aa_index_1": r1, "aa_index_2": r2, "ca_dist": float(d),
                    }
                )
    schema = {
        "structure_id": pl.Utf8, "structure_chain_1": pl.Utf8, "structure_chain_2": pl.Utf8,
        "aa_index_1": pl.Int64, "aa_index_2": pl.Int64, "ca_dist": pl.Float64,
    }
    return pl.DataFrame(rows, schema=schema)


def _element(atom_name: str) -> str:
    return atom_name[0] if atom_name else ""


def classify_contact(aa1: str, aa2: str, atom1: str, atom2: str, dist: float) -> str:
    """Heuristically classify a residue–residue contact from its closest atom pair.

    Returns one of ``salt_bridge``, ``hydrogen_bond``, ``aromatic``, ``hydrophobic``,
    ``polar``, ``other``. Cutoffs are pragmatic, not a force field — documented and kept
    in this pure function for easy tuning.
    """
    a1, a2 = (aa1, atom1), (aa2, atom2)
    if dist <= 4.0 and (
        (a1 in _ACIDIC and a2 in _BASIC) or (a1 in _BASIC and a2 in _ACIDIC)
    ):
        return "salt_bridge"
    e1, e2 = _element(atom1), _element(atom2)
    if dist <= 3.5 and e1 in ("N", "O") and e2 in ("N", "O"):
        return "hydrogen_bond"
    if aa1 in _AROMATIC_AA and aa2 in _AROMATIC_AA and e1 == "C" and e2 == "C":
        return "aromatic"
    if e1 == "C" and e2 == "C":
        return "hydrophobic"
    if "N" in (e1, e2) or "O" in (e1, e2):
        return "polar"
    return "other"


def contacts_table(structure: Structure, threshold: float = 5.0) -> pl.DataFrame:
    """Inter-chain residue contacts within ``threshold`` Å, classified by bond type.

    Args:
        structure: The structure.
        threshold: Distance cutoff in Å (``3 ≤ threshold ≤ 6``, default 5 — the original
            TCRen contact distance).

    Returns:
        Columns: ``structure_id, structure_chain_1, structure_chain_2, residue_index_1,
        residue_index_2, aa_index_1, aa_index_2, min_dist, contact_type, backbone_1,
        backbone_2``. One row per unordered inter-chain residue pair (TCRen-compatible).
    """
    if not (3.0 <= threshold <= 6.0):
        raise ValueError(f"threshold must be in [3, 6] Å, got {threshold}")
    raw = all_atom_contacts(structure, cutoff=threshold)

    # seq_index -> pdb_index lookup per chain (for the author-numbering display column).
    pdb_of = {
        (c.chain_id, r.seq_index): r.pdb_index for c in structure.chains for r in c.residues
    }

    rows = []
    for row in raw.iter_rows(named=True):
        c1, c2 = row["chain.id.from"], row["chain.id.to"]
        i1, i2 = row["residue.index.from"], row["residue.index.to"]
        rows.append(
            {
                "structure_id": structure.pdb_id,
                "structure_chain_1": c1,
                "structure_chain_2": c2,
                "residue_index_1": pdb_of.get((c1, i1)),
                "residue_index_2": pdb_of.get((c2, i2)),
                "aa_index_1": i1,
                "aa_index_2": i2,
                "min_dist": row["dist"],
                "contact_type": classify_contact(
                    row["residue.aa.from"], row["residue.aa.to"],
                    row["atom.from"], row["atom.to"], row["dist"],
                ),
                "backbone_1": row["atom.from"] in _BACKBONE,
                "backbone_2": row["atom.to"] in _BACKBONE,
            }
        )
    schema = {
        "structure_id": pl.Utf8, "structure_chain_1": pl.Utf8, "structure_chain_2": pl.Utf8,
        "residue_index_1": pl.Int64, "residue_index_2": pl.Int64,
        "aa_index_1": pl.Int64, "aa_index_2": pl.Int64, "min_dist": pl.Float64,
        "contact_type": pl.Utf8, "backbone_1": pl.Boolean, "backbone_2": pl.Boolean,
    }
    return pl.DataFrame(rows, schema=schema)
