"""Chemical typing of interface contacts — a DSSP-style annotation layer for TCR:pMHC contact maps.

Classifies each atom–atom contact in a :class:`~tcren.contactmap.ContactMap` interface into one chemical
type from heavy-atom geometry alone (no hydrogens, no external DSSP binary — models and many crystals lack
explicit H), by priority:

* ``salt_bridge``  — cationic (Lys/Arg/His N) ↔ anionic (Asp/Glu O), ≤ 4.0 Å
* ``hydrogen_bond``— two polar N/O atoms, ≤ 3.5 Å
* ``aromatic``     — two aromatic-ring atoms of aromatic residues (Phe/Tyr/Trp/His), ≤ 5.0 Å
* ``hydrophobic``  — two carbon atoms of apolar residues, ≤ 4.5 Å
* ``other``        — anything else within the contact-map cutoff

This replaces ad-hoc one-off H-bond counters: :func:`contact_type_counts` gives per-type contact and
residue-pair counts (the documented, reproducible source of an ``n_hbond``-style feature).
"""
from __future__ import annotations

import polars as pl

# residue single-letter sets
_CATIONIC_RES = {"K", "R", "H"}
_ANIONIC_RES = {"D", "E"}
_AROMATIC_RES = {"F", "Y", "W", "H"}
_APOLAR_RES = {"A", "V", "L", "I", "M", "F", "W", "P", "C"}

# side-chain atom names carrying formal charge
_CATIONIC_ATOMS = {"NZ", "NE", "NH1", "NH2", "ND1", "NE2"}
_ANIONIC_ATOMS = {"OD1", "OD2", "OE1", "OE2"}
# aromatic ring atoms (Phe/Tyr/Trp/His)
_RING_ATOMS = {"CG", "CD1", "CD2", "CE1", "CE2", "CZ", "NE1", "CE3", "CZ2", "CZ3", "CH2", "ND1", "NE2"}

_HBOND_MAX = 3.5
_SALT_MAX = 4.0
_AROMATIC_MAX = 5.0
_HYDROPHOBIC_MAX = 4.5

TYPES = ("salt_bridge", "hydrogen_bond", "aromatic", "hydrophobic", "other")


def _elem(atom: str) -> str:
    """Element from a PDB atom name (first alphabetic character)."""
    for ch in atom:
        if ch.isalpha():
            return ch
    return "?"


def _classify(aa_a: str, aa_b: str, atom_a: str, atom_b: str, dist: float) -> str:
    ea, eb = _elem(atom_a), _elem(atom_b)
    # salt bridge: one cationic side-chain N, one anionic side-chain O
    cat = (aa_a in _CATIONIC_RES and atom_a in _CATIONIC_ATOMS, aa_b in _CATIONIC_RES and atom_b in _CATIONIC_ATOMS)
    ani = (aa_a in _ANIONIC_RES and atom_a in _ANIONIC_ATOMS, aa_b in _ANIONIC_RES and atom_b in _ANIONIC_ATOMS)
    if dist <= _SALT_MAX and ((cat[0] and ani[1]) or (cat[1] and ani[0])):
        return "salt_bridge"
    # hydrogen bond: two polar (N/O) heavy atoms
    if dist <= _HBOND_MAX and ea in ("N", "O") and eb in ("N", "O"):
        return "hydrogen_bond"
    # aromatic: ring atoms of aromatic residues on both sides
    if (dist <= _AROMATIC_MAX and aa_a in _AROMATIC_RES and aa_b in _AROMATIC_RES
            and atom_a in _RING_ATOMS and atom_b in _RING_ATOMS):
        return "aromatic"
    # hydrophobic: carbon–carbon between apolar residues
    if dist <= _HYDROPHOBIC_MAX and ea == "C" and eb == "C" and aa_a in _APOLAR_RES and aa_b in _APOLAR_RES:
        return "hydrophobic"
    return "other"


def classify_contacts(interface_df: pl.DataFrame) -> pl.DataFrame:
    """Return ``interface_df`` with an added ``contact.type`` column (one of :data:`TYPES` per row).

    Args:
        interface_df: an interface frame from :meth:`ContactMap.interface`, carrying ``residue.aa.from/to``,
            ``atom.from/to`` and ``dist``.

    Returns:
        The same frame with a ``contact.type`` string column.
    """
    if interface_df.height == 0:
        return interface_df.with_columns(pl.lit(None, dtype=pl.Utf8).alias("contact.type"))
    types = [
        _classify(aa, ab, ta, tb, d)
        for aa, ab, ta, tb, d in zip(
            interface_df["residue.aa.from"].to_list(), interface_df["residue.aa.to"].to_list(),
            interface_df["atom.from"].to_list(), interface_df["atom.to"].to_list(),
            interface_df["dist"].to_list(),
        )
    ]
    return interface_df.with_columns(pl.Series("contact.type", types))


def contact_type_counts(cm, interface: str = "tcr_peptide", tcr_regions: str = "all") -> dict[str, int]:
    """Per-type contact counts + distinct residue-pair counts for one interface.

    Args:
        cm: a :class:`~tcren.contactmap.ContactMap`.
        interface: interface name (``"tcr_peptide"``, ``"tcr_mhc"``, ``"peptide_mhc"``).
        tcr_regions: passed through to :meth:`ContactMap.interface`.

    Returns:
        Mapping with ``n_<type>`` (atom-pair contacts of each type) and ``pairs_<type>`` (distinct
        residue–residue pairs having ≥1 contact of that type), e.g. ``pairs_hydrogen_bond`` is the
        documented ``n_hbond`` feature.
    """
    df = classify_contacts(cm.interface(interface, tcr_regions=tcr_regions))
    out = {f"n_{t}": 0 for t in TYPES} | {f"pairs_{t}": 0 for t in TYPES}
    if df.height == 0:
        return out
    for t in TYPES:
        sub = df.filter(pl.col("contact.type") == t)
        out[f"n_{t}"] = sub.height
        if sub.height:
            out[f"pairs_{t}"] = sub.select(
                ["chain.id.from", "residue.index.from", "chain.id.to", "residue.index.to"]
            ).unique().height
    return out
