"""Residue-pair opportunities for the recorded native contact model."""

from __future__ import annotations

from dataclasses import replace

import numpy as np
import polars as pl

from ..contactmap import ContactMap
from ..structure.model import MHC_TYPES, PEPTIDE_TYPE, TCR_TYPES, Structure
from .geometry import all_atom_contacts, representative_atom_contacts
from .table import directed_views, residue_annotation

_ROLES = (*TCR_TYPES, PEPTIDE_TYPE, *MHC_TYPES)
_INTERFACES = ("tcr_peptide", "tcr_mhc", "peptide_mhc")
_KEYS = [f"{field}.{side}" for side in ("from", "to") for field in ("chain.id", "residue.index")]
_JOIN_KEYS = [*_KEYS, "residue.aa.from", "residue.aa.to"]


def _native_structure(structure: Structure) -> Structure:
    chains = []
    for chain in structure.chains:
        if chain.chain_type not in _ROLES:
            continue
        residues = []
        for residue in chain.residues:
            atoms = tuple(a for a in residue.atoms if a.element.strip().upper() not in {"H", "D"})
            r = residue if len(atoms) == len(residue.atoms) else replace(residue, atoms=atoms)
            anchor = r.cb_or_ca
            if anchor is not None and np.isfinite(anchor).all():
                residues.append(r)
        chains.append(replace(chain, residues=residues))
    return replace(structure, chains=chains)


def native_atom_contacts(structure: Structure) -> pl.DataFrame:
    """Return all directed 5 Å heavy-atom contacts on resolved native-model endpoints."""
    contacts = ContactMap.from_structure(_native_structure(structure), atom_pairs=True)
    return pl.concat([
        contacts.interface(interface).with_columns(pl.lit(interface).alias("interface"))
        for interface in _INTERFACES
    ])


def native_opportunities(structure: Structure) -> pl.DataFrame:
    """Return the union of 12 Å Cβ/Cα neighbours and 5 Å heavy-atom contacts.

    Both residues must have a finite Cβ, or Cα where Cβ is absent. The table
    keeps one biologically directed row per pair, including non-contacting
    neighbours and contacts outside the neighbour shell.
    """
    selected = _native_structure(structure)
    observed = [r for c in selected.chains for r in c.residues]
    anchors = [r.cb_or_ca for r in observed]
    heavy_counts = [len(r.atoms) for r in observed]
    shell = representative_atom_contacts(selected, kind="cb", cutoff=12.0).select(
        *_JOIN_KEYS, pl.lit(True).alias("within_shell"))
    actual = all_atom_contacts(selected, cutoff=5.0, count_atoms=True).rename({
        "n_atom_contacts": "n_atom_pairs", "dist": "min_dist",
    }).select(*_JOIN_KEYS, "atom.from", "atom.to", "n_atom_pairs", "min_dist")
    edges = directed_views(shell.join(actual, on=_JOIN_KEYS, how="full", coalesce=True))
    residues = residue_annotation(selected).drop("residue.aa").with_row_index("anchor").with_columns(
        pl.Series("n_heavy", heavy_counts, dtype=pl.UInt32))
    for side in ("from", "to"):
        meta = residues.rename({c: f"{c}.{side}" for c in residues.columns})
        edges = edges.join(meta, on=[f"chain.id.{side}", f"residue.index.{side}"], how="inner")
    coords = np.asarray(anchors, dtype=float).reshape(-1, 3)
    distances = np.linalg.norm(coords[edges["anchor.from"].to_numpy()]
                               - coords[edges["anchor.to"].to_numpy()], axis=1)
    edges = edges.with_columns(
        pl.Series("cb_distance", distances),
        pl.lit(structure.pdb_id).alias("pdb.id"),
        pl.col("within_shell").fill_null(False),
        pl.col("min_dist").is_not_null().alias("contacted"),
        pl.col("n_atom_pairs").fill_null(0),
    )
    contacts = ContactMap(structure.pdb_id, edges)
    return pl.concat([
        contacts.interface(interface).select(
            "pdb.id", pl.lit(interface).alias("interface"), *_KEYS,
            "chain.type.from", "chain.type.to", "region.type.from", "region.type.to",
            pl.col("residue.aa.from").alias("aa.from"), pl.col("residue.aa.to").alias("aa.to"),
            "cb_distance", "within_shell", "contacted", "n_atom_pairs", "min_dist",
            "atom.from", "atom.to", "n_heavy.from", "n_heavy.to",
        ).rename({"n_heavy.from": "n_heavy_from", "n_heavy.to": "n_heavy_to"})
        for interface in _INTERFACES
    ])
