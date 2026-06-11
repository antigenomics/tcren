"""Contact-geometry parity against the legacy contact_maps_PDB.csv oracle.

Validates that the Python contact computation reproduces mir's TCR↔peptide contact set
(chain id, sequential residue index, and amino acid on both sides) for representative
structures, including the tricky edge cases:

* ``5m01`` — baseline mouse MHC-I complex.
* ``1ao7`` — baseline human MHC-I complex.
* ``5jhd`` — peptide with a non-standard N-terminal cap (AMN, kept as ``X`` at index 0).
* ``6v0y`` — peptide with internal citrulline (CIR) HETATM residues that mir skips.
* ``7qpj`` — structure with explicit hydrogens (contacts mediated by H atoms).
* ``9nmx`` — contact mediated by an alternate (altloc) conformer.

The full 312-structure sweep runs only under ``RUN_BENCHMARK=1``.
"""

from __future__ import annotations

import os
from pathlib import Path

import polars as pl
import pytest

from tcren.contacts import all_atom_contacts
from tcren.structure import parse_structure

REPO = Path(__file__).resolve().parents[2]
PDB_DIR = REPO / "data" / "PDB_structures"
CONTACT_MAPS = REPO / "data" / "contact_maps_PDB.csv"

_KEYS = [
    "chain.id.from",
    "residue.index.from",
    "residue.index.to",
    "residue.aa.from",
    "residue.aa.to",
]


def _tcr_peptide_contact_set(pdb_id: str, oracle: pl.DataFrame) -> set[tuple]:
    """Compute the oriented TCR→peptide contact set for a structure.

    Chain roles (which chains are TCR vs peptide) are taken from the oracle, since
    chain typing via arda arrives in a later milestone; this isolates the geometry.
    """
    orc = oracle.filter(pl.col("pdb.id") == pdb_id)
    tcr = set(orc["chain.id.from"].to_list())
    pep = set(orc["chain.id.to"].to_list())
    structure = parse_structure(PDB_DIR / f"{pdb_id}.pdb")
    con = all_atom_contacts(structure, cutoff=5.0)

    forward = con.filter(
        pl.col("chain.id.from").is_in(tcr) & pl.col("chain.id.to").is_in(pep)
    ).select(_KEYS)
    backward = con.filter(
        pl.col("chain.id.from").is_in(pep) & pl.col("chain.id.to").is_in(tcr)
    ).select(
        pl.col("chain.id.to").alias("chain.id.from"),
        pl.col("residue.index.to").alias("residue.index.from"),
        pl.col("residue.index.from").alias("residue.index.to"),
        pl.col("residue.aa.to").alias("residue.aa.from"),
        pl.col("residue.aa.from").alias("residue.aa.to"),
    )
    got = pl.concat([forward, backward]).unique()
    return set(map(tuple, got.rows()))


def _oracle_set(pdb_id: str, oracle: pl.DataFrame) -> set[tuple]:
    return set(
        map(tuple, oracle.filter(pl.col("pdb.id") == pdb_id).select(_KEYS).unique().rows())
    )


@pytest.mark.parametrize("pdb_id", ["5m01", "1ao7", "5jhd", "6v0y", "7qpj", "9nmx"])
def test_tcr_peptide_contacts_match_oracle(pdb_id):
    oracle = pl.read_csv(CONTACT_MAPS)
    assert _tcr_peptide_contact_set(pdb_id, oracle) == _oracle_set(pdb_id, oracle)


@pytest.mark.skipif(not os.getenv("RUN_BENCHMARK"), reason="set RUN_BENCHMARK=1 to run")
def test_all_structures_contacts_match_oracle():
    oracle = pl.read_csv(CONTACT_MAPS)
    mismatched = []
    for pdb_id in oracle["pdb.id"].unique().to_list():
        if not (PDB_DIR / f"{pdb_id}.pdb").exists():
            continue
        if _tcr_peptide_contact_set(pdb_id, oracle) != _oracle_set(pdb_id, oracle):
            mismatched.append(pdb_id)
    assert not mismatched, f"{len(mismatched)} structures mismatched: {mismatched}"
