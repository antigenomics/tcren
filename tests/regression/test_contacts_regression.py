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
PDB_DIR = REPO / "tests" / "assets" / "pdb"
CONTACT_MAPS = REPO / "tests" / "assets" / "oracle" / "data" / "contact_maps_PDB.csv"

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


# Structures deposited with explicit hydrogens. The legacy pipeline counted H-mediated pairs as
# residue contacts, so its contact set for these is a superset of ours: tcren now filters hydrogens
# in `_atom_arrays`, because otherwise the same complex scores differently depending only on whether
# the depositor modelled H (5jhd: 7 of 28 TCR:peptide contacts and -58.5% on Phi_tcr_pep; 7qpj: 8 of 33
# and +38.6%). Parity is asserted on the heavy-atom subset, which is what both pipelines mean by a
# 5 Å contact.
_HAS_HYDROGENS = {"5jhd", "7qpj"}


def _has_hydrogens(path) -> bool:
    return any(line.startswith(("ATOM", "HETATM")) and line[76:78].strip() == "H"
               for line in path.read_text().splitlines())


@pytest.mark.parametrize("pdb_id", ["5m01", "1ao7", "5jhd", "6v0y", "7qpj", "9nmx"])
def test_tcr_peptide_contacts_match_oracle(pdb_id):
    oracle = pl.read_csv(CONTACT_MAPS)
    got, want = _tcr_peptide_contact_set(pdb_id, oracle), _oracle_set(pdb_id, oracle)
    if pdb_id in _HAS_HYDROGENS:
        assert got <= want, "heavy-atom contacts must be a subset of the legacy H-inclusive set"
        assert len(got) >= 0.7 * len(want), "too many contacts lost; this is more than the H pairs"
    else:
        assert got == want


@pytest.mark.parametrize("pdb_id", sorted(_HAS_HYDROGENS))
def test_hydrogens_do_not_change_the_contact_set(pdb_id):
    """Stripping H from the file must give the same contacts as filtering them in the reader."""
    import dataclasses

    structure = parse_structure(PDB_DIR / f"{pdb_id}.pdb")
    stripped = parse_structure(PDB_DIR / f"{pdb_id}.pdb")
    for chain in stripped.chains:
        chain.residues = [dataclasses.replace(r, atoms=tuple(a for a in r.atoms if a.element != "H"))
                          for r in chain.residues]
    assert all_atom_contacts(structure, cutoff=5.0).equals(all_atom_contacts(stripped, cutoff=5.0))


@pytest.mark.skipif(not os.getenv("RUN_BENCHMARK"), reason="set RUN_BENCHMARK=1 to run")
def test_all_structures_contacts_match_oracle():
    oracle = pl.read_csv(CONTACT_MAPS)
    mismatched = []
    for pdb_id in oracle["pdb.id"].unique().to_list():
        path = PDB_DIR / f"{pdb_id}.pdb"
        if not path.exists():
            continue
        got, want = _tcr_peptide_contact_set(pdb_id, oracle), _oracle_set(pdb_id, oracle)
        # H-bearing depositions: the legacy set includes H-mediated pairs, so ours is a subset.
        ok = (got <= want) if _has_hydrogens(path) else (got == want)
        if not ok:
            mismatched.append(pdb_id)
    assert not mismatched, f"{len(mismatched)} structures mismatched: {mismatched}"
