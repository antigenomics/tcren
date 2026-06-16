"""Unit tests for the dataset analysis helpers (fast — uses the committed CSVs)."""

from __future__ import annotations

from pathlib import Path

import polars as pl

from tcren import analysis as an
from tcren.potential import derive_tcren, mj, tcren

REPO = Path(__file__).resolve().parents[2]
CONTACTS = REPO / "legacy" / "data" / "contact_maps_PDB.csv"
SUMMARY = REPO / "legacy" / "data" / "summary_PDB_structures.csv"


def _df():
    return an.load_interface_contacts(CONTACTS, SUMMARY)


def test_enriched_columns_present():
    df = _df()
    for col in ("peptide_pos", "peptide_len", "cdr3_len", "cdr3_rel_pos", "nonred"):
        assert col in df.columns
    # peptide position never exceeds peptide length
    assert df.filter(pl.col("peptide_pos") >= pl.col("peptide_len")).height == 0


def test_region_counts_dominated_by_cdr3():
    counts = an.region_contact_counts(_df())
    by_region = dict(
        counts.group_by("region.type.from").agg(pl.col("n_contacts").sum()).iter_rows()
    )
    assert by_region["CDR3"] > by_region.get("CDR1", 0)
    assert by_region["CDR3"] > by_region.get("CDR2", 0)
    assert by_region["CDR3"] > sum(v for k, v in by_region.items() if k.startswith("FR"))


def test_position_distributions():
    df = _df()
    pep = an.position_distribution(df, "peptide")
    assert set(pep.columns) == {"length", "position", "n_contacts"}
    assert (pep["position"] < pep["length"]).all()
    cdr3 = an.position_distribution(df, "cdr3b")
    assert cdr3.height > 0
    assert (cdr3["position"] >= 0).all()


def test_compare_potentials_derived_equals_published():
    nonred = pl.read_csv(SUMMARY).filter(pl.col("nonred"))["pdb.id"].to_list()
    derived = derive_tcren(pl.read_csv(CONTACTS), include=nonred)
    cmp = an.compare_potentials(derived, tcren())
    assert cmp["diff"].abs().max() < 1e-9  # reproduces the published potential


def test_potential_matrix_shape():
    m, froms, tos = an.potential_matrix(tcren())
    assert m.shape == (len(froms), len(tos))
    assert "C" not in froms  # classic TCRen drops Cys from the 'from' axis
    m_mj, _, _ = an.potential_matrix(mj())
    assert m_mj.shape[0] == m_mj.shape[1]
