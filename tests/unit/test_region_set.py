"""Opt-in TCR framework regions FR1-3 (S2): default 'all' unchanged; cdr/cdr+fr restrict."""

from __future__ import annotations

from pathlib import Path

import polars as pl
import pytest

from tcren.contactmap import TCR_REGIONS, ContactMap


def test_region_sets_definition():
    assert TCR_REGIONS["cdr"] == {"CDR1", "CDR2", "CDR3"}
    assert TCR_REGIONS["cdr+fr"] == {"CDR1", "CDR2", "CDR3", "FR1", "FR2", "FR3"}
    assert TCR_REGIONS["all"] is None  # 'all' = no filter = legacy behaviour


def _toy_cm() -> ContactMap:
    # Three TCR↔peptide contacts: CDR3, FR1, FR4. 'all' keeps all; 'cdr' keeps CDR3;
    # 'cdr+fr' keeps CDR3 + FR1 (FR4 excluded).
    contacts = pl.DataFrame(
        {
            "chain.type.from": ["TRA", "TRB", "TRB"],
            "chain.type.to": ["PEPTIDE", "PEPTIDE", "PEPTIDE"],
            "residue.aa.from": ["A", "L", "K"],
            "residue.aa.to": ["G", "G", "G"],
            "region.type.from": ["CDR3", "FR1", "FR4"],
            "region.type.to": ["PEPTIDE", "PEPTIDE", "PEPTIDE"],
            "residue.index.from": [10, 20, 30],
            "residue.index.to": [0, 1, 2],
            "region.start.from": [8, 18, 28],
            "region.start.to": [0, 0, 0],
            "pdb.id": ["toy", "toy", "toy"],
        }
    )
    return ContactMap(pdb_id="toy", contacts=contacts, peptide_length=3)


def test_toy_region_filter_counts_and_membership():
    cm = _toy_cm()
    n_all = cm.interface("tcr_peptide", tcr_regions="all").height
    n_cdr = cm.interface("tcr_peptide", tcr_regions="cdr").height
    n_cdrfr = cm.interface("tcr_peptide", tcr_regions="cdr+fr").height
    assert n_cdr <= n_cdrfr <= n_all
    assert n_cdr == 1 and n_cdrfr == 2 and n_all == 3  # CDR3 / +FR1 / +FR4
    # Default 'all' is byte-identical to no argument (legacy behaviour).
    assert cm.interface("tcr_peptide").equals(cm.interface("tcr_peptide", tcr_regions="all"))
    cdrfr = cm.interface("tcr_peptide", tcr_regions="cdr+fr")
    assert cdrfr.filter(pl.col("region.type.from").is_in(["FR1", "FR2", "FR3"])).height == 1


def test_unknown_region_set_raises():
    cm = _toy_cm()
    with pytest.raises(ValueError):
        cm.interface("tcr_peptide", tcr_regions="bogus")


# --- Real TCR-pMHC asset: requires arda annotation. ---

pytest.importorskip("arda")

from tcren.annotation import classify_chains  # noqa: E402
from tcren.structure import parse_structure  # noqa: E402

# 5jhd's TCR↔peptide interface has CDR + a single FR contact, so cdr < cdr+fr == all.
_FIXTURE = Path(__file__).resolve().parents[1] / "assets" / "pdb" / "5jhd.pdb"


def test_real_asset_region_ordering_and_fr_membership():
    s = parse_structure(_FIXTURE)
    classify_chains(s, organism="human")
    cm = ContactMap.from_structure(s)

    n_all = cm.interface("tcr_peptide", tcr_regions="all").height
    n_cdr = cm.interface("tcr_peptide", tcr_regions="cdr").height
    n_cdrfr = cm.interface("tcr_peptide", tcr_regions="cdr+fr").height
    assert n_cdr <= n_cdrfr <= n_all

    cdrfr = cm.interface("tcr_peptide", tcr_regions="cdr+fr")
    fr_rows = cdrfr.filter(pl.col("region.type.from").is_in(["FR1", "FR2", "FR3"]))
    assert fr_rows.height > 0  # cdr+fr includes framework FR1-3 rows

    # Default still byte-identical to no-filter on a real structure.
    assert cm.interface("tcr_peptide").equals(cm.interface("tcr_peptide", tcr_regions="all"))
