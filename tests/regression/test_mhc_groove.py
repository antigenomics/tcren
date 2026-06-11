"""Biological-topology checks for MHC groove partitioning.

The MHC fold dictates that the TCR docks on the groove helices (not the β-sheet floor)
while the peptide lies on the floor between the two helices. These tests assert that the
projected groove regions reproduce that topology. Requires arda + a built MHC reference.
"""

from __future__ import annotations

from pathlib import Path

import pytest

pytest.importorskip("arda")

pytestmark = pytest.mark.slow  # invokes arda / mmseqs per structure

from tcren import ContactMap, parse_structure
from tcren.annotation import classify_chains
from tcren.mhc import annotate_mhc

REPO = Path(__file__).resolve().parents[2]
PDB_DIR = REPO / "data" / "PDB_structures"
_HAVE_REF = (REPO / "database" / "mhc" / "alleles.aa.fasta").exists()
needs_ref = pytest.mark.skipif(not _HAVE_REF, reason="MHC reference not built")

_HELICES = {"HELIX_A1", "HELIX_A2", "HELIX_B1"}


def _region_counts(df, col="region.type.to"):
    return dict(df.group_by(col).len().iter_rows())


@needs_ref
@pytest.mark.parametrize("pdb_id,organism", [("5m01", "mouse"), ("4ozg", "human")])
def test_tcr_docks_on_helices_peptide_on_floor(pdb_id, organism):
    s = parse_structure(PDB_DIR / f"{pdb_id}.pdb")
    classify_chains(s, organism=organism)
    annotate_mhc(s)
    cm = ContactMap.from_structure(s)

    tcr_mhc = _region_counts(cm.interface("tcr_mhc"))
    pep_mhc = _region_counts(cm.interface("peptide_mhc"))

    tcr_total = sum(v for v in tcr_mhc.values() if v)
    tcr_helix = sum(v for k, v in tcr_mhc.items() if k in _HELICES)
    # The TCR contacts the MHC almost entirely through the groove helices.
    assert tcr_helix / tcr_total >= 0.85
    assert tcr_mhc.get("GROOVE_FLOOR", 0) <= 2

    # The peptide contacts both the floor and the helices that line the groove.
    assert pep_mhc.get("GROOVE_FLOOR", 0) > 0
    assert sum(v for k, v in pep_mhc.items() if k in _HELICES) > 0
