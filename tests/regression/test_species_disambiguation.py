"""Species autodetection by alignment score (human vs mouse).

arda always aligns to a fixed-organism germline reference, so tcren annotates each
structure against both human and mouse and keeps the higher total mmseqs alignment score
(``v_score``) over the receptor chains — the wrong species scores measurably lower. This
recovers the species even when a wrong reference would still return a plausible pairing
(e.g. BM3.3 mis-typed as αα under human; mouse TCRs that human refs call a clean αβ).
Adjudicated against TCR3D: BM3.3 (1fo0/1nam/2ol3) and the 2C TCR (1g6r/2ckb) are mouse αβ;
5xot is human βδ (Beta+Delta).
"""

from __future__ import annotations

from pathlib import Path

import pytest

from tcren import parse_structure
from tcren.annotation import classify_chains

arda = pytest.importorskip("arda")

pytestmark = pytest.mark.slow  # annotates each structure under both organisms

REPO = Path(__file__).resolve().parents[2]
PDB_DIR = REPO / "data" / "PDB_structures"


def _loci(structure) -> list[str]:
    return sorted(
        c.chain_type for c in structure.chains if c.chain_type in ("TRA", "TRB", "TRD", "TRG")
    )


# BM3.3 mis-types as αα under human refs; the 2C TCR (1g6r/2ckb) types as a *plausible* αβ
# under human refs — both must be recovered as mouse by the higher mouse alignment score.
@pytest.mark.parametrize("pdb_id", ["1fo0", "1nam", "2ol3", "1g6r", "2ckb", "5m01", "1mwa"])
def test_mouse_tcr_recovered_from_human_default(pdb_id):
    s = parse_structure(PDB_DIR / f"{pdb_id}.pdb")
    classify_chains(s, organism="human", autodetect_species=True)
    assert s.complex_species == "Mouse"
    assert _loci(s) == ["TRA", "TRB"]


@pytest.mark.parametrize("pdb_id", ["5xot", "6bj3", "6bj8"])
def test_exotic_beta_delta_stays_human(pdb_id):
    # Genuine βδ TCR: human scores higher than mouse, so the human βδ call is kept.
    s = parse_structure(PDB_DIR / f"{pdb_id}.pdb")
    classify_chains(s, organism="human", autodetect_species=True)
    assert s.complex_species == "Human"
    assert _loci(s) == ["TRB", "TRD"]


def test_score_overrides_wrong_explicit_organism():
    # A human TCR passed organism="mouse" is still corrected to Human by the score.
    s = parse_structure(PDB_DIR / "1ao7.pdb")
    classify_chains(s, organism="mouse", autodetect_species=True)
    assert s.complex_species == "Human"
    assert _loci(s) == ["TRA", "TRB"]


def test_autodetect_off_keeps_requested_organism():
    # Opting out must preserve the (wrong) requested-organism call verbatim.
    s = parse_structure(PDB_DIR / "1fo0.pdb")
    classify_chains(s, organism="human", autodetect_species=False)
    assert s.complex_species == "Human"
    assert _loci(s) == ["TRA", "TRA"]
