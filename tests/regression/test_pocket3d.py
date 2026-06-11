"""3D pocket/CDR view smoke test (slow — needs arda, MHC reference, native DB, py3Dmol)."""

from __future__ import annotations

from pathlib import Path

import pytest

pytest.importorskip("arda")
pytest.importorskip("py3Dmol")

pytestmark = pytest.mark.slow

from tcren.annotation import classify_chains
from tcren.mhc import annotate_mhc
from tcren.native import NativeDatabase
from tcren.structure import parse_structure
from tcren.viz import view_pocket_cdr

REPO = Path(__file__).resolve().parents[2]
_DB = NativeDatabase()
_HAVE_REF = (REPO / "database" / "mhc" / "alleles.aa.fasta").exists()


@pytest.mark.skipif(
    not (_DB.is_present() and _HAVE_REF), reason="native DB / MHC reference not built"
)
def test_view_pocket_cdr_builds():
    db = NativeDatabase()
    s = parse_structure(db.cif_for("1ao7"), pdb_id="1ao7")
    classify_chains(s, organism="human")
    annotate_mhc(s)
    view = view_pocket_cdr(s, db=db)
    html = view._make_html()
    assert len(html) > 1000
    assert "viewer" in html.lower() or "3dmol" in html.lower()
