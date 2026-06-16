"""3D pocket/CDR view smoke test (slow — needs arda, MHC reference, native DB, py3Dmol)."""

from __future__ import annotations

from pathlib import Path

import pytest

pytest.importorskip("arda")
pytest.importorskip("py3Dmol")

pytestmark = pytest.mark.slow

from tcren.annotation import classify_chains
from tcren.mhc import annotate_mhc
from tcren.paths import native_dir
from tcren.structure import parse_structure
from tcren.viz import view_pocket_cdr

REPO = Path(__file__).resolve().parents[2]
_HAVE_REF = (REPO / "database" / "mhc" / "alleles.aa.fasta").exists()
_HAVE_1AO7 = native_dir().is_dir() and any(native_dir().glob("1ao7.*"))


@pytest.mark.skipif(
    not (_HAVE_1AO7 and _HAVE_REF), reason="Native2026 references / MHC reference not built"
)
def test_view_pocket_cdr_builds():
    s = parse_structure(next(native_dir().glob("1ao7.*")), pdb_id="1ao7")
    classify_chains(s, organism="human")
    annotate_mhc(s)
    view = view_pocket_cdr(s)
    html = view._make_html()
    assert len(html) > 1000
    assert "viewer" in html.lower() or "3dmol" in html.lower()
