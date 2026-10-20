"""Adding an occluder cannot create area; separated partners bury no surface."""
import numpy as np
import pytest

from tcren.contacts.exposure import context_exposure
from tcren.structure.model import Atom, Chain, Residue, Structure


def test_exposure_preserves_context_direction_and_area():
    def chain(cid, kind, x):
        return Chain(cid, [Residue(0, 9, "A", "A", "ALA",
            (Atom("CA", "C", np.array([x, 0., 0.])),))], kind,
            "MHCI" if kind == "MHCa" else kind)
    s = Structure("toy", [chain("D", "TRA", 0.), chain("E", "TRB", -20.),
                          chain("P", "PEPTIDE", 3.), chain("M", "MHCa", 6.)])
    a = context_exposure(s)
    assert a.height == 4 and a["insertion_code"].to_list() == ["A"]*4
    p = a.filter(a["chain.id"] == "P").row(0, named=True)
    assert p["role"] == "target" and p["mhc_shielding_A2"] > 0
    assert p["peptide_alone_sasa_A2"] == pytest.approx(p["isolated_sasa_A2"] + p["mhc_shielding_A2"])
    assert np.all(a["partner_burial_A2"].to_numpy() >= 0)
    for c in s.by_type("TRA", "TRB"):
        c.residues[0].atoms[0].coord[:] += [0., 100., 0.]
    b = context_exposure(s)
    assert b["partner_burial_A2"].sum() == 0
    assert np.array_equal(a["isolated_sasa_A2"], b["isolated_sasa_A2"])
    from tcren.exports import export_structure
    assert export_structure(s, "contact-exposure", None, None).equals(b)
    s.chain("P").residues[0].atoms[0].coord[0] = np.nan
    with pytest.raises(ValueError, match="finite"):
        context_exposure(s)


def test_exposure_cli_needs_no_potential(monkeypatch, tmp_path):
    import polars as pl
    from typer.testing import CliRunner
    from tcren import exports
    from tcren.cli import app

    def export(sources, output, **kwargs):
        assert kwargs["quantity"] == "contact-exposure"
        assert kwargs["peptide_potential"] is None and kwargs["mhc_potential"] is None
        return pl.DataFrame({"status": ["ok"]})

    monkeypatch.setattr(exports, "export_batch", export)
    result = CliRunner().invoke(app, ["export", "-s", "raw.pdb", "--quantity",
                                     "contact-exposure", "-o", str(tmp_path/"atoms.tsv")])
    assert result.exit_code == 0, result.output
