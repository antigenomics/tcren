"""CLI regression for ``tcren energy`` and ``tcren mechanics`` (slow — arda/mmseqs).

Guards the wiring of the two manuscript-facing feature extractors onto the CLI: the DOPE
interface-energy scorer (``interface_energy`` / the ``_relax`` kernel) and the interface
mechanics / koff-proxy battery (``stiffness_tensor`` + ``rupture`` + ``coupling_residues``).
The underlying functions are unit-tested in ``test_relax.py`` / ``test_mechanics.py``; this
pins the CLI column contract on a real structure.
"""

from __future__ import annotations

from pathlib import Path

import polars as pl
import pytest

pytestmark = pytest.mark.slow  # invokes arda / mmseqs to classify chains

ASSET = Path(__file__).resolve().parents[1] / "assets" / "pdb" / "1ao7.pdb"


def test_energy_cli_reports_native_relax_and_gap(tmp_path):
    from typer.testing import CliRunner

    from tcren.cli import app

    out = tmp_path / "energy.csv"
    result = CliRunner().invoke(app, ["energy", "-s", str(ASSET), "-o", str(out), "--relax"])
    assert result.exit_code == 0, result.output

    df = pl.read_csv(out)
    assert df.height == 1 and df["pdb.id"][0] == "1ao7"
    # Native interface is favourable (DOPE energy negative); relaxation only lowers it,
    # so the gap = e_native - e_relax is >= 0.
    assert df["e_native"][0] < 0
    assert df["e_relax"][0] <= df["e_native"][0] + 1e-6
    assert df["gap"][0] >= -1e-6


def test_mechanics_cli_reports_stiffness_rupture_coupling(tmp_path):
    from typer.testing import CliRunner

    from tcren.cli import app

    out = tmp_path / "mechanics.csv"
    result = CliRunner().invoke(app, ["mechanics", "-s", str(ASSET), "-o", str(out)])
    assert result.exit_code == 0, result.output

    df = pl.read_csv(out)
    assert df.height == 1 and df["pdb.id"][0] == "1ao7"
    # A real TCR-pMHC interface has many springs; the tensile stiffness and rupture force
    # are finite and positive; K_tens is part of S_tot.
    assert df["n_spring"][0] >= 3
    assert df["K_tens"][0] > 0 and df["rupture_force"][0] > 0
    assert df["K_tens"][0] <= df["S_tot"][0] + 1e-6
    for col in ("couple_pep", "couple_total"):
        assert col in df.columns
