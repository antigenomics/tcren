"""Smoke tests for the user-facing commands, one per command.

``cli.py`` was the least-covered module in the package (17% — 375 of 452 statements) while being
the only surface most users touch. The library underneath is well covered, so what was untested was
specifically the *wiring*: an option renamed, an output column dropped, a command that raises before
it writes anything. Those are release regressions, and they are cheap to catch.

Each test asserts the two things a caller depends on and nothing more — the command exits 0, and the
table it writes carries its key column plus whatever that command exists to produce. Numerical
correctness belongs to the unit tests of the functions underneath; duplicating it here would only
mean two places to update when a number legitimately moves.

``rank``, ``energy`` and ``mechanics`` already have their own regression files and are not repeated.
"""

from __future__ import annotations

from pathlib import Path

import polars as pl
import pytest
from typer.testing import CliRunner

ASSETS = Path(__file__).resolve().parents[1] / "assets" / "pdb"
ASSET = ASSETS / "1ao7.pdb"
NATIVE_PEPTIDE = "LLFGYPVYV"  # 1ao7 — HTLV-1 Tax on HLA-A*02:01


def run(*args) -> object:
    """Invoke the CLI and fail with its own output, which is where the useful message is."""
    from tcren.cli import app

    result = CliRunner().invoke(app, [str(a) for a in args])
    assert result.exit_code == 0, f"{args!r} exited {result.exit_code}\n{result.output}"
    return result


def test_info_reports_version_and_potential():
    """`tcren info` must not need arda, a network or a reference DB — it is the first thing a user
    runs when something is wrong, so it has to work in the environment where something is wrong."""
    import tcren

    out = run("info").output
    assert tcren.__version__ in out
    assert "potential" in out.lower()


# Everything below annotates chains, so it needs arda + mmseqs.
arda = pytest.mark.slow


@arda
def test_annotate_emits_region_markup(tmp_path):
    pytest.importorskip("arda")
    out = tmp_path / "ann.csv"
    run("annotate", "-s", ASSET, "-o", out)
    df = pl.read_csv(out)
    assert df.height > 0
    assert {"chain.id", "chain.type"} <= set(df.columns), df.columns


@arda
def test_contacts_emits_an_annotated_interface_table(tmp_path):
    pytest.importorskip("arda")
    out = tmp_path / "contacts.csv"
    run("contacts", "-s", ASSET, "-o", out)
    df = pl.read_csv(out)
    assert df.height > 0, "a real TCR:pMHC crystal has interface contacts"
    assert {"chain.id.from", "chain.id.to"} <= set(df.columns), df.columns


@arda
def test_scoring_reports_three_interfaces_and_delta(tmp_path):
    pytest.importorskip("arda")
    out = tmp_path / "scores.csv"
    run("scoring", "-s", ASSET, "-o", out, "--delta", "--no-superimpose")
    df = pl.read_csv(out)
    assert df.height == 1 and df["pdb.id"][0] == "1ao7"
    # Phi is decomposed over the three interfaces and they sum to the total.
    parts = ["F_tcr_pep", "F_tcr_mhc", "F_pep_mhc"]
    assert set(parts) <= set(df.columns), df.columns
    assert df["F_total"][0] == pytest.approx(sum(df[c][0] for c in parts), abs=1e-6)
    assert "dF_pep_mhc" in df.columns, "--delta must add the poly-alanine reference"


@arda
def test_recognize_reports_descriptors_scores_and_mechanics(tmp_path):
    """The manuscript's path: one command, one table, all three column families in it."""
    pytest.importorskip("arda")
    out = tmp_path / "rec.tsv"
    run("recognize", "-s", ASSET, "-o", out, "--scores", "--mechanics")
    df = pl.read_csv(out, separator="\t")
    assert df.height == 1 and df["complex.id"][0] == "1ao7"
    assert "error" not in df.columns, df.to_dicts()
    assert {"burial", "chain_balance", "n_hbond"} <= set(df.columns)      # descriptors
    assert {"p_real", "q_bind", "s_strain"} <= set(df.columns)            # scores
    assert {"K_tens", "rupture_work", "couple_total"} <= set(df.columns)  # --mechanics
    assert 0.0 < df["p_real"][0] < 1.0


@arda
def test_ddg_alanine_scan_covers_every_peptide_position(tmp_path):
    pytest.importorskip("arda")
    out = tmp_path / "ddg.csv"
    run("ddg", "-s", ASSET, "-o", out, "--native", NATIVE_PEPTIDE, "--alanine-scan")
    df = pl.read_csv(out)
    assert {"pos", "wt_aa", "ddG"} <= set(df.columns), df.columns
    assert df.height == len(NATIVE_PEPTIDE), "one row per peptide position"
    assert df["wt_aa"].to_list() == list(NATIVE_PEPTIDE), "wt_aa must track the native sequence"
    assert df["ddG"].is_finite().all()


@arda
def test_ddg_requires_a_mode(tmp_path):
    """Neither mode given is a usage error, not a silent empty table."""
    pytest.importorskip("arda")
    from tcren.cli import app

    r = CliRunner().invoke(app, ["ddg", "-s", str(ASSET), "-o", str(tmp_path / "d.csv"),
                                 "--native", NATIVE_PEPTIDE])
    assert r.exit_code != 0
    assert "alanine-scan" in r.output


@arda
def test_score_ranks_candidate_epitopes(tmp_path):
    pytest.importorskip("arda")
    cands = tmp_path / "cands.txt"
    cands.write_text(f"{NATIVE_PEPTIDE}\nSLYNTVATL\nGILGFVFTL\n")
    out = tmp_path / "score.csv"
    run("score", "-s", ASSET, "-c", cands, "-o", out)
    df = pl.read_csv(out)
    assert df.height >= 3, "one row per candidate at least"


@arda
def test_cpl_predicts_a_full_response_matrix(tmp_path):
    pytest.importorskip("arda")
    out = tmp_path / "cpl.csv"
    run("cpl", "-s", ASSET, "-o", out)
    df = pl.read_csv(out)
    # A CPL matrix is 20 amino acids at each peptide position.
    assert df.height == 20 * len(NATIVE_PEPTIDE) or df.width >= 20, (df.shape, df.columns)


def test_surface_writes_the_featureless_scalars_and_the_extra_outputs(tmp_path):
    pytest.importorskip("arda")
    out, svg, dist = tmp_path / "surface.csv", tmp_path / "svg", tmp_path / "d.csv"
    run("surface", "-s", ASSET, "-o", out, "--svg", svg, "--compare", dist)
    df = pl.read_csv(out)
    assert df.height == 1
    assert {"relief", "peak_to_valley", "frac_above_ridge", "phobic_centre"} <= set(df.columns)
    assert df["peptide"][0] == NATIVE_PEPTIDE
    assert list(svg.glob("*.svg")), "--svg wrote no map"
    assert dist.exists()


def test_a_missing_structure_is_one_line_not_a_traceback(monkeypatch, capsys, tmp_path):
    """``-s`` takes files, globs, manifests and archives, so Typer cannot check it — the CLI must."""
    import sys

    from tcren.cli import main

    monkeypatch.setattr(sys, "argv",
                        ["tcren", "contacts", "-s", "/no/such.pdb", "-o", str(tmp_path / "c.csv")])
    with pytest.raises(SystemExit) as exc:
        main()
    assert exc.value.code == 1
    err = capsys.readouterr().err.strip()
    assert err.startswith("Error:") and "\n" not in err, err
