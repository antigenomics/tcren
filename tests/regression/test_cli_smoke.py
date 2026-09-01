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
    # Phi is decomposed over the three interfaces; the total is their Native2026-normalised sum,
    # since the three are scored with potentials that are not on a common scale.
    parts = ["Phi_tcr_pep", "Phi_tcr_mhc", "Phi_pep_mhc"]
    assert set(parts) <= set(df.columns), df.columns
    from tcren.pipeline import _phi_scale, _resolve_potentials
    pots = _resolve_potentials(None)
    want = sum(df[c][0] / _phi_scale(i, pots[i])
               for c, i in zip(parts, ("tcr_peptide", "tcr_mhc", "peptide_mhc")))
    assert df["Phi_total"][0] == pytest.approx(want, abs=1e-6)
    assert "dPhi_pep_mhc" in df.columns, "--delta must add the poly-alanine reference"


@arda
def test_recognize_reports_descriptors_and_mechanics(tmp_path):
    """The manuscript's path: one command, one table, descriptors and mechanics in it.

    The fitted composites (`p_real`, `p_bind`, `p_forced`) and the cohort-relative `q_bind` /
    `s_strain` that this command used to append were removed in 2.26.0 -- their coefficients were
    frozen against training sets that no longer exist. Scoring is `tcren recognize --features`.
    """
    pytest.importorskip("arda")
    out = tmp_path / "rec.tsv"
    run("recognize", "-s", ASSET, "-o", out, "--mechanics")
    df = pl.read_csv(out, separator="\t")
    assert df.height == 1 and df["complex.id"][0] == "1ao7"
    assert "error" not in df.columns, df.to_dicts()
    assert {"burial", "chain_balance", "n_hbond"} <= set(df.columns)      # descriptors
    assert {"K_tens", "rupture_work", "couple_total"} <= set(df.columns)  # --mechanics
    assert not {"p_real", "p_bind", "p_forced", "q_bind", "s_strain"} & set(df.columns)


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


@arda
def test_potts_score_emits_the_energy_block_that_s_score_consumes(tmp_path):
    """The column `reliability.s_score` names as its Pi block has to be the column this writes.

    It was not: 2.15.0 emitted only `energy = E(sigma)`, the opposite sign, so the three-block
    `S` was unreachable from the shipped package and every caller silently fell back to two
    blocks. The contract is a sign, so assert the sign.
    """
    pytest.importorskip("arda")
    from tcren.reliability import PI_FROZEN

    out = tmp_path / "potts.tsv"
    run("potts", "score", "-s", ASSET, "-o", out)
    df = pl.read_csv(out, separator="\t")
    assert df.height == 1 and PI_FROZEN in df.columns, df.columns
    assert df[PI_FROZEN][0] == pytest.approx(-df["energy"][0])


@arda
def test_assess_writes_the_three_blocks_a_caller_decides_on(tmp_path):
    """`tcren assess` is the one command that turns a folder of models into a decision."""
    pytest.importorskip("arda")
    feats, out = tmp_path / "feats.tsv", tmp_path / "assessed.tsv"
    run("features", "-s", ASSET, "-i", "placement,interface,topology,energetics", "-o", feats)
    run("assess", "--features", feats, "-o", out)
    df = pl.read_csv(out, separator="\t")
    assert df.height == 1 and df["complex.id"][0] == "1ao7"
    assert {"S", "p_binder", "rank", "percentile"} <= set(df.columns), df.columns
    assert 0.0 < df["p_binder"][0] < 1.0


@arda
def test_diagnose_corrects_the_generator_confidence_and_shows_its_parts(tmp_path):
    """`tcren diagnose` answers "it says it is confident -- what should I believe instead"."""
    import numpy as np

    pytest.importorskip("arda")
    feats, out = tmp_path / "feats.tsv", tmp_path / "diagnosed.tsv"
    run("features", "-s", ASSET, "-i", "placement,interface,topology,energetics,potts", "-o", feats)
    # the assets are crystals and carry no generator output, so the confidence is supplied here
    pl.read_csv(feats, separator="\t").with_columns(
        pl.lit(0.88).alias("iptm")).write_csv(feats, separator="\t")
    run("diagnose", "--features", feats, "--confidence", "iptm", "-o", out)
    df = pl.read_csv(out, separator="\t")
    assert df.height == 1 and df["complex.id"][0] == "1ao7"
    assert {"p_confidence", "delta_logit", "p_corrected", "S"} <= set(df.columns), df.columns
    assert 0.0 < df["p_corrected"][0] < 1.0
    # the decomposition is the point: the two probabilities differ by exactly delta_logit
    lo = lambda p: float(np.log(p / (1 - p)))  # noqa: E731
    assert lo(df["p_corrected"][0]) - lo(df["p_confidence"][0]) == pytest.approx(
        df["delta_logit"][0], abs=1e-9)


def test_diagnose_lists_its_frozen_corrections_without_a_features_table():
    """The listing flag must not require the table it is helping the caller prepare."""
    assert "tcrvdb|ipTM" in run("diagnose", "--list-references").output


def test_diagnose_refuses_a_table_with_no_confidence_column(tmp_path):
    """A missing confidence is a clear error naming the column, never a silent skip."""
    from tcren.cli import app

    feats = tmp_path / "feats.tsv"
    pl.DataFrame({"complex.id": ["x"], **{c: [0.0] for c in
                  ("burial", "n_pep_contacted", "chain_balance", "n_hbond",
                   "D2_pep24", "fp_b0_frac_r7", "H_cell", "L_canon", "ab_imb")}}
                 ).write_csv(feats, separator="\t")
    r = CliRunner().invoke(app, ["diagnose", "--features", str(feats), "--confidence", "iptm"])
    assert r.exit_code != 0
    assert "confidence" in r.output


@arda
def test_potts_map_closes_the_pairs_onto_the_two_grids_a_caller_reads(tmp_path):
    """`--by loop` is the contact-frequency map, `--by position` the residue-importance profile.

    Both are aggregations of `--by pair`, so the invariant worth asserting at the CLI boundary is
    that they agree with it: the loop grid cannot have more rows than the pair table, the position
    grid cannot have more than the loop grid, and every frequency is a probability.
    """
    pytest.importorskip("arda")
    got = {}
    for by in ("pair", "loop", "position"):
        out = tmp_path / f"map_{by}.tsv"
        run("potts", "map", "-s", ASSET, "--by", by, "-o", out, "-w", "1")
        got[by] = pl.read_csv(out, separator="\t")

    assert got["pair"].height >= got["loop"].height >= got["position"].height > 0
    for by in ("loop", "position"):
        d = got[by]
        assert {"pdb.id", "pos.par", "aa.par", "p_any", "p_expected", "n_pairs", "n_observed",
                "observed"} <= set(d.columns), d.columns
        assert d["p_any"].min() >= 0.0 and d["p_any"].max() <= 1.0
        # p_any is P(at least one) and p_expected the expected count, so the count dominates
        assert (d["p_expected"] >= d["p_any"] - 1e-9).all()
        assert set(d["observed"].unique()) <= {0, 1}
    assert "region.rec" in got["loop"].columns and "region.rec" not in got["position"].columns
    # the peptide of 1ao7 is a 9-mer, and every position that has an available pair appears once
    assert got["position"].height == got["position"].select("pos.par").n_unique()


def test_potts_map_rejects_an_unknown_grouping(tmp_path):
    """A bad --by must fail before any structure is parsed, and name the valid choices."""
    from tcren.cli import app

    res = CliRunner().invoke(app, ["potts", "map", "-s", str(ASSET), "--by", "residue",
                                   "-o", str(tmp_path / "x.tsv")])
    assert res.exit_code != 0
    assert "loop|position|pair" in res.output


def test_potts_scan_emits_an_equimolar_referenced_energy_per_substitution(tmp_path):
    """`scan` must cover every position `map --by position` does, twenty residues each.

    The invariant worth asserting at the CLI boundary is the reference: `dF` is taken against the
    mean over the twenty residues at a position, so it sums to zero down each position. If that
    ever stops holding the table is referenced against something else and every downstream number
    carries the offset.
    """
    pytest.importorskip("arda")
    scan = tmp_path / "scan.tsv"
    run("potts", "scan", "-s", ASSET, "-o", scan)
    d = pl.read_csv(scan, separator="\t")
    assert {"pdb.id", "pos.par", "aa.par", "log_z0", "dF", "n_pairs",
            "is_observed"} <= set(d.columns), d.columns

    pos = tmp_path / "pos.tsv"
    run("potts", "map", "-s", ASSET, "--by", "position", "-o", pos, "-w", "1")
    covered = set(pl.read_csv(pos, separator="\t")["pos.par"])
    assert set(d["pos.par"]) == covered
    assert d.height == 20 * len(covered)

    per_pos = d.group_by("pdb.id", "pos.par").agg(pl.col("dF").sum().alias("s"),
                                                  pl.col("is_observed").sum().alias("obs"))
    assert per_pos["s"].abs().max() < 1e-8
    assert set(per_pos["obs"].unique()) == {1}
