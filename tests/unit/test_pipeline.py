"""End-to-end pipeline: annotate → superimpose → contacts → per-interface scores."""

from __future__ import annotations

from pathlib import Path

import pytest

pytest.importorskip("arda")

from tcren.pipeline import run, score_row  # noqa: E402

_FIXTURE = Path(__file__).resolve().parents[1] / "assets" / "pdb" / "1ao7.pdb"


def test_pipeline_no_superimpose_scores_three_interfaces():
    res = run(_FIXTURE, superimpose=False)
    assert set(res.scores) == {"tcr_peptide", "tcr_mhc", "peptide_mhc", "total"}
    assert res.oriented is None and res.rmsd is None
    assert res.markup.height > 0 and res.contacts.height > 0
    # Phi_TCRpMHC = c_TP Phi_TP + c_TM Phi_TM + c_PM Phi_PM. The three interfaces are scored with
    # different potentials, whose matrices are not on a common scale, so the total is the
    # Native2026-normalised sum rather than the raw one -- an unweighted sum would be 2.6x more
    # sensitive to a presentation contact than to a recognition one.
    from tcren.pipeline import _INTERFACE_POTENTIAL, _phi_scale, _resolve_potentials

    pots = _resolve_potentials(None)
    want = sum(res.scores[i] / _phi_scale(i, pots[i]) for i in _INTERFACE_POTENTIAL)
    assert res.scores["total"] == pytest.approx(want)
    assert res.scores["total"] != pytest.approx(
        res.scores["tcr_peptide"] + res.scores["tcr_mhc"] + res.scores["peptide_mhc"])


def test_pipeline_superimpose_sets_canonical_frame(tmp_path):
    from tcren.docking import run_folder

    canonical = tmp_path / "canonical"
    run_folder(_FIXTURE, canonical)
    res = run(_FIXTURE, superimpose=True, db_dir=canonical)
    assert res.oriented is not None and res.rmsd is not None
    assert {c.chain_id for c in res.oriented.chains} <= set("ABCDE")
    row = score_row(res)
    assert row["pdb.id"] == "1ao7" and row["mhc.class"] == "MHCI"
    assert row["Phi_tcr_pep"] == res.scores["tcr_peptide"]


def test_pipeline_reference_aa_adds_delta_f():
    """--delta path: ΔΦ per interface, ΔΦ_tcr_mhc ≡ 0, ΔΦ_total the weighted TP + PM sum."""
    from tcren.pipeline import _phi_scale, _resolve_potentials

    res = run(_FIXTURE, superimpose=False, reference_aa="A")
    assert res.scores["delta_tcr_mhc"] == 0.0          # the peptide is not in that interface
    pots = _resolve_potentials(None)
    assert res.scores["delta_total"] == pytest.approx(
        res.scores["delta_tcr_peptide"] / _phi_scale("tcr_peptide", pots["tcr_peptide"])
        + res.scores["delta_peptide_mhc"] / _phi_scale("peptide_mhc", pots["peptide_mhc"])
    )
    # ΔΦ = Φ(peptide) − Φ(poly-Ala) is a genuine difference, not a copy of Φ
    assert res.scores["delta_tcr_peptide"] != res.scores["tcr_peptide"]
    row = score_row(res)
    assert row["dPhi_total"] == res.scores["delta_total"]
    assert "dPhi_total" not in score_row(run(_FIXTURE, superimpose=False))


def test_score_row_exposes_raw_peptide_sum_without_changing_historical_totals():
    """Export the existing TP+PM quantity without borrowing the historical scale weights."""
    import polars as pl
    from tcren.pipeline import PipelineResult

    scores = {"tcr_peptide": 2.0, "tcr_mhc": 3.0, "peptide_mhc": 4.0, "total": 71.0,
              "delta_tcr_peptide": 1.25, "delta_tcr_mhc": 0.0,
              "delta_peptide_mhc": -4.5, "delta_total": 23.0}
    result = PipelineResult("toy", [], pl.DataFrame(), pl.DataFrame(), scores)
    row = score_row(result)
    assert row["dPhi_peptide_sum"] == -3.25
    assert row["Phi_total"] == 71.0 and row["dPhi_total"] == 23.0
    del scores["delta_total"]
    assert "dPhi_peptide_sum" not in score_row(result)


@pytest.mark.parametrize("contact_weight", ["residue", "atomic"])
@pytest.mark.parametrize("regions", ["all", "cdr", "cdr+fr"])
def test_structural_reference_reuses_one_map(monkeypatch, contact_weight, regions):
    from tcren.annotation import classify_chains
    from tcren.contactmap import ContactMap
    from tcren.energetics.mutation import reference_delta
    from tcren.mhc import annotate_mhc
    from tcren.pipeline import _resolve_potentials
    from tcren.structure.io import import_structure
    from tcren.structure.model import PEPTIDE_TYPE
    import tcren.energetics.mutation as mutation

    structure = import_structure(_FIXTURE)
    classify_chains(structure)
    structure.mhc_calls = annotate_mhc(structure)
    peptide = next(c.sequence() for c in structure.chains if c.chain_type == PEPTIDE_TYPE)
    cm = ContactMap.from_structure(structure, count_atoms=contact_weight == "atomic")
    potentials = _resolve_potentials(None)
    expected = {i: reference_delta(cm, peptide, potentials[i], interface=i,
                                  structure=structure, tcr_regions=regions,
                                  contact_weight=contact_weight)
                for i in ("tcr_peptide", "tcr_mhc", "peptide_mhc")}
    original = mutation._mutant_map
    rebuilds = []

    def counted(*args, **kwargs):
        rebuilds.append(1)
        return original(*args, **kwargs)

    monkeypatch.setattr(mutation, "_mutant_map", counted)
    result = run(structure, typed=True, superimpose=False, reference_aa="A",
                 reference_mode="structural", contact_weight=contact_weight,
                 tcr_regions=regions)
    assert len(rebuilds) == 1
    assert {i: result.scores[f"delta_{i}"] for i in expected} == expected


@pytest.mark.parametrize("contact_weight", ["residue", "atomic"])
@pytest.mark.parametrize("mode", ["virtual", "structural"])
def test_receptor_reference_scores_both_interfaces_once(monkeypatch, contact_weight, mode):
    import polars as pl
    from test_footprint import _full_complex

    from tcren.contactmap import ContactMap
    from tcren.energetics.scoring import _interface_energy
    from tcren.pipeline import _resolve_potentials
    from tcren.refine import substitute
    from tcren.structure.model import TCR_TYPES, RegionMarkup

    structure = _full_complex()
    structure.mhc_calls = []
    peptide_chain, = structure.by_type("PEPTIDE")
    peptide_chain.regions = [
        RegionMarkup("PEPTIDE", 0, len(peptide_chain.residues) - 1,
                     peptide_chain.sequence(), peptide_chain.residues)
    ]
    native_map = ContactMap.from_structure(structure, count_atoms=contact_weight == "atomic")
    pots = _resolve_potentials(None)
    aa = "A" if mode == "structural" else "V"
    targets = {(chain.chain_id, res.seq_index): aa
               for chain in structure.by_type(*TCR_TYPES) for res in chain.residues}
    original = substitute.substitute_residues
    expected_calls = [targets] if mode == "structural" else []
    if mode == "structural":
        reference_map = ContactMap.from_structure(
            original(structure, targets), count_atoms=contact_weight == "atomic",
        )
    calls = []

    def counted(source, changes):
        calls.append(changes)
        return original(source, changes)

    monkeypatch.setattr(substitute, "substitute_residues", counted)
    kwargs = dict(typed=True, superimpose=False, reference_aa=aa, reference_mode=mode,
                  contact_weight=contact_weight, tcr_regions="cdr")
    result = score_row(run(structure, reference_side="tcr", **kwargs))
    assert calls == expected_calls
    assert "dPhi_peptide_sum" not in result
    assert result["dPhi_tcr_ref_pep_mhc"] == 0.0
    for interface, suffix in (("tcr_peptide", "tcr_pep"), ("tcr_mhc", "tcr_mhc")):
        contacts = native_map.interface(interface, tcr_regions="cdr")
        reference = (reference_map.interface(interface, tcr_regions="cdr")
                     if mode == "structural" else
                     contacts.with_columns(pl.lit(aa).alias("residue.aa.from")))
        assert result[f"dPhi_tcr_ref_{suffix}"] == pytest.approx(
            _interface_energy(contacts, pots[interface], contact_weight)
            - _interface_energy(reference, pots[interface], contact_weight)
        )
    assert result["dPhi_tcr_sum"] == pytest.approx(
        result["dPhi_tcr_ref_tcr_pep"] + result["dPhi_tcr_ref_tcr_mhc"]
    )
    peptide = score_row(run(structure, **kwargs))
    both = score_row(run(structure, reference_side="both", **kwargs))
    assert {key: both[key] for key in peptide} == peptide
    assert both["dPhi_tcr_sum"] == result["dPhi_tcr_sum"]


def test_scoring_cli_accepts_reference_side_and_records_it(tmp_path, monkeypatch):
    from click import unstyle
    import polars as pl
    from test_footprint import _full_complex
    from typer.testing import CliRunner

    from tcren import cli

    structure = _full_complex()
    structure.mhc_calls = []
    monkeypatch.setattr(cli, "_iter_typed", lambda *a, **kw: [structure])
    monkeypatch.setattr("tcren.mhc.annotate_mhc_batch", lambda *a, **kw: None)
    out = tmp_path / "scores.csv"
    args = ["scoring", "-s", "test.pdb", "--no-superimpose", "-t", "1", "-o", str(out)]
    runner = CliRunner()
    result = runner.invoke(cli.app, args + ["--delta", "--reference-side", "tcr"])
    assert result.exit_code == 0, result.output
    assert "dPhi_tcr_sum" in pl.read_csv(out).columns
    from tcren.provenance import read
    assert read(out)["reference_side"] == "tcr"
    result = runner.invoke(cli.app, args + ["--reference-side", "tcr"])
    assert result.exit_code == 2
    assert "--delta" in unstyle(result.output)
