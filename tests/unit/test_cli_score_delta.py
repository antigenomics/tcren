"""`tcren score --delta`: the two poly-alanine references, and the per-template candidate manifest.

The referenced energy is what makes two ranking units commensurate, and it comes in two forms that
are NOT interchangeable -- the reason this file exists is that the benchmark computed the
structural one by hand because no CLI emitted it.
"""

from __future__ import annotations

from pathlib import Path

import polars as pl
import pytest

pytest.importorskip("arda")

from typer.testing import CliRunner  # noqa: E402

from tcren.cli import _read_candidates, app  # noqa: E402

_FIXTURE = Path(__file__).resolve().parents[1] / "assets" / "pdb" / "1ao7.pdb"
_PEPTIDE = "LLFGYPVYV"


def _run(tmp_path, *args, candidates=_PEPTIDE):
    cand = tmp_path / "cand.txt"
    cand.write_text(candidates if candidates.endswith("\n") else candidates + "\n")
    out = tmp_path / "out.csv"
    res = CliRunner().invoke(app, ["score", "-s", str(_FIXTURE), "-c", str(cand),
                                   "-o", str(out), *args])
    assert res.exit_code == 0, res.output
    return pl.read_csv(out)


def test_delta_off_by_default(tmp_path):
    assert "dphi" not in _run(tmp_path).columns


def test_structural_and_virtual_references_differ(tmp_path):
    """Both are Phi(candidate) - Phi(polyAla); they disagree on what polyAla is allowed to touch.

    Virtually, the baseline is charged for every contact the real side chains make. Structurally,
    the peptide is truncated and the map recomputed, so the baseline is what backbone plus Cbeta
    alone reach -- fewer contacts, and a more negative reference.
    """
    from tcren.annotation import classify_chains
    from tcren.contactmap import ContactMap
    from tcren.ddg import reference_delta
    from tcren.potential import tcren2
    from tcren.structure.io import import_structure

    st = _run(tmp_path, "--delta", "--reference", "structural")
    vi = _run(tmp_path, "--delta", "--reference", "virtual")
    assert st["score"][0] == vi["score"][0]          # the candidate's own energy is untouched
    assert st["dphi"][0] != vi["dphi"][0]

    # Which way they differ is a property of the structure, not an invariant: on 1ao7 the rebuilt
    # reference sits nearer zero and on 3qib further from it. What IS the contract is that each
    # mode reproduces `reference_delta` with and without a structure to rebuild from.
    struct = import_structure(_FIXTURE)
    classify_chains(struct)   # in place; the CLI does this in one batched pass
    cm = ContactMap.from_structure(struct, cutoff=5.0)
    pot = tcren2()
    assert st["dphi"][0] == pytest.approx(
        reference_delta(cm, _PEPTIDE, pot, structure=struct), rel=1e-9)
    assert vi["dphi"][0] == pytest.approx(
        reference_delta(cm, _PEPTIDE, pot), rel=1e-9)


def test_delta_on_one_template_is_score_minus_a_constant(tmp_path):
    """Phi(reference) does not depend on the candidate, so it is computed once per template."""
    two = _run(tmp_path, "--delta", candidates=f"{_PEPTIDE}\nLLFGYPVYA\n")
    assert two.height == 2
    offsets = (two["score"] - two["dphi"]).to_list()
    assert offsets[0] == pytest.approx(offsets[1])


def test_bad_reference_mode_is_rejected(tmp_path):
    cand = tmp_path / "c.txt"
    cand.write_text(_PEPTIDE + "\n")
    res = CliRunner().invoke(app, ["score", "-s", str(_FIXTURE), "-c", str(cand),
                                   "-o", str(tmp_path / "o.csv"), "--reference", "nonsense"])
    assert res.exit_code != 0


def test_manifest_selects_a_pool_per_template(tmp_path):
    """A benchmark of many templates, each with its own pool, is ONE call and one typing pass."""
    man = tmp_path / "pool.tsv"
    man.write_text("complex.id\tpeptide\n1ao7\tLLFGYPVYV\n1ao7\tLLFGYPVYA\n6zzz\tLLFGYPVYV\n")
    out = tmp_path / "m.csv"
    res = CliRunner().invoke(app, ["score", "-s", str(_FIXTURE), "-c", str(man), "-o", str(out)])
    assert res.exit_code == 0, res.output
    table = pl.read_csv(out)
    # Only 1ao7's own two rows: the pool named for a template that is not in -s scores nothing.
    assert table.height == 2
    assert set(table["complex.id"]) == {"1ao7"}


def test_plain_candidate_file_still_reads_as_a_flat_list(tmp_path):
    plain = tmp_path / "p.txt"
    plain.write_text("peptide\nLLFGYPVYV\nLLFGYPVYA\n")
    assert _read_candidates(plain) == ["LLFGYPVYV", "LLFGYPVYA"]


def test_manifest_is_detected_by_its_header(tmp_path):
    man = tmp_path / "m.tsv"
    man.write_text("complex.id\tpeptide\n1AO7\tLLFGYPVYV\n")
    assert _read_candidates(man) == {"1ao7": ["LLFGYPVYV"]}


def test_bundled_potts_model_resolves_by_name():
    """`potts contacts -m` took a path only, so a shipped model had no name on the command line."""
    from tcren.cli import _load_potts

    default, sc = _load_potts(None), _load_potts("potts_tcr_peptide_sc")
    assert sc.coupling_pinned and sc.beta_matrix != default.beta_matrix
    with pytest.raises(Exception):
        _load_potts("no_such_model")


def test_batched_ids_are_resolved_for_the_set_not_per_file(tmp_path):
    """Seven files sharing a prefix are seven structures, not one.

    `structure_id_from_path` takes the part before the first underscore, which is right for
    `4x5w_renumbered.cif` and wrong for a cohort whose names carry metadata there --
    `302TIL_A0206_9_TCR-pMHC` and its six siblings all came back as `302TIL`. `iter_structures`
    has always decided per SET (`resolve_structure_ids`); the batched annotation path did not, so
    `tcren scoring` and `tcren score` disagreed about what a structure is called.
    """
    import shutil

    from tcren.annotation.batch import iter_typed

    for tag in ("A0206_9", "B3501_10", "C0303_9"):
        shutil.copy(_FIXTURE, tmp_path / f"302TIL_{tag}_TCR-pMHC.pdb")
    got = [s.pdb_id for s in iter_typed(tmp_path, "human")]
    assert len(set(got)) == 3, got
    assert all(g.startswith("302TIL_") for g in got), got


def test_potential_table_writes_the_shipped_matrix(tmp_path):
    """`potential-table` resolves the same names `score -p` does, and computes nothing."""
    import polars as pl
    from typer.testing import CliRunner

    from tcren.cli import app
    from tcren.potential import tcren2_em

    out = tmp_path / "em.csv"
    r = CliRunner().invoke(app, ["potential-table", "-p", "tcren2_em", "-o", str(out)])
    assert r.exit_code == 0, r.output
    got = pl.read_csv(out)
    want = tcren2_em().matrix
    assert got.height == want.height == 380
    assert set(got["residue.aa.from"]) == set(want["residue.aa.from"])
    # cysteine is a peptide-side residue only, so 20 of the 400 cells never exist
    assert "C" not in set(got["residue.aa.from"])
    assert "C" in set(got["residue.aa.to"])


def test_potential_table_rejects_an_unknown_name(tmp_path):
    from typer.testing import CliRunner

    from tcren.cli import app

    r = CliRunner().invoke(app, ["potential-table", "-p", "nope", "-o", str(tmp_path / "x.csv")])
    assert r.exit_code != 0
