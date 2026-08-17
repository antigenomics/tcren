"""Unit tests for peptide-position annotation, weighting, and the per-position energy profile."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import polars as pl
import pytest

from tcren.scoring import (
    POSITION_SCHEMES,
    central_strain,
    peptide_positions,
    position_profile,
    position_weights,
    score_peptides,
)

PDB_DIR = Path(__file__).resolve().parents[1] / "assets" / "pdb"


@pytest.fixture(scope="module")
def annotated():
    """1ao7 (class I) and 4ozg (class II), chain-typed and MHC-annotated."""
    pytest.importorskip("arda")
    from tcren.annotation import classify_chains
    from tcren.contactmap import ContactMap
    from tcren.mhc import annotate_mhc
    from tcren.structure import parse_structure

    out = {}
    for pid in ("1ao7", "4ozg"):
        s = parse_structure(PDB_DIR / f"{pid}.pdb")
        classify_chains(s, organism="human", autodetect_species=True)
        annotate_mhc(s)
        s.pdb_id = pid
        out[pid] = (s, ContactMap.from_structure(s))
    return out


# --- annotation -----------------------------------------------------------------------------------
@pytest.mark.slow
def test_positions_are_one_based_and_within_the_peptide(annotated):
    s, cm = annotated["1ao7"]
    ann = peptide_positions(cm, s)
    pos = ann["peptide.pos"].to_numpy()
    assert pos.min() >= 1
    assert pos.max() <= cm.peptide_length
    assert set(ann["peptide.role"].unique()) <= {"anchor", "tcr_facing"}


@pytest.mark.slow
def test_class_i_anchors_are_p2_and_the_c_terminus(annotated):
    s, cm = annotated["1ao7"]
    ann = peptide_positions(cm, s)
    anchors = set(ann.filter(pl.col("peptide.role") == "anchor")["peptide.pos"].to_list())
    assert anchors <= {2, cm.peptide_length}
    assert 2 in anchors


@pytest.mark.slow
def test_class_ii_register_uses_the_real_sequence_not_the_contacts(annotated):
    """4ozg is HLA-DQ2.5 + gliadin APQPELPYPQPG; its core starts at P2, so P1/P4/P6/P9 of the core
    are peptide positions 2/5/7/10. Reassembling the peptide from contacts alone gets this wrong,
    because the class-II heuristic slides a 9-mer window over a sequence full of gaps."""
    from tcren.refine.anchors import native_peptide, predict_anchors

    s, _cm = annotated["4ozg"]
    pep = native_peptide(s)
    assert pep == "APQPELPYPQPG"
    d = predict_anchors(pep, s)
    assert d.mhc_class == "MHCII"
    assert [a + 1 for a in d.anchors] == [2, 5, 7, 10]


# --- weights --------------------------------------------------------------------------------------
@pytest.mark.slow
def test_uniform_weights_leave_the_score_untouched(annotated):
    """The default must be bit-identical, or every existing number moves."""
    from tcren.potential import tcren

    s, cm = annotated["1ao7"]
    ann = peptide_positions(cm, s)
    base = score_peptides(cm, ["LLFGYPVYV"], tcren())["score"][0]
    weighted = score_peptides(cm, ["LLFGYPVYV"], tcren(),
                              weights=position_weights(ann, "uniform"))["score"][0]
    assert weighted == base


@pytest.mark.slow
@pytest.mark.parametrize("scheme", POSITION_SCHEMES)
def test_weights_are_bounded_and_one_per_contact(annotated, scheme):
    s, cm = annotated["1ao7"]
    ann = peptide_positions(cm, s)
    w = position_weights(ann, scheme)
    assert w.shape == (ann.height,)
    assert (w >= 0).all() and (w <= 1).all()


@pytest.mark.slow
def test_central_weighting_peaks_in_the_middle(annotated):
    s, cm = annotated["1ao7"]
    ann = peptide_positions(cm, s).with_columns(
        pl.Series("w", position_weights(peptide_positions(cm, s), "central")))
    by_pos = ann.group_by("peptide.pos").agg(pl.col("w").first()).sort("peptide.pos")
    w = by_pos["w"].to_numpy()
    assert w.argmax() not in (0, len(w) - 1)
    assert w[0] < w.max() and w[-1] < w.max()


@pytest.mark.slow
def test_tcr_facing_weighting_zeroes_exactly_the_anchors(annotated):
    s, cm = annotated["1ao7"]
    ann = peptide_positions(cm, s)
    w = position_weights(ann, "tcr_facing")
    is_anchor = np.asarray(ann["peptide.role"].to_list()) == "anchor"
    assert (w[is_anchor] == 0).all()
    assert (w[~is_anchor] == 1).all()


def test_unknown_scheme_raises():
    with pytest.raises(ValueError, match="scheme must be"):
        position_weights(pl.DataFrame({"peptide.pos": [1]}), "middle")


# --- the profile ----------------------------------------------------------------------------------
@pytest.mark.slow
@pytest.mark.parametrize("pdb_id", ["1ao7", "4ozg"])
def test_profile_sums_to_the_interface_energy(annotated, pdb_id):
    """A decomposition that does not add up is not a decomposition."""
    from tcren.pipeline import _interface_energy
    from tcren.potential import tcren

    s, cm = annotated[pdb_id]
    prof = position_profile(cm, tcren(), s)
    total = _interface_energy(cm.interface("tcr_peptide"), tcren())
    assert prof["phi"].sum() == pytest.approx(total, abs=1e-9)


@pytest.mark.slow
def test_profile_contact_counts_sum_to_the_interface_size(annotated):
    from tcren.potential import tcren

    s, cm = annotated["1ao7"]
    prof = position_profile(cm, tcren(), s)
    assert prof["n_contacts"].sum() == cm.interface("tcr_peptide").height


@pytest.mark.slow
def test_central_strain_is_a_subset_of_the_profile(annotated):
    from tcren.potential import tcren

    s, cm = annotated["1ao7"]
    prof = position_profile(cm, tcren(), s)
    strain = central_strain(prof)
    assert np.isfinite(strain)
    assert abs(strain) <= abs(prof["phi"]).sum() + 1e-9


def test_central_strain_of_an_empty_profile_is_nan():
    assert np.isnan(central_strain(pl.DataFrame({"peptide.pos": [], "phi": []})))


# --- the contact-type filter ----------------------------------------------------------------------
@pytest.mark.slow
def test_type_filter_drops_only_proximity_contacts(annotated):
    from tcren.contact_types import UNTYPED, residue_pair_types, type_weights

    s, _cm = annotated["1ao7"]
    typed = residue_pair_types(s, "tcr_peptide")
    w = type_weights(typed)
    dropped = np.asarray(typed["contact.type"].to_list())[w == 0]
    assert set(dropped) <= set(UNTYPED)
    assert 0 < w.sum() < len(w)


def test_type_filter_needs_the_v2_booleans():
    from tcren.contact_types import type_weights

    with pytest.raises(ValueError, match="is_<type>"):
        type_weights(pl.DataFrame({"contact.type": ["polar"]}))
