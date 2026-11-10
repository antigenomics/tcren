import numpy as np
import pytest

from tcren.complementarity.prepare import ranking_frame, reference_frame
from tcren.structure.model import Atom, Chain, RegionMarkup, Residue, Structure


def chain(identity, role, points, region):
    residues = [Residue(i, i + 1, "", "A", "ALA",
                        (Atom("CA", "C", np.array(point, dtype=float)),))
                for i, point in enumerate(points)]
    return Chain(identity, residues, role, regions=[
        RegionMarkup(region, 0, len(residues), "A" * len(residues), residues)])


def structure():
    return Structure("reference", [
        chain("P", "PEPTIDE", [[-2, 0, 2], [0, 0, 3], [2, 0, 2]], "PEPTIDE"),
        chain("M", "MHCa", [[-3, -2, 0], [3, -2, 0], [3, 2, 0], [-3, 2, 0]], "GROOVE_FLOOR"),
        chain("A", "TRA", [[-3, -1, 4], [-2, 2, 5], [0, 0, 6]], "FR1"),
        chain("B", "TRB", [[3, -1, 4], [2, 2, 5], [1, 0, 6]], "FR1"),
    ])


@pytest.mark.parametrize("task,expected", [("tcr", 4), ("epitope", 6)])
def test_fixed_components_are_selected_by_task(task, expected):
    fixed, mobile = structure(), structure()
    frame = reference_frame(fixed)
    aligned = ranking_frame(mobile, fixed, frame, task=task)
    assert aligned.anchor_count == expected
    assert aligned.rmsd < 1e-12


def test_epitope_frame_preserves_mobile_peptide_displacement():
    fixed, mobile = structure(), structure()
    frame = reference_frame(fixed)
    for res in mobile.chain("P").residues:
        res.atoms[0].coord[:] += [0, 0, 2]
    aligned = ranking_frame(mobile, fixed, frame, task="epitope")
    np.testing.assert_allclose(aligned.transform(mobile.chain("P").residues[0].ca),
                               frame.transform(fixed.chain("P").residues[0].ca) + [0, 0, 2], atol=1e-12)


def test_epitope_origin_is_fixed_tcr_mass_center():
    fixed = structure()
    frame = reference_frame(fixed, task="epitope")
    points = [r.ca for c in fixed.by_type("TRA", "TRB") for r in c.residues]
    np.testing.assert_allclose(frame.transform(np.mean(points, axis=0)), 0, atol=1e-12)
    assert frame.origin_component == "tcr-variable-domain"


def test_epitope_apex_uses_groove_facing_cdr_tip():
    fixed = structure()
    for chain_id in ("A", "B"):
        selected = fixed.chain(chain_id).residues[:1]
        fixed.chain(chain_id).regions.append(RegionMarkup("CDR3", 0, 1, "A", selected))
    frame = reference_frame(fixed, task="epitope", origin="apex")
    np.testing.assert_allclose(frame.transform(fixed.chain("A").residues[0].ca), 0, atol=1e-12)
    assert frame.origin_component == "tcr-variable-domain"


def test_missing_framework_and_mixed_construct_fail():
    fixed, mobile = structure(), structure()
    frame = reference_frame(fixed)
    mobile.chain("A").regions.clear()
    mobile.chain("B").regions.clear()
    with pytest.raises(ValueError):
        ranking_frame(mobile, fixed, frame, task="epitope")
    mobile.chains.append(chain("M2", "MHCa", [[0, 0, 0]], "GROOVE_FLOOR"))
    with pytest.raises(ValueError, match="multiple"):
        ranking_frame(mobile, fixed, frame)


def test_peptide_required_even_for_epitope_reference():
    source = structure()
    source.chains = [c for c in source.chains if c.chain_type != "PEPTIDE"]
    with pytest.raises(ValueError, match="exactly one peptide"):
        reference_frame(source)
