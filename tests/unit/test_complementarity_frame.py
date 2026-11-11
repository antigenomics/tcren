"""Rigid frame contracts, including degeneracy and shared-reference overlays."""
import numpy as np
import pytest

from tcren.complementarity import canonical_frame, fit_frame

PEPTIDE = np.array([[-2., 0, 2], [0, 0, 3], [2, 0, 2]])
GROOVE = np.array([[-3., -2, 0], [3, -2, 0], [3, 2, 0], [-3, 2, 0]])


def test_apex_axes_and_chirality():
    frame = canonical_frame(PEPTIDE, GROOVE, reference_id="ref", origin="apex")
    np.testing.assert_allclose(frame.transform(PEPTIDE)[1], 0)
    assert frame.transform(PEPTIDE)[-1, 0] > frame.transform(PEPTIDE)[0, 0]
    assert np.linalg.det(frame.rotation) == pytest.approx(1)
    assert frame.transform(GROOVE)[:, 2].max() < 0


def test_mass_origin_is_not_centroid():
    frame = canonical_frame(PEPTIDE, GROOVE, reference_id="ref",
                            peptide_atoms=PEPTIDE, peptide_masses=[1, 10, 1])
    np.testing.assert_allclose(frame.transform(np.average(PEPTIDE, axis=0,
                                                        weights=[1, 10, 1])), 0)
    with pytest.raises(ValueError, match="atoms and masses"):
        canonical_frame(PEPTIDE, GROOVE, reference_id="ref")


def test_rotated_mobile_inherits_reference_without_recentering():
    reference = canonical_frame(PEPTIDE, GROOVE, reference_id="ref", origin="apex")
    rotation = np.array([[0., 1, 0], [-1, 0, 0], [0, 0, 1]])
    translation = np.array([9., -4, 7])
    mobile = GROOVE @ rotation + translation
    fit = fit_frame(mobile, GROOVE, reference)
    np.testing.assert_allclose(fit.transform(PEPTIDE @ rotation + translation),
                               reference.transform(PEPTIDE), atol=1e-12)
    assert fit.rmsd < 1e-12
    assert fit.anchor_count == 4
    moved_peptide = (PEPTIDE + [0, 0, 1]) @ rotation + translation
    np.testing.assert_allclose(fit.transform(moved_peptide),
                               reference.transform(PEPTIDE) + [0, 0, 1], atol=1e-12)


@pytest.mark.parametrize("bad", [np.zeros((3, 3)), [[0, 0, 0], [1, 0, 0]],
                                     [[0, 0, 0], [1, 0, 0], [np.nan, 1, 0]]])
def test_degenerate_groove_fails(bad):
    with pytest.raises(ValueError):
        canonical_frame(PEPTIDE, bad, reference_id="ref", origin="apex")


def test_collinear_fit_fails():
    frame = canonical_frame(PEPTIDE, GROOVE, reference_id="ref", origin="apex")
    with pytest.raises(ValueError, match="collinear"):
        fit_frame(np.zeros((3, 3)), np.zeros((3, 3)), frame)
