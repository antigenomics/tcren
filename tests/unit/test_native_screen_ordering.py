"""Pose failures and unavailable cases cannot outrank a scored interface."""
import polars as pl

from tcren.native_ranking import rank_native_screen


def test_native_order_uses_only_scored_contact_bearing_poses():
    data = pl.DataFrame({
        "unit": ["a"] * 5 + ["b"] * 2,
        "iptm": [0.1, 0.2, 0.8, 0.8, 0.9, 0.3, 0.4],
        "dPhi_peptide_sum": [1.0] * 7,
        "d_charge": [0.1] * 7,
        "coverage": [0.5] * 7,
        "n_atom_contacts_tp": [0, 10, 20, 30, None, 2, 0],
    })
    result = rank_native_screen(data, task="epitope", ranking_unit="unit", min_cohort=2)
    assert result["screen_score"][0] == 0.0
    assert result["screen_score"][1] < 0.0
    assert result["screen_rank"].to_list() == [None, 2, 1, 1, None, None, None]
