"""Boundary regressions identified by the release review; no benchmark fitting."""
from copy import deepcopy
from dataclasses import replace

import numpy as np
import polars as pl
import pytest
from polars.testing import assert_frame_equal

from test_contact_observations import complex_fixture
from test_motif_screen import table
from tcren.contacts.additive import contact_counts, oriented_contacts
from tcren.contacts.glm import atom_totals
from tcren.contacts.observations import contact_observations, select_conformer
from tcren.screen import CohortKind, bundled_reference, screen_scores


def test_contact_counts_accept_both_physical_endpoint_orders():
    structure = complex_fixture()
    reverse = contact_observations(structure, interface="tcr_peptide")
    assert reverse["role.from"].unique().to_list() == ["peptide"]
    oriented = oriented_contacts(reverse)
    assert oriented["theta.from_rad"].equals(reverse["theta.to_rad"].rename("theta.from_rad"))
    forward_structure = replace(structure, chains=[
        replace(c, chain_id={"T": "A", "P": "Z"}.get(c.chain_id, c.chain_id))
        for c in structure.chains])
    forward = contact_observations(forward_structure, interface="tcr_peptide")
    assert forward["role.from"].unique().to_list() == ["tcr"]
    ids = pl.DataFrame({"pdb.id": [structure.pdb_id, "zero"]})
    for function in (contact_counts, atom_totals):
        assert_frame_equal(function(reverse, ids), function(forward, ids))
        both = pl.concat([forward.with_columns(pl.lit("forward").alias("pdb.id")),
                          reverse.with_columns(pl.lit("reverse").alias("pdb.id"))])
        result = function(both, pl.DataFrame({"pdb.id": ["forward", "reverse"]}))
        assert result.drop("pdb.id").row(0) == result.drop("pdb.id").row(1)


def test_numeric_topology_identifiers_survive_additive_join(tmp_path, monkeypatch):
    from tcren.contacts import additive

    ids = ["0", "01"]
    pl.DataFrame({"id": ids, "dataset": ["tcrvdb"] * 2, "epitope": ["E"] * 2,
                  "y": [0, 1], "iptm": [.7, .8], "prepared": [.5, .6]}).write_parquet(tmp_path / "base.parquet")
    contacts = contact_observations(complex_fixture(), interface="tcr_peptide")
    pl.concat([contacts.with_columns(pl.lit(i).alias("pdb.id")) for i in ids]).write_parquet(tmp_path / "contacts.parquet")
    pl.DataFrame({"pdb.id": ids, "graph": ["tp_pm_tm"] * 2, "available": [True] * 2,
                  "cyclo_frac": [.3, .4]}).write_csv(tmp_path / "topology.tsv", separator="\t")
    pl.DataFrame({"pdb.id": [i for i in ids for _ in range(2)],
                  "chain.type": ["TRA", "TRB"] * 2, "region.type": ["CDR3"] * 4,
                  "residue.index": [0] * 4, "residue.aa": ["A", "R"] * 2}).write_parquet(tmp_path / "atoms.parquet")
    def stop_after_join(models, *args):
        assert models["id"].to_list() == ids
        raise RuntimeError("validated joins; no fitting requested")
    monkeypatch.setattr(additive, "splits", stop_after_join)
    monkeypatch.setattr(additive, "candidate_spec", lambda: [])
    with pytest.raises(RuntimeError, match="validated joins"):
        additive.command(tmp_path / "contacts.parquet", tmp_path / "atoms.parquet",
                         tmp_path / "topology.tsv", tmp_path / "base.parquet", tmp_path / "out", 1)


def test_confidence_only_resource_override_and_undefined_structural_ratio():
    resource = deepcopy(bundled_reference(CohortKind.TCR_RANKING))
    resource["motif_confidence"]["parameters"] = {"b": 0., "m": 0.}
    with np.errstate(all="raise"):
        out = screen_scores(table(), resource, min_cohort=1, weight=1.)
        np.testing.assert_array_equal(out["screen_score"], table()["iptm"])
        assert out["topology_weight"].to_list() == [0., 0.]
        assert out["motif_weight"].to_list() == [0., 0.]
        with pytest.raises(ValueError, match="defined Betti:motif ratio"):
            screen_scores(table(), resource, min_cohort=1, weight=.5)


def test_neighborhood_contacts_and_stacking_share_selected_conformer(monkeypatch):
    from tcren.topology import neighborhood
    from tcren.contact_types import stacked_pairs

    original = complex_fixture()
    selected = select_conformer(original)
    assert any(a.altloc == "B" for c in original.chains for r in c.residues for a in r.atoms)
    seen = []
    def selected_stacking(structure):
        assert not any(a.altloc == "B" for c in structure.chains for r in c.residues for a in r.atoms)
        # The whole MHC (including non-helix floor) must remain present.
        assert len(structure.chain("M").residues) == 2
        seen.append(structure)
        return stacked_pairs(structure)
    monkeypatch.setattr(neighborhood, "stacked_pairs", selected_stacking)
    assert neighborhood.neighborhood_features(original) == neighborhood.neighborhood_features(selected)
    assert_frame_equal(neighborhood.neighborhood_contacts(original),
                       neighborhood.neighborhood_contacts(selected))
    assert len(seen) == 4
