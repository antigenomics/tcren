import polars as pl

from tcren.contacts.mixture import atom_contact_mixture


def test_atom_mixture_averages_within_residue_pairs_and_accounts_for_unknown_cells():
    rows = pl.DataFrame({
        "pdb.id": ["mixed"] * 3 + ["positive"],
        "interface": ["tcr_peptide"] * 4,
        "chain.id.from": ["T"] * 4,
        "chain.id.to": ["P"] * 4,
        "residue.index.from": [0, 0, 1, 0],
        "residue.index.to": [0, 0, 0, 0],
        "atom.from": ["CA", "CB", "CA", "CB"],
        "atom.to": ["O", "O", "N", "O"],
    })
    effects = {"tcr_peptide/CA|O": -2., "tcr_peptide/CB|O": 2.}
    result = atom_contact_mixture(rows, effects, ["empty", "positive", "mixed"])
    assert result["id"].to_list() == ["empty", "positive", "mixed"]
    assert result["pi_atom"][0] is None
    assert result["pi_atom"][1] > 0.99
    assert result["pi_atom"][2] == 0.5
    assert result["n_atom_pairs_unknown"].to_list() == [0, 0, 1]
    assert result["n_residue_pairs"].to_list() == [0, 1, 2]
