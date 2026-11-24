import polars as pl

from tcren.contacts.veto import contact_veto_features


def test_veto_uses_contacted_edges_and_preserves_full_graph():
    rows = pl.DataFrame({
        "pdb.id": ["pose"] * 5,
        "interface": ["tcr_peptide"] * 5,
        "chain.id.from": ["R"] * 5,
        "chain.id.to": ["P"] * 5,
        "residue.index.from": [0, 0, 1, 1, 2],
        "residue.index.to": [0, 1, 0, 1, 2],
        "aa.from": ["A", "A", "C", "C", "A"],
        "aa.to": ["D", "E", "D", "E", "D"],
        "contacted": [True, True, True, True, False],
    })
    effects = {"tcr_peptide/A|D": 1., "tcr_peptide/A|E": 1., "tcr_peptide/C|D": 1.}
    scored = contact_veto_features(rows, effects).row(0, named=True)
    assert scored["id"] == "pose"
    assert scored["full_n_edges"] == 4
    assert scored["full_b1"] == 1
    assert scored["kept_n_edges"] == 3
    assert scored["kept_b1"] == scored["vetoed_b1"] == 0
    assert scored["kept_edge_frac"] == .75
    assert scored["n_unknown_edges"] == scored["vetoed_n_edges"] == 1
    assert scored["core_size"] == 1.
