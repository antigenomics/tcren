"""Unique atom/edge counts and unpenalized coefficient reconstruction."""

import numpy as np
import polars as pl
from scipy.special import expit
from tcren.contacts.glm import atom_totals, fit_glm


def test_totals_and_unpenalized_fit():
    c = pl.DataFrame(
        {
            "pdb.id": ["x"] * 3,
            "interface": ["tcr_peptide"] * 3,
            "role.from": ["tcr"] * 3,
            "role.to": ["peptide"] * 3,
            "dist": [3.0] * 3,
            "chain.id.from": ["A", "A", "B"],
            "chain.type.from": ["TRA", "TRA", "TRB"],
            "residue.index.from": [1, 1, 1],
            "atom.from": ["CA", "CB", "CA"],
            "chain.id.to": ["P"] * 3,
            "residue.index.to": [1] * 3,
            "atom.to": ["O", "O", "O"],
        }
    )
    r = atom_totals(c, pl.DataFrame({"pdb.id": ["x", "zero"]}))
    assert r.select("peptide_atoms", "alpha_atoms", "beta_atoms", "tcr_peptide_edges").rows() == [
        (1, 2, 1, 2),
        (0, 0, 0, 0),
    ]
    x = np.array([-2, -1, 0, 0, 1, 2.0])[:, None]
    y = np.array([0, 1, 0, 1, 0, 1])
    p, fit = fit_glm(x, y)
    np.testing.assert_allclose(p, expit(fit["intercept"] + x @ np.array(fit["coefficients"])))
    assert p[0] < p[-1] and fit["max_abs_gradient"] < 1e-6


def test_cli_outputs_both_baselines_and_fits(tmp_path):
    from tcren.contacts.glm import command

    ids = [str(i) for i in range(6)]
    iptm = np.array([0.1, 0.3, 0.4, 0.5, 0.7, 0.9])
    pl.DataFrame(
        {
            "id": ids,
            "dataset": ["tcrvdb"] * 6,
            "epitope": ["E"] * 6,
            "y": [0, 1, 0, 1, 0, 1],
            "iptm": iptm,
            "prepared": 0.55 * iptm + 0.45 * 0.4,
        }
    ).write_parquet(tmp_path / "baseline.parquet")
    pl.DataFrame(
        {"pdb.id": ids, "graph": ["tp_pm_tm"] * 6, "cyclo_frac": [0.4] * 6, "available": [True] * 6}
    ).write_csv(tmp_path / "topology.tsv", separator="\t")
    pl.DataFrame(
        {
            "pdb.id": ids,
            "interface": ["tcr_peptide"] * 6,
            "role.from": ["tcr"] * 6,
            "role.to": ["peptide"] * 6,
            "dist": [3.0] * 6,
            "chain.id.from": ["A"] * 6,
            "chain.type.from": ["TRA"] * 6,
            "residue.index.from": [1] * 6,
            "atom.from": ["CA"] * 6,
            "chain.id.to": ["P"] * 6,
            "residue.index.to": [1] * 6,
            "atom.to": ["O"] * 6,
        }
    ).write_parquet(tmp_path / "contacts.parquet")
    command(
        tmp_path / "contacts.parquet",
        tmp_path / "topology.tsv",
        tmp_path / "baseline.parquet",
        tmp_path / "out",
    )
    result = pl.read_parquet(tmp_path / "out/predictions.parquet")
    assert result.height == 24
    assert set(result["arm"]) == {"iptm", "betti", "raw_glm", "normalized_glm"}
