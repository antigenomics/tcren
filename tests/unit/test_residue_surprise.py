"""Small exact examples for molecular contact cycles and mock-tail surprise."""

import numpy as np
import polars as pl
import pytest

from tcren.contacts.residue_surprise import residue_counts, fit_reference, score_residues


def test_covalent_residues_make_repeated_contacts_cycles():
    atoms = pl.DataFrame(
        {
            "pdb.id": ["x"] * 15,
            "chain.id": ["A"] * 5 + ["P"] * 10,
            "residue.index": [1] * 5 + [1] * 5 + [2] * 5,
            "residue.aa": ["A"] * 15,
            "chain.type": ["TRA"] * 5 + ["PEPTIDE"] * 10,
            "region.type": ["CDR1"] * 5 + ["PEPTIDE"] * 10,
            "atom.name": ["N", "CA", "C", "O", "CB"] * 3,
        }
    )
    contacts = pl.DataFrame(
        {
            "pdb.id": ["x"] * 3,
            "chain.id.from": ["A"] * 3,
            "residue.index.from": [1] * 3,
            "atom.from": ["CB"] * 3,
            "chain.id.to": ["P"] * 3,
            "residue.index.to": [1, 1, 2],
            "atom.to": ["N", "CA", "CB"],
            "dist": [4.0] * 3,
        }
    )
    rows = residue_counts(atoms, contacts)
    a = rows.filter(pl.col("chain") == "TRA").row(0, named=True)
    assert (a["occupied"], a["pairs"], a["partners"], a["cycles"]) == (1, 3, 2, 1)
    # The same physical atom is counted repeatedly; internal aromatic cycles
    # never enter m-q. Reciprocal peptide counts retain the distinct atom pairs.
    assert rows.filter(pl.col("chain") == "PEPTIDE")["cycles"].to_list() == [1, 0]
    side = residue_counts(atoms, contacts, "sidechain")
    assert side["available"].to_list() == [1, 1, 1]
    assert side.filter(pl.col("chain") == "TRA")["cycles"].item() == 1
    assert side.filter(pl.col("chain") == "PEPTIDE")["cycles"].to_list() == [0, 0]
    with pytest.raises(ValueError, match="Duplicate"):
        residue_counts(atoms, pl.concat([contacts, contacts.head(1)]))


def test_empirical_tail_keeps_ties_zeros_and_missing_support():
    reference = pl.DataFrame(
        {
            "id": ["a", "b", "c"],
            "chain": ["TRA"] * 3,
            "region": ["CDR1"] * 3,
            "aa": ["A"] * 3,
            "available": [5] * 3,
            "complete": [True] * 3,
            "occupied": [0, 1, 2],
            "pairs": [0, 2, 3],
            "cycles": [0, 1, 2],
        }
    )
    resource = fit_reference(reference)
    result = score_residues(
        reference.with_columns(pl.col("available").cast(pl.UInt32)), resource, "cycles"
    )
    np.testing.assert_allclose(result["tail_probability"], [1, 3 / 4, 2 / 4])
    np.testing.assert_allclose(result["bits"], [0, -np.log2(3 / 4), 1])
    weighted = fit_reference(reference.with_columns(pl.Series("weight", [1.0, 1.0, 2.0])))
    weighted_rows = score_residues(reference, weighted, "cycles")
    np.testing.assert_allclose(weighted_rows["tail_probability"], [1.0, 0.8, 0.6])
    assert weighted_rows["reference_residues"].to_list() == [3, 3, 3]
    assert weighted_rows["reference_mass"].to_list() == [4.0, 4.0, 4.0]
    positive = fit_reference(reference.with_columns(pl.lit(2).alias("cycles")))
    relative = score_residues(reference, resource, "cycles", positive)
    np.testing.assert_allclose(relative["log_tail_ratio"], [0, np.log2(4 / 3), 1])
    absent = reference.with_columns(pl.lit("R").alias("aa"))
    assert score_residues(absent, resource, "cycles")["bits"].null_count() == 3
    with pytest.raises(ValueError, match="nonnegative integers"):
        fit_reference(reference.with_columns(pl.lit(-1).alias("cycles")))


def test_evaluation_keeps_labelled_rows_outside_old_mask(tmp_path):
    from tcren.contacts.residue_surprise import evaluate_command

    baseline = pl.DataFrame(
        dict(
            id=["a", "b", "c", "d"],
            dataset=["d"] * 4,
            stratum=["s"] * 4,
            epitope=["e"] * 4,
            y=[0, 1, 0, 1],
            iptm=[0.1, 0.2, 0.3, 0.4],
            prepared=[0.1, 0.2, 0.3, 0.4],
            matched=[True, True, False, False],
        )
    )
    baseline.write_parquet(tmp_path / "baseline.parquet")
    folder = tmp_path / "cycles"
    folder.mkdir()
    pl.DataFrame(
        dict(id=["a", "b", "c", "d"], percentile_geometric=[0.0, 1.0, 0.0, 1.0])
    ).write_csv(folder / "scores.tsv", separator="\t")
    evaluate_command(
        baseline=tmp_path / "baseline.parquet",
        scores=[folder / "scores.tsv"],
        output=tmp_path / "out",
        comparators=None,
        weights=[0.25],
        enrichment=False,
    )
    result = pl.read_csv(tmp_path / "out/metrics.tsv", separator="\t")
    assert set(result.filter(pl.col("population") == "all_labelled")["n"]) == {4}
    assert set(result.filter(pl.col("population") == "observed_baselines")["n"]) == {2}


def test_cycle_increment_matches_explicit_covalent_ring_graph():
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import connected_components

    # Two internally cyclic residues and one linear residue; all atoms are
    # linked by covalent edges before the four observed inter-residue contacts.
    covalent = [(0, 1), (1, 2), (2, 0), (3, 4), (4, 5), (5, 3), (6, 7)]
    contact = [(0, 3), (0, 4), (2, 5), (1, 6)]

    def beta(edges):
        graph = coo_matrix(([1] * len(edges), np.array(edges).T), shape=(8, 8))
        return len(edges) - 8 + connected_components(graph, directed=False, return_labels=False)

    assert beta(covalent + contact) - beta(covalent) == len(contact) - 2 == 2


def test_mock_fit_excludes_benchmark_pairs_and_retains_all_controls(tmp_path):
    import json
    from tcren.contacts.residue_surprise import fit_command
    import hashlib

    def write_identity(path, specs):
        pl.DataFrame(
            [
                dict(id=i, TRA=a, TRB=b, PEPTIDE=p, sha256=hashlib.sha256(i.encode()).hexdigest())
                for i, a, b, p in specs
            ]
        ).write_csv(path.with_suffix(".identities.tsv"), separator="\t")

    def table(specs):
        rows = []
        for model, a, b, p in specs:
            for chain, region, aa in [
                ("TRA", "CDR3", a),
                ("TRB", "CDR3", b),
                ("PEPTIDE", "PEPTIDE", p),
            ]:
                rows.append(
                    dict(
                        id=model,
                        chain_id=chain,
                        residue_index=1,
                        chain=chain,
                        region=region,
                        aa=aa,
                        available=5,
                        complete=True,
                        occupied=0,
                        pairs=0,
                        cycles=0,
                    )
                )
        return pl.DataFrame(rows)

    ref = table(
        [("a", "A", "R", "G"), ("b", "A", "A", "G"), ("c", "A", "A", "F"), ("d", "R", "R", "F")]
    )
    ref.write_parquet(tmp_path / "ref.parquet")
    write_identity(
        tmp_path / "ref.parquet",
        [("a", "A", "R", "G"), ("b", "A", "A", "G"), ("c", "A", "A", "F"), ("d", "R", "R", "F")],
    )
    write_identity(tmp_path / "query.parquet", [("q", "A", "A", "G")])
    table([("q", "A", "A", "G")]).write_parquet(tmp_path / "query.parquet")
    pl.DataFrame({"pdb.id": ["a", "b", "c", "d"], "y": [0] * 4}).write_csv(
        tmp_path / "meta.tsv", separator="\t"
    )
    fit_command(
        features=tmp_path / "ref.parquet",
        query_features=tmp_path / "query.parquet",
        metadata=tmp_path / "meta.tsv",
        output=tmp_path / "ref.json",
        reference_label=0,
        balance_epitopes=False,
    )
    result = json.loads((tmp_path / "ref.json").read_text())
    assert result["excluded_models"] == 1
    assert result["models"] == 3
    assert result["downsampled_models"] == 0
    identities = pl.read_csv(tmp_path / "ref.identities.tsv", separator="\t").with_columns(
        pl.when(pl.col("id") == "b").then(pl.lit("VA")).otherwise(pl.col("TRA")).alias("TRA")
    )
    identities.write_csv(tmp_path / "ref.identities.tsv", separator="\t")
    fit_command(
        features=tmp_path / "ref.parquet",
        query_features=tmp_path / "query.parquet",
        metadata=tmp_path / "meta.tsv",
        output=tmp_path / "different-full-TCR.json",
        reference_label=0,
        balance_epitopes=False,
    )
    assert json.loads((tmp_path / "different-full-TCR.json").read_text())["models"] == 4
    pl.DataFrame({"pdb.id": ["a", "b", "c", "d"], "y": [0, None, 0, 0]}).write_csv(
        tmp_path / "unknown.tsv", separator="\t"
    )
    with pytest.raises(ValueError, match="Metadata"):
        fit_command(
            features=tmp_path / "ref.parquet",
            query_features=tmp_path / "query.parquet",
            metadata=tmp_path / "unknown.tsv",
            output=tmp_path / "unknown.json",
            reference_label=0,
            balance_epitopes=False,
        )

    # Neither rare epitopes nor repeated CDR3 pairs cause control removal.
    specs = [
        ("a", "A", "R", "G"),
        ("b", "R", "R", "F"),
        ("c", "N", "R", "F"),
        ("d", "D", "R", "F"),
        ("e", "C", "R", "Y"),
        ("f", "E", "R", "Y"),
        ("g", "Q", "R", "Y"),
        ("h", "Q", "R", "Y"),
    ]
    table(specs).write_parquet(tmp_path / "ref2.parquet")
    write_identity(tmp_path / "ref2.parquet", specs)
    pl.DataFrame({"pdb.id": [r[0] for r in specs], "y": [0] * 8}).write_csv(
        tmp_path / "meta2.tsv", separator="\t"
    )
    fit_command(
        features=tmp_path / "ref2.parquet",
        query_features=tmp_path / "query.parquet",
        metadata=tmp_path / "meta2.tsv",
        output=tmp_path / "ref2.json",
        reference_label=0,
        balance_epitopes=False,
    )
    second = json.loads((tmp_path / "ref2.json").read_text())
    assert second["models"] == 8
    assert second["downsampled_models"] == 0
