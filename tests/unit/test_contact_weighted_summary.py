"""Signed weighted loss identities retain zero-contact models and denominators."""
import gzip

import polars as pl


def test_weighted_summary_signed_identity_and_zero_models(tmp_path):
    from tcren.contacts.weighted_summary import summarize_weighted
    from tcren.contacts.plausibility import comparison_command
    run, accepted, previous, out = [tmp_path/p for p in ["run", "accepted", "previous", "out"]]
    for root in [run, accepted]:
        (root/"scoring").mkdir(parents=True)
    previous.mkdir()
    labels = pl.DataFrame(dict(id=["a", "b", "c"], dataset=["test"]*3, task=["ranking"]*3,
        cohort=["case"]*3, mhc_class=["I"]*3, template_coverage=["unknown"]*3,
        y=[1, 0, 1], iptm=[.9, .2, .5]))
    scores = pl.DataFrame(dict(id=["a", "b", "c"], potential=["TCRen2"]*3, mhc_class=["I"]*3,
        original=[1., 0., 0.], weighted=[0., 0., 0.], n_contacts=[2, 0, 0], missing_potential_entries=[0]*3))
    edges = pl.DataFrame({"pdb.id": ["a"]*2, "potential": ["TCRen2"]*2,
        "chain.id.from": ["A"]*2, "chain.id.to": ["P"]*2,
        "residue.index.from": [0, 1], "residue.index.to": [2, 2],
        "original_contribution": [-2., 3.], "contact_weight": [.75, .5],
        "weighted_contribution": [-1.5, 1.5], "weighted_increment": [.5, -1.5],
        "geometry_status": ["supported_cdr_geometry"]*2,
        "calibration_source": ["native_tcr_peptide", "external_antibody_antigen"],
        "n_atom_pairs": [2, 3], "chain.type.from": ["TRA"]*2,
        "region.type.from": ["CDR3"]*2, "pos.to": [1, 1],
        "residue.aa.from": ["S", "R"], "residue.aa.to": ["G"]*2})
    edges = edges.with_columns(pl.lit(2).alias("n_distance"),
        *[pl.lit(1).alias("n_"+family) for family in ["omega", "theta", "phi"]],
        pl.lit(5).alias("calibration_n_edges"), pl.lit(2).alias("calibration_n_structures"),
        pl.lit(.3).alias("minimum_tail"))
    ledger = pl.DataFrame(dict(id=["a", "b", "c"], sha256=["raw-a", "raw-b", "raw-c"], status=["ok"]*3))
    for root in [run, accepted]:
        labels.write_csv(root/"labels.tsv", separator="\t")
        scores.write_csv(root/"scoring/scores.tsv", separator="\t")
        ledger.write_csv(root/"scoring/accounting.tsv", separator="\t")
        with gzip.open(root/"scoring/contacts.tsv.gz", "wt") as stream:
            stream.write(edges.write_csv(separator="\t"))
    scores.with_columns(pl.Series("chemical_real", [1, 0, 0]), pl.Series("chemical_proximity", [1, 0, 0])).write_csv(previous/"models.tsv", separator="\t")
    edges.with_columns(pl.Series("union_type", ["polar", "vdw"]),
        pl.Series("chemical_real", [True, False]), pl.Series("chemical_proximity", [False, True])).rename({"pdb.id": "id"}).write_parquet(previous/"contacts.parquet")
    comparison_command(run/"labels.tsv", run/"scoring/scores.tsv", run/"comparison")
    summarize_weighted(run, accepted, previous, out)
    models = pl.read_csv(out/"models.tsv", separator="\t")
    a = models.filter(pl.col("id")=="a").row(0, named=True)
    assert a["positive_magnitude_lost"]==1.5
    assert a["favourable_magnitude_lost"]==.5
    assert a["weighted_increment"]==a["weighted"]-a["original"]==-1.
    assert models.filter(pl.col("id")=="c")["contact_weight"].item()==0
    motifs = pl.read_csv(out/"motifs.tsv", separator="\t")
    assert motifs["n_models"].to_list()==[2, 2]
    ranks = pl.read_csv(out/"rank_changes.tsv", separator="\t")
    assert ranks.height==2
    assert ranks.filter(pl.col("arm")=="weighted")["delta_auc"].item()==.25
