import polars as pl
import pytest
from tcren.validation.prediction_comparison import compare_peptide_contact_sets


def test_empty_complete_complex_differs_from_missing_component():
    chains = pl.DataFrame(
        [
            {"complex.id": model, "chain.type": role, "sequence": seq}
            for model in ["reference", "candidate", "empty", "missing"]
            for role, seq in [("TRA", "AC"), ("TRB", "DE"), ("PEPTIDE", "FG")]
            if not (model == "missing" and role == "TRB")
        ]
    )
    rows = []
    for model, pairs in [
        ("reference", [("TRA", 0, 0), ("TRB", 1, 1)]),
        ("candidate", [("TRA", 0, 0), ("TRA", 1, 1)]),
    ]:
        for role, i, j in pairs:
            rows.append(
                {
                    "complex.id": model,
                    "chain.type.from": role,
                    "chain.type.to": "PEPTIDE",
                    "residue.index.from": i,
                    "residue.index.to": j,
                    "residue.aa.from": {"TRA": "AC", "TRB": "DE"}[role][i],
                    "residue.aa.to": "FG"[j],
                }
            )
    contacts = pl.DataFrame(rows)
    output = compare_peptide_contact_sets(contacts, chains, reference_id="reference")
    lookup = {r["complex.id"]: r for r in output.iter_rows(named=True)}
    assert lookup["candidate"]["contact_jaccard"] == pytest.approx(1 / 3)
    assert lookup["candidate"]["shared_residue_pairs"] == 1
    assert lookup["empty"]["contact_jaccard"] == 0
    assert lookup["empty"]["candidate_residue_pairs"] == 0
    assert lookup["missing"]["contact_jaccard"] is None
    assert lookup["missing"]["candidate_residue_pairs"] is None
    with pytest.raises(ValueError, match="Repeated directed residue contact"):
        compare_peptide_contact_sets(
            pl.concat([contacts, contacts.head(1)]), chains, reference_id="reference"
        )
    with pytest.raises(ValueError, match="amino acids differ"):
        compare_peptide_contact_sets(
            contacts.with_columns(pl.lit("W").alias("residue.aa.to")),
            chains,
            reference_id="reference",
        )
