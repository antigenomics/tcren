"""Read and cross-check exported table summaries without rebuilding contacts."""
from __future__ import annotations

from collections import Counter
import json
from pathlib import Path

import polars as pl


def summarize_exports(directory: str | Path) -> dict:
    """Count exported rows and verify graph incidence totals against summary rows.

    Read one partition at a time. This verifies serialization/count consistency;
    it does not certify a structural pose or reproduce an experimental endpoint.
    """
    directory = Path(directory)
    summary = pl.read_parquet(directory / "summary.parquet")
    if summary["structure_id"].n_unique() != summary.height:
        raise ValueError("duplicate structure IDs in summary")
    expected = {row["structure_id"]: row for row in summary.to_dicts()}
    fields = {"atoms": "atoms", "residues": "residues", "bonds": "covalent_bonds",
              "atom_contacts": "atom_contacts", "residue_contacts": "residue_contacts"}
    for table, quantity in fields.items():
        counts, incidence = Counter(), Counter()
        extra = "contact_count" if table in {"atoms", "residues", "residue_contacts"} else None
        columns = ["structure_id"] + ([extra] if extra else [])
        for path in sorted((directory / "tables" / table).glob("*.parquet")):
            data = pl.read_parquet(path, columns=columns)
            for row in data.group_by("structure_id").agg(pl.len().alias("rows")).iter_rows(named=True):
                counts[row["structure_id"]] += row["rows"]
            if extra:
                for row in data.group_by("structure_id").agg(pl.col(extra).sum().alias("incidence")).iter_rows(named=True):
                    incidence[row["structure_id"]] += row["incidence"]
        if counts.keys() - expected.keys():
            raise ValueError(f"{table}: unknown structure IDs")
        for identity, row in expected.items():
            if counts[identity] != row[quantity]:
                raise ValueError(f"{identity}: {table} row count differs from summary")
            if extra:
                factor = 1 if table == "residue_contacts" else 2
                if incidence[identity] != factor * row["atom_contacts"]:
                    raise ValueError(f"{identity}: {table} contact incidence differs")
    peptide_contacts = Counter()
    for path in sorted((directory / "tables/residue_contacts").glob("*.parquet")):
        roles = pl.read_parquet(directory / "tables/residues" / path.name,
                                columns=["structure_id", "id", "role"])
        edges = pl.read_parquet(path)
        paired = edges.join(roles.rename({"id": "source", "role": "source_role"}),
                            on=["structure_id", "source"], how="left", validate="m:1").join(
            roles.rename({"id": "target", "role": "target_role"}),
            on=["structure_id", "target"], how="left", validate="m:1")
        if paired["source_role"].null_count() or paired["target_role"].null_count():
            raise ValueError("residue contacts contain unknown residue IDs")
        selected = paired.filter(
            (pl.col("source_role").is_in(["alpha", "beta"]) & (pl.col("target_role") == "peptide")) |
            (pl.col("target_role").is_in(["alpha", "beta"]) & (pl.col("source_role") == "peptide")))
        for row in selected.group_by("structure_id").agg(pl.col("contact_count").sum()).iter_rows(named=True):
            peptide_contacts[row["structure_id"]] += row["contact_count"]
    for identity, row in expected.items():
        if peptide_contacts[identity] != row["tcr_peptide_atom_contacts"]:
            raise ValueError(f"{identity}: TCR–peptide contacts differ from summary")
        if row["no_tcr_peptide_interface"] != (peptide_contacts[identity] == 0):
            raise ValueError(f"{identity}: no-interface flag differs from contacts")
    totals = {quantity: int(summary[quantity].sum()) for quantity in
              [*fields.values(), "tcr_peptide_atom_contacts"]}
    errors = json.loads((directory / "errors.json").read_text()) if (directory / "errors.json").exists() else []
    excluded = json.loads((directory / "excluded.json").read_text()) if (directory / "excluded.json").exists() else []
    return dict(structures=summary.height, reference_groups=summary["reference"].n_unique(),
                status="partial" if errors else "complete",
                totals=totals, export_errors=errors, excluded=excluded, row_and_incidence_checks="passed",
                zero_tcr_peptide_interfaces=int(summary["no_tcr_peptide_interface"].sum()),
                by_reference=summary.group_by("reference", maintain_order=True).agg(
                    pl.len().alias("structures"), pl.col("tcr_peptide_atom_contacts").sum(),
                    pl.col("no_tcr_peptide_interface").sum().alias("zero_tcr_peptide_interfaces")).to_dicts(),
                interpretation="Observed molecular graph counts; no binding or pose-accuracy claim")


def summary_command(directory: Path):
    """Summarize and check complementarity exports without rerunning annotation."""
    import typer

    try:
        report = summarize_exports(directory)
    except (ValueError, FileNotFoundError) as exc:
        raise typer.BadParameter(str(exc)) from exc
    typer.echo(json.dumps(report, indent=2))
