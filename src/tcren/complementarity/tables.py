"""Bounded typed Parquet exports of the owning molecular scene tables."""
from __future__ import annotations

import json
import gzip
from pathlib import Path

import polars as pl

INT, STRING, FLOAT, BOOL = pl.Int64, pl.String, pl.Float64, pl.Boolean
_SCHEMAS = {
    "structures": dict(reference=STRING, task=STRING, origin=STRING, origin_component=STRING,
                       rmsd=FLOAT, anchors=INT, atoms=INT, residues=INT, covalent_bonds=INT,
                       atom_contacts=INT, residue_contacts=INT, tcr_peptide_atom_contacts=INT,
                       no_tcr_peptide_interface=BOOL, warnings=INT),
    "atoms": dict(id=INT, key=STRING, ordinal=INT, residue=INT, name=STRING, element=STRING, altloc=STRING,
                  backbone=BOOL, x=FLOAT, y=FLOAT, z=FLOAT, covalent_degree=INT, contact_count=INT,
                  contact_degree=INT, partner_residue_count=INT, partner_chain_count=INT),
    "residues": dict(id=INT, key=STRING, chain=STRING, role=STRING, region=STRING, seq_index=INT, number=INT,
                     insertion=STRING, aa=STRING, resname=STRING, ca_atom=INT, atom_count=INT,
                     contacting_atom_count=INT, contact_count=INT, contact_degree=INT,
                     covalent_degree=INT, partner_chain_count=INT),
    "bonds": dict(source=INT, target=INT, kind=STRING, distance=FLOAT, method=STRING),
    "atom_contacts": dict(source=INT, target=INT, distance=FLOAT),
    "residue_contacts": dict(source=INT, target=INT, contact_count=INT),
}


def structure_record(scene):
    """Extract the owning scene's flat summary without rebuilding its graph."""
    values = {**scene["frame"], **scene["summary"], "task": scene.get("task"),
              "warnings": len(scene["warnings"])}
    return {"structure_id": scene["id"], **{key: values.get(key) for key in _SCHEMAS["structures"]}}


def write_table_partition(table, directory, filename):
    """Write one typed table in the two supported compressed serializations."""
    directory = Path(directory)
    table.write_parquet(directory / filename, compression="zstd")
    tsv = str(Path(filename).with_suffix(".tsv.gz"))
    with gzip.GzipFile(filename=str(directory / tsv), mode="wb", mtime=0) as stream:
        table.write_csv(stream, separator="\t", quote_style="non_numeric")
    return tsv


def write_summaries(records, root):
    """Write JSON and flat typed summaries, returning the owning typed table."""
    root = Path(root)
    summary = pl.DataFrame(records, schema={"structure_id": STRING, **_SCHEMAS["structures"]})
    write_table_partition(summary, root, "summary.parquet")
    totals = {name: int(summary[name].sum()) for name in
              ("atoms", "residues", "covalent_bonds", "atom_contacts", "residue_contacts",
               "tcr_peptide_atom_contacts", "warnings")}
    (root / "summary.json").write_text(json.dumps(dict(
        structures=summary.height, totals=totals, records=summary.to_dicts()),
        allow_nan=False, indent=2), encoding="utf-8")
    return summary


class SceneTables:
    """Serialize at most eight scenes per partition, including typed empty tables.

    IDs are local to a structure: every table's composite namespace includes
    ``structure_id``. Residue coordinates resolve ``ca_atom`` in the atoms table.
    No contacts, degrees or coordinates are recalculated here.
    """

    def __init__(self, directory: Path):
        self.directory = directory
        directory.mkdir()
        self.rows = {name: [] for name in _SCHEMAS}
        self.count = self.partition = 0
        self.manifest = []
        self.summaries = []

    def add(self, scene):
        for name, schema in _SCHEMAS.items():
            if name == "structures":
                source = [structure_record(scene)]
            else:
                source = scene[name]
            for original in source:
                row = {key: original.get(key) for key in schema}
                if name == "atoms":
                    row.update(zip(("x", "y", "z"), original["xyz"]))
                row["structure_id"] = scene["id"]
                self.rows[name].append(row)
                if name == "structures":
                    self.summaries.append(row)
        self.count += 1
        if self.count == 8:
            self.flush()

    def flush(self):
        if not self.count:
            return
        for name, schema in _SCHEMAS.items():
            folder = self.directory / name
            folder.mkdir(exist_ok=True)
            filename = f"{name}/part-{self.partition:05d}.parquet"
            table = pl.DataFrame(self.rows[name], schema={"structure_id": STRING, **schema})
            tsv = write_table_partition(table, self.directory, filename)
            key = (["structure_id"] if name == "structures" else
                   ["structure_id", "id"] if name in {"atoms", "residues"} else
                   ["structure_id", "source", "target"])
            self.manifest.append(dict(table=name, file=filename, rows=table.height,
                                      tsv_file=tsv, key=key,
                                      schema={k: str(v) for k, v in table.schema.items()}))
        self.rows = {name: [] for name in _SCHEMAS}
        self.count = 0
        self.partition += 1

    def close(self):
        self.flush()
        (self.directory / "manifest.json").write_text(json.dumps(dict(
            schema="tcren.complementarity/2",
            residue_coordinates="atoms joined on (structure_id, ca_atom=id)",
            files=self.manifest), indent=2), encoding="utf-8")
        write_summaries(self.summaries, self.directory.parent)
