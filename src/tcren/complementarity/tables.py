"""Bounded typed Parquet exports of the owning molecular scene tables."""
from __future__ import annotations

import json
from pathlib import Path

import polars as pl

INT, STRING, FLOAT, BOOL = pl.Int64, pl.String, pl.Float64, pl.Boolean
_SCHEMAS = {
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

    def add(self, scene):
        for name, schema in _SCHEMAS.items():
            for original in scene[name]:
                row = {key: original.get(key) for key in schema}
                if name == "atoms":
                    row.update(zip(("x", "y", "z"), original["xyz"]))
                row["structure_id"] = scene["id"]
                self.rows[name].append(row)
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
            table.write_parquet(self.directory / filename, compression="zstd")
            self.manifest.append(dict(table=name, file=filename, rows=table.height,
                                      schema={k: str(v) for k, v in table.schema.items()}))
        self.rows = {name: [] for name in _SCHEMAS}
        self.count = 0
        self.partition += 1

    def close(self):
        self.flush()
        (self.directory / "manifest.json").write_text(json.dumps(dict(
            schema="tcren.complementarity/2", key=["structure_id", "id"],
            residue_coordinates="atoms joined on (structure_id, ca_atom=id)",
            files=self.manifest), indent=2), encoding="utf-8")
