"""Recover external contact frames only through exact float32 atom correspondence.

The deposited antibody contacts use IMGT numbering. Current author numbering is
never assumed equal. Exact chemical identity and deposited float32 coordinates
must identify one selected raw atom; ambiguous/missing matches remain missing.
"""
from __future__ import annotations

import gzip
import json
from pathlib import Path
import time

import numpy as np
import polars as pl
import typer
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.Data.PDBData import protein_letters_3to1

from .observations import _angle, _dihedral

ATOM_KEY = ["resname", "atom", "x", "y", "z"]
RES_KEY = ["chain", "resid", "icode", "resname"]
FRAMES = [prefix + axis for prefix in ("n.", "ca.", "cb.") for axis in "xyz"]


def endpoint_frames(endpoints: pl.DataFrame, raw: pl.DataFrame) -> pl.DataFrame:
    """Map endpoints without tolerance, numbering assumptions or ambiguous copies."""
    unique = raw.filter(~raw.select(ATOM_KEY).is_duplicated())
    frames = raw.group_by(RES_KEY).agg([
        pl.col(axis).filter(pl.col("atom") == name).first().alias(prefix + axis)
        for name, prefix in (("N", "n."), ("CA", "ca."), ("CB", "cb."))
        for axis in "xyz"])
    unique = unique.join(frames, on=RES_KEY, validate="m:1").with_columns(pl.lit(True).alias("matched"))
    return endpoints.join(unique, on=ATOM_KEY, how="left", validate="m:1", maintain_order="left").with_columns(
        pl.col("matched").fill_null(False))


def raw_atoms(path: Path) -> pl.DataFrame:
    """First model, one lexicographic conformer per residue, observed atoms only."""
    with gzip.open(path, "rt") as handle:
        data = MMCIF2Dict(handle)
    fields = {"resname": "auth_comp_id", "atom": "auth_atom_id", "chain": "auth_asym_id",
              "resid": "auth_seq_id", "icode": "pdbx_PDB_ins_code", "altloc": "label_alt_id",
              "element": "type_symbol", "model": "pdbx_PDB_model_num",
              **dict(zip("xyz", ("Cartn_x", "Cartn_y", "Cartn_z")))}
    table = pl.DataFrame({key: data["_atom_site." + value] for key, value in fields.items()})
    table = table.filter((pl.col("model") == table["model"][0]) & ~pl.col("element").is_in(["H", "D"]))
    table = table.with_columns(pl.col("x", "y", "z").cast(pl.Float32),
                              pl.col("icode", "altloc").replace({".": "", "?": ""}))
    labels = table.group_by(RES_KEY).agg(pl.col("altloc").filter(pl.col("altloc") != "").min().fill_null("").alias("selected"))
    table = table.join(labels, on=RES_KEY, validate="m:1").filter(
        (pl.col("altloc") == "") | (pl.col("altloc") == pl.col("selected")))
    return table.sort("altloc").unique(subset=RES_KEY + ["atom"], keep="last", maintain_order=True).drop("selected", "model")


def add_orientations(table: pl.DataFrame) -> pl.DataFrame:
    """Use the same observed-frame geometry as the native contact export."""
    def xyz(prefix, side):
        return table.select([prefix + a + "." + side for a in "xyz"]).to_numpy().astype(float)
    ca, cb, n = [{side: xyz(prefix, side) for side in ("from", "to")}
                 for prefix in ("ca.", "cb.", "n.")]
    values = {"omega_rad": _dihedral(ca["from"], cb["from"], cb["to"], ca["to"]),
              "cb_distance_A": np.linalg.norm(cb["to"] - cb["from"], axis=1)}
    for side, other in (("from", "to"), ("to", "from")):
        values["theta." + side + "_rad"] = _dihedral(n[side], ca[side], cb[side], cb[other])
        values["phi." + side + "_rad"] = _angle(ca[side], cb[side], cb[other])
    return table.with_columns([pl.Series(k, v).fill_nan(None) for k, v in values.items()])


def external_geometry(root: Path, output: Path, excluded: set[str], *, limit: int | None = None):
    """Audit/recover raw geometry one structure at a time, with compact final output.

    Only projected immutable contact columns are held (~856k rows); one mmCIF
    dictionary and one structure's joins coexist. No sequence annotation/search.
    """
    output.mkdir(exist_ok=False, parents=True)
    keys = ["PDB"] + [f"{s} {field}" for s in ("Ab", "Ag") for field in ("chain ID", "resi ID", "atom")]
    columns = keys + [f"{s} {field}" for s in ("Ab", "Ag") for field in ("resi aa", "atom coord")] + ["Contact distance"]
    table = pl.read_parquet(root / "Output_contact.parquet", columns=columns).filter(~pl.col("PDB").is_in(sorted(excluded)))
    # Exact duplicates are known in the deposit. Conflicting observations fail.
    distinct = table.unique(maintain_order=True)
    if distinct.select(keys).is_duplicated().any():
        raise ValueError("conflicting deposited physical atom-pair observations")
    table = distinct.sort("PDB")
    pdbs = table["PDB"].unique(maintain_order=True).to_list()
    if limit is not None:
        pdbs = pdbs[:limit]
    counts, start = [], time.monotonic()
    with gzip.open(output / "external_geometry.tsv.gz", "wt") as stream:
        first = True
        for index, pdb in enumerate(pdbs):
            block = table.filter(pl.col("PDB") == pdb)
            raw = raw_atoms(root / "rcsb-20261018" / (pdb + ".cif.gz"))
            parts = []
            for label, side in (("Ab", "from"), ("Ag", "to")):
                coords = block[label + " atom coord"]
                points = block.select(pl.col(label + " resi aa").alias("resname"), pl.col(label + " atom").alias("atom"))
                points = points.with_columns([coords.list.get(i).alias(axis) for i, axis in enumerate("xyz")])
                mapped = endpoint_frames(points, raw)
                parts.append(mapped.rename({c: c + "." + side for c in mapped.columns}))
            result = pl.concat([block.select(pl.col("PDB").alias("pdb.id"), pl.col("Contact distance").alias("dist"),
                *[pl.col(s + " " + f).alias(prefix + "." + side) for s, side in (("Ab", "from"), ("Ag", "to"))
                  for f, prefix in (("chain ID", "source_chain"), ("resi ID", "source_resid"))]), *parts], how="horizontal")
            result = result.with_columns([
                (pl.col("resname." + side).replace_strict(protein_letters_3to1, default=None) + ":" + pl.col("atom." + side)).alias("chem." + side)
                for side in ("from", "to")])
            result = add_orientations(result)
            stream.write(result.write_csv(separator="\t", include_header=first))
            first = False
            counts.append({"pdb.id": pdb, "atom_pairs": result.height,
                "matched_both": result.select((pl.col("matched.from") & pl.col("matched.to")).sum()).item(),
                "both_frames": result.select(pl.all_horizontal(pl.col([p + "." + s for s in ("from", "to") for p in FRAMES]).is_not_null()).sum()).item(),
                "residue_pairs": block.select(keys[:1] + [k for k in keys[1:] if not k.endswith("atom")]).unique().height})
            if index % 50 == 0:
                print(json.dumps({"structures": index + 1, "total": len(pdbs), "seconds": time.monotonic() - start}), flush=True)
    pl.DataFrame(counts).write_csv(output / "correspondence.tsv", separator="\t")
    report = {"structures": len(counts), "excluded_pdb_ids": sorted(excluded), "diagnostic_limit": limit,
        "deposited_atom_rows_after_exclusion": distinct.height, "elapsed_s": time.monotonic()-start,
        "atom_pairs": sum(r["atom_pairs"] for r in counts), "matched_both": sum(r["matched_both"] for r in counts),
        "interpretation": "Exact float32 chemical/coordinate correspondence; unmatched or ambiguous endpoints have no recovered frame"}
    (output / "summary.json").write_text(json.dumps(report, indent=2) + "\n")
    return report


def geometry_command(
    source: Path = typer.Option(..., "--source", exists=True, file_okay=False),
    out: Path = typer.Option(..., "-o", "--out"),
    exclude: Path = typer.Option(..., "--exclude", exists=True, dir_okay=False, help="Native PDB IDs, one per line"),
    limit: int | None = typer.Option(None, "--limit", min=1, help="Diagnostic subset only"),
):
    """Verify deposited external contacts against selected current raw structures."""
    excluded = set(exclude.read_text().lower().split())
    print(json.dumps(external_geometry(source, out, excluded, limit=limit)))
