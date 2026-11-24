#!/usr/bin/env python3
"""Build the frozen native reference from verified saved analysis tables."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path

import tcren  # noqa: F401 -- initialize Polars limits before importing it

import numpy as np
import polars as pl

from tcren.provenance import read
from tcren.structure.io import structure_paths, structure_stem

ARMS = ("native", "real", "mock", "tcrvdb", "cpl")
GEOMETRY = (
    "gap_mean", "gap_sd", "gap_vol", "gap_index", "gap_asym", "gap_depth", "gap_height",
    "d_h", "d_charge", "d_phobic", "interlock", "interlock_frac", "shape_r", "charge_r",
    "phobic_r", "charge_product", "phobic_product", "coverage",
)


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def summary(values: np.ndarray) -> dict:
    values = values[np.isfinite(values)]
    if not values.size:
        raise ValueError("native summary has no finite observations")
    centre = float(np.median(values))
    spread = float(np.median(np.abs(values - centre)) * 1.4826)
    return {"n": len(values), "centre": centre, "spread": spread or 1.0}


def observe(path: Path, alphabet: str) -> pl.LazyFrame:
    frame = pl.scan_parquet(path)
    if alphabet == "residue":
        frame = frame.filter(pl.col("contacted"))
        endpoints = ("aa.from", "aa.to")
    else:
        endpoints = ("atom.from", "atom.to")
    return frame.select(
        "pdb.id", "interface", pl.concat_str(endpoints, separator="|").alias("pair"),
    )


def independence_prior(cells: list[str], weight: np.ndarray) -> np.ndarray:
    # Factorize endpoint marginals within each interface over the frozen support.
    parts = [key.rpartition("/") for key in cells]
    interfaces = np.asarray([part[0] for part in parts])
    prior = np.zeros(len(cells))
    for interface in np.unique(interfaces):
        indices = np.flatnonzero(interfaces == interface)
        w = weight[indices]
        left, right = zip(*(parts[i][2].split("|") for i in indices))
        _, left_idx = np.unique(left, return_inverse=True)
        _, right_idx = np.unique(right, return_inverse=True)
        left_weight = np.bincount(left_idx, weights=w)[left_idx]
        right_weight = np.bincount(right_idx, weights=w)[right_idx]
        product = left_weight * right_weight
        prior[indices] = product / product.sum() * w.sum()
    return prior / prior.sum()


def contact_profiles(paths: list[Path], alphabet: str, native_ids: set[str]) -> dict:
    sources = [observe(path, alphabet) for path in paths]
    native = sources[0].collect()
    if set(native["pdb.id"]) - native_ids:
        raise ValueError("native contacts contain structures outside the inventory")
    support = pl.concat([source.select("interface", "pair").unique() for source in sources])
    support = support.unique().collect()
    profiles = {}
    for scope in ("complex", "tcr_peptide"):
        if scope == "complex":
            cell = pl.concat_str("interface", "pair", separator="/").alias("cell")
            observed = native.select("pdb.id", cell)
            vocabulary = support.select(cell)
        else:
            observed = native.filter(pl.col("interface") == scope).select(
                "pdb.id", pl.col("pair").alias("cell"),
            )
            vocabulary = support.filter(pl.col("interface") == scope).select(
                pl.col("pair").alias("cell"),
            )
        cells = sorted(vocabulary["cell"])
        ids = observed["pdb.id"].unique(maintain_order=True).to_list()
        by_id = {value: i for i, value in enumerate(ids)}
        by_cell = {value: i for i, value in enumerate(cells)}
        counts = np.zeros((len(ids), len(cells)), dtype=np.int64)
        for pdb_id, key, n in observed.group_by("pdb.id", "cell").len().iter_rows():
            counts[by_id[pdb_id], by_cell[key]] = n
        totals = counts.sum(axis=1)
        observations = int(totals.sum())
        if observations == 0:
            raise ValueError(f"{alphabet}/{scope} has no native observations")
        structure_mean = (counts / totals[:, None]).mean(axis=0)
        prior = independence_prior(cells, structure_mean * observations + 1e-9)
        q = (observations * structure_mean + 500 * prior) / (observations + 500)
        profiles[f"{alphabet}/{scope}"] = {
            "n_native_structures": len(ids), "n_native_observations": observations,
            "n_native_cells": int(np.count_nonzero(counts.sum(axis=0))),
            "probability": dict(zip(cells, (q / q.sum()).tolist())),
        }
    return profiles


def build(args: argparse.Namespace) -> dict:
    coordinates = structure_paths(args.native_coordinates)
    coordinate_manifest = {structure_stem(path): digest(path) for path in coordinates}
    ids = set(coordinate_manifest)
    if len(ids) != 387:
        raise ValueError("expected all 387 native coordinate files")

    markup = pl.read_csv(args.markup).filter(
        (pl.col("chain.type") == "PEPTIDE") & (pl.col("region.type") == "PEPTIDE")
    )
    peptides = (markup.filter(pl.col("chain.id") == pl.col("chain.id").first().over("pdb.id"))
                .sort("pdb.id", "residue.index").group_by("pdb.id")
                .agg(pl.col("residue.aa").str.join("").alias("peptide")))
    scores = pl.read_csv(args.scoring).select("pdb.id", pl.col("dPhi_peptide_sum").alias("energy"))
    if set(peptides["pdb.id"]) != ids or set(scores["pdb.id"]) != ids:
        raise ValueError("native peptide and energy tables must cover the coordinate inventory")
    energies = peptides.join(scores, on="pdb.id", validate="1:1")
    if not energies["energy"].is_finite().fill_null(False).all():
        raise ValueError("native energies must be finite")
    anchors = energies.sort("pdb.id").group_by("peptide").agg(
        pl.len().alias("n"), pl.col("energy").median().alias("centre"),
        pl.col("energy").std().alias("spread"),
    )
    if anchors.height != 239:
        raise ValueError("unexpected native peptide count")

    geometry = pl.read_csv(args.complementarity, infer_schema_length=5000)
    geometry = geometry.unique(subset="structure.id", keep="first")
    if set(geometry["structure.id"]) != ids:
        raise ValueError("native geometry must cover the coordinate inventory")
    counts = (observe(args.atom[0], "atom").filter(pl.col("interface") == "tcr_peptide")
              .group_by("pdb.id").len().collect())
    provenance = read(args.scoring)
    if not provenance or provenance["reference_mode"] != "structural":
        raise ValueError("need the native structural reference scoring receipt")
    indexes = provenance["searched_indexes"]
    return {
        "format": "tcren-native-reference-1",
        "version": "native387-v1",
        "native_structures": len(ids),
        "cutoff_A": 5.0,
        "energy_reference": {
            "side": "peptide", "mode": "structural", "aa": "A", "contact_weight": "residue",
            "tcr_regions": "all", "potentials": {"tcr_peptide": "tcren2", "peptide_mhc": "mj"},
        },
        "size": summary(counts["len"].to_numpy()),
        "energy_by_peptide": {
            peptide: {"n": n, "centre": centre, "spread": spread}
            for peptide, n, centre, spread in anchors.iter_rows()
        },
        "geometry": {term: summary(geometry[term].cast(pl.Float64).to_numpy())
                     for term in GEOMETRY},
        "contact_profiles": {
            **contact_profiles(args.atom, "atom", ids),
            **contact_profiles(args.residue, "residue", ids),
        },
        "contact_smoothing": {"pseudocount": 500, "marginal_floor": 1e-9,
                              "vocabulary_arms": ARMS, "counts_arm": "native"},
        "provenance": {
            "native_coordinate_manifest_sha256": sha256(
                json.dumps(coordinate_manifest, sort_keys=True).encode(),
            ).hexdigest(),
            "native_tables_sha256": {name: digest(getattr(args, name))
                                     for name in ("markup", "scoring", "complementarity")},
            "support_tables_sha256": {
                "atom": dict(zip(ARMS, map(digest, args.atom))),
                "residue": dict(zip(ARMS, map(digest, args.residue))),
            },
            "scoring_tcren_version": provenance["tcren"],
            "scoring_registry_digest": provenance["registry_digest"],
            "builder_versions": {"numpy": np.__version__, "polars": pl.__version__},
            "searched_index_files_sha256": {
                **{f"arda_{k}": v["files_sha256"] for k, v in indexes["arda"].items()},
                "mhc": indexes["mhc"]["files_sha256"],
            },
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    for flag in ("native-coordinates", "markup", "scoring", "complementarity", "out"):
        parser.add_argument("--" + flag, required=True, type=Path)
    for flag in ("atom", "residue"):
        parser.add_argument("--" + flag, nargs=len(ARMS), required=True, type=Path,
                            metavar="TABLE", help="tables in native, real, mock, TCRvdb, CPL order")
    args = parser.parse_args()
    payload = json.dumps(build(args), sort_keys=True, separators=(",", ":"), allow_nan=False) + "\n"
    raw = payload.encode("utf-8")
    args.out.write_bytes(raw)
    print(sha256(raw).hexdigest())


if __name__ == "__main__":
    main()
