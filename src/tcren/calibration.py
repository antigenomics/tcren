#!/usr/bin/env python3
"""Class-separated cutoff calibration from fresh CLI contacts and full residue identities.

Keep one 8-A-contacting cohort per class. Generate twenty full-peptide shuffles
once per complex; preserve class-I P2/Pomega or structurally resolved class-II P1/P4/P6/P9.
Reuse the same decoys at every cutoff and potential. Export individual energies;
R owns complex-paired, epitope-cluster bootstrap summaries and publication macros.
This calibrates empirical potentials in sample, not held-out binding predictions.
"""

from __future__ import annotations

import shutil
from types import SimpleNamespace
from importlib.resources import files
import typer
import hashlib
import json
from pathlib import Path
import sys
import subprocess
from tempfile import TemporaryDirectory

import tcren
import numpy as np
import polars as pl
from tcren import analysis
from tcren.potential import Potential
from tcren.refine.anchors import anchor_indices, _core_anchor_score


GRID = [3.5 + 0.5 * i for i in range(10)]
SEED, N_DECOY = 20260725, 20
AA = "ACDEFGHIKLMNPQRSTVWY"
KEY = ["pdb.id", "chain.id.to", "residue.index.to"]


RES_KEY = KEY

def positional_anchors(tp: pl.DataFrame, metadata: pl.DataFrame) -> pl.DataFrame:
    """Map class-specific anchors to zero-based indices of the COMPLETE observed peptide.

    ContactMap residue.index is a chain sequence index, never a PDB residue number
    or an offset from the first contacting residue. Class II uses the existing
    TCRen sequence-register predictor; these are inferred, not measured registers.
    """
    from tcren.refine.anchors import anchor_indices
    per = tp.select("pdb.id", "epitope").unique().join(
        metadata.select(pl.col("id").alias("pdb.id"), "mhc_class"),
        on="pdb.id", how="left", validate="1:1")
    if per["mhc_class"].null_count() or not set(per["mhc_class"]) <= {"MHCI", "MHCII"}:
        raise ValueError("Every peptide requires an explicit MHCI/MHCII annotation")
    keep = [{"pdb.id": r["pdb.id"], "residue.index.to": j}
            for r in per.iter_rows(named=True)
            for j in anchor_indices(r["epitope"], r["mhc_class"])]
    return tp.select(*RES_KEY).unique().join(pl.DataFrame(keep),
        on=["pdb.id", "residue.index.to"], how="semi")


def epitope_weights(tp: pl.DataFrame) -> dict[str, float]:
    """1/n per complex, n = number of complexes sharing that epitope."""
    per = tp.select("pdb.id", "epitope").unique()
    n = per.group_by("epitope").agg(pl.len().alias("n"))
    w = per.join(n, on="epitope").with_columns((1.0 / pl.col("n")).alias("w"))
    return dict(zip(w["pdb.id"].to_list(), w["w"].to_list()))


def derive_cli(contacts, output):
    """Apply the exact historical recipe through the owning published CLI."""
    votes = epitope_weights(contacts)
    payload = contacts.select("pdb.id", "residue.aa.from", "residue.aa.to").with_columns(
        contacts["pdb.id"].replace_strict(votes, return_dtype=pl.Float64).alias("vote")
    ).write_csv()
    command = [sys.executable, "-m", "tcren.cli", "derive-potential",
               "-i", "/dev/stdin", "--variant", "classic", "--pseudocount", "1",
               "--weight-column", "vote", "--keep-cys", "-o", str(output)]
    # Polars' path reader requires a seekable file; this is transient CLI transport,
    # deleted immediately, never a retained or reusable computed intermediate.
    with TemporaryDirectory(prefix="tcren-derive-") as tmp:
        source = Path(tmp) / "contacts.csv"
        source.write_text(payload)
        command[command.index("/dev/stdin")] = str(source)
        result = subprocess.run(command, text=True, capture_output=True, timeout=60)
        Path(str(output)+".log").write_text(result.stdout+result.stderr)
        result.check_returncode()
    Path(str(output)+".json").write_text(json.dumps({
        "argv": command, "stdin_sha256": hashlib.sha256(payload.encode()).hexdigest(),
        "weighting": "1/n structures sharing the epitope, within this class/cutoff",
        "drop_cys": False, "pseudocount": 1, "tcren": tcren.__version__,
        "output_sha256": hashlib.sha256(output.read_bytes()).hexdigest(),
    }, indent=2)+"\n")
    return Potential.from_csv(output)


def prepare(contacts, residues, metadata, registers=None, unavailable_chemistry=None):
    """Join on physical residue identity; no contact-derived sequence offsets."""
    pep = residues.filter(pl.col("chain.type") == "PEPTIDE").with_columns(
        pl.col("peptide.pos").cast(pl.Int64)
    )
    seq = (
        pep.sort(["pdb.id", "peptide.pos"])
        .group_by("pdb.id")
        .agg(
            pl.col("residue.aa").str.join("").alias("epitope"),
            pl.col("chain.id").n_unique().alias("n_chains"),
        )
    )
    assert seq.filter(pl.col("n_chains") != 1).is_empty(), "Multiple peptide chains"
    cohort = (
        contacts.select("pdb.id")
        .unique()
        .join(seq, on="pdb.id", validate="1:1")
        .join(
            metadata.select(pl.col("id").alias("pdb.id"), "mhc_class"),
            on="pdb.id",
            how="left",
            validate="1:1",
        )
        .sort("pdb.id")
    )
    assert cohort["mhc_class"].null_count() == 0
    assert set(cohort["mhc_class"]) <= {"MHCI", "MHCII"}
    if unavailable_chemistry is None:
        raise ValueError("Source-verified chemical scope table is required")
    unavailable = seq.join(metadata.select(pl.col("id").alias("pdb.id"), "mhc_class"), on="pdb.id").join(
        unavailable_chemistry.select("pdb.id", "reason"), on="pdb.id")
    cohort = cohort.join(unavailable.select("pdb.id"), on="pdb.id", how="anti")
    assert all(set(s) <= set(AA) for s in cohort["epitope"]), "Nonstandard peptide"
    if registers is None:
        raise ValueError("Structural class-II register export is required")
    best = {r["pdb.id"]: r for r in registers.filter(pl.col("rank") == 1).iter_rows(named=True)}
    assert set(cohort.filter(pl.col("mhc_class") == "MHCII")["pdb.id"]) <= best.keys()
    identity = pep.select(
        "pdb.id",
        pl.col("chain.id").alias("chain.id.to"),
        pl.col("residue.index").alias("residue.index.to"),
        pl.col("peptide.pos").alias("peptide_pos"),
        pl.col("residue.aa").alias("expected_aa"),
    )
    assert identity.select(KEY).is_duplicated().sum() == 0
    keys = [
        "pdb.id",
        "chain.id.from",
        "residue.index.from",
        "chain.id.to",
        "residue.index.to",
        "residue.aa.from",
        "residue.aa.to",
    ]
    cd = (
        contacts.group_by(keys)
        .agg(pl.col("dist").min().alias("dmin"))
        .join(identity, on=KEY, how="left", validate="m:1")
        .join(
            cohort.select("pdb.id", "epitope", "mhc_class"), on="pdb.id", validate="m:1"
        )
    )
    assert cd["peptide_pos"].null_count() == 0
    assert cd.filter(pl.col("residue.aa.to") != pl.col("expected_aa")).is_empty()
    # Export index and anchor assertions explicitly, including noncontacting residues.
    records = []
    for r in cohort.iter_rows(named=True):
        s, cls = r["epitope"], r["mhc_class"]
        reg = best[r["pdb.id"]] if cls == "MHCII" else None
        anchors = tuple(reg[f"P{pos}.peptide_pos"] - 1 for pos in (1, 4, 6, 9)) if reg else anchor_indices(s, cls)
        assert anchors, f"No class-specific anchors: {r['pdb.id']}"
        starts = (
            [_core_anchor_score(s[i : i + 9]) for i in range(len(s) - 8)]
            if cls == "MHCII"
            else []
        )
        for j in range(len(s)):
            records.append(
                {
                    "pdb.id": r["pdb.id"],
                    "peptide_pos": j + 1,
                    "mhc_class": cls,
                    "epitope": s,
                    "residue.aa": s[j],
                    "is_anchor": j in anchors,
                    "register_source": "MHC-only structural correspondence"
                    if cls == "MHCII"
                    else "P2/Pomega convention",
                    "tied_best_registers": sum(v == max(starts) for v in starts)
                    if starts
                    else 1,
                }
            )
    anchor_rows = pl.DataFrame(records).join(
        pep.select(
            "pdb.id",
            "peptide.pos",
            "chain.id",
            "residue.index",
            "residue.index.pdb",
            "insertion_code",
        ),
        left_on=["pdb.id", "peptide_pos"],
        right_on=["pdb.id", "peptide.pos"],
        validate="1:1",
    )
    # Cross-check the shared helper used by the other historical derivation caller.
    expected = anchor_rows.filter(pl.col("is_anchor")).select(
        "pdb.id",
        pl.col("chain.id").alias("chain.id.to"),
        pl.col("residue.index").alias("residue.index.to"),
    )
    expected = cd.select(KEY).unique().join(expected, on=KEY, how="semi").sort(KEY)
    class_i = cohort.filter(pl.col("mhc_class") == "MHCI").select("pdb.id")
    assert positional_anchors(cd.join(class_i, on="pdb.id", how="semi"), metadata).sort(KEY).equals(
        expected.join(class_i, on="pdb.id", how="semi").sort(KEY))
    return cd.sort(keys), cohort, anchor_rows, unavailable


def command(
    contacts: Path = typer.Option(..., exists=True, dir_okay=False, help="Fresh 8-A CLI contacts with residue/register sidecars"),
    metadata: Path = typer.Option(..., exists=True, dir_okay=False),
    chemistry: Path = typer.Option(..., exists=True, dir_okay=False, help="Explicit source-verified chemical exclusions"),
    output: Path = typer.Option(..., "--output", "-o"),
    comparison_energies: Path | None = typer.Option(None, exists=True, dir_okay=False),
    previous_pooled: Path = typer.Option(..., exists=True, dir_okay=False),
    previous_replay: Path = typer.Option(..., exists=True, dir_okay=False),
):
    """Reproduce the fixed 20-decoy/class-separated cutoff protocol; requires Rscript/data.table.

    This explicit historical calibration uses the established 10 cutoffs and 2000
    epitope-cluster draws. It does not change any scoring default or model resource.
    """
    rscript = shutil.which("Rscript")
    if rscript is None:
        raise ValueError("cutoff-calibration requires optional R with Rscript and data.table on PATH")
    subprocess.run([rscript, "-e", "stopifnot(requireNamespace('data.table', quietly=TRUE))"], check=True, timeout=30)
    a = SimpleNamespace(contacts=contacts, metadata=metadata, output=output)
    wheel = {"version": tcren.__version__, "producer_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    a.output.mkdir(parents=True, exist_ok=False)
    (a.output / "potentials").mkdir()
    (a.output / "wheel.json").write_text(json.dumps(wheel, indent=2)+"\n")
    c = pl.read_csv(a.contacts, separator="\t")
    r = pl.read_csv(str(a.contacts) + ".residues.tsv", separator="\t")
    md = pl.read_csv(a.metadata, separator="\t")
    reg_path = Path(str(a.contacts) + ".registers.tsv")
    registers = pl.read_csv(reg_path, separator="\t")
    chemistry_path = chemistry
    cd, cohort, anchors, unavailable = prepare(c, r, md, registers,
        pl.read_csv(chemistry_path, separator="\t"))
    registers.write_csv(a.output / "structural_registers.tsv", separator="\t")
    anchor_map = {pdb: set(rows.filter(pl.col("is_anchor"))["peptide_pos"] - 1)
                  for (pdb,), rows in anchors.group_by("pdb.id")}
    unavailable.write_csv(a.output / "unavailable_register.tsv", separator="\t")
    anchors.write_csv(a.output / "peptide_anchors.tsv", separator="\t")
    cohort.write_csv(a.output / "cohort.tsv", separator="\t")
    all_ids = r.select("pdb.id").unique()
    all_ids.join(cohort.select("pdb.id"), on="pdb.id", how="anti").with_columns(
        pl.lit(
            "No exported TCR:peptide contact within 8 A; inspect CLI accounting"
        ).alias("reason")
    ).write_csv(a.output / "not_in_contact_cohort.tsv", separator="\t")
    rows, drift, matrices, membership, oracles = [], [], [], [], []
    for cls in ["MHCI", "MHCII"]:
        members = cohort.filter(pl.col("mhc_class") == cls)
        data = cd.filter(pl.col("mhc_class") == cls)
        ids = members["pdb.id"].to_list()
        n = len(ids)
        model = {v: i for i, v in enumerate(ids)}
        rng = np.random.default_rng(SEED)
        mutations = {}
        # Keep historical RNG order for unchanged antigens; newly corrected antigens follow.
        previous = dict(zip(md["id"], md["epitope"]))
        mutation_order = sorted(members.iter_rows(named=True),
                                key=lambda v: (previous[v["pdb.id"]] != v["epitope"], v["pdb.id"]))
        for v in mutation_order:
            sequence = np.array([AA.index(x) for x in v["epitope"]])
            free = np.array(
                [
                    j
                    for j in range(len(sequence))
                    if j not in anchor_map[v["pdb.id"]]
                ]
            )
            seqs = np.tile(sequence, (N_DECOY + 1, 1))
            for rep in range(1, N_DECOY + 1):
                seqs[rep, free] = rng.permutation(sequence[free])
            mutations[v["pdb.id"]] = seqs
            membership.append(
                {
                    "pdb.id": v["pdb.id"],
                    "mhc_class": cls,
                    "epitope": v["epitope"],
                    "n_free": len(free),
                    "decoy_sequences": ";".join(
                        "".join(AA[j] for j in row) for row in seqs[1:]
                    ),
                }
            )
        fixed = derive_cli(data.filter(pl.col("dmin") <= 5),
                           a.output / "potentials" / (cls + "_fixed5.csv"))
        ref = None
        for cutoff in GRID:
            tp = data.filter(pl.col("dmin") <= cutoff)
            derived = derive_cli(tp, a.output / "potentials" / f"{cls}_{cutoff:.1f}.csv")
            M, rr, cc = analysis.potential_matrix(derived)
            matrix = np.array(
                [[M[list(rr).index(x), list(cc).index(y)] for y in AA] for x in AA]
            )
            matrices.append((cls, cutoff, matrix))
            if cutoff == 5:
                ref = matrix
            mi = np.array([model[v] for v in tp["pdb.id"]], dtype=int)
            ri = np.array([AA.index(v) for v in tp["residue.aa.from"]], dtype=int)
            pj = np.stack(
                [
                    mutations[v][:, int(pos) - 1]
                    for v, pos in tp.select("pdb.id", "peptide_pos").iter_rows()
                ]
            )
            # Fixed AA-pair counts make identical contact compositions exactly tied.
            bins = (
                (mi[:, None] * (N_DECOY + 1) + np.arange(N_DECOY + 1)) * 400
                + ri[:, None] * 20
                + pj
            )
            counts = np.bincount(
                bins.ravel(), minlength=n * (N_DECOY + 1) * 400
            ).reshape(n, N_DECOY + 1, 400)
            for name, potential in [
                ("Fixed 5 Å", fixed),
                ("Re-derived", derived),
                ("Published", tcren.potential.tcren()),
            ]:
                F = analysis.potential_long(potential)
                lookup = {
                    (x, y): z
                    for x, y, z in F.select(
                        "residue.aa.from", "residue.aa.to", "value"
                    ).iter_rows()
                }
                mat = np.array([[lookup.get((x, y), np.nan) for y in AA] for x in AA])
                # Match the library's _interface_energy: absent entries contribute no term.
                defined = np.isfinite(mat)
                mat = np.nan_to_num(mat)
                from tcren.energetics.scoring import _interface_energy

                oracle = _interface_energy(
                    tp.filter(pl.col("pdb.id") == ids[0]), potential
                )
                energies = counts @ mat.ravel()
                assert abs(energies[0, 0] - oracle) < 1e-10
                oracles.append(
                    {
                        "mhc_class": cls,
                        "cutoff": cutoff,
                        "potential": name,
                        "pdb.id": ids[0],
                        "library_energy": oracle,
                        "table_energy": energies[0, 0],
                    }
                )
                missing = counts @ (~defined).ravel()
                for i, pdb in enumerate(ids):
                    for rep in range(N_DECOY + 1):
                        rows.append(
                            {
                                "pdb.id": pdb,
                                "epitope": members["epitope"][i],
                                "mhc_class": cls,
                                "cutoff": cutoff,
                                "potential": name,
                                "replicate": rep,
                                "energy": energies[i, rep],
                                "n_contacts": int(np.sum(mi == i)),
                                "n_unparameterized_contacts": int(missing[i, rep]),
                            }
                        )
        for mc, cutoff, matrix in matrices:
            if mc != cls:
                continue
            drift.append(
                {
                    "mhc_class": cls,
                    "cutoff": cutoff,
                    "matrix_R": float(np.corrcoef(matrix.ravel(), ref.ravel())[0, 1]),
                    "matrix_MAD": float(np.mean(np.abs(matrix - ref))),
                }
            )
    pl.DataFrame(rows).write_csv(a.output / "energies.tsv", separator="\t")
    pl.DataFrame(oracles).write_csv(
        a.output / "library_correspondence.tsv", separator="\t"
    )
    pl.DataFrame(drift).write_csv(a.output / "matrix_drift.tsv", separator="\t")
    pl.DataFrame(membership).write_csv(a.output / "decoys.tsv", separator="\t")
    (a.output / "inputs.json").write_text(
        json.dumps(
            {
                "seed": SEED,
                "n_decoys": N_DECOY,
                "tcren": tcren.__version__,
                "inputs": {
                    str(v): hashlib.sha256(v.read_bytes()).hexdigest()
                    for v in [
                        a.contacts,
                        Path(str(a.contacts) + ".residues.tsv"),
                        a.metadata,
                        reg_path,
                        chemistry_path,
                        previous_pooled,
                        previous_replay,
                    ]
                },
                "argv": sys.argv,
            },
            indent=2,
        )
        + "\n"
    )
    print(
        f"Exported {len(rows)} native/decoy energies, {len(cohort)} complexes; classes separate."
    )

    script = files("tcren.data").joinpath("cutoff_summary.R")
    if comparison_energies is not None:
        comparison_script = Path(str(files("tcren.data").joinpath("calibration_compare.R")))
        subprocess.run([rscript, str(comparison_script), str(comparison_energies),
                        str(output / "energies.tsv"), str(output)], check=True)
        (output / "comparison_sources.json").write_text(json.dumps({
            "inputs": {str(comparison_energies): hashlib.sha256(comparison_energies.read_bytes()).hexdigest()},
            "resource_sha256": hashlib.sha256(comparison_script.read_bytes()).hexdigest()}, indent=2)+"\n")
    command_line = [rscript, str(script), str(output), str(previous_pooled), str(previous_replay)]
    with (output / "summary.log").open("w") as log:
        subprocess.run(command_line, stdout=log, stderr=subprocess.STDOUT, check=True, timeout=600)
    (output / "summary_sources.json").write_text(json.dumps({"argv": command_line,
        "R_script_sha256": hashlib.sha256(script.read_bytes()).hexdigest(),
        "comparison_inputs": {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in [previous_pooled, previous_replay]},
        "seed": SEED, "bootstrap_draws": 2000, "optional_runtime": "Rscript + data.table"}, indent=2)+"\n")
