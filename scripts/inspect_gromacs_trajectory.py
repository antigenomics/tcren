"""Export a coherent continuous protein trajectory and inspect it with tcren."""

import os

for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "POLARS_MAX_THREADS"):
    os.environ[key] = "1"

import argparse
import hashlib
import json
import subprocess
from pathlib import Path

import polars as pl

from tcren.validation.md import inspect_md_trajectory


def sha256(path):
    result = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            result.update(block)
    return result.hexdigest()


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--gmx", type=Path, required=True)
    p.add_argument("--tpr", type=Path, required=True)
    p.add_argument("--trajectory", type=Path, required=True)
    p.add_argument(
        "--reference", type=Path, help="RMSD reference; defaults to the compact first saved frame"
    )
    p.add_argument(
        "--potential",
        type=Path,
        required=True,
        help="Two-column time (ps) / whole-system potential (kJ/mol) XVG",
    )
    p.add_argument("--components", required=True, help="JSON component name -> chain IDs")
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--cutoff", type=float, default=5.0)
    p.add_argument("--docking", action="store_true")
    p.add_argument("--threads", type=int, default=1)
    args = p.parse_args()
    output = args.output.resolve()
    output.mkdir(exist_ok=False)
    commands = []
    whole = output / "protein-whole.xtc"
    anchor = output / "protein-first-compact.pdb"
    continuous = output / "protein-continuous.pdb"

    def convert(name, options, selection):
        command = [str(args.gmx.resolve()), "trjconv", *map(str, options)]
        commands.append({"command": command, "selection": selection})
        with (output / (name + ".log")).open("w") as log:
            subprocess.run(
                command,
                input=selection,
                text=True,
                stdout=log,
                stderr=subprocess.STDOUT,
                check=True,
            )

    convert(
        "whole",
        ["-s", args.tpr.resolve(), "-f", args.trajectory.resolve(), "-pbc", "whole", "-o", whole],
        "Protein\n",
    )
    convert(
        "compact-anchor",
        ["-s", args.tpr.resolve(), "-f", whole, "-pbc", "cluster", "-dump", "0", "-o", anchor],
        "Protein\nProtein\n",
    )
    convert(
        "continuous", ["-s", anchor, "-f", whole, "-pbc", "nojump", "-o", continuous], "Protein\n"
    )
    reference = args.reference if args.reference is not None else anchor
    result = inspect_md_trajectory(
        continuous,
        reference,
        json.loads(args.components),
        args.potential,
        contact_cutoff=args.cutoff,
        docking=args.docking,
        threads=args.threads,
    )
    for name in ["rmsd", "contact_counts", "contact_occupancy", "docking"]:
        if result[name]:
            pl.DataFrame(result[name]).write_csv(output / (name + ".tsv"), separator="\t")
    (output / "minimum-potential-frame.pdb").write_text(result["selected_frame_pdb"])
    (output / "first-frame.pdb").write_text(result["first_frame_pdb"])
    (output / "last-frame.pdb").write_text(result["last_frame_pdb"])
    (output / "saved-frame-times.json").write_text(
        json.dumps(result["saved_frame_times_ps"]) + "\n"
    )
    (output / "minimum-potential-frame.json").write_text(
        json.dumps(result["minimum_potential"], indent=2) + "\n"
    )
    sources = [args.tpr, args.trajectory, reference, args.potential]
    provenance = {
        "state": "complete",
        "sources": {str(p.resolve()): sha256(p) for p in sources},
        "commands": commands,
        "components": json.loads(args.components),
        "cutoff_A": args.cutoff,
        "docking": args.docking,
        "threads": args.threads,
        "occupancy_scope": "fraction of saved frames; observed residue pairs only",
        "pbc": "molecules whole, compact first-frame anchor, continuous nojump; later frames are not clustered",
    }
    (output / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    print(
        json.dumps(
            {
                "output": str(output),
                "rmsd_rows": len(result["rmsd"]),
                "contact_pairs": len(result["contact_occupancy"]),
                "docking_frames": len(result["docking"]),
            }
        )
    )
