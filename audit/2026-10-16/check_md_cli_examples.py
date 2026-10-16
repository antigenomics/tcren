"""Verify installed MD CLI against APIs on public class-separated coordinate fixtures.

Two rigidly translated crystal frames and synthetic energies exercise software only,
not molecular dynamics or a scientific control. Requires the owning repo environment.
Run with --structures CANONICAL_DIR --output NEW_DIR. No result reuse/resume.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time

import tcren
import numpy as np
import polars as pl
import psutil
from threadpoolctl import threadpool_info

from tcren.structure import parse_structure
from tcren.structure.io import write_pdb
from tcren.validation.md import inspect_md_trajectory
from tcren.validation.structures import close_atom_pairs


def save(path, value):
    path.write_text(json.dumps(value, indent=2) + "\n")


def run(command, root):
    env = dict(os.environ, SLURM_CPUS_PER_TASK="2")
    for key in ("OMP_NUM_THREADS", "OMP_THREAD_LIMIT", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                "POLARS_MAX_THREADS", "POLARS_ASYNC_THREAD_COUNT", "POLARS_MAX_BLOCKING_THREAD_COUNT",
                "RAYON_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS", "PYTHONPATH"):
        env.pop(key, None)
    seen = {}
    peak_rss = peak_threads = missed = 0
    started = time.monotonic()
    with (root / f"command-{len(commands)}.log").open("w") as log:
        proc = subprocess.Popen(command, env=env, stdout=log, stderr=subprocess.STDOUT)
        parent = psutil.Process(proc.pid)
        while proc.poll() is None:
            try:
                tree = [parent, *parent.children(recursive=True)]
            except psutil.NoSuchProcess:
                tree = []
            rss = threads = 0
            for item in tree:
                try:
                    with item.oneshot():
                        n = item.num_threads()
                        rss += item.memory_info().rss
                        threads += n
                        cmd = item.cmdline()
                        cpu = item.cpu_times().user + item.cpu_times().system
                    old = seen.get(str(item.pid), {})
                    seen[str(item.pid)] = {"role": "mmseqs" if any("mmseqs" in v for v in cmd) else "python",
                                           "max_threads": max(old.get("max_threads", 0), n),
                                           "cpu_seconds": cpu,
                                           "search_budget": [cmd[i+1] for i, v in enumerate(cmd[:-1]) if v == "--threads"]}
                except (psutil.NoSuchProcess, psutil.AccessDenied):
                    missed += 1
            peak_rss, peak_threads = max(peak_rss, rss), max(peak_threads, threads)
            if time.monotonic() - started > 90 or rss > 4 * 1024**3:
                for child in parent.children(recursive=True):
                    child.kill()
                proc.kill()
                proc.wait()
                raise RuntimeError("CLI example exceeded 90 seconds/4 GiB")
            time.sleep(0.02)
    record = {"command": command, "wall_seconds": time.monotonic() - started,
              "peak_tree_rss_bytes": peak_rss, "peak_live_threads": peak_threads,
              "processes": seen, "unreadable_process_samples": missed, "exit_code": proc.returncode}
    commands.append(record)
    save(root / "commands-resources.json", commands)
    if proc.returncode:
        raise RuntimeError(f"Command failed: {command}; inspect command log")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--structures", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    root = args.output.resolve()
    root.mkdir(parents=True, exist_ok=False)
    commands, checks = [], []
    for pdb_id, mhc_class in [("1oga", "I"), ("1fyt", "II")]:
        folder = root / pdb_id
        folder.mkdir()
        source = args.structures / f"{pdb_id}.pdb.gz"
        reference = folder / "reference.pdb"
        structure = parse_structure(source, keep_hydrogens=False)
        # The native PDB parser stores float32 coordinates: translated read-back
        # is not bit-exact geometry. Bound that fixture error from coordinate scale;
        # CLI/API comparisons below remain exact, without a numerical tolerance.
        coordinate_scale = max(np.linalg.norm(a.coord) for c in structure.chains for r in c.residues for a in r.atoms)
        coordinate_bound = 4 * np.finfo(np.float32).eps * coordinate_scale
        write_pdb(structure, reference, keep_hydrogens=False)
        translated = folder / "translated.pdb"
        write_pdb(structure, translated, transform=lambda x: x + [1, 2, 3], keep_hydrogens=False)
        trajectory = folder / "trajectory.pdb"
        trajectory.write_text("".join(
            f"TITLE     t= {time_ps}\nMODEL        {index}\n" +
            "\n".join(line for line in path.read_text().splitlines() if line != "END") + "\nENDMDL\n"
            for index, (time_ps, path) in enumerate([(0, reference), (10, translated)], 1)))
        potential = folder / "potential.xvg"
        potential.write_text("0 -10\n10 -20\n")
        components = {"TCR": ["A", "B"], "peptide": ["C"], "MHC": ["D", "E"] if mhc_class == "II" else ["D"]}
        component_file = folder / "components.json"
        save(component_file, components)
        base = [sys.executable, "-m", "tcren.cli", "md", "inspect", "--trajectory", str(trajectory),
                "--reference", str(reference), "--potential", str(potential),
                "--components", str(component_file), "--mhc-class", mhc_class, "--pbc-whole", "--docking"]
        budgets = [1, 2, 0] if mhc_class == "I" else [1]
        for budget in budgets:
            out = folder / f"inspect-{budget}"
            run([*base, "-t", str(budget), "-o", str(out)], root)
            expected = inspect_md_trajectory(trajectory, reference, components, potential,
                                             contact_cutoff=5, docking=True, threads=1)
            for name in ["rmsd", "contact_counts", "contact_occupancy", "docking"]:
                actual = pl.read_csv(out / f"{name}.tsv", separator="\t").to_dicts()
                assert actual == expected[name], name
            assert json.loads((out / "minimum-potential-frame.json").read_text()) == expected["minimum_potential"]
            assert (out / "minimum-potential-frame.pdb").read_text() == expected["selected_frame_pdb"]
            assert max(row["rmsd_A"] for row in expected["rmsd"]) < coordinate_bound
            checks.append({"id": pdb_id, "mhc_class": mhc_class, "budget": budget,
                           "rmsd_rows": len(expected["rmsd"]), "contact_count_rows": len(expected["contact_counts"]),
                           "occupancy_rows": len(expected["contact_occupancy"]), "docking_rows": len(expected["docking"]),
                           "maximum_rmsd_A": max(row["rmsd_A"] for row in expected["rmsd"]),
                           "float32_coordinate_bound_A": float(coordinate_bound),
                           "source_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
                           "components": components, "exact_api_equality": True})
        quality = folder / "quality"
        run([sys.executable, "-m", "tcren.cli", "md", "quality", "-s", str(reference),
             "--mhc-class", mhc_class, "-o", str(quality)], root)
        assert pl.read_csv(quality / "close_atom_pairs.tsv", separator="\t").to_dicts() == close_atom_pairs(reference)
    save(root / "checks.json", {"purpose": "software fixtures; two translated crystal frames and synthetic energies, not MD",
                                "checks": checks, "tcren": tcren.__version__, "module": tcren.__file__,
                                "numpy": np.__version__, "polars_threads": pl.thread_pool_size(),
                                "native_pools": threadpool_info()})
    print(json.dumps(checks, indent=2))
