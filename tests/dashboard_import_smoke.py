"""Smoke-test both dashboard CLI modes on bundled public PDB fixtures."""

import csv
import json
import subprocess
import sys
from pathlib import Path
from tempfile import TemporaryDirectory

FIXTURES = Path(__file__).parent / "assets" / "pdb"


def main():
    with TemporaryDirectory() as scratch:
        root = Path(scratch)
        for task, label_kind, structures, labels in (
            ("tcr", "binary", ("1ao7", "1bd2"), ("1", "0")),
            ("epitope", "continuous", ("1fo0", "1nam"), ("2.5", "-0.2")),
        ):
            metadata = root / f"{task}.tsv"
            with metadata.open("w", newline="", encoding="utf-8") as stream:
                writer = csv.writer(stream, delimiter="\t")
                writer.writerow(("structure_id", "structure_path", "iptm", "plddt", "label"))
                for identity, label in zip(structures, labels, strict=True):
                    writer.writerow((identity, FIXTURES / f"{identity}.pdb", "0.8", "90", label))
            bundle = root / task
            subprocess.run(
                [sys.executable, "-m", "tcren", "dashboard", "import", str(metadata),
                 "--task", task, "--label-kind", label_kind, "-o", str(bundle), "-t", "2"],
                check=True,
            )
            with (bundle / "scores.tsv").open(encoding="utf-8") as stream:
                scores = list(csv.DictReader(stream, delimiter="\t"))
            assert [row["structure_id"] for row in scores] == list(structures)
            assert [float(row["label"]) for row in scores] == list(map(float, labels))
            assert all(row["score_unavailable"] for row in scores)
            for view in ("maps", "views"):
                assert (bundle / view / "index.html").is_file()
                index = json.loads((bundle / view / "index.json").read_text(encoding="utf-8"))
                assert [row["id"] for row in index] == list(structures)
            print(f"Verified {task} dashboard with {len(scores)} fixture models", flush=True)


if __name__ == "__main__":
    main()
