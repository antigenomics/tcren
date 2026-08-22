"""Every shipped potential must be reproducible by the documented CLI invocation.

The package ships more than one matrix. Without a machine-checkable link from each file
back to the command that produced it, provenance drifts: a matrix gets regenerated with
different flags, or an unrelated file is dropped into ``data/`` and nothing notices.
``src/tcren/data/potentials.json`` records the recipe; this module re-runs it.

Entries whose ``source`` is ``structure-dir`` need the reference structures, which are
fetched separately (``tcren fetch-data``); those skip when the folder is absent. An entry
names its own CLI subcommand via ``command`` (default ``derive-potential``), because the
DFIRE-derived matrices come out of ``derive-dfire``.
"""

from __future__ import annotations

import json
import shutil
import subprocess
import sys
from pathlib import Path

import polars as pl
import pytest

REPO = Path(__file__).resolve().parents[2]
DATA = REPO / "src" / "tcren" / "data"
MANIFEST = DATA / "potentials.json"


def _manifest() -> dict:
    return json.loads(MANIFEST.read_text())["potentials"]


def _derivable() -> list[str]:
    return sorted(k for k, e in _manifest().items()
                  if e["source"] in ("contact-maps", "structure-dir"))


def _cli() -> list[str]:
    exe = shutil.which("tcren")
    return [exe] if exe else [sys.executable, "-m", "tcren"]


def _structure_dir(name: str) -> Path | None:
    from tcren.paths import data_dir
    for base in (Path(data_dir()), Path.home() / "hf" / "tcren_structures"):
        cand = base / name
        if cand.is_dir():
            return cand
    return None


def _read(path: Path) -> pl.DataFrame:
    d = pl.read_csv(path)
    val = "TCRen" if "TCRen" in d.columns else "value"
    return d.select("residue.aa.from", "residue.aa.to", pl.col(val).alias("v"))


def test_every_shipped_matrix_has_a_recipe():
    """A CSV in data/ that looks like a potential must appear in the manifest."""
    declared = {e["file"] for e in _manifest().values()}
    shipped = {
        p.name for p in DATA.glob("*.csv")
        if {"residue.aa.from", "residue.aa.to"} <= set(pl.read_csv(p, n_rows=1).columns)
    }
    assert shipped <= declared, (
        f"shipped potentials with no derivation recipe in {MANIFEST.name}: "
        f"{sorted(shipped - declared)}"
    )


def test_every_recipe_points_at_a_shipped_file():
    for name, e in _manifest().items():
        assert (DATA / e["file"]).is_file(), f"{name}: {e['file']} is missing from data/"


def test_no_new_matrix_ships_without_provenance():
    """Our own matrices must be reproducible; a published third-party table need not be.

    ``known_unresolved`` records historical files whose provenance cannot be recovered.
    Anything else with ``source: unknown`` is a new matrix shipped without a recipe, which
    is how three mutually inconsistent potentials came to be in circulation at once.
    """
    doc = json.loads(MANIFEST.read_text())
    allowed = set(doc.get("known_unresolved", []))
    unknown = {n for n, e in doc["potentials"].items() if e["source"] == "unknown"}
    assert unknown - allowed == set(), (
        "matrices shipped without a derivation recipe: " + ", ".join(sorted(unknown - allowed)))
    assert allowed - unknown == set(), (
        "known_unresolved lists entries that are no longer unresolved (or no longer exist): "
        + ", ".join(sorted(allowed - unknown)))


@pytest.mark.parametrize("name", _derivable())
def test_recipe_reproduces_the_shipped_matrix(name, tmp_path):
    e = _manifest()[name]
    shipped = DATA / e["file"]
    out = tmp_path / "derived.csv"

    cmd = [*_cli(), e.get("command", "derive-potential"), "-o", str(out)]
    for opt in ("variant", "pseudocount"):
        if opt in e:
            cmd += [f"--{opt}", str(e[opt])]
    cmd += e["flags"]

    if e["source"] == "contact-maps":
        cmd += ["-i", str(REPO / e["contacts"])]
        if "summary" in e:
            cmd += ["--summary", str(REPO / e["summary"])]
    else:
        sdir = _structure_dir(e["structure_dir"])
        if sdir is None:
            pytest.skip(f"{e['structure_dir']} not present; run `tcren fetch-data`")
        pytest.importorskip("arda")
        cmd += ["--structure-dir", str(sdir)]

    r = subprocess.run(cmd, capture_output=True, text=True)
    assert r.returncode == 0, f"derive-potential failed for {name}:\n{r.stderr[-2000:]}"

    got, want = _read(out), _read(shipped)
    j = got.join(want, on=["residue.aa.from", "residue.aa.to"], how="full",
                 coalesce=True, suffix="_want")
    assert j.filter(pl.col("v").is_null() | pl.col("v_want").is_null()).height == 0, (
        f"{name}: derived and shipped matrices cover different residue pairs")
    max_abs = j.select((pl.col("v") - pl.col("v_want")).abs().max()).item()
    assert max_abs == pytest.approx(0.0, abs=1e-9), (
        f"{name}: shipped {e['file']} is not what its recipe produces (max |d| = {max_abs})")
