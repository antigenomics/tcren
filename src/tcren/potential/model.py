"""Pairwise residue-level statistical potentials.

A :class:`Potential` is a long-form table of pairwise amino-acid energies keyed on
``(residue.aa.from, residue.aa.to)``. The "from" side is conventionally the TCR
residue and the "to" side the antigen (peptide) residue, matching the orientation of
the legacy R pipeline. Potentials can be loaded from the two CSV layouts shipped with
the project (wide and long) and exported to a dense matrix for fast scoring.
"""

from __future__ import annotations

from dataclasses import dataclass
from importlib import resources
from pathlib import Path

import numpy as np
import polars as pl

#: 20 standard amino acids (one-letter), TCRen ordering used in the paper.
AA20: tuple[str, ...] = (
    "L", "F", "I", "M", "V", "W", "Y", "C", "H", "A",
    "G", "P", "T", "S", "Q", "N", "D", "E", "R", "K",
)

#: Alphabet of the alignment-matrix variant: 21 amino acids plus the gap symbol.
AA21: tuple[str, ...] = (
    "A", "I", "L", "V", "R", "H", "K", "C", "M", "S", "T",
    "D", "E", "N", "Q", "G", "P", "Y", "F", "W", "-",
)

_LONG_COLUMNS = {"residue.aa.from", "residue.aa.to", "potential", "value"}


@dataclass(slots=True)
class Potential:
    """A pairwise amino-acid potential in long form.

    Attributes:
        name: Identifier of the potential (e.g. ``"TCRen"``, ``"MJ"``, ``"Keskin"``).
        matrix: Long-form table with columns ``residue.aa.from``, ``residue.aa.to``,
            ``value``.
        alphabet: Amino-acid symbols present on each axis.
    """

    name: str
    matrix: pl.DataFrame
    alphabet: tuple[str, ...]

    def value(self, aa_from: str, aa_to: str) -> float:
        """Return the energy for an ordered residue pair.

        Args:
            aa_from: One-letter code of the "from" (TCR) residue.
            aa_to: One-letter code of the "to" (antigen) residue.

        Returns:
            The pairwise energy.

        Raises:
            KeyError: If the pair is absent from the potential.
        """
        hit = self.matrix.filter(
            (pl.col("residue.aa.from") == aa_from) & (pl.col("residue.aa.to") == aa_to)
        )
        if hit.height == 0:
            raise KeyError(f"pair ({aa_from!r}, {aa_to!r}) not in potential {self.name!r}")
        return float(hit["value"][0])

    def as_matrix(self) -> tuple[np.ndarray, dict[str, int]]:
        """Return a dense ``(n, n)`` matrix and an amino-acid → index map.

        Rows are indexed by ``residue.aa.from``, columns by ``residue.aa.to``. Missing
        pairs are filled with ``nan``.
        """
        index = {aa: i for i, aa in enumerate(self.alphabet)}
        n = len(self.alphabet)
        dense = np.full((n, n), np.nan, dtype=np.float64)
        for row in self.matrix.iter_rows(named=True):
            fr, to = row["residue.aa.from"], row["residue.aa.to"]
            if fr in index and to in index:
                dense[index[fr], index[to]] = row["value"]
        return dense, index

    def to_csv(self, path: str | Path) -> None:
        """Write the potential to a long-form CSV (``from, to, value``)."""
        self.matrix.write_csv(str(path))

    @classmethod
    def from_csv(
        cls,
        path: str | Path,
        name: str | None = None,
        value_col: str | None = None,
    ) -> "Potential":
        """Load a potential from a CSV, auto-detecting wide vs long layout.

        Two layouts are supported:

        * **wide** — ``residue.aa.from, residue.aa.to, <name>`` (e.g.
          ``TCRen_potential.csv`` with a ``TCRen`` value column).
        * **long** — ``residue.aa.from, residue.aa.to, potential, value`` (e.g.
          ``MJ_Keskin_potentials.csv``); load a single named potential from it.

        Args:
            path: Path to the CSV file.
            name: Which potential to select (long layout) or the name to assign
                (wide layout). Defaults to the value-column name (wide) and is
                required when a long file holds more than one potential.
            value_col: Override the value column name for the wide layout.

        Returns:
            The loaded :class:`Potential`.
        """
        df = pl.read_csv(str(path))
        cols = set(df.columns)

        if _LONG_COLUMNS.issubset(cols):
            potentials = df["potential"].unique().to_list()
            if name is None:
                if len(potentials) != 1:
                    raise ValueError(
                        f"{path} holds multiple potentials {potentials!r}; pass name="
                    )
                name = potentials[0]
            sel = df.filter(pl.col("potential") == name).select(
                "residue.aa.from", "residue.aa.to", "value"
            )
            if sel.height == 0:
                raise ValueError(f"potential {name!r} not found in {path}")
            alphabet = _infer_alphabet(sel)
            return cls(name=name, matrix=sel, alphabet=alphabet)

        # Wide layout: the third column carries the values.
        key_cols = ["residue.aa.from", "residue.aa.to"]
        candidates = [c for c in df.columns if c not in key_cols]
        if value_col is None:
            if len(candidates) != 1:
                raise ValueError(
                    f"cannot infer value column in {path}; candidates={candidates!r}"
                )
            value_col = candidates[0]
        long = df.select(*key_cols, pl.col(value_col).alias("value"))
        alphabet = _infer_alphabet(long)
        return cls(name=name or value_col, matrix=long, alphabet=alphabet)


def _infer_alphabet(long: pl.DataFrame) -> tuple[str, ...]:
    """Union of symbols on both axes, ordered against the known alphabets."""
    seen = set(long["residue.aa.from"].to_list()) | set(long["residue.aa.to"].to_list())
    for known in (AA20, AA21):
        if seen <= set(known):
            return tuple(a for a in known if a in seen)
    return tuple(sorted(seen))


def _bundled(filename: str) -> Path:
    """Resolve a CSV shipped under ``tcren/data``."""
    return resources.files("tcren.data").joinpath(filename)


def tcren() -> Potential:
    """Load the bundled classic TCRen potential."""
    return Potential.from_csv(_bundled("TCRen_potential.csv"), name="TCRen")


def mj() -> Potential:
    """Load the bundled Miyazawa–Jernigan potential."""
    return Potential.from_csv(_bundled("MJ_Keskin_potentials.csv"), name="MJ")


def keskin() -> Potential:
    """Load the bundled Keskin contact potential."""
    return Potential.from_csv(_bundled("MJ_Keskin_potentials.csv"), name="Keskin")
