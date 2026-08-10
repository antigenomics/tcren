"""Pairwise residue-level statistical potentials.

A :class:`Potential` is a long-form table of pairwise amino-acid energies keyed on
``(residue.aa.from, residue.aa.to)``. The "from" side is conventionally the TCR
residue and the "to" side the antigen (peptide) residue, matching the orientation of
the legacy R pipeline. Potentials can be loaded from the two CSV layouts shipped with
the project (wide and long) and exported to a dense matrix for fast scoring.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from functools import lru_cache
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
class PotentialDecomposition:
    """A potential split as ``e(a, b) = mean + H(a) + H(b) + J(a, b)``.

    Attributes:
        name: Name of the potential this came from.
        mean: The grand mean of the matrix.
        one_body: ``H``, indexed like ``index``; the per-residue part.
        pair: ``J``, double-centred, so every row and column sums to zero.
        index: Amino-acid → row/column index.
    """

    name: str
    mean: float
    one_body: np.ndarray
    pair: np.ndarray
    index: dict[str, int]

    def h(self, aa: str) -> float:
        """One-body term of a residue."""
        return float(self.one_body[self.index[aa]])

    def j(self, aa: str, bb: str) -> float:
        """Pair-specific term of a residue pair, with the one-body parts removed."""
        return float(self.pair[self.index[aa], self.index[bb]])

    def energy(self, aa: str, bb: str) -> float:
        """Reassemble the original contact energy; equals the potential's own value."""
        return self.mean + self.h(aa) + self.h(bb) + self.j(aa, bb)


@dataclass(slots=True)
class HydrophobicityFit:
    """A potential approximated as ``C0 + C1 (q_a + q_b) + C2 q_a q_b``.

    Attributes:
        name: Name of the potential this came from.
        c0, c1, c2: Fitted coefficients.
        q: One value per residue, from the leading eigenvector; orders by hydrophobicity.
        index: Amino-acid → index into ``q``.
        r2: Fraction of the matrix variance the three-parameter form reproduces.
        eigenvalue_share: ``|lambda_1| / sum |lambda|``, i.e. how nearly rank-one the
            matrix is to begin with.
    """

    name: str
    c0: float
    c1: float
    c2: float
    q: np.ndarray
    index: dict[str, int]
    r2: float
    eigenvalue_share: float

    def value(self, aa: str, bb: str) -> float:
        """The fitted contact energy for a residue pair."""
        qa, qb = self.q[self.index[aa]], self.q[self.index[bb]]
        return self.c0 + self.c1 * (qa + qb) + self.c2 * qa * qb

    def one_body(self, aa: str) -> float:
        """``C1 q_a`` -- the per-residue term, which is what ``H(a)`` refers to."""
        return self.c1 * float(self.q[self.index[aa]])


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
    # Lazily-built dense form (see as_matrix); not part of identity/repr.
    _matrix_cache: tuple[np.ndarray, dict[str, int]] | None = field(
        default=None, init=False, repr=False, compare=False
    )

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
        pairs are filled with ``nan``. The dense form is cached (the table is immutable), so
        repeated scoring/energy calls over one potential rebuild it once. Callers treat the
        returned array as read-only.
        """
        if self._matrix_cache is not None:
            return self._matrix_cache
        index = {aa: i for i, aa in enumerate(self.alphabet)}
        n = len(self.alphabet)
        dense = np.full((n, n), np.nan, dtype=np.float64)
        for row in self.matrix.iter_rows(named=True):
            fr, to = row["residue.aa.from"], row["residue.aa.to"]
            if fr in index and to in index:
                dense[index[fr], index[to]] = row["value"]
        self._matrix_cache = (dense, index)
        return self._matrix_cache

    def decompose(self) -> "PotentialDecomposition":
        """Split the potential into a one-body part and a genuinely pairwise part.

        A contact energy is not purely an interaction. Burying any residue against any
        partner costs or gains something that depends only on that residue -- its transfer
        propensity -- and only what is left after removing those one-body terms is an
        interaction between the two identities. Miyazawa and Jernigan make this split
        explicitly; here it is taken directly off the matrix, which needs no solvent
        reference and works for any potential:

            e(a, b) = mean + H(a) + H(b) + J(a, b)

        with ``H(a)`` the row mean of ``a`` less the grand mean, and ``J`` the double-centred
        remainder, whose every row and column sums to zero. The split is exact and unique.

        Why it matters for scoring: an additive per-position model can already absorb
        ``mean`` and both ``H`` terms, because they depend on one residue each. ``J`` is the
        only part that cannot be written as a sum over positions, so it is the only part a
        per-position model is actually missing.

        Returns:
            A :class:`PotentialDecomposition`.

        Raises:
            ValueError: If the dense matrix is not square and symmetric, since the split is
                only defined for an undirected potential (TCRen is *directed* and must not
                be decomposed this way).
        """
        dense, index = self.as_matrix()
        if dense.shape[0] != dense.shape[1] or not np.allclose(dense, dense.T, equal_nan=True):
            raise ValueError(
                f"potential {self.name!r} is not symmetric; the one-body/pair split is "
                "only defined for an undirected potential"
            )
        grand = float(np.nanmean(dense))
        row = np.nanmean(dense, axis=1) - grand
        pair = dense - grand - row[:, None] - row[None, :]
        return PotentialDecomposition(
            name=self.name, mean=grand, one_body=row, pair=pair, index=index
        )

    def hydrophobicity_fit(self) -> "HydrophobicityFit":
        """Fit ``e(a,b) = C0 + C1 (q_a + q_b) + C2 q_a q_b`` -- one number per residue.

        Where the one-body term comes from, for a matrix that does not ship one. Miyazawa
        and Jernigan derive their own one-body terms from residue--solvent contact energies,
        which the bundled matrices do not carry, so that route is unavailable here. Li,
        Tang and Wingreen showed it is not needed: the MJ matrix is dominated by a single
        eigenvalue, and reconstructing it from the leading eigenvector ``q`` gives the form
        above, with ``q`` ordering the residues by hydrophobicity.

        The consequence is worth stating plainly, because it limits what any MJ-based score
        can express. Not only is the one-body part a function of ``q``; so is the
        interaction, which is just ``C2 q_a q_b``. A potential of that shape knows how
        hydrophobic each residue is and nothing else -- it cannot prefer one specific pair
        of side chains over another pair of equal hydrophobicity.

        Reference: Li H, Tang C, Wingreen NS. Nature of driving force for protein folding:
        a result from analyzing the statistical potential. Phys Rev Lett. 1997;79:765.
        arXiv:cond-mat/9512111.

        Returns:
            A :class:`HydrophobicityFit`. Check its ``r2`` before relying on it; the form is
            an approximation, not an identity, unlike :meth:`decompose`.

        Raises:
            ValueError: If the matrix is not square and symmetric.
        """
        dense, index = self.as_matrix()
        if dense.shape[0] != dense.shape[1] or not np.allclose(dense, dense.T, equal_nan=True):
            raise ValueError(
                f"potential {self.name!r} is not symmetric; the hydrophobicity fit is only "
                "defined for an undirected potential"
            )
        eigenvalues, eigenvectors = np.linalg.eigh(dense)
        lead = int(np.argmax(np.abs(eigenvalues)))
        q = eigenvectors[:, lead]
        if q.mean() < 0:                      # sign of an eigenvector is arbitrary
            q = -q
        n = dense.shape[0]
        design = np.empty((n * n, 3))
        design[:, 0] = 1.0
        design[:, 1] = (q[:, None] + q[None, :]).reshape(-1)
        design[:, 2] = (q[:, None] * q[None, :]).reshape(-1)
        target = dense.reshape(-1)
        coef, *_ = np.linalg.lstsq(design, target, rcond=None)
        fitted = (design @ coef).reshape(n, n)
        ss_res = float(((dense - fitted) ** 2).sum())
        ss_tot = float(((dense - dense.mean()) ** 2).sum())
        return HydrophobicityFit(
            name=self.name, c0=float(coef[0]), c1=float(coef[1]), c2=float(coef[2]),
            q=q, index=index, r2=1.0 - ss_res / ss_tot,
            eigenvalue_share=float(abs(eigenvalues[lead]) / np.abs(eigenvalues).sum()),
        )

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


@lru_cache(maxsize=None)
def tcren() -> Potential:
    """Load the bundled classic TCRen potential (cached; treat as read-only)."""
    return Potential.from_csv(_bundled("TCRen_potential.csv"), name="TCRen")


@lru_cache(maxsize=None)
def mj() -> Potential:
    """Load the bundled Miyazawa–Jernigan potential (cached; treat as read-only)."""
    return Potential.from_csv(_bundled("MJ_Keskin_potentials.csv"), name="MJ")


@lru_cache(maxsize=None)
def keskin() -> Potential:
    """Load the bundled Keskin contact potential (cached; treat as read-only)."""
    return Potential.from_csv(_bundled("MJ_Keskin_potentials.csv"), name="Keskin")
