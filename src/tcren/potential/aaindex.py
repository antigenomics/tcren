"""Every published residue--residue contact matrix in AAindex3, as :class:`Potential` objects.

AAindex3 is GenomeNet's section of *statistical protein contact potentials*: 47 matrices over the
20 amino acids, each transcribed by its curators from a published table. The whole flat file is
bundled (``data/aaindex3.txt``, 80 kB) rather than a hand-picked subset converted to our own
format, for three reasons: the provenance is then the upstream record itself, a reader can diff the
bundled file against a fresh download, and adding a matrix to a comparison costs a string rather
than a transcription.

Two of tcren's own bundled potentials were **identified** against this file rather than guessed:
``mj()`` is ``MIYS990106`` and ``keskin()`` is ``KESO980101``, both matching 400 of 400 cells
exactly (see :func:`identify`). The MJ one had carried an "upstream table unknown" warning since
2026-08-11 and is Miyazawa--Jernigan **1999**, not 1985 and not 1996.

Not every entry is a pairwise energy. :func:`catalogue` reports ``kind`` for each:

``energy``
    A pairwise contact energy. The 42 entries a scoring pipeline can use.
``count``
    Observed contact *counts*, not energies (``TANS760102``, ``MIYS960103``).
``distance``
    Side-chain centre separations in angstroms (``BONM030104``--``BONM030106``).

and ``symmetric`` for whether the matrix equals its transpose; the six ``ZHAC*`` entries are
environment-dependent (row secondary structure vs column secondary structure) and three of those
are asymmetric by construction, so they are directed potentials and must not be decomposed.

Example:
    >>> from tcren.potential import aaindex, catalogue
    >>> catalogue().filter(pl.col("kind") == "energy").height          # doctest: +SKIP
    42
    >>> aaindex("MOOG990101").name                                    # doctest: +SKIP
    'MOOG990101'
"""
from __future__ import annotations

import re
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

import numpy as np
import polars as pl

from .model import Potential

#: Entries whose cells are not an energy, so a scoring call must not reach for them.
NON_ENERGY: dict[str, str] = {
    "TANS760102": "count", "MIYS960103": "count",
    "BONM030104": "distance", "BONM030105": "distance", "BONM030106": "distance",
}


@dataclass(frozen=True, slots=True)
class AAindexEntry:
    """One parsed AAindex3 record.

    Attributes:
        accession: The ``H`` field, e.g. ``"MIYS990106"``.
        description: The ``D`` field, one line of prose.
        authors: The ``A`` field, verbatim.
        title: The ``T`` field, verbatim (AAindex truncates long titles).
        journal: The ``J`` field, verbatim.
        pmid: The ``R`` field with the ``PMID:`` prefix stripped, or ``""``. **AAindex sometimes
            cites the paper that tabulated a matrix rather than the one that derived it** --
            ``MIYS850102`` carries Bastolla 2001 -- so verify before citing.
        rows: Amino-acid symbols down the rows.
        cols: Amino-acid symbols across the columns.
        matrix: Dense ``(20, 20)`` array, missing cells ``nan``.
        kind: ``"energy"``, ``"count"`` or ``"distance"``.
    """

    accession: str
    description: str
    authors: str
    title: str
    journal: str
    pmid: str
    rows: str
    cols: str
    matrix: np.ndarray
    kind: str

    @property
    def symmetric(self) -> bool:
        """Whether the matrix equals its transpose (so the one-body split is defined)."""
        return bool(self.rows == self.cols and np.allclose(self.matrix, self.matrix.T,
                                                           equal_nan=True))

    @property
    def n_missing(self) -> int:
        """Cells AAindex records as ``-`` or ``NA``; four entries have 39 of them."""
        return int(np.isnan(self.matrix).sum())

    def to_potential(self) -> Potential:
        """This entry as a :class:`Potential`, named for its accession. Missing cells are dropped."""
        rows = [{"residue.aa.from": self.rows[i], "residue.aa.to": self.cols[j],
                 "value": float(self.matrix[i, j])}
                for i in range(len(self.rows)) for j in range(len(self.cols))
                if not np.isnan(self.matrix[i, j])]
        alphabet = tuple(sorted(set(self.rows) | set(self.cols)))
        return Potential(name=self.accession, matrix=pl.DataFrame(rows), alphabet=alphabet)


def _bundled_aaindex3() -> Path:
    return Path(__file__).resolve().parent.parent / "data" / "aaindex3.txt"


def parse_aaindex3(text: str) -> dict[str, AAindexEntry]:
    """Parse an AAindex3 flat file into one :class:`AAindexEntry` per accession.

    The ``M`` block is **lower-triangular including the diagonal** for a symmetric entry and a full
    rectangle for a directed one; both forms appear in the file and both are handled. A cell of
    ``-`` or ``NA`` becomes ``nan`` rather than an error, because four entries genuinely omit one
    residue's whole row and column.

    Args:
        text: The contents of an ``aaindex3`` flat file.

    Returns:
        Accession -> entry.

    Raises:
        ValueError: If a record's ``M`` block matches neither shape, which would mean the upstream
            format changed and silently mis-parsing it would corrupt every downstream score.
    """
    out: dict[str, AAindexEntry] = {}
    for block in text.split("\n//\n"):
        if not block.strip():
            continue
        acc = re.search(r"^H (\S+)", block, re.M).group(1)
        m = re.search(r"^M rows = (\S+?),?\s+cols = (\S+)", block, re.M)
        if m is None:
            continue
        rows, cols = m.group(1).rstrip(","), m.group(2)
        body = [ln for ln in block[block.index(m.group(0)) + len(m.group(0)):]
                .strip().splitlines() if ln.strip()]
        vals = [[_num(v) for v in ln.split()] for ln in body]
        A = np.full((len(rows), len(cols)), np.nan)
        if rows == cols and all(len(v) == i + 1 for i, v in enumerate(vals)):
            for i, v in enumerate(vals):                       # lower triangle, mirrored
                for j, x in enumerate(v):
                    A[i, j] = A[j, i] = x
        elif len(vals) == len(rows) and all(len(v) == len(cols) for v in vals):
            A = np.array(vals, float)                          # full rectangle, as given
        else:
            raise ValueError(f"{acc}: M block is neither a lower triangle nor a full rectangle")
        out[acc] = AAindexEntry(
            accession=acc, description=_field(block, "D"), authors=_field(block, "A"),
            title=_field(block, "T"), journal=_field(block, "J"),
            pmid=_field(block, "R").replace("PMID:", "").strip(),
            rows=rows, cols=cols, matrix=A, kind=NON_ENERGY.get(acc, "energy"))
    return out


def _num(v: str) -> float:
    try:
        return float(v)
    except ValueError:
        return float("nan")


def _field(block: str, tag: str) -> str:
    m = re.search(rf"^{tag} (.*)$", block, re.M)
    return m.group(1).strip() if m else ""


@lru_cache(maxsize=None)
def _entries() -> dict[str, AAindexEntry]:
    return parse_aaindex3(_bundled_aaindex3().read_text())


def aaindex(accession: str) -> Potential:
    """Load one AAindex3 matrix as a :class:`Potential`.

    Args:
        accession: e.g. ``"MOOG990101"``. Case-insensitive.

    Returns:
        The potential, named for the accession.

    Raises:
        KeyError: If the accession is not in AAindex3.
        ValueError: If the entry is a contact *count* or a *distance* table, which are in the file
            but are not energies -- scoring a contact map with one is a silent category error, so
            it is refused rather than allowed through.
    """
    e = _entries().get(accession.upper())
    if e is None:
        raise KeyError(f"{accession!r} is not an AAindex3 accession; see catalogue()")
    if e.kind != "energy":
        raise ValueError(f"{e.accession} holds {e.kind} values, not contact energies "
                         f"({e.description!r}); load it with entry() if that is what you want")
    return e.to_potential()


def entry(accession: str) -> AAindexEntry:
    """The parsed record itself, including the non-energy tables :func:`aaindex` refuses."""
    e = _entries().get(accession.upper())
    if e is None:
        raise KeyError(f"{accession!r} is not an AAindex3 accession; see catalogue()")
    return e


def catalogue() -> pl.DataFrame:
    """Every bundled AAindex3 entry, one row each, with what a caller needs to choose between them.

    Columns: ``accession``, ``kind``, ``symmetric``, ``n_missing``, ``mean``, ``min``, ``max``,
    ``description``, ``authors``, ``journal``, ``pmid``.

    ``mean`` is the column to read for reference state: a matrix with mean near zero is a
    *pair-contact* form with the one-body transfer term removed, one with a large negative mean is
    a *raw contact energy* that still carries it. Comparing across the two answers a different
    question from comparing within (see :meth:`Potential.components`).
    """
    rows = []
    for e in _entries().values():
        rows.append({"accession": e.accession, "kind": e.kind, "symmetric": e.symmetric,
                     "n_missing": e.n_missing, "mean": float(np.nanmean(e.matrix)),
                     "min": float(np.nanmin(e.matrix)), "max": float(np.nanmax(e.matrix)),
                     "description": e.description, "authors": e.authors,
                     "journal": e.journal, "pmid": e.pmid})
    return pl.DataFrame(rows).sort("accession")


def identify(potential: Potential, tol: float = 1e-9) -> list[tuple[str, float]]:
    """Which AAindex3 entries a potential matches, best first, as ``(accession, max |delta|)``.

    Written for a matrix whose upstream table was never recorded: run it and the answer is either
    an exact match or a shortlist. ``mj()`` and ``keskin()`` were identified this way, at max
    ``|delta| = 0`` over all 400 cells, with every other candidate off by at least 0.66.

    Args:
        potential: The matrix to identify.
        tol: Report matches at or below this max absolute difference first; everything else
            follows in order, so a near-miss is visible rather than silently dropped.

    Returns:
        ``(accession, max_abs_delta)`` for every entry sharing the alphabet, ascending by delta.
    """
    ours = {(r["residue.aa.from"], r["residue.aa.to"]): r["value"]
            for r in potential.matrix.iter_rows(named=True)}
    hits = []
    for e in _entries().values():
        theirs = {(e.rows[i], e.cols[j]): e.matrix[i, j]
                  for i in range(len(e.rows)) for j in range(len(e.cols))}
        shared = [k for k in ours if k in theirs and not np.isnan(theirs[k])]
        if len(shared) < 100:
            continue
        hits.append((e.accession, float(max(abs(ours[k] - theirs[k]) for k in shared))))
    hits.sort(key=lambda h: h[1])
    del tol                                                   # ordering already puts matches first
    return hits
