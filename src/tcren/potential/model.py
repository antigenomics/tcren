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
from .._provenance import not_in_tcren2

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

    def scale(self) -> float:
        r"""The potential's own energy scale: the standard deviation over its defined pairs.

        Two potentials Boltzmann-inverted from different contact statistics are not on a common
        scale, so summing energies read off them weights whichever has the wider matrix. Measured
        over the shipped matrices: TCRen2 0.4880, MJ 0.3270, **Keskin 1.3181** -- so an unweighted
        :math:`\Phi_{\mathrm{TCR:pep}} + \Phi_{\mathrm{TCR:MHC}} + \Phi_{\mathrm{pep:MHC}}` is
        2.70x more sensitive to a presentation contact than to a recognition one when the
        presentation interfaces are scored with Keskin.

        Dividing each interface energy by its potential's scale makes the three terms commensurate.
        The coefficient is a property of the matrix alone -- no cohort, no label, no fit.

        Diagonal and off-diagonal entries are pooled and each unordered pair counted once, since
        the matrix is symmetric in use.
        """
        m, _ = self.as_matrix()
        v = np.asarray(m, float)
        iu = np.triu_indices_from(v)
        w = v[iu]
        w = w[np.isfinite(w)]
        if w.size < 2:
            raise ValueError(f"potential {self.name!r} has {w.size} defined pairs; no scale")
        return float(w.std(ddof=0))

    def offset(self) -> float:
        """The potential's mean over its defined pairs (see :meth:`scale`).

        An additive offset multiplied by a contact count is a contact count, not an energy, so a
        potential with a large one -- Keskin's mean is -3.5630 and every entry is negative -- makes
        its interface energy read mostly as interface size. Subtracting it leaves the identity
        preference, which is what the other channels do not already carry.
        """
        m, _ = self.as_matrix()
        v = np.asarray(m, float)
        iu = np.triu_indices_from(v)
        w = v[iu]
        w = w[np.isfinite(w)]
        if not w.size:
            raise ValueError(f"potential {self.name!r} has no defined pairs")
        return float(w.mean())

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

    def components(self) -> dict[str, "Potential"]:
        """The three additive parts of :meth:`decompose`, each as a scorable :class:`Potential`.

        :meth:`decompose` splits the matrix as ``e(a,b) = mean + H(a) + H(b) + J(a,b)``. Because an
        interface score is a *sum over contacts*, that split carries straight through to the score:

        ============  ==========================  ====================================
        component     matrix                      what its interface sum equals
        ============  ==========================  ====================================
        ``"size"``    the grand mean everywhere   ``mean x (number of contacts)``
        ``"comp"``    ``H(a) + H(b)``             a degree-weighted composition term
        ``"pair"``    ``J(a, b)``                 the interaction, one-body parts gone
        ============  ==========================  ====================================

        So scoring a structure with each part in turn says which of three very different things a
        potential is reading on that interface: how *big* it is, what it is *made of*, or which
        residue *faces which*. That distinction is not cosmetic -- a matrix with no positive entries
        has a large negative mean, so its interface sum is dominated by the contact count, and a
        result obtained with one can be an interface-area effect wearing a chemical name. The three
        parts sum back to the original exactly, which the unit tests assert.

        Returns:
            ``{"size": ..., "comp": ..., "pair": ...}``, each named ``<this potential>_<part>``.

        Raises:
            ValueError: If the potential is not symmetric (see :meth:`decompose`).
        """
        d = self.decompose()
        aas = sorted(d.index, key=lambda a: d.index[a])
        H, J = d.one_body, d.pair
        parts = {
            "size": lambda a, b: d.mean,
            "comp": lambda a, b: float(H[d.index[a]] + H[d.index[b]]),
            "pair": lambda a, b: float(J[d.index[a], d.index[b]]),
        }
        return {
            tag: Potential(
                name=f"{self.name}_{tag}",
                matrix=pl.DataFrame([
                    {"residue.aa.from": a, "residue.aa.to": b, "value": float(fn(a, b))}
                    for a in aas for b in aas
                ]),
                alphabet=tuple(aas),
            )
            for tag, fn in parts.items()
        }

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
@not_in_tcren2('The 2022 matrix, kept for reproducing published results. TCRen2 is tcren.potential.tcren2() and is the default since 2.11.0; the two correlate at r = 0.867 with a maximum absolute difference of 0.943 and are not interchangeable.')
def tcren() -> Potential:
    """Load the bundled classic TCRen potential (cached; treat as read-only)."""
    return Potential.from_csv(_bundled("TCRen_potential.csv"), name="TCRen")


@lru_cache(maxsize=None)
def tcren2() -> Potential:
    """Load the bundled TCRen2 potential (cached; treat as read-only).

    The redundancy-balanced derivation over the **362 fully annotated αβ** ``Native2026``
    crystals, the default TCR:peptide potential since 2.11.0, and the matrix the TCRen2
    manuscript reports. It is **not** interchangeable with :func:`tcren`: the two correlate at
    Pearson *r* = 0.867 with a maximum absolute difference of 0.943 over a range of 2.95, so
    scores computed under one cannot be compared with scores under the other.
    """
    return Potential.from_csv(_bundled("TCRen2_potential.csv"), name="TCRen2")




@lru_cache(maxsize=None)
def mj() -> Potential:
    """Load the bundled Miyazawa–Jernigan potential (cached; treat as read-only).

    **Identified 2026-08-29: this is AAindex3 ``MIYS990106``, Miyazawa & Jernigan 1999** -- not
    1985 and not 1996, which is what the "upstream table not recorded" warning here used to say.
    All 400 cells match the AAindex record exactly (``identify(mj())`` returns
    ``("MIYS990106", 0.0)``), and the next-closest entry in the whole of AAindex3 is off by 0.65,
    so the identification is unique. Every score in the package is built on this file and it is
    left byte-for-byte untouched; what changed is that it can now be cited.

    It takes both signs with a mean of −0.079, so it is a contact-*pair* matrix with the one-body
    transfer term removed; :func:`mj1996` and :func:`keskin` are raw contact energies (mean ≈ −3.3)
    and :func:`betancourt` is the other pair-form matrix. Compare like with like, and see
    :meth:`Potential.components` for why the distinction changes what a comparison measures.

    Reference: Miyazawa S, Jernigan RL. Self-consistent estimation of inter-residue protein contact
    energies based on an equilibrium mixture approximation of residues. Proteins. 1999;34(1):49-68.
    doi:10.1002/(SICI)1097-0134(19990101)34:1<49::AID-PROT5>3.0.CO;2-L.
    """
    return Potential.from_csv(_bundled("MJ_Keskin_potentials.csv"), name="MJ")


@lru_cache(maxsize=None)
def mj1996() -> Potential:
    """Miyazawa--Jernigan 1996 inter-residue contact energies, ``e_ij``, in RT units.

    The 20x20 attractive contact energies of Table 3, re-evaluated by the authors on 1168
    structures. Every entry is negative, from ``-7.37`` (Leu--Leu) to ``-0.12``, and
    Ala--Ala is ``-2.72``; a raw contact matrix looks like this, and the bundled
    :func:`mj` matrix does not, which is how the two are told apart.

    Provenance is recorded because the older bundled matrix has none: the numbers here were
    transcribed from a published copy of Table 3 (AAindex accession MIYS960101) and checked
    against a second independent copy, agreeing on the alphabet order ``CMFILVWYAGTSNQDEHRKP``,
    on Ala--Ala, and on the full range. They correlate with the bundled ``MJ`` matrix at
    ``r = 0.89``, so the two are related but not the same quantity, and the bundled one is
    *not* the double-centred pair part of this one (``r = 0.51``). What the bundled matrix
    actually is remains unresolved.

    The companion repulsive packing-density term of the same paper is **not** included; it
    is a function of coordination number rather than of a residue pair, so it does not fit
    the :class:`Potential` shape and nothing here uses it.

    Reference: Miyazawa S, Jernigan RL. Residue-residue potentials with a favorable contact
    pair term and an unfavorable high packing density term, for simulation and threading.
    J Mol Biol. 1996;256(3):623-644. doi:10.1006/jmbi.1996.0114.
    """
    return Potential.from_csv(_bundled("MJ1996_contact_energies.csv"), name="MJ1996")


@lru_cache(maxsize=None)
def mj_partition_energy() -> dict[str, float]:
    """Miyazawa--Jernigan effective partition energies, one value per residue.

    The one-body term of the MJ framework: the energy of transferring a residue from water
    into the protein interior, which is what a contact energy carries in addition to any
    interaction between two identities. A pairwise matrix cannot supply this on its own, so
    it is bundled separately rather than derived.

    Larger is more hydrophobic: Phe 4.37, Met 4.22, Ile 4.17 at one end, Lys 1.23, Asp 1.67,
    Asn 1.70 at the other. Note the sign convention is opposite to a contact energy, where
    lower is more favourable.

    Provenance: AAindex accession MIYS850101, retrieved from two endpoints of the AAindex
    database that returned identical values. As an independent check, this scale correlates
    at ``r = +0.98`` with the hydrophobicity axis recovered by
    :meth:`Potential.hydrophobicity_fit` from :func:`mj1996`, which was transcribed from a
    different source entirely.

    Reference: Miyazawa S, Jernigan RL. Estimation of effective interresidue contact
    energies from protein crystal structures: quasi-chemical approximation. Macromolecules.
    1985;18:534-552.

    Returns:
        Amino acid one-letter code → partition energy. The mapping is cached; copy it before
        mutating.
    """
    table = pl.read_csv(_bundled("MJ1985_partition_energies.csv"))
    return {row["residue.aa"]: float(row["value"]) for row in table.iter_rows(named=True)}


@lru_cache(maxsize=None)
def keskin() -> Potential:
    """Load the bundled Keskin contact potential (cached; treat as read-only).

    **Identified 2026-08-29 as AAindex3 ``KESO980101``**, "Quasichemical transfer energy derived
    from interfacial regions", matching all 400 cells exactly with the next-closest AAindex3 entry
    off by 2.77. That is the *solvent-mediated* form; the companion ``KESO980102`` is the
    residue-mediated one, also available through :func:`aaindex`.

    Every entry is negative, from ``-7.23`` to ``-0.50``, so this is a raw contact matrix in
    the same reference state as :func:`mj1996` and **not** in the pair-contact reference state
    of the bundled :func:`mj` (mixed sign, mean ``-0.08``). Compare it against ``mj1996``;
    comparing it against ``mj`` compares two different reference states as well as two
    different derivations.

    Reference: Keskin O, Bahar I, Badretdinov AY, Ptitsyn OB, Jernigan RL. Empirical
    solvent-mediated potentials hold for both intra-molecular and inter-molecular
    inter-residue interactions. Protein Sci. 1998;7(12):2578-2586. doi:10.1002/pro.5560071211.
    """
    return Potential.from_csv(_bundled("MJ_Keskin_potentials.csv"), name="Keskin")


@lru_cache(maxsize=None)
def betancourt() -> Potential:
    """Betancourt--Thirumalai contact energies, the ``B`` matrix, in RT units.

    Miyazawa--Jernigan re-referenced with **Thr as the reference solvent**, which is why every
    Thr entry is exactly ``0.00``; the remaining 190 cross terms and 19 self terms run ``-1.34``
    (Cys--Cys) to ``+0.66``. Mixed sign with a mean near zero, so it is a pair-contact matrix in
    the same reference state as the bundled :func:`mj`, and that is the matrix to compare it
    against. The authors report it gives "hydrophobicities that are in very good agreement with
    experiment", and it is the potential Schueler-Furman et al. found generalises an
    MJ-based peptide--MHC groove score across alleles where MJ itself worked only for
    hydrophobic-pocket alleles.

    Provenance: parsed from AAindex3 accession ``BETM990101`` ("Modified version of the
    Miyazawa-Jernigan transfer energy"), lower-triangular over ``ARNDCQEGHILKMFPSTWYV``, never
    retyped. Three properties are asserted at build time: the Thr row is zero, the matrix is
    symmetric, and all 400 cells are present.

    Reference: Betancourt MR, Thirumalai D. Pair potentials for protein folding: choice of
    reference states and sensitivity of predicted native states to variations in the interaction
    schemes. Protein Sci. 1999;8(2):361-369. doi:10.1110/ps.8.2.361.
    """
    return Potential.from_csv(_bundled("BT1999_contact_energies.csv"), name="BT1999")
