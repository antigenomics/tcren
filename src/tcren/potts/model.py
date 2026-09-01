"""The coupled Potts model over a TCR:pMHC contact map — parameters and serialisation.

A **site** ``a = (i, j)`` is an *available* residue pair: a receptor residue ``i`` and a partner
residue ``j`` whose Cα atoms lie within :attr:`PottsModel.radius`. ``sigma_a = 1`` iff the pair
formed a heavy-atom contact within :attr:`PottsModel.cutoff`. A whole contact map is the
configuration ``sigma``, and the model is a Boltzmann distribution over it::

    E(sigma) = - sum_a eta_a sigma_a - 1/2 sum_{a,b} A_ab sigma_a sigma_b
    P(sigma) = exp(-E(sigma)) / Z ,   Z = sum over all 2^n configurations

``eta_a`` is the one-body log-odds of a contact at that site and ``A`` the coupling matrix built
from the kernel coefficients. Without ``A`` the partition function factorises site by site and
``log Z = sum_a log(1 + exp(eta_a))`` in closed form; that factorised model is the ``beta = 0``
reference the partition function of the coupled one is estimated against
(:func:`tcren.potts.ais_log_z`).

The one-body term is additive over categorical blocks::

    eta_a = alpha + h_rec(aa_i) + h_par(aa_j) + J(aa_i, aa_j)
            + g_dist(bin) + g_region(region_i) + g_role(role_j) + g_class(partner class)

so the fields carry single-residue propensity, ``J`` the pair chemistry, and ``g`` the geometry and
annotation the fields must be adjusted for. Every block is in the **zero-sum (Ising) gauge** — each
one sums to zero and ``J`` is double-centred — which puts as much as possible in the fields and as
little as necessary in the couplings, and makes ``J`` comparable to a double-centred
:class:`~tcren.potential.Potential`.

``J`` is either free (a full asymmetric 20x20, ``coupling_matrix_name`` is ``None``) or **one scale
on a bundled potential**: ``J = beta_matrix * (-M_dc)`` with ``M_dc`` the named potential
double-centred. The second form is 1 parameter against 400 and is how a potential is tested on an
interface.
"""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field
from importlib import resources
from pathlib import Path

import numpy as np

#: Residue alphabet for the fields and the coupling. Alphabetical, and written into every
#: serialised model, so a stored ``h``/``J`` never depends on an import-time ordering.
AA: tuple[str, ...] = tuple("ACDEFGHIKLMNPQRSTVWY")

#: Receptor region levels for the ``g_region`` block. ``other`` absorbs anything unannotated.
REGIONS: tuple[str, ...] = tuple(
    f"{c}:{r}" for c in ("TRA", "TRB") for r in ("CDR1", "CDR2", "CDR3", "FR1", "FR2", "FR3", "FR4")
) + ("other",)

#: Partner-role levels: the peptide's anchor call, then the four MHC groove regions.
ROLES: tuple[str, ...] = ("tcr_facing", "anchor", "HELIX_A1", "HELIX_A2", "HELIX_B1", "GROOVE_FLOOR")

#: Partner classes. A joint model carries both; a single-interface model still stores the block.
CLASSES: tuple[str, ...] = ("peptide", "mhc")

#: Width of a Cα-distance bin, Å. The distance profile is nonparametric rather than linear so
#: that every block stays categorical and the design stays a sum of one-hots.
DBIN: float = 0.5

#: Within-loop coupling classes, one representative per ``{d, -d}`` pair over ``|di| <= 2``,
#: ``|dj| <= 2`` minus ``(0, 0)``. ``di`` is a receptor sequence offset **inside one loop**,
#: ``dj`` a partner offset inside one chain.
OFFSETS: tuple[tuple[int, int], ...] = (
    (0, 1), (0, 2), (1, -2), (1, -1), (1, 0), (1, 1), (1, 2),
    (2, -2), (2, -1), (2, 0), (2, 1), (2, 2),
)

#: Cross-loop coupling classes: two sites in *different* hypervariable loops at partner offset
#: ``|dj|``, crossed with whether the two loops sit on the same receptor chain.
CROSS_DJ: tuple[int, ...] = (0, 1, 2)

#: Cross-loop couplings are restricted to the hypervariable loops. A framework residue at the
#: periphery is not gripping anything, and admitting every loop pair makes the coupling graph
#: dense enough to cost the sampler far more than it buys.
CDR_LOOPS: tuple[str, ...] = tuple(
    f"{c}:{r}" for c in ("TRA", "TRB") for r in ("CDR1", "CDR2", "CDR3")
)


def kernel_names(joint: bool) -> list[str]:
    """Names of the coupling classes, in the order :func:`tcren.potts.edges` returns them."""
    names = [f"K({di:+d},{dj:+d})" for di, dj in OFFSETS]
    names += [f"L(|dj|={dj},{'same' if s else 'cross'}chain)" for dj in CROSS_DJ for s in (1, 0)]
    if joint:
        names += ["M(same receptor residue, peptide vs MHC)"]
    return names


@dataclass
class PottsModel:
    """Fitted parameters of the coupled contact-map model, plus everything needed to re-apply them.

    Attributes:
        alpha: Intercept — the log-odds of a contact at the mean of every block.
        h_rec: Receptor field, one entry per :data:`AA`. Positive = engages an available partner
            more often than average; negative = declines it.
        h_par: Partner field, same convention and alphabet.
        coupling: Free 20x20 ``J(receptor aa, partner aa)``, double-centred, or ``None`` when
            ``coupling_matrix_name`` is set.
        beta_matrix: Scale on the named potential when ``J`` is not free.
        coupling_matrix_name: Bundled potential the coupling is fixed to (``tcren2``, ``mj``,
            ``mj1996``, ``keskin``, ``tcren``), or ``None`` for a free ``J``.
        g_dist: Cα-distance-bin coefficients, bin ``b`` covering ``[b*DBIN, (b+1)*DBIN)`` Å.
        g_region: Receptor-region coefficients over :data:`REGIONS`.
        g_role: Partner-role coefficients over :data:`ROLES`.
        g_class: Partner-class coefficients over :data:`CLASSES`.
        kernel: Coupling coefficients in the order of :func:`kernel_names`.
        kernel_se: Cluster-robust standard errors for ``kernel``, clustered on the structure.
        radius: Availability radius, Å (Cα–Cα). A pair beyond it is not a site.
        cutoff: Contact definition, Å (closest heavy atom).
        joint: Whether the model was fitted over both partners at once.
        alphabet: The residue alphabet ``h`` and ``J`` are indexed by.
        n_structures / n_sites / n_contacts: What the fit was estimated on.
        pseudo_loglik: Penalised pseudo-log-likelihood at the optimum.
        notes: Free text, e.g. the structure set.
    """

    alpha: float
    h_rec: list[float]
    h_par: list[float]
    g_dist: list[float]
    g_region: list[float]
    g_role: list[float]
    g_class: list[float]
    kernel: list[float]
    coupling: list[list[float]] | None = None
    beta_matrix: float | None = None
    coupling_matrix_name: str | None = None
    kernel_se: list[float] = field(default_factory=list)
    radius: float = 15.0
    cutoff: float = 5.0
    joint: bool = False
    alphabet: str = "".join(AA)
    regions: list[str] = field(default_factory=lambda: list(REGIONS))
    roles: list[str] = field(default_factory=lambda: list(ROLES))
    classes: list[str] = field(default_factory=lambda: list(CLASSES))
    dbin: float = DBIN
    n_structures: int = 0
    n_sites: int = 0
    n_contacts: int = 0
    pseudo_loglik: float = float("nan")
    notes: str = ""

    # -- coupling -----------------------------------------------------------------------------

    def coupling_array(self) -> np.ndarray:
        """``J`` as a 20x20 array, whether it was fitted freely or fixed to a named potential."""
        if self.coupling is not None:
            return np.asarray(self.coupling, dtype=float)
        if self.coupling_matrix_name is None:
            raise ValueError("model has neither a free coupling nor a coupling_matrix_name")
        return float(self.beta_matrix) * centred_potential(self.coupling_matrix_name,
                                                           tuple(self.alphabet))

    def n_parameters(self) -> int:
        """Free parameters, counting a fixed-matrix coupling as the single scale it is."""
        n = 1 + len(self.h_rec) + len(self.h_par) + len(self.g_dist) + len(self.g_region)
        n += len(self.g_role) + len(self.g_class) + len(self.kernel)
        return n + (400 if self.coupling is not None else 1)

    # -- serialisation ------------------------------------------------------------------------

    def to_json(self, path: str | Path) -> None:
        """Write the model as JSON (round-trips through :meth:`from_json`)."""
        Path(path).write_text(json.dumps(asdict(self), indent=1) + "\n")

    @classmethod
    def from_json(cls, path: str | Path) -> "PottsModel":
        """Read a model written by :meth:`to_json`."""
        return cls(**json.loads(Path(path).read_text()))

    @classmethod
    def bundled(cls, name: str = "potts_tcr_peptide") -> "PottsModel":
        """Load a model shipped under ``tcren/data`` (default: the TCR:peptide interface)."""
        text = resources.files("tcren.data").joinpath(f"{name}.json").read_text()
        return cls(**json.loads(text))


def centred_potential(name: str, alphabet: tuple[str, ...] = AA, *,
                      centre: bool = True) -> np.ndarray:
    """``-M`` onto ``alphabet``, double-centred by default — positive means more contact.

    A :class:`~tcren.potential.Potential` is signed *negative is favourable*; a coupling here is a
    log-odds where positive is *more likely*, hence the sign flip.

    The two settings answer two different questions and are not interchangeable.

    ``centre=True`` (default) removes every one-body term, so a matrix used as a fixed coupling
    contributes **nothing** to the single-residue marginals — the fields carry those, refitted
    freely, and competing matrices are compared on pair structure alone. This is the setting for
    *ranking potentials against each other*.

    ``centre=False`` keeps the raw matrix. Use it to *reproduce a referenced contact-map score*: the
    peptide-referenced energy of :func:`tcren.ddg.reference_delta` is a difference of one-body sums,
    so double-centring re-injects a burial-scaled composition term ``n_i·c(a)`` — the position's
    contact count times the potential's partner-residue column mean — and the identity fails. Pinned
    uncentred, the model's referenced energy equals the potential's own up to the fitted scale.

    Args:
        name: One of ``tcren2``, ``tcren``, ``mj``, ``mj1996``, ``keskin``.
        alphabet: Residue order of the returned array.
        centre: Double-centre the matrix. See above for which setting a task wants.

    Returns:
        A ``(len(alphabet), len(alphabet))`` array; cells absent from the potential are ``0``.
    """
    from ..topology.pose import _double_centred
    from ..potential import model as _m

    loader = {"tcren2": _m.tcren2, "tcren": _m.tcren, "mj": _m.mj,
              "mj1996": _m.mj1996, "keskin": _m.keskin}.get(name)
    if loader is None:
        raise ValueError(f"unknown potential {name!r}; expected one of tcren2, tcren, mj, "
                         f"mj1996, keskin")
    pot = loader()
    M, index = _double_centred(pot) if centre else pot.as_matrix()
    out = np.zeros((len(alphabet), len(alphabet)))
    for i, a in enumerate(alphabet):
        for j, b in enumerate(alphabet):
            if a in index and b in index:
                v = M[index[a], index[b]]
                out[i, j] = 0.0 if np.isnan(v) else -float(v)
    return out
