"""Per-structure pose consistency: do the tight contacts carry the favourable chemistry?

:func:`tcren.cohort.coupling` measures the forced-pose signature **across a cohort** as
``C* = corr(Q, dPhi)``: in a genuine complex a better interface holds more favourable contacts, so
the two channels rise together, while a generator that manufactures a pose optimises contacts
without the interface and breaks that tie. It is the right diagnostic and the wrong estimator for a
user with two or three models --- at ``n = 2`` the sample correlation is +-1 by construction, and its
sign is wrong in roughly a third to a half of draws (``bench/scripts/coupling_smalln.py``).

The tie it measures also holds **within one structure**, over that structure's own contacts. In a
crystal the residue pairs that sit tightest are the ones whose identities are complementary, because
that is what selected the pose; a pose built to satisfy a contact-density prior has no such
alignment. Correlating contact tightness against contact favourability *inside* a single complex
therefore reads the same physics from ``n = 1`` structure, over its ~20--120 interface pairs.

Three superimposable maps over one residue-pair index carry it, and
:func:`tcren.contacts.multi_contacts` already returns all three stacked (``layer`` column):

* ``d1`` --- closest heavy-atom distance (the 5 A contact definition used everywhere else);
* ``d2`` --- Cbeta distance (Calpha for glycine), i.e. where the side chains point;
* ``d3`` --- Calpha distance, i.e. where the backbones sit.

The chemistry axis is ``J``, not the raw potential entry. A contact energy splits as
``e(a,b) = mean + H_tcr(a) + H_pep(b) + J(a,b)``, where the two ``H`` terms depend on one residue
each and ``J`` is the double-centred remainder. Correlating distance against raw ``e`` would partly
measure *which* residues happen to sit at the interface rather than whether they suit each other;
``J`` is the pair-specific part, and complementarity lives there.

:meth:`tcren.Potential.decompose` performs that split but only for a *symmetric* matrix, and TCRen2
is deliberately directional (TCR residue by peptide residue; symmetrising it costs measurable
accuracy). :func:`_double_centred` therefore applies the same two-way centring to the matrix as
given, which is well defined whether or not it is symmetric and reduces to ``decompose`` when it is.

Every descriptor is oriented **higher = more crystal-like**, so they compose with
:func:`tcren.cohort.q_score` under the same all-descriptors-higher-is-better convention.

Evaluation (ROC/PR/CI) belongs downstream in the benchmark repo, not here.
"""

from __future__ import annotations

import numpy as np
import polars as pl


from ..contacts.definitions import ContactDefinition, multi_contacts
# Both moved to the layer that owns them on 2026-09-01: `_interface_layers` builds the
# d1/d2/d3 contact layers (contacts), `_double_centred` operates on a potential matrix
# (hamiltonian). `potts` was reaching UP into this module for each of them.
from ..contacts.definitions import _KEY, _REP_BUILD_CUTOFF, _interface_layers  # noqa: F401
from ..potential.model import _double_centred  # noqa: F401
from ..structure.model import MHC_TYPES, PEPTIDE_TYPE, RECEPTOR_TYPES, Structure

__all__ = ["pose_consistency", "POSE_FEATURES", "POSE_FEATURES_CONTACT",
           "POSE_FEATURES_SHELL", "POSE_FEATURES_DEGREE"]

#: The cross-map descriptors :func:`pose_consistency` returns, each oriented positive-is-crystal-like.
#: These are the ``k`` terms a pose score standardizes against a native-crystal reference.
POSE_FEATURES = (
    # --- read over the realized 5 A contacts (n ~ 20-30) ---
    "c_local",
    "e_tight_minus_loose",
    "frac_close_favourable",
    "frac_cb_close_engaged",
    "m_face_tp",
    "margin_energy_slope",
    # --- read over the Calpha approach shell, which is ~10x larger and includes the pairs that
    #     are close but form nothing (suffix _tp = TCR:peptide, _tm = TCR:MHC) ---
    "ca_energy_coupling_tp",
    "ca_energy_slope_tp",
    "frac_ca_close_engaged_tp",
    "ca_cb_agreement_tp",
    "ca_energy_coupling_tm",
    "ca_energy_slope_tm",
    "frac_ca_close_engaged_tm",
    "ca_cb_agreement_tm",
    # --- is the contact budget spread over residues, or hoarded by a few over-reaching ones? ---
    "degree_evenness_tp",
    "frac_well_coordinated_tp",
)

#: The six read over the **realized 5 A contacts**. Grouped by what they read, not by what they
#: score. This is the subset to use for **provenance** (crystal against generated): on 374 crystals
#: against 2,000 AlphaFold models it gives ROC 0.706 [0.675, 0.736] against 0.629 for the full set,
#: because the shell and degree terms carry little provenance signal and dilute the whitened sum.
POSE_FEATURES_CONTACT = POSE_FEATURES[:6]

#: The eight read over the **Calpha approach shell** --- an order of magnitude more pairs than the
#: contact set, and the only ones that see residues that are close but form nothing.
POSE_FEATURES_SHELL = POSE_FEATURES[6:14]

#: The two describing how the contact budget is **distributed over receptor residues**.
POSE_FEATURES_DEGREE = POSE_FEATURES[14:]

# The Cbeta layer threshold for "side chains are near enough that an interaction is expected".
# Distinct from the layer build cutoff below, which is deliberately generous so that every d1
# contact also has a Cbeta and a Calpha distance available to pair against.
_CB_CLOSE = 8.0
# "close" on the Calpha axis: the standard 8 A Calpha proxy for a residue-residue contact, so
# `frac_ca_close_engaged` reads "backbones near enough to touch, but do they actually?".
_CA_CLOSE = 8.0
# The Calpha approach shell. `tcren.interface_graph` reads the same shell for `ca_cb_agreement_*`,
# so the radius has one home rather than a default on each caller.
_CA_SHELL_RADIUS = 12.0
# A residue reaching more partners than this is over-coordinated: in crystals a side chain typically
# contacts one to three residues across the interface, and a long one parked against five or six is
# the shape a contact-density objective produces.
_MAX_TYPICAL_DEGREE = 3





def _pair_j(aa_from, aa_to, jmat, index) -> np.ndarray:
    """Vectorised gather of ``J(a, b)``; NaN outside the alphabet or for an unobserved cell."""
    i = np.array([index.get(a, -1) for a in aa_from], dtype=np.int64)
    j = np.array([index.get(b, -1) for b in aa_to], dtype=np.int64)
    out = np.full(len(i), np.nan)
    ok = (i >= 0) & (j >= 0)
    out[ok] = jmat[i[ok], j[ok]]
    return out


def _spearman(x: np.ndarray, y: np.ndarray) -> float:
    """Spearman rho over the finite pairs; NaN when fewer than 3 remain or either side is constant."""
    from scipy.stats import spearmanr

    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() < 3 or np.std(x[ok]) < 1e-12 or np.std(y[ok]) < 1e-12:
        return float("nan")  # guard before the call: scipy warns rather than returning quietly
    return float(spearmanr(x[ok], y[ok]).statistic)




def _degree_descriptors(contacts: pl.DataFrame, suffix: str) -> dict[str, float]:
    """Contact-degree structure: is the interface spread over residues, or hoarded by a few?

    A real interface distributes its contacts --- a given receptor side chain typically reaches one
    to three peptide residues. A pose optimised for contact count can instead park one long residue
    (Arg, Lys, Trp) against five or six at once, which is cheap for a density objective and rare in
    a crystal. Read as the participation ratio of the degree distribution (1 = perfectly even) plus
    the fraction of receptor residues that are *not* over-coordinated.

    Counted on the **receptor side only**. A peptide residue lies inside the groove ringed by CDR
    loops and carries a high partner count in every real complex, so pooling both sides would
    measure peptide burial rather than receptor over-reach.
    """
    out = {f"degree_evenness{suffix}": float("nan"),
           f"frac_well_coordinated{suffix}": float("nan"),
           f"max_degree{suffix}": float("nan")}
    if contacts.is_empty():
        return out
    deg = (contacts.group_by(["key.tcr.chain", "key.tcr.res"]).len()["len"]
           .to_numpy().astype(float))
    if not len(deg):
        return out
    out[f"max_degree{suffix}"] = float(deg.max())
    out[f"frac_well_coordinated{suffix}"] = float((deg <= _MAX_TYPICAL_DEGREE).mean())
    # participation ratio, normalised to [0, 1]: 1 when every residue carries the same degree,
    # ~1/n when a single residue hoards every contact
    out[f"degree_evenness{suffix}"] = float(deg.sum() ** 2 / (len(deg) * (deg ** 2).sum()))
    return out


def _ca_map_descriptors(wide: pl.DataFrame, jmat, jindex, cutoff: float,
                        ca_radius: float, suffix: str) -> dict[str, float]:
    """Descriptors of the Calpha interface *neighbourhood*, not just the realized contacts.

    The realized 5 A contacts of one TCR:peptide interface number ~20-30, which is thin for a
    correlation. The Calpha neighbourhood within ``ca_radius`` is an order of magnitude larger and
    contains the pairs that matter most here: those whose backbones are close but which form **no**
    contact. A crystal packs complementary residues into its close-approach shell; a pose built to
    satisfy a contact-density prior fills that shell without regard to which identities are there.
    """
    out = {f"ca_energy_coupling{suffix}": float("nan"),
           f"ca_energy_slope{suffix}": float("nan"),
           f"frac_ca_close_engaged{suffix}": float("nan"),
           f"ca_cb_agreement{suffix}": float("nan"),
           f"n_ca_near{suffix}": float("nan")}
    near = wide.filter(pl.col("d3").is_not_null() & (pl.col("d3") <= ca_radius))
    out[f"n_ca_near{suffix}"] = float(near.height)
    if near.height < 3:
        return out

    d_ca = near["d3"].to_numpy()
    j = _pair_j(near["aa.tcr"].to_list(), near["aa.pep"].to_list(), jmat, jindex)
    fav = -j
    # closer Calpha should mean more complementary chemistry, over the WHOLE approach shell
    out[f"ca_energy_coupling{suffix}"] = _spearman(-d_ca, fav)
    ok = np.isfinite(fav)
    if ok.sum() >= 3 and np.std(d_ca[ok]) > 1e-12:
        out[f"ca_energy_slope{suffix}"] = float(-np.polyfit(d_ca[ok], fav[ok], 1)[0])

    # the user's "close Calpha that do not form good contacts", read directly
    close = near.filter(pl.col("d3") <= _CA_CLOSE)
    if close.height:
        out[f"frac_ca_close_engaged{suffix}"] = float(close["d1"].is_not_null().mean())

    both = near.filter(pl.col("d2").is_not_null())
    if both.height >= 3:
        # do the side chains track the backbone? they do in a real interface
        out[f"ca_cb_agreement{suffix}"] = _spearman(both["d3"].to_numpy(), both["d2"].to_numpy())
    return out


def pose_consistency(
    structure: Structure, potential=None, cutoff: float = 5.0,
    ca_radius: float = _CA_SHELL_RADIUS
) -> dict[str, float]:
    """Cross-map consistency descriptors of one TCR:peptide interface.

    Reads whether the structure's *tight* contacts are its *complementary* ones --- the
    within-structure analogue of :func:`tcren.cohort.coupling`, and unlike it defined for a single
    complex. Every value is oriented so that **higher is more crystal-like**.

    Args:
        structure: a chain-typed complex (:func:`tcren.annotation.classify_chains` run) with a
            peptide chain and at least one receptor chain.
        potential: the residue-pair potential whose double-centred ``J`` supplies the chemistry
            axis; defaults to the bundled TCRen2 matrix.
        cutoff: the heavy-atom contact cutoff (A) defining the d1 layer and the contact margin.

    Returns:
        A dict with :data:`POSE_FEATURES` plus ``n_contacts`` (the pair count every value rests on)
        and ``n_cb_close``. Descriptors that cannot be estimated --- fewer than three contacts, a
        constant axis, no Cbeta-close pairs --- come back as ``nan`` rather than a made-up number.

    Note:
        This is a *pose* readout, not a binder score: it says whether the geometry and the chemistry
        of one model agree, not whether the receptor binds.
    """
    if potential is None:
        from ..potential import tcren2

        potential = tcren2()
    jmat, jindex = _double_centred(potential)

    wide = _interface_layers(structure, cutoff)
    contacts = wide.filter(pl.col("d1").is_not_null())
    n = contacts.height
    out: dict[str, float] = {k: float("nan") for k in POSE_FEATURES}
    out["n_contacts"] = float(n)
    out["n_cb_close"] = float("nan")

    # --- the Cbeta-engagement descriptor lives on the d2 layer, not on the contacts ---------------
    cb_close = wide.filter(pl.col("d2").is_not_null() & (pl.col("d2") <= _CB_CLOSE))
    out["n_cb_close"] = float(cb_close.height)
    if cb_close.height:
        out["frac_cb_close_engaged"] = float(cb_close["d1"].is_not_null().mean())

    # --- do the side chains lean in? Cbeta closer than Calpha means they point at each other ------
    # A mean over pairs, so unlike the correlations below it is defined for a single contact.
    both = contacts.filter(pl.col("d2").is_not_null() & pl.col("d3").is_not_null())
    if both.height:
        out["m_face_tp"] = float(
            (both["d3"].to_numpy() - both["d2"].to_numpy()).mean()
        )

    # The degree structure and the Calpha-shell descriptors carry their own guards and do NOT need
    # three *contacts* -- the shell is an order of magnitude larger than the contact set -- so they
    # are computed before the contact-correlation block bails out.
    out.update(_degree_descriptors(contacts, "_tp"))
    out.update(_ca_map_descriptors(wide, jmat, jindex, cutoff, ca_radius, "_tp"))
    # The TCR:MHC interface is the same physics on a much larger shell, and it is the half a
    # generator has the most freedom to invent: the peptide is short and anchored, the MHC is not.
    out.update(_ca_map_descriptors(
        _interface_layers(structure, cutoff, partner=MHC_TYPES),
        jmat, jindex, cutoff, ca_radius, "_tm"))

    # Everything below is a correlation, a slope or a tercile split over the realized contacts,
    # and needs at least three of them.
    if n < 3:
        return out

    d1 = contacts["d1"].to_numpy()
    j = _pair_j(contacts["aa.tcr"].to_list(), contacts["aa.pep"].to_list(), jmat, jindex)
    margin = cutoff - d1          # positional slack: higher = the pair sits deeper than the cutoff
    fav = -j                      # favourability: J is an energy, so lower J is better

    out["c_local"] = _spearman(margin, fav)

    ok = np.isfinite(fav)
    if ok.sum() >= 3 and np.std(margin[ok]) > 1e-12:
        # per-Angstrom slope of favourability on slack; the signed, physically-scaled companion
        out["margin_energy_slope"] = float(np.polyfit(margin[ok], fav[ok], 1)[0])

    if ok.sum() >= 3:
        lo, hi = np.quantile(d1[ok], [1 / 3, 2 / 3])
        tight, loose = fav[ok][d1[ok] <= lo], fav[ok][d1[ok] >= hi]
        if len(tight) and len(loose):
            out["e_tight_minus_loose"] = float(tight.mean() - loose.mean())
        med = np.median(d1[ok])
        close = j[ok][d1[ok] <= med]
        if len(close):
            out["frac_close_favourable"] = float((close < 0).mean())
    return out


def _selfcheck() -> None:
    """Assert the descriptors read the sign they claim, on a hand-built two-pair interface."""
    index = {"A": 0, "B": 1, "C": 2}
    jm = np.array([[0.0, -1.0, 1.0], [-1.0, 0.0, 0.0], [1.0, 0.0, 0.0]])

    # three pairs spanning J = -1 (complementary), 0, +1 (repulsive)
    j = _pair_j(["A", "B", "A"], ["B", "C", "C"], jm, index)
    assert np.allclose(j, [-1.0, 0.0, 1.0]), j
    # crystal-like: the complementary pair sits tightest, the repulsive one loosest
    rho = _spearman(5.0 - np.array([3.0, 4.0, 4.9]), -j)
    assert rho == 1.0, f"tight+favourable must give c_local = +1, got {rho}"
    # forced: the same chemistry with the distance order inverted
    rho_forced = _spearman(5.0 - np.array([4.9, 4.0, 3.0]), -j)
    assert rho_forced == -1.0, f"inverted pose must give c_local = -1, got {rho_forced}"
    # an unknown residue must not silently score as zero
    assert np.isnan(_pair_j(["A"], ["X"], jm, index)[0])
    assert np.isnan(_spearman(np.array([1.0, 2.0]), np.array([1.0, 2.0])))  # n < 3
    assert np.isnan(_spearman(np.ones(5), np.arange(5.0)))                  # constant axis

    # double-centring: rows and columns of J sum to zero, for an ASYMMETRIC matrix too
    class _P:
        _m = np.array([[1.0, 2.0, 9.0], [3.0, 0.0, 1.0], [5.0, 4.0, 2.0]])
        def as_matrix(self):
            return self._m, {"A": 0, "B": 1, "C": 2}

    J, _ = _double_centred(_P())
    assert np.allclose(J.sum(axis=0), 0) and np.allclose(J.sum(axis=1), 0), J
    # a nan cell must not poison the whole row
    class _Pn(_P):
        _m = np.array([[1.0, 2.0, np.nan], [3.0, 0.0, 1.0], [5.0, 4.0, 2.0]])

    Jn, _ = _double_centred(_Pn())
    assert np.isnan(Jn[0, 2]) and np.isfinite(Jn[0, 0])
    print("pose selfcheck ok")


if __name__ == "__main__":
    _selfcheck()
