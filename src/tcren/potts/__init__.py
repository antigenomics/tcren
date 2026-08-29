"""A coupled Potts model over the TCR:pMHC contact map.

The configuration is the contact map itself. Sites are the residue pairs that *could* have
contacted — Cα within ``radius`` — and ``sigma_a = 1`` iff a heavy-atom contact formed::

    E(sigma) = - sum_a eta_a sigma_a - 1/2 sum_{a,b} A_ab sigma_a sigma_b
    P(sigma) = exp(-E(sigma)) / Z

``eta`` is an additive one-body model (fields for each side, a pair coupling ``J``, and the
geometry the fields must be adjusted for); ``A`` couples neighbouring cells of the map. Fitting is
penalised pseudolikelihood and needs no ``Z``; scoring gets ``Z`` by annealed importance sampling
from the uncoupled model, whose partition function is exact.

Typical use::

    from tcren.potts import PottsModel, available_pairs, score_sites, contact_probabilities

    pairs = available_pairs(structure)                     # the sites
    model = PottsModel.bundled()                           # or fit_potts(pairs)
    score_sites(pairs, model)                              # energy, log Z, likelihood
    contact_probabilities(pairs, model)                    # per-site contact probability
    contact_map(pairs, model, by="loop")                   # loop x peptide-position frequency map

The CLI mirrors it: ``tcren potts fit``, ``tcren potts score``, ``tcren potts contacts``,
``tcren potts map``.
"""

from .fit import DEFAULT_RIDGE, cluster_se, design, fit_potts, gauge, irls, kernel_table
from .kernel import (bucket_edges, colour, coupling_matrix, edges, neighbour_counts)
from .model import (AA, CDR_LOOPS, CLASSES, CROSS_DJ, DBIN, OFFSETS, REGIONS, ROLES,
                    PottsModel, centred_potential, kernel_names)
from .sample import (ais_log_z, count_free_energy, delta_f_empty, delta_f_threshold, energy,
                     exact_log_z, factorised_log_z, gibbs, mu_star, tilt_mean)
from .score import (bound_unbound, connected_correlations, contact_map, contact_probabilities,
                    peptide_free_energy,
                    count_profile, sample_maps, score_sites, score_structure)
from .sites import GROOVE_REGIONS, MHC_PARTNER, available_pairs, eta, site_codes

__all__ = [
    "AA", "CDR_LOOPS", "CLASSES", "CROSS_DJ", "DBIN", "OFFSETS", "REGIONS", "ROLES",
    "GROOVE_REGIONS", "MHC_PARTNER", "DEFAULT_RIDGE",
    "PottsModel", "centred_potential", "kernel_names",
    "available_pairs", "site_codes", "eta",
    "edges", "neighbour_counts", "bucket_edges", "colour", "coupling_matrix",
    "design", "irls", "gauge", "cluster_se", "fit_potts", "kernel_table",
    "gibbs", "ais_log_z", "exact_log_z", "factorised_log_z", "energy",
    "tilt_mean", "mu_star", "count_free_energy", "delta_f_empty", "delta_f_threshold",
    "bound_unbound", "count_profile",
    "score_sites", "contact_probabilities", "contact_map",
    "peptide_free_energy", "connected_correlations",
    "sample_maps",
    "score_structure",
]
