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
    energy_ensemble(pairs, model)                          # <Phi>, Var Phi, free energy, entropy
    map_ensemble(pairs, model)                             # the configurations, with probabilities
    threaded_ensemble(pairs, model, candidates)            # the same per candidate sequence

`available_pairs` needs a structure. :mod:`tcren.potts.geom` models the Ca distance matrix instead,
in length-normalised coordinates, and emits the same schema — so a receptor:peptide pair can be
placed from sequence and lengths alone, with no template of the right shape::

    prior = fit_geom_prior(structures)                     # or GeomPrior.from_json(...)
    predicted_sites(prior, loops, peptide)                 # the same columns available_pairs emits

That prior sees one geometry per crystal. :mod:`tcren.potts.cross` removes the scarcity: every
complex in ``Canonical2026`` sits in one frame, so a TCR from one crystal against a pMHC from
another is a pseudo-complex needing no transform, and within MHC class there are 79,061 of them
against 359 solved::

    rec, par = residue_tables(structures)                  # co-framed, chain-typed, MHC-annotated
    cross_pairs(rec, par)                                  # every pairing's Ca distances

The marginals settle the *mean* interaction energy and nothing else: `Phi` is linear in `sigma`, so
`<Phi> = sum_a p_a phi_a` exactly — but only with `p_model`, the coupled marginal, and only read at
the sequence being scored. :mod:`tcren.potts.ensemble` carries the read-outs a set of marginals
cannot determine, and re-reads the marginals per candidate.

The CLI mirrors it: ``tcren potts fit``, ``tcren potts score``, ``tcren potts contacts``,
``tcren potts map``, ``tcren potts ensemble``, ``tcren potts cross``.
"""

from .fit import DEFAULT_RIDGE, cluster_se, design, fit_potts, gauge, irls, kernel_table
from .kernel import (bucket_edges, colour, coupling_matrix, edges, neighbour_counts)
from .model import (AA, CDR_LOOPS, CLASSES, CROSS_DJ, DBIN, OFFSETS, REGIONS, ROLES,
                    PottsModel, centred_potential, kernel_names)
from .em import LATENT_MARGIN, consistency, fit_potts_em, impute, latent_mask
from .cross import (CROSS_COLUMNS, MIN_CA, PAR_COLUMNS, REC_COLUMNS, cross_from_structures, cross_pairs,
                    partner_residues, receptor_residues, residue_tables, write_cross_pairs)
from .geom import (BUNDLED_ENSEMBLE, BUNDLED_PRIOR, REGIONS_MODELLED, GeomPrior, feasible,
                   fit_geom_prior, loops_from_structure, predict_contact_map,
                   predicted_distances, predicted_ensemble, predicted_sites)
from .ensemble import energy_ensemble, map_ensemble, site_energies, threaded_ensemble
from .sample import (ais_log_z, count_free_energy, delta_f_empty, delta_f_threshold, energy,
                     exact_log_z, factorised_log_z, gibbs, mu_star, tilt_mean)
from .score import (bound_unbound, connected_correlations, contact_map, contact_probabilities,
                    map_agreement, peptide_free_energy,
                    count_profile, sample_maps, score_sites, score_structure)
from .selfcon import (LOGIT_CLIP, contact_posterior, fit_selfconsistent, observation_field,
                      weighted_contacts)
from .sites import (GROOVE_REGIONS, MHC_PARTNER, MHC_RECEPTOR_REGIONS, available_pairs, eta,
                    site_codes)

__all__ = [
    "AA", "CDR_LOOPS", "CLASSES", "CROSS_DJ", "DBIN", "OFFSETS", "REGIONS", "ROLES",
    "GROOVE_REGIONS", "MHC_PARTNER", "MHC_RECEPTOR_REGIONS", "DEFAULT_RIDGE",
    "PottsModel", "centred_potential", "kernel_names",
    "available_pairs", "site_codes", "eta",
    "edges", "neighbour_counts", "bucket_edges", "colour", "coupling_matrix",
    "design", "irls", "gauge", "cluster_se", "fit_potts", "kernel_table",
    "gibbs", "ais_log_z", "exact_log_z", "factorised_log_z", "energy",
    "tilt_mean", "mu_star", "count_free_energy", "delta_f_empty", "delta_f_threshold",
    "bound_unbound", "count_profile",
    "score_sites", "contact_probabilities", "map_agreement", "contact_map",
    "peptide_free_energy", "connected_correlations",
    "sample_maps",
    "score_structure",
    "energy_ensemble", "map_ensemble", "threaded_ensemble", "site_energies",
    "LOGIT_CLIP", "observation_field", "contact_posterior", "weighted_contacts",
    "fit_selfconsistent",
    "LATENT_MARGIN", "latent_mask", "impute", "fit_potts_em", "consistency",
    "BUNDLED_ENSEMBLE", "BUNDLED_PRIOR", "GeomPrior", "REGIONS_MODELLED", "feasible",
    "fit_geom_prior", "loops_from_structure", "predict_contact_map", "predicted_distances",
    "predicted_ensemble", "predicted_sites",
    "CROSS_COLUMNS", "MIN_CA", "PAR_COLUMNS", "REC_COLUMNS", "cross_from_structures", "cross_pairs",
    "partner_residues", "receptor_residues", "residue_tables", "write_cross_pairs",
]
