"""Descriptors: the catalogue, the computation, and the batch dispatch, in three layers.

============================  ==========================================================
:mod:`~tcren.descriptors.catalogue`  what every column is -- names, families, invariance
                              classes, units, definitions and known defects. Pure data.
:mod:`~tcren.descriptors.compute`    structure -> values for the interface block, and the
                              calls out to the energetics, topology, potts and kinetics
                              modules that own the rest.
:mod:`~tcren.descriptors.table`      a whole structure set -> one row each, with the single
                              batched annotation pass and the process pool.
============================  ==========================================================

:mod:`tcren.recognition` re-exports all three under the name every caller already uses.
"""
from .catalogue import *  # noqa: F401,F403
from .catalogue import (  # noqa: F401
    CDR3_FRAME_FEATURES, DESCRIPTORS, DETAIL, FAMILIES, FULL_FEATURES,
    INTERFACE_SYMMETRY_FEATURES, INVARIANCE, INVARIANCE_CLASSES, PEPTIDE_INTERNAL_FEATURES,
    RECOGNITION_FEATURES, STATUS, TCR_PLACEMENT_FEATURES, descriptors,
)
from .compute import (  # noqa: F401
    _burial, _cdr3_frame_features, _chain_balance, _extent, _interface_symmetry,
    _peptide_internal_columns, _placement_columns, _stability_clash_columns, _symmetry_columns,
    recognition_features,
)
from .table import _featurise_families, _featurise_one, recognition_table  # noqa: F401
