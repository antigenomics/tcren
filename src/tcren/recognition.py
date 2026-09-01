"""Structure -> descriptors. The public name; the implementation is :mod:`tcren.descriptors`.

Split into three layers on 2026-09-01 because this module had grown to 1,151 lines doing three
different jobs -- holding the catalogue, computing the interface block, and dispatching a batch --
and a caller that only wanted to ask what a column means was importing all of it. The names are
unchanged and every existing import keeps working:

* :mod:`tcren.descriptors.catalogue` -- ``DESCRIPTORS``, ``INVARIANCE``, ``DETAIL``, ``STATUS``,
  ``FAMILIES``, :func:`descriptors`. Data and selection, no arithmetic.
* :mod:`tcren.descriptors.compute` -- :func:`recognition_features` and the interface terms this
  package computes itself; it *calls* :mod:`tcren.pipeline`, :mod:`tcren.footprint`,
  :mod:`tcren.potts` and :mod:`tcren.mechanics` for the other families rather than reimplementing
  them.
* :mod:`tcren.descriptors.table` -- :func:`recognition_table`, the batched set-level path.

New code may import from the three directly; nothing is deprecated here.
"""
from __future__ import annotations

from .descriptors import *  # noqa: F401,F403
from .descriptors import (  # noqa: F401  - the private names the package's own modules reach for
    _burial,
    _cdr3_frame_features,
    _chain_balance,
    _extent,
    _featurise_families,
    _featurise_one,
    _interface_symmetry,
    _peptide_internal_columns,
    _placement_columns,
    _stability_clash_columns,
    _symmetry_columns,
)
from .descriptors.catalogue import (  # noqa: F401
    _CDR3_FRAME_KEYS,
    _CT_TYPES,
    _EPS,
    _FAMILY_ALIASES,
    _TCR_TYPES,
)
