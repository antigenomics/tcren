"""Deprecated location; moved to :mod:`tcren.topology.footprint` on 2026-09-01.

Re-exported so every existing import keeps working. See :mod:`tcren.topology` for why the module now
lives where it does.
"""
from __future__ import annotations

from .topology import footprint as _moved

# Everything, not `import *`. Two things `import *` would drop, and both were reached for within
# an hour of the move: underscore names, which the sibling modules and the numeric tests use, and
# public names left out of the module's own `__all__`. A shim that is not transparent is worse
# than no shim, because the failure surfaces far from the cause.
globals().update({k: v for k, v in vars(_moved).items() if not k.startswith("__")})
del _moved
