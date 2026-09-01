"""Forces, stiffness and the kinetics proxies: the contact map as something that can break.

Everything here is mechanical rather than thermodynamic -- a stiffness in N/m, a rupture force in N,
a work in J, a margin in Angstrom. No potential enters `springs` or `stability`; the network is
built from the geometry alone, which is what makes the off-rate proxy independent of the energy
channel it is compared against."""
from __future__ import annotations

# `import *` skips underscore names, and two of them are part of what callers reach for: the
# coupling tally the descriptor row needs and the pull direction the rupture test pins.

# The submodules by name, plus their private helpers -- see the note in any of the deprecated
# top-level shims for why the underscore names travel too.

# The submodules by name, and every name they define -- see the note in any of the deprecated
# top-level shims for why `import *` is not enough.
from . import springs, stability, dynamics  # noqa: E402,F401

for _m in (springs, stability, dynamics):
    globals().update({k: v for k, v in vars(_m).items() if not k.startswith("__")})
del _m
