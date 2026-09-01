"""The shape of the contact set, and of the surface it sits on. No energy anywhere.

Every quantity here is a count, a share, an entropy, a Betti number or an Angstrom -- nothing in this
package loads a potential, which is what makes these descriptors independent of the energy channel
they are reported beside. `footprint` partitions the contacts and measures how evenly they spread;
`graph` reads the same contact map as a bipartite graph and the Calpha/Cbeta maps as matrices;
`surface` is the pMHC face the receptor actually meets; `pose` asks whether the tight contacts and
the favourable chemistry are the same contacts."""
from __future__ import annotations



# The submodules by name, plus their private helpers -- see the note in any of the deprecated
# top-level shims for why the underscore names travel too.

# The submodules by name, and every name they define -- see the note in any of the deprecated
# top-level shims for why `import *` is not enough.
from . import footprint, graph, surface, pose  # noqa: E402,F401

for _m in (footprint, graph, surface, pose):
    globals().update({k: v for k, v in vars(_m).items() if not k.startswith("__")})
del _m
