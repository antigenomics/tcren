"""Everything measured in kT: sums of a pair potential over a contact set, and differences of them.

`scoring` owns the interface energy sum and the per-contact weighting the others use; `ddg` is the
change on mutation, virtually or on rebuilt coordinates; `rotamers` averages the contact map over
side-chain states so one modelled rotamer does not decide the answer. The potentials themselves are
:mod:`tcren.potential` and :mod:`tcren.potts`, one layer down -- this package APPLIES a Hamiltonian,
it does not define one."""
from __future__ import annotations



from .mutation import *  # noqa: F401,F403


# The submodules by name, plus their private helpers -- see the note in any of the deprecated
# top-level shims for why the underscore names travel too.

# The submodules by name, and every name they define -- see the note in any of the deprecated
# top-level shims for why `import *` is not enough.
from . import scoring, mutation, rotamers  # noqa: E402,F401

for _m in (scoring, mutation, rotamers):
    globals().update({k: v for k, v in vars(_m).items() if not k.startswith("__")})
del _m
