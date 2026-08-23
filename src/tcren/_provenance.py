"""Marking what is, and is not, part of the shipped TCRen2 derivation.

The package accumulated alternatives faster than it used them: three redundancy-weighting schemes,
a leave-one-out derivation, a second matrix from 2022. Exactly one combination produces the matrix
the paper reports, and a reader of the source cannot tell which. That is how three mutually
inconsistent potentials came to be in circulation at once.

:func:`not_in_tcren2` marks the rest. It is not ``@deprecated`` -- most of these are correct,
tested and worth keeping, and some were measured and rejected, which is a result rather than a
defect. It says only: **this is not how the shipped matrix is made**, with the reason. The recipe
itself is machine-readable in ``data/potentials.json`` and enforced by
``tests/regression/test_shipped_potentials.py``.

Grep ``not_in_tcren2`` for the full list.
"""

from __future__ import annotations

import inspect

from collections.abc import Callable
from typing import TypeVar

T = TypeVar("T", bound=Callable)

#: The exact recipe behind ``TCRen2_potential.csv``, as recorded in ``data/potentials.json``.
TCREN2_RECIPE = "tcren derive-potential --structure-dir Native2026 --balance both"


def not_in_tcren2(reason: str) -> Callable[[T], T]:
    """Record that a callable takes no part in producing the shipped TCRen2 matrix.

    Args:
        reason: One line saying what this is for instead, and -- where it was measured against
            TCRen2 and not adopted -- what the measurement showed.

    Returns:
        A decorator that prepends a Sphinx admonition to the docstring and sets
        ``__not_in_tcren2__``. It does not warn at runtime and does not wrap the call: these are
        supported entry points, they are just not the paper's.

    Example:
        >>> @not_in_tcren2("an alternative nobody ships")
        ... def f():
        ...     '''Do a thing.'''
        >>> f.__not_in_tcren2__
        'an alternative nobody ships'
    """
    def decorate(obj: T) -> T:
        note = (f".. admonition:: Not part of the TCRen2 derivation\n"
                f"   :class: caution\n\n"
                f"   {reason}\n\n"
                f"   The shipped matrix comes from ``{TCREN2_RECIPE}`` and nothing else.\n\n")
        # cleandoc, not a bare concat: the note is unindented and a function docstring's body is
        # indented by four, so prepending one to the other makes RST read the body as a block quote
        # -- which turns `**kwargs` into an unterminated strong start-string and fails the -W docs
        # build. cleandoc strips the common indent first, so the result is one flat block.
        obj.__doc__ = note + inspect.cleandoc(obj.__doc__ or "")
        obj.__not_in_tcren2__ = reason
        return obj
    return decorate


def marked() -> list[tuple[str, str]]:
    """Every callable carrying :func:`not_in_tcren2`, as ``(dotted name, reason)``.

    Walks the package, so a new marking appears here without anyone maintaining a list. This is
    what ``OBSOLETE.md`` is generated from and what the test that keeps it honest compares against.

    Example:
        >>> names = [n for n, _ in marked()]
        >>> "tcren.potential.redundancy.cluster_weights" in names
        True
    """
    import importlib
    import pkgutil

    import tcren

    found: dict[str, str] = {}
    for mod_info in pkgutil.walk_packages(tcren.__path__, "tcren."):
        try:
            mod = importlib.import_module(mod_info.name)
        except Exception:  # noqa: BLE001 -- an optional extra missing is not this function's problem
            continue
        for name, obj in vars(mod).items():
            if getattr(obj, "__module__", None) != mod_info.name:
                continue
            reason = getattr(obj, "__not_in_tcren2__", None)
            if reason:
                found[f"{mod_info.name}.{name}"] = reason
    return sorted(found.items())


if __name__ == "__main__":  # pragma: no cover -- regenerates OBSOLETE.md
    for _name, _reason in marked():
        print(f"- `{_name}` --- {_reason}")
