"""Deprecated name for :mod:`tcren.docking`; kept so existing imports keep working.

Renamed on 2026-09-01. The package holds a native reimplementation of TCRdock's rigid-body
parameterisation alongside the crossing and incident angles, the canonical frame and the
superposition helpers -- it is the docking geometry, and ``orient`` said only that something gets
moved. ``orient.docking`` is now :mod:`tcren.docking.angles`, because ``docking.docking`` stutters.
"""
from __future__ import annotations

from ..docking import *  # noqa: F401,F403
