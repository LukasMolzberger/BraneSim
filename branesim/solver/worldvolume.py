"""World-volume data container — the canonical solved/seeded field stack.

A ``WorldVolume`` is the (N+1, n_nodes, m_ambient) stack of node positions over
the time loop, plus the action/lattice parameters and a solver report.  It is
produced by the block solver (``solver/bvp.py``), the breather eigen-solver
(``solver/breather.py``), or written directly by the initialization layer
(a seed worldvolume).

This module is intentionally solver-agnostic: it defines the data structure
only, so any layer can depend on it without pulling in a particular solver.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from branesim.core.conventions import ActionParams, LatticeParams


@dataclass
class WorldVolume:
    """World-volume produced by a solve (or written as a seed).

    Attributes
    ----------
    slices : ndarray, shape (N+1, n_nodes, m_ambient)
        World-volume stack; slices[l] is slice l.
    params : ActionParams
    lattice_params : LatticeParams
    solver_report : dict
        Telemetry from the solver (residual norms, timing, etc.).
    """

    slices: np.ndarray
    params: ActionParams
    lattice_params: LatticeParams
    solver_report: dict = field(default_factory=dict)

    @property
    def n_slices(self) -> int:
        return self.slices.shape[0] - 1  # number of time steps

    @property
    def n_nodes(self) -> int:
        return self.slices.shape[1]

    @property
    def m_ambient(self) -> int:
        return self.slices.shape[2]