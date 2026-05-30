"""Forward Verlet initial-value problem (IVP) solver — the regression baseline.

Implements the Störmer–Verlet stencil (ARCHITECTURE.md §1.4 IVP mode):

    R^{l+1} = 2 R^l - R^{l-1} + (dt^2 / m) * F(R^l)

This is the exact discrete Euler–Lagrange equation of the brane action
(discrete_4d_brane_action.md §3) — Verlet IS the discrete variational
integrator of this action.

The IVP mode prescribes two past slices (l=0, l=1) as boundary data and
marches forward.  It is used only for regression; the block BVP solver
(solver/bvp.py, increment 2) is the foundational mode.

Boundary conditions for the spacelike slice are handled by the lattice
neighbor table (periodic or open per-axis); no additional spatial BC logic is
needed here.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from branesim.core.action import spacelike_force
from branesim.core.conventions import ActionParams, LatticeParams
from branesim.core.lattice import SpacelikeLattice


# ---------------------------------------------------------------------------
# Data containers
# ---------------------------------------------------------------------------


@dataclass
class IVPProblem:
    """Specification of a forward-IVP march.

    Parameters
    ----------
    lattice : SpacelikeLattice
        Precomputed spacelike neighbor topology.
    params : ActionParams
        Action parameters (k_s, alpha, dt, n_slices).
    mass : float
        Node mass m = rho * a^dim.
    R0 : ndarray, shape (n_nodes, m_ambient)
        Node positions on slice l=0.
    R1 : ndarray, shape (n_nodes, m_ambient)
        Node positions on slice l=1.
    """

    lattice: SpacelikeLattice
    params: ActionParams
    mass: float
    R0: np.ndarray
    R1: np.ndarray

    def __post_init__(self) -> None:
        if self.R0.shape != self.R1.shape:
            raise ValueError("R0 and R1 must have the same shape")
        n_nodes = self.lattice.n_nodes
        if self.R0.shape[0] != n_nodes:
            raise ValueError(
                f"R0 has {self.R0.shape[0]} rows but lattice has {n_nodes} nodes"
            )


@dataclass
class WorldVolume:
    """World-volume produced by a march or solve.

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


# ---------------------------------------------------------------------------
# Forward Verlet march
# ---------------------------------------------------------------------------


def march(problem: IVPProblem) -> WorldVolume:
    """Forward Verlet march: produce the world-volume from IVP boundary data.

    Uses the Störmer–Verlet stencil:
        R^{l+1} = 2 R^l - R^{l-1} + (dt^2 / m) * F(R^l)

    The march is equivalent to stepping the discrete EL equations forward in
    time — producing slices l=2..N from the two past slices l=0,1.

    Parameters
    ----------
    problem : IVPProblem
        Boundary data (R0, R1) and action/lattice parameters.

    Returns
    -------
    WorldVolume
        Full world-volume array of shape (N+1, n_nodes, m_ambient).
        slices[0] = R0, slices[1] = R1, slices[2..N] from Verlet.
    """
    import time as _time

    lattice = problem.lattice
    params = problem.params
    mass = problem.mass
    dt = params.dt
    dt2_over_m = (dt * dt) / mass
    N = params.n_slices  # total time steps; world-volume has N+1 slices

    n_nodes, m_ambient = problem.R0.shape
    world = np.empty((N + 1, n_nodes, m_ambient), dtype=np.float64)
    world[0] = problem.R0
    world[1] = problem.R1

    t0 = _time.perf_counter()

    R_prev = world[0]
    R_curr = world[1]

    for l in range(1, N):
        F = spacelike_force(R_curr, lattice, params)
        R_next = 2.0 * R_curr - R_prev + dt2_over_m * F
        world[l + 1] = R_next
        R_prev = R_curr
        R_curr = R_next

    elapsed = _time.perf_counter() - t0

    solver_report = {
        "mode": "ivp",
        "n_slices": N,
        "walltime_s": elapsed,
    }

    return WorldVolume(
        slices=world,
        params=params,
        lattice_params=lattice.params,
        solver_report=solver_report,
    )
