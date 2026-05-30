"""Discrete brane residual: the shared physics primitive.

    R_p^l = m * (R_p^{l+1} - 2 R_p^l + R_p^{l-1}) / dt^2  -  F_p^l
           = -grad_{R_p^l} S

The residual is zero at every interior node of a valid world-volume
(ARCHITECTURE.md §1.3, discrete_4d_brane_action.md §3).

Key uses:
  1. Solver quality check: ‖R‖ ≈ 0 after IVP march (acceptance criterion 1).
  2. BVP root-find target: drive R → 0 over interior nodes (increment 2).
  3. Diagnostics: report ‖R‖ as solve quality metric.

The residual is computed matrix-free — no Jacobian is assembled.

Conventions:
  - ``world`` has shape ``(L+1, n_nodes, m_ambient)``; index 0 is slice l=0.
  - Interior slices: l = 1 .. L-1.
  - Boundary slices (l=0 and l=L): residual is set to zero (they are fixed).
  - The temporal difference uses model (a): zero-rest-length quadratic
    increment (default; model (b) extension is left for the BVP increment).
"""

from __future__ import annotations

import numpy as np

from branesim.core.action import spacelike_force
from branesim.core.conventions import ActionParams
from branesim.core.lattice import SpacelikeLattice


def residual(
    world: np.ndarray,
    lattice: SpacelikeLattice,
    params: ActionParams,
    mass: float,
) -> np.ndarray:
    """Compute the residual R at every node of the world-volume.

    R_p^l = m * (R^{l+1} - 2R^l + R^{l-1}) / dt^2  -  F_p^l

    Boundary slices (l=0 and l=L) are set to zero — they are not interior
    nodes and carry prescribed data.

    Parameters
    ----------
    world : ndarray, shape (L+1, n_nodes, m_ambient)
        World-volume; world[l] is the position array on slice l.
    lattice : SpacelikeLattice
        Spacelike neighbor topology.
    params : ActionParams
        Action parameters (k_s, alpha, dt).
    mass : float
        Node mass m = rho * a^dim.

    Returns
    -------
    res : ndarray, shape (L+1, n_nodes, m_ambient)
        Residual array.  res[0] = res[-1] = 0 (boundary slices).
        res[l] = 0 iff slice l satisfies the discrete EL equation.
    """
    n_slices_plus_1, n_nodes, m_ambient = world.shape
    L = n_slices_plus_1 - 1  # number of time steps (L+1 slices, indices 0..L)
    dt = params.dt
    dt2 = dt * dt

    res = np.zeros_like(world)

    for l in range(1, L):  # interior slices only
        # Discrete second temporal derivative (model a)
        accel = (world[l + 1] - 2.0 * world[l] + world[l - 1]) / dt2
        # Spacelike force on this slice
        F = spacelike_force(world[l], lattice, params)
        res[l] = mass * accel - F

    # Boundary slices stay zero — they carry fixed data.
    return res


def residual_norm(
    world: np.ndarray,
    lattice: SpacelikeLattice,
    params: ActionParams,
    mass: float,
) -> float:
    """Scalar L2 norm of the interior residual.

    Excludes boundary slices (l=0 and l=L) which are prescribed.
    Returns sqrt(sum of squared residual components over all interior nodes).
    """
    res = residual(world, lattice, params, mass)
    # res[0] and res[-1] are zero by construction; sum includes them harmlessly
    return float(np.sqrt(np.sum(res ** 2)))
