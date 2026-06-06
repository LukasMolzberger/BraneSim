"""Discrete brane residual: the shared physics primitive.

    𝓡_p^l = temporal_force_p^l  −  F_spatial_p^l

The temporal force is derived from the same central-force spring law as the
spacelike links (½k_t(|ΔR|−r_t)²), parameterized by the temporal rest
length r_t:

Linear limit (r_t = 0):

    temporal_force = m * (R^{l+1} - 2R^l + R^{l-1}) / dt^2

    This is the exact zero-rest-length stencil; the residual becomes

    𝓡_p^l = m * (R_p^{l+1} - 2 R_p^l + R_p^{l-1}) / dt^2  -  F_p^l

    Fast explicit march path.  Used for small-amplitude wave/dispersion
    validation.  NOT a separate model — the r_t→0 limit of the spring law.

Temporal spring (r_t > 0, canonical substrate):

    With ΔR^+ = R_p^{l+1} - R_p^l,  ΔR^- = R_p^{l-1} - R_p^l,

    temporal_force = k_t * [
        (|ΔR^+| - r_t) * ΔR^+ / |ΔR^+|
      + (|ΔR^-| - r_t) * ΔR^- / |ΔR^-|
    ]

    Full m_ambient-component Euclidean norm (including the timelike component).
    No minimum-image on the temporal axis. No epsilon guard, no clamp —
    matches the spacelike_force policy (principles §3.2, non-negotiable #4).

    At r_t→0: (|ΔR| - 0) * ΔR / |ΔR| = ΔR, so

        temporal_force|_{r_t=0} = k_t * (ΔR^+ + ΔR^-)
                                 = k_t * (R^{l+1} - 2R^l + R^{l-1})

    With k_t = m/dt^2, this is term-for-term the linear-limit stencil.
    The regression gate requires ‖𝓡_spring(r_t=0) - 𝓡_linear‖ < 1e-12.

The residual is zero at every interior node of a valid world-volume
(ARCHITECTURE.md §1.3, discrete_4d_brane_action.md §3).

Key uses:
  1. Solver quality check: ‖𝓡‖ ≈ 0 after IVP march (acceptance criterion 1).
  2. BVP root-find target: drive 𝓡 → 0 over interior nodes (increment 2).
  3. Diagnostics: report ‖𝓡‖ as solve quality metric.

The residual is computed matrix-free — no Jacobian is assembled.

Conventions:
  - ``world`` has shape ``(L+1, n_nodes, m_ambient)``; index 0 is slice l=0.
  - Interior slices: l = 1 .. L-1.
  - Boundary slices (l=0 and l=L): residual is set to zero (they are fixed).
"""

from __future__ import annotations

import numpy as np

from branesim.core.action import spacelike_force
from branesim.core.conventions import ActionParams
from branesim.core.lattice import SpacelikeLattice


# ---------------------------------------------------------------------------
# Temporal force helpers
# ---------------------------------------------------------------------------


def _temporal_force_linear(
    world: np.ndarray,
    l: int,
    mass: float,
    dt: float,
) -> np.ndarray:
    """Temporal force at interior slice l — linear/Verlet limit (r_t = 0).

    temporal_force = m * (R^{l+1} - 2R^l + R^{l-1}) / dt^2

    This is the exact zero-rest-length stencil: the r_t→0 limit of the
    central-force temporal spring.  Bit-identical to the r_t=0 linear/Verlet path.

    Parameters
    ----------
    world : ndarray, shape (L+1, n_nodes, m_ambient)
    l : int
        Interior slice index (1 <= l <= L-1).
    mass : float
    dt : float

    Returns
    -------
    ndarray, shape (n_nodes, m_ambient)
    """
    dt2 = dt * dt
    return mass * (world[l + 1] - 2.0 * world[l] + world[l - 1]) / dt2


def _temporal_force_spring(
    world: np.ndarray,
    l: int,
    k_t: float,
    r_t: float,
) -> np.ndarray:
    """Temporal force at interior slice l — central-force spring (r_t > 0).

    With ΔR^+ = R_p^{l+1} - R_p^l  and  ΔR^- = R_p^{l-1} - R_p^l,

        temporal_force = k_t * [
            (|ΔR^+| - r_t) * ΔR^+ / |ΔR^+|
          + (|ΔR^-| - r_t) * ΔR^- / |ΔR^-|
        ]

    Full m_ambient-component Euclidean norm — includes the timelike component.
    No minimum-image on the temporal axis.
    No epsilon guard, no clamp (principles §3.2, non-negotiable #4).

    At r_t→0 this collapses term-for-term to the linear stencil with k_t=m/dt².

    Parameters
    ----------
    world : ndarray, shape (L+1, n_nodes, m_ambient)
    l : int
        Interior slice index (1 <= l <= L-1).
    k_t : float
        Temporal spring constant.
    r_t : float
        Temporal rest length (> 0 for the canonical prestressed substrate;
        the r_t=0 caller should use _temporal_force_linear instead).

    Returns
    -------
    ndarray, shape (n_nodes, m_ambient)
    """
    dR_plus = world[l + 1] - world[l]   # (n_nodes, m_ambient)
    dR_minus = world[l - 1] - world[l]  # (n_nodes, m_ambient)

    # Full ambient-norm (keepdims for broadcasting into the force sum)
    dist_plus = np.linalg.norm(dR_plus, axis=1, keepdims=True)   # (n_nodes, 1)
    dist_minus = np.linalg.norm(dR_minus, axis=1, keepdims=True)  # (n_nodes, 1)

    # Central-force spring: (|ΔR| - r_t) * ΔR / |ΔR|.
    # No eps guard — same deliberate policy as spacelike_force (principles §3.2).
    # In the canonical substrate the vacuum itself drifts (|ΔR| = r_t > 0), so a
    # temporal bond is never zero; a degenerate zero-drift input (R^{l+1}=R^l) is
    # a caller error, not silently guarded.
    force_plus = k_t * (dist_plus - r_t) * dR_plus / dist_plus
    force_minus = k_t * (dist_minus - r_t) * dR_minus / dist_minus

    return force_plus + force_minus


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def residual(
    world: np.ndarray,
    lattice: SpacelikeLattice,
    params: ActionParams,
    mass: float,
) -> np.ndarray:
    """Compute the residual 𝓡 at every node of the world-volume.

    Routes on ``params.r_t``:

    r_t == 0 (linear/Verlet limit — fast exact path):
        𝓡_p^l = m * (R^{l+1} - 2R^l + R^{l-1}) / dt^2  -  F_p^l

    r_t > 0 (canonical prestressed substrate):
        𝓡_p^l = k_t * [(|ΔR^+| - r_t) * ΔR^+/|ΔR^+|
                       + (|ΔR^-| - r_t) * ΔR^-/|ΔR^-|]  -  F_p^l

    Boundary slices (l=0 and l=L) are set to zero — they are not interior
    nodes and carry prescribed data.

    The spacelike force F_p^l is the same in both cases — it uses the
    existing spacelike_force() function from action.py.

    Parameters
    ----------
    world : ndarray, shape (L+1, n_nodes, m_ambient)
        World-volume; world[l] is the position array on slice l.
    lattice : SpacelikeLattice
        Spacelike neighbor topology.
    params : ActionParams
        Action parameters (k_s, alpha, dt, r_t, k_t).
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

    use_spring = (params.r_t > 0.0)
    if use_spring:
        k_t = params.resolved_k_t(mass)
        r_t = params.r_t

    res = np.zeros_like(world)

    for l in range(1, L):  # interior slices only
        # Temporal force — routes on r_t
        if use_spring:
            temp_force = _temporal_force_spring(world, l, k_t, r_t)
        else:
            temp_force = _temporal_force_linear(world, l, mass, dt)

        # Spacelike force on this slice
        F = spacelike_force(world[l], lattice, params)

        res[l] = temp_force - F

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
