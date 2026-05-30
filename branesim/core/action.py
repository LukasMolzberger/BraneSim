"""Discrete brane action: spacelike potential V, kinetic T, total S, and force F.

Settled physics (ARCHITECTURE.md §1.2, discrete_4d_brane_action.md §2):

    V^l = (k_s / 4) * sum_p sum_{delta in N_s} (|R_{p+delta}^l - R_p^l| - alpha*a)^2

    T^{l+1/2} = sum_p (m/2) * ((R_p^{l+1} - R_p^l) / dt)^2

    S[R] = sum_l dt * (T^{l+1/2} - V^l)

    F_p^l = -dV^l/dR_p^l
           = k_s * sum_{delta in N_s} (|R_{p+delta} - R_p| - alpha*a) * hat(R_{p+delta} - R_p)

The full Euclidean NORM of the link vector is used (geometric nonlinearity).

Minimum-image convention for periodic boundaries
-------------------------------------------------
The position array stores the full (unwrapped) node positions.  For nodes
at periodic boundaries, the raw ``pos[q] - pos[p]`` difference may be
``(N-1)*a`` in magnitude instead of the correct ``-a``.  We apply the
minimum-image convention on the ``dim`` spacelike axes: for each periodic
axis with box length ``L_i = N_i * a``, shift the delta by ``-L_i`` if
``delta > L_i/2`` and by ``+L_i`` if ``delta < -L_i/2``.  The ``m_ambient``
components beyond axis ``dim-1`` (the temporal/extra-ambient ones) are left
unchanged.

Non-negotiables enforced here:
- No clamp/cutoff on the norm (principles §3.2, §4, non-negotiable #4).
- No back-reaction from diagnostics (action/force depend only on R).
"""

from __future__ import annotations

import numpy as np

from branesim.core.conventions import ActionParams
from branesim.core.lattice import SpacelikeLattice


# ---------------------------------------------------------------------------
# Minimum-image delta helper
# ---------------------------------------------------------------------------


def _minimum_image_delta(
    raw_delta: np.ndarray,
    dim: int,
    periodic_axes: tuple[bool, ...],
    box_lengths: np.ndarray,
) -> np.ndarray:
    """Apply minimum-image convention on the first ``dim`` components.

    Parameters
    ----------
    raw_delta : ndarray, shape (n_valid, m_ambient)
        Raw ``pos[q] - pos[p]`` vectors.
    dim : int
        Number of spacelike axes.
    periodic_axes : tuple of bool, length ``dim``
    box_lengths : ndarray, shape (dim,)
        Box length ``N_i * a`` along each axis.

    Returns
    -------
    delta : ndarray, shape (n_valid, m_ambient)
        Delta with minimum-image shift applied to periodic axes.
    """
    delta = raw_delta.copy()
    for axis in range(dim):
        if not periodic_axes[axis]:
            continue
        L = box_lengths[axis]
        # shift: if delta > L/2 subtract L; if delta < -L/2 add L
        delta[:, axis] -= L * np.round(delta[:, axis] / L)
    return delta


# ---------------------------------------------------------------------------
# Spacelike potential V^l and force F_p^l on a single slice
# ---------------------------------------------------------------------------


def spacelike_potential(
    positions: np.ndarray,
    lattice: SpacelikeLattice,
    params: ActionParams,
) -> float:
    """Spacelike potential energy V^l for a single spacelike slice.

    V^l = (k_s / 4) * sum_p sum_{delta} (|R_{p+delta} - R_p| - alpha*a)^2

    The factor k_s/4 avoids double-counting each bond: each directed link
    (p→q) and its reverse (q→p) both appear in the sum, but represent the
    same spring.  Combined with the two endpoints that share each spring, the
    prefactor 1/4 is the standard convention matching the derivation in
    discrete_4d_brane_action.md §2.

    Parameters
    ----------
    positions : ndarray, shape (n_nodes, m_ambient)
        Node positions on the slice.
    lattice : SpacelikeLattice
        Precomputed neighbor table.
    params : ActionParams
        Action parameters (k_s, alpha).

    Returns
    -------
    float
        Total spacelike potential energy V^l.
    """
    k_s = params.k_s
    alpha_a = params.alpha * lattice.params.spacing
    dim = lattice.dim
    periodic_axes = lattice.params.periodic_axes
    box_lengths = np.array(
        [n * lattice.params.spacing for n in lattice.params.grid_shape],
        dtype=np.float64,
    )

    V = 0.0
    for nb_idx in range(lattice.n_neighbors):
        nb_ids = lattice.neighbors[:, nb_idx]
        valid = nb_ids >= 0
        if not np.any(valid):
            continue

        p_pos = positions[valid]
        q_pos = positions[nb_ids[valid]]
        raw_delta = q_pos - p_pos

        delta = _minimum_image_delta(raw_delta, dim, periodic_axes, box_lengths)
        dist = np.linalg.norm(delta, axis=1)

        strain = dist - alpha_a
        V += np.sum(strain ** 2)

    return (k_s / 4.0) * V


def spacelike_force(
    positions: np.ndarray,
    lattice: SpacelikeLattice,
    params: ActionParams,
) -> np.ndarray:
    """Spacelike force F_p^l = -dV^l/dR_p^l on a single slice.

    F_p = k_s * sum_{delta in N_s} (|R_{p+delta} - R_p| - alpha*a)
              * (R_{p+delta} - R_p) / |R_{p+delta} - R_p|

    Uses the minimum-image convention for periodic boundaries.

    The unit-vector computation ``delta / dist`` performs no epsilon guard
    by design (ARCHITECTURE §3.2): the central-force law is non-coercive and
    callers are responsible for ensuring dist > 0.  Adding an epsilon guard
    would silently alter forces at near-zero separations and mask bugs.

    Parameters
    ----------
    positions : ndarray, shape (n_nodes, m_ambient)
    lattice : SpacelikeLattice
    params : ActionParams

    Returns
    -------
    forces : ndarray, shape (n_nodes, m_ambient)
        Net force on each node from all spacelike springs.
    """
    k_s = params.k_s
    alpha_a = params.alpha * lattice.params.spacing
    dim = lattice.dim
    periodic_axes = lattice.params.periodic_axes
    box_lengths = np.array(
        [n * lattice.params.spacing for n in lattice.params.grid_shape],
        dtype=np.float64,
    )

    forces = np.zeros_like(positions)

    for nb_idx in range(lattice.n_neighbors):
        nb_ids = lattice.neighbors[:, nb_idx]
        valid = nb_ids >= 0
        if not np.any(valid):
            continue

        valid_ids = np.where(valid)[0]
        p_pos = positions[valid_ids]
        q_pos = positions[nb_ids[valid_ids]]

        raw_delta = q_pos - p_pos
        delta = _minimum_image_delta(raw_delta, dim, periodic_axes, box_lengths)

        # Full Euclidean norm — the geometric nonlinearity
        dist = np.linalg.norm(delta, axis=1, keepdims=True)

        # No clamp: trust the physics (principles §3.2, non-negotiable #4)
        strain = dist - alpha_a
        # ARCHITECTURE §3.2: central-force law is non-coercive; no eps guard
        # by design (callers ensure dist>0).
        unit = delta / dist

        link_force = k_s * strain * unit
        np.add.at(forces, valid_ids, link_force)

    return forces


# ---------------------------------------------------------------------------
# Kinetic energy T^{l+1/2}
# ---------------------------------------------------------------------------


def kinetic_energy(
    R_l: np.ndarray,
    R_l1: np.ndarray,
    mass: float,
    dt: float,
) -> float:
    """Kinetic energy T^{l+1/2} between slices l and l+1 (model a).

    T^{l+1/2} = sum_p (m/2) * ((R_p^{l+1} - R_p^l) / dt)^2

    The temporal link has zero rest length (model a).
    The squared norm is the full ambient dot product.

    Parameters
    ----------
    R_l : ndarray, shape (n_nodes, m_ambient)  — slice l
    R_l1 : ndarray, shape (n_nodes, m_ambient) — slice l+1
    mass : float  — node mass m
    dt : float    — temporal step Delta_t

    Returns
    -------
    float
        Scalar kinetic energy T^{l+1/2}.
    """
    vel = (R_l1 - R_l) / dt
    return float(0.5 * mass * np.sum(vel ** 2))


# ---------------------------------------------------------------------------
# Total action S[R] over a world-volume
# ---------------------------------------------------------------------------


def action(
    world: np.ndarray,
    lattice: SpacelikeLattice,
    params: ActionParams,
    mass: float,
) -> float:
    """Discrete Lorentzian brane action S[R] over the full world-volume.

    S[R] = sum_l dt * (T^{l+1/2} - V^l)

    IMPORTANT: S is Lorentzian and therefore a SADDLE — unbounded below.
    The solver must root-find grad S = 0, not minimize S.
    (ARCHITECTURE.md §1.3, discrete_4d_brane_action.md §5)

    Parameters
    ----------
    world : ndarray, shape (N+1, n_nodes, m_ambient)
        World-volume stack; world[l] is slice l.
    lattice : SpacelikeLattice
    params : ActionParams
    mass : float

    Returns
    -------
    float
        Total discrete action S.
    """
    N_plus_1 = world.shape[0]
    dt = params.dt

    S = 0.0
    for l in range(N_plus_1 - 1):
        T = kinetic_energy(world[l], world[l + 1], mass, dt)
        V = spacelike_potential(world[l], lattice, params)
        S += dt * (T - V)
    return S
