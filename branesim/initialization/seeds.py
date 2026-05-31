"""Axis-aligned colour-seed initial conditions for VSH soliton candidates.

Background
----------
These seeds implement the Phase 2 baryon-candidate menu from
``paper/baryon_simulation_roadmap.md`` (BAYRON_SIMULATION_ROADMAP.md).
The defining idea is **colour = axis-alignment**: for the locked seeds the
lateral displacement component ``i`` is locked to the spatial radial unit
vector ``x_hat_i``, which is the VSH-L=1 / hedgehog SU(3) colour structure
(backbone #20).

All seeds operate on the spacelike-slice coordinate frame described by
``SpacelikeLattice.reference_positions(m)``:

  ref   -- shape (n_nodes, m_ambient)      flat lattice positions
  coords -- ref[:, :dim]                    spatial (in-brane) coordinates
  centre -- box centre                      = mean(coords)
  dx    -- coords - centre                  displacement from centre
  r     -- ||dx||                           radial distance
  x_hat -- dx / max(r, eps)                 unit radial vector (zero at r=0)

Lateral (colour) components: indices 0..dim-1 of the ambient vector.
Amplitude (X4) channel:      index ``dim`` (requires m >= dim+1).

Return value
------------
Every public function returns ``R0 = ref + displacement`` with shape
``(n_nodes, m_ambient)`` and dtype float64.  No energy term is evaluated,
no force is applied — these are pure initial displacements.

Principles compliance (principles.md §2, §3, §7.4, §7.6)
---------------------------------------------------------
- No solver state is touched (Layer C: initialization only).
- No extra energy / force / damping is introduced.
- No clamps or saturation rules.
- All core arithmetic is dimension-agnostic: ``dim = len(grid_shape)``; loops
  over ``range(dim)``; no hard-coded 3-D indices.
- dim=1 edge case: seeds that use the amplitude channel (component ``dim``)
  require ``m >= dim+1 = 2``; hedgehog and axis_triplet work at m=dim as well
  (amplitude channel displacement is zero).
"""

from __future__ import annotations

import math
from typing import Sequence

import numpy as np

from branesim.core.lattice import SpacelikeLattice


# ---------------------------------------------------------------------------
# Radial profile helpers (module-private)
# ---------------------------------------------------------------------------

_EPS = 1e-30   # guard against exact r=0 division


def _profile_gaussian(r: np.ndarray, w: float) -> np.ndarray:
    """f(r) = exp(-(r/w)^2).  f(0)=1, f(inf)=0, finite everywhere."""
    return np.exp(-((r / w) ** 2))


def _profile_sech(r: np.ndarray, w: float) -> np.ndarray:
    """f(r) = 1/cosh(r/w).  f(0)=1, f(inf)=0, finite everywhere."""
    return 1.0 / np.cosh(r / w)


def _profile_power2(r: np.ndarray, w: float) -> np.ndarray:
    """f(r) = 1/(1+(r/w)^2).  f(0)=1, f(inf)=0, finite everywhere."""
    return 1.0 / (1.0 + (r / w) ** 2)


_RADIAL_PROFILES: dict[str, object] = {
    "gaussian": _profile_gaussian,
    "sech": _profile_sech,
    "power2": _profile_power2,
}

_SKYRME_PROFILES_F: dict[str, object] = {
    "power2": lambda r, w, **_: math.pi / (1.0 + (r / w) ** 2),
    "tanh": lambda r, w, s=3.0, **_: math.pi * (1.0 - np.tanh(s * (r / w - 1.0))) / 2.0,
}


def _radial_geometry(
    lattice: SpacelikeLattice, m: int
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return (ref, dx, r, x_hat) for the given lattice and ambient dimension.

    Parameters
    ----------
    lattice : SpacelikeLattice
    m : int
        Ambient dimension m.

    Returns
    -------
    ref   : (n_nodes, m) float64 — flat reference positions
    dx    : (n_nodes, dim) float64 — displacement from box centre
    r     : (n_nodes,) float64 — Euclidean distance from centre
    x_hat : (n_nodes, dim) float64 — unit radial vector (zero at r=0 node)
    """
    dim = lattice.params.dim
    ref = lattice.reference_positions(m)                  # (n_nodes, m)
    coords = ref[:, :dim]                                 # (n_nodes, dim)
    centre = coords.mean(axis=0)                          # (dim,)
    dx = coords - centre                                  # (n_nodes, dim)
    r = np.linalg.norm(dx, axis=1)                        # (n_nodes,)
    # Zero at exact-centre node (odd grid); safe division elsewhere.
    r_safe = np.where(r > _EPS, r, 1.0)                  # avoid /0
    x_hat = dx / r_safe[:, None]                          # (n_nodes, dim)
    x_hat[r <= _EPS] = 0.0                                # zero at centre
    return ref, dx, r, x_hat


# ---------------------------------------------------------------------------
# Metadata helpers
# ---------------------------------------------------------------------------

def _metadata_hedgehog(profile_shape: str) -> dict:
    return {
        "ansatz": "hedgehog",
        "J": 0,
        "L": 1,
        "B_winding": 0,
        "colour_structure": "locked_to_x_hat",
        "profile_shape": profile_shape,
    }


def _metadata_skyrme_twisted(profile_shape: str) -> dict:
    return {
        "ansatz": "skyrme_twisted_hedgehog",
        "J": 0,
        "L": 1,
        "B_winding": 1,
        "colour_structure": "locked_to_x_hat",
        "profile_shape": profile_shape,
    }


def _metadata_axis_triplet(profile_shape: str) -> dict:
    return {
        "ansatz": "axis_triplet",
        "J": 1,
        "L": 0,
        "B_winding": 0,
        "colour_structure": "unlocked_axis_scalars",
        "role": "negative_control",
        "profile_shape": profile_shape,
    }


# ---------------------------------------------------------------------------
# Public seed functions
# ---------------------------------------------------------------------------


def hedgehog(
    lattice: SpacelikeLattice,
    m: int,
    u0: float,
    w: float,
    profile_shape: str = "gaussian",
) -> tuple[np.ndarray, dict]:
    """Hedgehog seed — colour locked to radial direction (J=0, L=1, B=0).

    Displacement:
        xi^i = u0 * f(r) * x_hat^i    for i in 0..dim-1
        xi^j = 0                       for j >= dim   (amplitude channel zero)

    ``f`` is the selected radial profile (``gaussian``, ``sech``, ``power2``).
    The unit radial vector ``x_hat`` locks colour index ``i`` to spatial
    direction ``x^i``; this is the canonical SU(3) hedgehog structure
    (backbone #20, VSH J=0, L=1 partial wave).

    At the exact-centre node (only on odd grids) ``x_hat = 0`` so the
    displacement is zero there.

    Parameters
    ----------
    lattice : SpacelikeLattice
    m : int
        Ambient dimension.  Must be >= dim.
    u0 : float
        Peak amplitude at r=0.
    w : float
        Profile half-width in lattice units.
    profile_shape : str
        One of ``"gaussian"``, ``"sech"``, ``"power2"``.

    Returns
    -------
    R0 : ndarray, shape (n_nodes, m), float64
        Full spacelike-slice configuration = ref + displacement.
    meta : dict
        Metadata: ansatz, J, L, B_winding, colour_structure, profile_shape.
    """
    if profile_shape not in _RADIAL_PROFILES:
        raise ValueError(
            f"Unknown profile_shape {profile_shape!r}; "
            f"choose from {sorted(_RADIAL_PROFILES)}"
        )
    if m < lattice.params.dim:
        raise ValueError(
            f"m ({m}) must be >= dim ({lattice.params.dim})"
        )
    dim = lattice.params.dim
    ref, dx, r, x_hat = _radial_geometry(lattice, m)     # (n_nodes, dim)
    f = _RADIAL_PROFILES[profile_shape](r, w)             # (n_nodes,)  scalar

    disp = np.zeros_like(ref)
    # Lateral channels: colour locked to x_hat
    disp[:, :dim] = (u0 * f[:, None]) * x_hat            # (n_nodes, dim)
    # Amplitude channel (index dim and above) stays zero.

    return ref + disp, _metadata_hedgehog(profile_shape)


def skyrme_twisted_hedgehog(
    lattice: SpacelikeLattice,
    m: int,
    u0: float,
    w: float,
    profile_shape: str = "power2",
    tanh_steepness: float = 3.0,
) -> tuple[np.ndarray, dict]:
    """Skyrme-twisted hedgehog — pi_3 winding B=1 (J=0, L=1, B=1).

    Displacement (Skyrme Ansatz on S^3):
        xi^i   = u0 * sin(F(r)) * x_hat^i     i in 0..dim-1
        xi^dim = u0 * cos(F(r))                amplitude (X4) channel

    Boundary conditions on F:
        F(0) = pi    -> sin=0, cos=-1  (south pole of S^3)
        F(inf) = 0   -> sin=0, cos=+1  (north pole of S^3)
        F(w) = pi/2

    As r scans 0 -> inf, the map (sin(F)*x_hat, cos(F)) wraps the spatial
    R^3 once around S^3, realising a topological winding number B=1 (the
    pi_3(S^3) = Z generator; Skyrme/Hopf-class stabilisation).

    Requires ``m >= dim + 1`` so the amplitude (X4) channel exists.

    Profile shapes for F:
        ``power2``  F(r) = pi / (1 + (r/w)^2)
        ``tanh``    F(r) = pi * (1 - tanh(steepness*(r/w - 1))) / 2

    Parameters
    ----------
    lattice : SpacelikeLattice
    m : int
        Ambient dimension.  Must be >= dim+1.
    u0 : float
        S^3 radius (amplitude of the full (xi, dX4) vector at every point).
    w : float
        Profile half-width; F(w) = pi/2 by construction.
    profile_shape : str
        ``"power2"`` or ``"tanh"``.
    tanh_steepness : float
        Steepness parameter ``s`` for the tanh profile (default 3).

    Returns
    -------
    R0 : ndarray, shape (n_nodes, m), float64
    meta : dict
    """
    if profile_shape not in _SKYRME_PROFILES_F:
        raise ValueError(
            f"Unknown profile_shape {profile_shape!r} for Skyrme-twisted; "
            f"choose from {sorted(_SKYRME_PROFILES_F)}"
        )
    if m < lattice.params.dim + 1:
        raise ValueError(
            f"Skyrme-twisted hedgehog requires m >= dim+1 = {lattice.params.dim+1}; "
            f"got m={m}"
        )
    dim = lattice.params.dim
    ref, dx, r, x_hat = _radial_geometry(lattice, m)

    F_fn = _SKYRME_PROFILES_F[profile_shape]
    F = F_fn(r, w, s=tanh_steepness)                     # (n_nodes,)

    disp = np.zeros_like(ref)
    # Lateral channels: sin(F) * x_hat  (colour locked)
    disp[:, :dim] = (u0 * np.sin(F))[:, None] * x_hat
    # Amplitude channel: cos(F)
    disp[:, dim] = u0 * np.cos(F)

    return ref + disp, _metadata_skyrme_twisted(profile_shape)


def axis_triplet(
    lattice: SpacelikeLattice,
    m: int,
    u0: float,
    w: float,
    weights: Sequence[float] | None = None,
    profile_shape: str = "gaussian",
) -> tuple[np.ndarray, dict]:
    """Axis-triplet seed — NEGATIVE CONTROL (J=1, L=0, B=0).

    Each lateral component ``i`` is an independent isotropic scalar bump,
    NOT locked to the radial direction:

        xi^i = u0 * weights[i] * g(r)    for i in 0..dim-1
        xi^j = 0                         for j >= dim

    Because the colour index is *not* locked to ``x_hat``, this seed has
    a different radial-lock metric than the hedgehog.  The expected value
    for the radial lock
        radial_lock = sum_nodes (xi . x_hat)^2 / sum_nodes |xi|^2
    is 1/dim for an isotropic weights vector (each component contributes
    equally but only its projection onto x_hat is counted), compared to 1.0
    for perfectly locked seeds.  This is the falsification discriminator.

    Parameters
    ----------
    lattice : SpacelikeLattice
    m : int
        Ambient dimension.  Must be >= dim.
    u0 : float
        Overall amplitude scale.
    w : float
        Profile half-width.
    weights : sequence of float, length dim, optional
        Per-axis amplitude multipliers.  Defaults to [1.0, ...] (isotropic).
    profile_shape : str
        One of ``"gaussian"``, ``"sech"``, ``"power2"``.

    Returns
    -------
    R0 : ndarray, shape (n_nodes, m), float64
    meta : dict
    """
    if profile_shape not in _RADIAL_PROFILES:
        raise ValueError(
            f"Unknown profile_shape {profile_shape!r}; "
            f"choose from {sorted(_RADIAL_PROFILES)}"
        )
    if m < lattice.params.dim:
        raise ValueError(
            f"m ({m}) must be >= dim ({lattice.params.dim})"
        )
    dim = lattice.params.dim

    if weights is None:
        weights = [1.0] * dim
    weights = np.asarray(weights, dtype=np.float64)
    if weights.shape != (dim,):
        raise ValueError(
            f"weights must have length dim={dim}; got shape {weights.shape}"
        )

    ref, dx, r, x_hat = _radial_geometry(lattice, m)
    g = _RADIAL_PROFILES[profile_shape](r, w)             # (n_nodes,)

    disp = np.zeros_like(ref)
    # Each lateral component is an independent scalar: weight[i] * g(r)
    # NOT multiplied by x_hat[i] — that is the key structural difference.
    for i in range(dim):
        disp[:, i] = u0 * weights[i] * g                  # isotropic, not locked

    return ref + disp, _metadata_axis_triplet(profile_shape)
