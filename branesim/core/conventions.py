"""Physical conventions, units, and closed-form dispersion relations.

All canonical constants and the closed-form dynamical-matrix eigenvalues for
the ``2*dim``-neighbor axial-only lattice. Dimension-agnostic: works for
``dim`` in {1, 2, 3} (and beyond).

Settled physics (ARCHITECTURE.md §3.3, discrete_4d_brane_action.md §4):
  - alpha := rest_length / spacing   (default 0.2; alpha=0 max prestress)
  - dimensionless units: k_s = a = rho = 1  →  m = rho * a**dim = 1
  - c_L^2 = k_s * a^2 / m
  - c_T^2 = (1 - alpha) * k_s * a^2 / m
  - Dynamical-matrix eigenvalues (spacelike, axial-only stencil):
      omega_a^2(k) = (2 k_s / rho) * [alpha * h_a + (1-alpha) * sum_b h_b]
      h_i := 1 - cos(k_i * a)
    Eigenvectors are the Cartesian unit vectors at every k and alpha
    (backbone #15, Sprint-2 #9).

No torch/numpy import at module level — pure arithmetic helpers that can be
called from tests or diagnostics without a GPU context.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Sequence

import numpy as np


# ---------------------------------------------------------------------------
# Parameter containers
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class LatticeParams:
    """Spacelike lattice geometry.

    Parameters
    ----------
    grid_shape : tuple of int
        Number of nodes along each spacelike axis.  len(grid_shape) == dim.
    spacing : float
        Lattice spacing ``a`` (dimensionless units: default 1.0).
    periodic_axes : tuple of bool
        Whether each axis wraps periodically.  Must match len(grid_shape).
    axial_weight : float
        Uniform per-link weight for the 2*dim axial bonds (default 1.0).
    """

    grid_shape: tuple[int, ...]
    spacing: float = 1.0
    periodic_axes: tuple[bool, ...] = field(default_factory=tuple)
    axial_weight: float = 1.0

    def __post_init__(self) -> None:
        if not self.grid_shape:
            raise ValueError("grid_shape must be non-empty")
        if any(n < 1 for n in self.grid_shape):
            raise ValueError("All grid dimensions must be >= 1")
        if self.spacing <= 0.0:
            raise ValueError("spacing must be positive")
        if self.axial_weight <= 0.0:
            raise ValueError("axial_weight must be positive")
        # Fill periodic_axes with False if not supplied
        if not self.periodic_axes:
            object.__setattr__(
                self, "periodic_axes", tuple(False for _ in self.grid_shape)
            )
        if len(self.periodic_axes) != len(self.grid_shape):
            raise ValueError("periodic_axes length must equal len(grid_shape)")

    @property
    def dim(self) -> int:
        return len(self.grid_shape)

    @property
    def n_nodes(self) -> int:
        n = 1
        for s in self.grid_shape:
            n *= s
        return n


@dataclass(frozen=True)
class ActionParams:
    """Parameters for the discrete brane action.

    Parameters
    ----------
    k_s : float
        Spring constant (dimensionless units: default 1.0).
    alpha : float
        Rest-length ratio ``alpha = rest_length / spacing`` (default 0.2).
    rho : float
        Mass density (dimensionless units: default 1.0).
    dt : float
        Temporal step size ``Delta_t``.
    n_slices : int
        Number of timelike slices ``N`` (world-volume has slices 0..N).
    m_ambient : int or None
        Number of ambient components (default dim+1; canonical d=3 → m=4).
    temporal_model : str
        ``"a"`` (zero-rest-length kinetic) or ``"b"`` (central-force temporal
        spring).  Default ``"a"`` — the validated IVP model.
    r_t : float
        Temporal rest length for model ``"b"``; must be 0.0 for model ``"a"``.
    """

    k_s: float = 1.0
    alpha: float = 0.2
    rho: float = 1.0
    dt: float = 0.1
    n_slices: int = 10
    m_ambient: int | None = None
    temporal_model: str = "a"
    r_t: float = 0.0

    def __post_init__(self) -> None:
        if not (0.0 <= self.alpha <= 1.0):
            raise ValueError(f"alpha must be in [0, 1]; got {self.alpha}")
        if self.dt <= 0.0:
            raise ValueError("dt must be positive")
        if self.n_slices < 1:
            raise ValueError("n_slices must be >= 1")
        if self.temporal_model not in ("a", "b"):
            raise ValueError("temporal_model must be 'a' or 'b'")
        if self.temporal_model == "a" and self.r_t != 0.0:
            raise ValueError("r_t must be 0.0 for temporal_model='a'")

    def mass(self, lattice: LatticeParams) -> float:
        return self.rho * lattice.spacing ** lattice.dim

    def ambient_dim(self, spatial_dim: int) -> int:
        """Ambient component count m (default d+1)."""
        return self.m_ambient if self.m_ambient is not None else spatial_dim + 1


# ---------------------------------------------------------------------------
# Closed-form speeds
# ---------------------------------------------------------------------------


def c_longitudinal(k_s: float = 1.0, a: float = 1.0, m: float = 1.0) -> float:
    """Longitudinal (L) light-cone speed: c_L = sqrt(k_s * a^2 / m).

    In dimensionless units k_s=a=m=1 → c_L=1.
    """
    return math.sqrt(k_s * a * a / m)


def c_transverse(
    k_s: float = 1.0, a: float = 1.0, m: float = 1.0, alpha: float = 0.2
) -> float:
    """Transverse (T) light-cone speed: c_T = sqrt((1-alpha) * k_s * a^2 / m).

    At alpha=0.2, c_T/c_L = sqrt(0.8) = 0.8944.
    """
    return math.sqrt((1.0 - alpha) * k_s * a * a / m)


def speed_ratio(alpha: float) -> float:
    """c_T / c_L = sqrt(1 - alpha).  Regression target at alpha=0.2: 0.8944."""
    if not (0.0 <= alpha <= 1.0):
        raise ValueError(f"alpha must be in [0, 1]; got {alpha}")
    return math.sqrt(1.0 - alpha)


# ---------------------------------------------------------------------------
# Dynamical-matrix eigenvalues (closed form, dimension-agnostic)
# ---------------------------------------------------------------------------


def d_of_k_eigenvalues(
    k: Sequence[float] | np.ndarray,
    alpha: float,
    k_s: float = 1.0,
    rho: float = 1.0,
    a: float = 1.0,
) -> np.ndarray:
    """Closed-form spacelike dynamical-matrix eigenvalues for the axial-only stencil.

    For a ``dim``-dimensional cubic lattice with only 2*dim axial bonds the
    stiffness matrix is diagonal in the Cartesian basis at every ``k`` and
    ``alpha`` (backbone #15, Sprint-2 #9).  The eigenvalues are:

        omega_a^2(k) = (2 k_s / rho) * [alpha * h_a + (1-alpha) * sum_b h_b]
        h_i = 1 - cos(k_i * a)

    Parameters
    ----------
    k : array_like, shape (dim,)
        Wavevector in inverse-length units.
    alpha : float
        Prestretch ratio (0 = max prestress, 1 = no prestress).
    k_s, rho, a : float
        Spring constant, mass density, lattice spacing.

    Returns
    -------
    omega_sq : ndarray, shape (dim,)
        Squared eigenfrequencies in axis order.  Eigenvectors are the
        Cartesian unit vectors (independent of k and alpha).
    """
    k_arr = np.asarray(k, dtype=float).ravel()
    dim = len(k_arr)
    h = 1.0 - np.cos(k_arr * a)  # shape (dim,)
    h_sum = float(h.sum())
    prefactor = 2.0 * k_s / rho
    return prefactor * (alpha * h + (1.0 - alpha) * h_sum)
