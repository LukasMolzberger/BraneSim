"""branesim.electron.electron_initialization

Electron initialization for a 3D brane embedded in 4D.

Design goal
-----------
Seed an *electron-like* localized excitation as a **rotating standing wave**
defined by a spherical-harmonic angular pattern and a spherical-Bessel radial
profile in *reference/material coordinates*.

We represent the (complex) spherical harmonic mode using **two real quadratures**
in a single embedding component (typically X^4 = positions[:, 3]):

    u(x, t) = Re( R(r) Y_lm(θ,φ) e^{-i ω t} )
            = Re(RY) cos(ωt) + Im(RY) sin(ωt)

Therefore, at t=0:
    u(x,0) = Re(RY)
    v(x,0) = ∂_t u(x,0) = ω Im(RY)

This yields a real field pattern u(x,t)=A(r,θ) cos(m φ - ω t) for m≠0,
which corresponds to a *phase pattern rotating* around the symmetry axis.

Containment
-----------
We optionally add a static X^4 deformation “well” w(r) (also written into the
same embedding component) intended as a first-pass containment scaffold.
The actual self-consistent equilibrium between wave energy and deformation is a
nonlinear coupled problem and generally needs numerical relaxation; this module
only provides an initializer.

Notes
-----
* This initializer is intended for **rest-length/reference-metric** simulations
  (fixed neighbor graph / fixed material Laplacian). The spherical construction
  happens in reference space.
* Uses SciPy for spherical harmonics and spherical Bessel functions.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Tuple, Dict, Any

import numpy as np

from scipy.special import sph_harm, spherical_jn
from scipy.optimize import brentq

import torch

from branesim.core.state import BraneState, Dimensionality
from branesim.core.grid import BraneGrid


@dataclass(frozen=True)
class ElectronModeSpec:
    """Specification for an electron spherical-harmonic mode (3D brane only).

    Parameters are expressed in **simulation units** unless noted.

    Attributes
    ----------
    l, m : int
        Spherical harmonic degree and order.
    n : int
        Radial mode index (nth positive zero of j_l). n=1 is the first zero.
    radius : float
        Effective cavity radius a (sim units). The radial wave number is
        k = α_{l,n} / a where α_{l,n} is the nth zero of j_l.
    amplitude : float
        Peak wave displacement amplitude (sim length units).
    center : (float, float, float)
        Center of the electron mode in reference coordinates (sim units).
    field_component : int
        Embedding component used for the electron field (default: 3 => X^4).
    wave_speed : float
        Wave speed used to set ω = c * k (default: 1.0 in sim units).

    containment_depth : float
        Depth of an added static containment well in the same component.
        Set to 0 to disable.
    containment_sigma : float
        Gaussian sigma for containment well. If None, uses radius/2.

    smooth_edge : float
        Smooth edge thickness for the radial cutoff window. If <=0, uses
        a hard mask r<=radius.
    """

    l: int = 1
    m: int = 1
    n: int = 1
    radius: float = 20.0
    amplitude: float = 0.5
    center: Tuple[float, float, float] = (0.0, 0.0, 0.0)
    field_component: int = 3
    wave_speed: float = 1.0

    containment_depth: float = 0.2
    containment_sigma: Optional[float] = None
    smooth_edge: float = 2.0


def _spherical_bessel_zero(l: int, n: int) -> float:
    """Return the nth positive root of the spherical Bessel j_l(x).

    Uses robust bracketing + Brent root-finding.
    """
    if l < 0:
        raise ValueError("l must be >= 0")
    if n <= 0:
        raise ValueError("n must be >= 1")

    # Heuristic: zeros are roughly spaced ~π. Start scanning intervals.
    # For l>0, the first zero is > 0.
    roots = []
    x_left = 1e-6
    step = np.pi
    max_scan = 20000  # plenty for small l,n

    def f(x: float) -> float:
        return spherical_jn(l, x)

    x = x_left
    scanned = 0
    while len(roots) < n and scanned < max_scan:
        x_next = x + step
        # Ensure we bracket a sign change; if not, keep scanning.
        f1 = f(x)
        f2 = f(x_next)
        if np.isnan(f1) or np.isnan(f2):
            x = x_next
            scanned += 1
            continue
        if f1 == 0.0:
            roots.append(x)
        elif f1 * f2 < 0:
            root = float(brentq(f, x, x_next, maxiter=200, xtol=1e-12))
            roots.append(root)
        x = x_next
        scanned += 1

    if len(roots) < n:
        raise RuntimeError(f"Failed to find {n} roots for j_{l}. Found {len(roots)}.")
    return roots[n - 1]


def _compute_spherical_coordinates(
    coords_xyz: np.ndarray,
    center: Tuple[float, float, float],
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Compute (r, theta, phi) from Euclidean xyz coordinates.

    theta is polar angle in [0, π], phi is azimuth in [0, 2π).
    """
    x = coords_xyz[:, 0] - center[0]
    y = coords_xyz[:, 1] - center[1]
    z = coords_xyz[:, 2] - center[2]

    r = np.sqrt(x * x + y * y + z * z)
    # Avoid division by zero
    r_safe = np.maximum(r, 1e-12)
    theta = np.arccos(np.clip(z / r_safe, -1.0, 1.0))
    phi = np.arctan2(y, x)
    phi = np.mod(phi, 2.0 * np.pi)
    return r, theta, phi


def initialize_electron_mode_3d(
    state: BraneState,
    grid: BraneGrid,
    spec: ElectronModeSpec,
    *,
    set_velocities: bool = True,
    normalize: bool = True,
    return_debug: bool = False,
) -> Optional[Dict[str, Any]]:
    """Initialize an electron-like spherical harmonic mode on a 3D brane.

    This function writes to `state` in-place.

    Parameters
    ----------
    state : BraneState
        Must be a 3D brane state (Dimensionality.THREE_D).
        `state.initialize_flat_configuration(...)` must have been called.
    grid : BraneGrid
        Must correspond to the state grid (same shape and spacing).
    spec : ElectronModeSpec
        Mode and geometry parameters in simulation units.
    set_velocities : bool
        If True, sets v(x,0)=ω Im(RY) to seed a rotating phase pattern.
    normalize : bool
        If True, scales the wave component so that max|u_wave| == spec.amplitude
        within the cutoff window.
    return_debug : bool
        If True, return a dictionary with computed scalars (k, ω, etc.).

    Returns
    -------
    Optional[dict]
        Debug info if requested.
    """
    if state.dimension != Dimensionality.THREE_D:
        raise ValueError("initialize_electron_mode_3d requires a 3D brane state")
    if grid.dimension != Dimensionality.THREE_D:
        raise ValueError("initialize_electron_mode_3d requires a 3D brane grid")
    if state.grid_shape != grid.grid_shape:
        raise ValueError("state.grid_shape and grid.grid_shape must match")

    # Reference-space coordinates in simulation units.
    # grid.get_spatial_coordinates() is already scaled by grid.spacing.
    coords_t = grid.get_spatial_coordinates().detach().cpu().to(torch.float64)
    coords = coords_t.numpy()  # (N, 3)

    r, theta, phi = _compute_spherical_coordinates(coords, spec.center)

    # Radial wave number from spherical Bessel boundary condition.
    alpha_ln = _spherical_bessel_zero(spec.l, spec.n)
    k = alpha_ln / max(spec.radius, 1e-12)
    omega = spec.wave_speed * k

    # Radial profile
    R = spherical_jn(spec.l, k * r)  # regular at r=0

    # Angular profile (complex spherical harmonic)
    Y = sph_harm(spec.m, spec.l, phi, theta)  # complex64/128

    # Combine
    mode_real = R * np.real(Y)
    mode_imag = R * np.imag(Y)

    # Cutoff window
    if spec.smooth_edge is not None and spec.smooth_edge > 0:
        # Smooth Heaviside-like window ~1 inside, ~0 outside
        delta = float(spec.smooth_edge)
        window = 0.5 * (1.0 - np.tanh((r - spec.radius) / max(delta, 1e-12)))
    else:
        window = (r <= spec.radius).astype(np.float64)

    mode_real *= window
    mode_imag *= window

    # Optional normalization so amplitude has intuitive meaning.
    scale = 1.0
    if normalize:
        max_abs = float(np.max(np.abs(mode_real)))
        if max_abs < 1e-12:
            max_abs = 1.0
        scale = spec.amplitude / max_abs
        mode_real *= scale
        mode_imag *= scale

    # Containment well (static deformation) in same component
    containment = 0.0
    if spec.containment_depth is not None and abs(spec.containment_depth) > 0:
        sigma = spec.containment_sigma
        if sigma is None:
            sigma = 0.5 * spec.radius
        sigma = max(float(sigma), 1e-12)
        containment = -float(spec.containment_depth) * np.exp(-(r * r) / (2.0 * sigma * sigma))
        containment *= window

    # Build displacement and velocity fields
    N = state.num_points
    u = torch.zeros((N, 4), device=state.device, dtype=state.dtype)
    v = torch.zeros((N, 4), device=state.device, dtype=state.dtype)

    u_comp = torch.from_numpy((containment + mode_real).astype(np.float32)).to(state.device, dtype=state.dtype)
    u[:, spec.field_component] = u_comp

    if set_velocities:
        v_comp = torch.from_numpy((omega * mode_imag).astype(np.float32)).to(state.device, dtype=state.dtype)
        v[:, spec.field_component] = v_comp

    # Apply via init API (keeps rest_positions intact)
    state.set_kinematics(u=u, v=v)

    # Respect fixed boundary conditions if present
    if hasattr(state, "fixed_mask") and state.fixed_mask is not None:
        state.apply_fixed_boundaries()

    if not return_debug:
        return None

    return {
        "alpha_ln": alpha_ln,
        "k": k,
        "omega": omega,
        "scale": scale,
        "radius": spec.radius,
        "l": spec.l,
        "m": spec.m,
        "n": spec.n,
    }
