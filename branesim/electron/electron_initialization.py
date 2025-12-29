"""electron.electron_initialization

Electron initialization for a 3D brane embedded in 4D.

This initializer is built for the **reference-metric / rest-length elasticity**
family used in this project:

- The spherical construction happens in **reference/material coordinates**
  (fixed grid, fixed neighbor graph / fixed material Laplacian).
- We seed a *localized rotating standing wave* using a spherical-harmonic angular
  pattern and a spherical-Bessel radial profile.

Rotating phase field via spatial polarization
--------------------------------------------
A real-valued simulator needs a quadrature pair to represent a rotating phase.
We create a complex mode

    ψ(x) = R(r) Y_{l m}(θ, φ)

and map it into a **circularly polarized vector displacement** on the embedding
space, using two orthonormal polarization vectors p1, p2 (a plane):

    u(x,0) = scale * ( Re(ψ) p1 + Im(ψ) p2 )
    v(x,0) = scale * ω * ( Im(ψ) p1 - Re(ψ) p2 )

This corresponds to a local rotation in the (p1,p2)-plane with angular frequency
ω = c * k, and for m≠0 it produces the desired azimuthal phase winding.

Containment scaffold (X^4)
--------------------------
Separately, we can add a *static* Gaussian "well" in a chosen embedding
component (default: X^4) as a first-pass containment scaffold:

    w(r) = -depth * exp( -r^2 / (2 σ^2) )

A fully self-consistent equilibrium between wave leakage and geometric
containment is a coupled nonlinear problem (field + deformation) and will
generally require numerical relaxation/optimization; this module only provides
an initializer.

Dependencies
------------
Uses SciPy:
- scipy.special.sph_harm
- scipy.special.spherical_jn
- scipy.optimize.brentq
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
    """Specification for an electron spherical-harmonic seed (3D brane only).

    All lengths are in **simulation units**.

    Parameters
    ----------
    l, m : int
        Spherical harmonic degree and order.
    n : int
        Radial mode index: nth positive zero of j_l. n=1 is the first zero.
    radius : float
        Effective cavity radius a (sim units). Radial wave number:
        k = α_{l,n} / a.
    amplitude : float
        Target peak magnitude of the oscillatory displacement vector (sim units).
    center : (float, float, float)
        Center of the mode in reference coordinates.
    wave_speed : float
        Wave speed used to set ω = c * k in sim units.

    polarization : str
        Which embedding-plane to use for the rotating mode.
        - "xy", "xz", "yz": simple spatial planes
        - "spatial": a plane spanning all spatial dims (x,y,z)
        - "spatial_x4": a plane spanning (x,y,z) AND X^4 (both quadratures touch X^4)
        - "all": a plane that also includes X^4

    polarization_p1 / polarization_p2 : Optional tuples
        Explicit polarization vectors (length = embedding dimension). If set,
        they override `polarization`.

    containment_component : int
        Embedding component receiving the static containment well (default: 3 => X^4).

    containment_depth : float
        Depth of the containment well (0 disables).
    containment_sigma : Optional[float]
        Gaussian sigma; if None uses 0.5*radius.

    smooth_edge : float
        Smooth edge thickness for radial cutoff. If <=0, use hard mask r<=radius.

    Backward compatibility
    ----------------------
    field_component : Optional[int]
        Deprecated alias used in the earlier patch; treated as containment_component.
    """

    l: int = 1
    m: int = 1
    n: int = 1

    radius: float = 20.0
    amplitude: float = 0.5
    center: Tuple[float, float, float] = (0.0, 0.0, 0.0)
    wave_speed: float = 1.0

    polarization: str = "spatial"
    polarization_p1: Optional[Tuple[float, ...]] = None
    polarization_p2: Optional[Tuple[float, ...]] = None

    containment_component: int = 3
    containment_depth: float = 0.2
    containment_sigma: Optional[float] = None

    smooth_edge: float = 2.0

    # Deprecated: do not use for new code
    field_component: Optional[int] = None


def _orthonormalize_plane(p1: np.ndarray, p2: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Return orthonormal (p1,p2) via Gram–Schmidt."""
    p1 = np.asarray(p1, dtype=np.float64)
    p2 = np.asarray(p2, dtype=np.float64)

    n1 = float(np.linalg.norm(p1))
    if n1 == 0.0:
        raise ValueError("polarization_p1 must not be zero")
    p1 = p1 / n1

    p2 = p2 - float(np.dot(p2, p1)) * p1
    n2 = float(np.linalg.norm(p2))
    if n2 == 0.0:
        raise ValueError("polarization_p2 must not be collinear with polarization_p1")
    p2 = p2 / n2

    return p1, p2


def _default_plane(embed_dim: int, mode: str) -> Tuple[np.ndarray, np.ndarray]:
    p1 = np.zeros(embed_dim, dtype=np.float64)
    p2 = np.zeros(embed_dim, dtype=np.float64)

    if mode == "xy":
        p1[0] = 1.0
        p2[1] = 1.0
    elif mode == "xz":
        p1[0] = 1.0
        p2[2] = 1.0
    elif mode == "yz":
        p1[1] = 1.0
        p2[2] = 1.0
    elif mode == "spatial":
        if embed_dim < 3:
            raise ValueError("polarization='spatial' requires embedding dim >= 3")
        # A plane that touches all spatial axes
        p1[0:3] = 1.0
        p2[0] = 1.0
        p2[1] = -1.0
    elif mode == "spatial_x4":
        # Like 'spatial', but explicitly couples into X^4 at t=0.
        # This is useful if you want X^4 to carry a *wave* (not just a static well)
        # and to have nonzero momentum right away.
        if embed_dim < 4:
            return _default_plane(embed_dim, "spatial")
        p1[0:3] = 1.0
        p1[3] = 1.0
        p2[0] = 1.0
        p2[1] = -1.0
        p2[3] = 1.0
    elif mode == "all":
        if embed_dim < 4:
            return _default_plane(embed_dim, "spatial")
        p1[:] = 1.0
        p2[0] = 1.0
        p2[1] = -1.0
    else:
        raise ValueError(f"Unknown polarization mode: {mode!r}")

    return _orthonormalize_plane(p1, p2)


def _spherical_bessel_zero(l: int, n: int) -> float:
    """Return the nth positive root of the spherical Bessel j_l(x)."""
    if l < 0:
        raise ValueError("l must be >= 0")
    if n <= 0:
        raise ValueError("n must be >= 1")

    # Heuristic: zeros are spaced ~π. Scan intervals until we find n roots.
    roots = []
    x = 1e-6
    step = np.pi

    def f(t: float) -> float:
        return float(spherical_jn(l, t))

    # Safety cap
    for _ in range(40000):
        if len(roots) >= n:
            break
        x_next = x + step
        f1 = f(x)
        f2 = f(x_next)
        if f1 == 0.0:
            roots.append(x)
        elif f1 * f2 < 0.0:
            roots.append(float(brentq(f, x, x_next, maxiter=200, xtol=1e-12)))
        x = x_next

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

    Writes to `state` in-place via `state.set_kinematics(u=..., v=...)`.

    The oscillatory part is seeded as a *rotating vector mode* (polarization plane)
    so that the electron is not confined to X^4 at t=0.

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

    embed_dim = int(state.positions.shape[1])

    # Back-compat: treat field_component as containment_component
    containment_component = int(spec.containment_component)
    if spec.field_component is not None:
        containment_component = int(spec.field_component)

    # Reference-space coordinates in simulation units.
    coords_t = grid.get_spatial_coordinates().detach().cpu().to(torch.float64)
    coords = coords_t.numpy()  # (N, 3)

    r, theta, phi = _compute_spherical_coordinates(coords, spec.center)

    # Radial wave number from spherical Bessel boundary condition.
    alpha_ln = _spherical_bessel_zero(spec.l, spec.n)
    k = alpha_ln / max(float(spec.radius), 1e-12)
    omega = float(spec.wave_speed) * k

    # Radial profile (regular at r=0)
    R = spherical_jn(spec.l, k * r)

    # Angular profile (complex spherical harmonic)
    Y = sph_harm(spec.m, spec.l, phi, theta)

    # Combine complex scalar mode
    psi = R * Y
    Re = np.real(psi)
    Im = np.imag(psi)

    # Cutoff window
    if spec.smooth_edge is not None and spec.smooth_edge > 0:
        delta = float(spec.smooth_edge)
        window = 0.5 * (1.0 - np.tanh((r - float(spec.radius)) / max(delta, 1e-12)))
    else:
        window = (r <= float(spec.radius)).astype(np.float64)

    Re *= window
    Im *= window

    # Choose polarization plane
    if spec.polarization_p1 is not None and spec.polarization_p2 is not None:
        p1 = np.asarray(spec.polarization_p1, dtype=np.float64)
        p2 = np.asarray(spec.polarization_p2, dtype=np.float64)
        if p1.shape[0] != embed_dim or p2.shape[0] != embed_dim:
            raise ValueError("polarization_p1/p2 must have length == embedding dimension")
        p1, p2 = _orthonormalize_plane(p1, p2)
    else:
        p1, p2 = _default_plane(embed_dim, spec.polarization)

    # Normalize by peak magnitude |u| = sqrt(Re^2 + Im^2) (because p1,p2 are orthonormal)
    scale = 1.0
    if normalize:
        mag = np.sqrt(Re * Re + Im * Im)
        max_mag = float(np.max(np.abs(mag)))
        if max_mag < 1e-12:
            max_mag = 1.0
        scale = float(spec.amplitude) / max_mag
        Re *= scale
        Im *= scale

    # Containment well (static deformation) in containment_component (default X^4)
    containment = None
    if spec.containment_depth is not None and abs(spec.containment_depth) > 0:
        sigma = spec.containment_sigma
        if sigma is None:
            sigma = 0.5 * float(spec.radius)
        sigma = max(float(sigma), 1e-12)
        containment = -float(spec.containment_depth) * np.exp(-(r * r) / (2.0 * sigma * sigma))
        containment *= window

    # Build displacement and velocity fields (u,v) in embedding space
    # u = Re*p1 + Im*p2
    # v = ω*(Im*p1 - Re*p2)
    N = state.num_points
    u_np = (Re[:, None] * p1[None, :]) + (Im[:, None] * p2[None, :])
    v_np = omega * ((Im[:, None] * p1[None, :]) - (Re[:, None] * p2[None, :]))

    u = torch.zeros((N, embed_dim), device=state.device, dtype=state.dtype)
    v = torch.zeros((N, embed_dim), device=state.device, dtype=state.dtype)

    u += torch.from_numpy(u_np.astype(np.float32)).to(state.device, dtype=state.dtype)

    if set_velocities:
        v += torch.from_numpy(v_np.astype(np.float32)).to(state.device, dtype=state.dtype)

    if containment is not None:
        c = torch.from_numpy(containment.astype(np.float32)).to(state.device, dtype=state.dtype)
        if not (0 <= containment_component < embed_dim):
            raise ValueError(f"containment_component={containment_component} out of range for embed_dim={embed_dim}")
        u[:, containment_component] = u[:, containment_component] + c

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
        "radius": float(spec.radius),
        "l": int(spec.l),
        "m": int(spec.m),
        "n": int(spec.n),
        "polarization": spec.polarization,
        "p1": p1,
        "p2": p2,
        "containment_component": containment_component,
    }
