"""Electron initialization: rotating spherical-harmonic cavity seed (3D brane embedded in 4D).

This module provides a deterministic initializer that writes an *electron-like* rotating
standing-wave mode into a BraneState, without any additional "containment scaffold".

Key construction (real-valued simulation)
----------------------------------------
We start from a complex scalar cavity mode

    ψ(x) = R(r) Y_{ℓm}(θ,φ),

where Y_{ℓm} is a complex spherical harmonic and R(r) is a spherical-Bessel radial profile
with a smooth cutoff at r≈a. We realize this complex phase as a circularly polarized
vector oscillation in embedding space using two orthonormal polarization vectors p1,p2:

    u(x,0) = Re(ψ) p1 + Im(ψ) p2
    v(x,0) = ω ( Im(ψ) p1 - Re(ψ) p2 )

This guarantees a rotating local phase (for m≠0) and *immediately* excites multiple
embedding components (including X^4 if the polarization plane includes it).

Important design decision
-------------------------
In the current rest-length spring force model, embedding components are coupled through
the Euclidean spring length, but *forces act component-wise*. If X^4 starts identically
zero everywhere, its springs have no X^4-extension and X^4 remains exactly zero.
Therefore, if we want X^4 to participate, we must include it in the polarization plane
(e.g. polarization='spatial_x4' or custom p1,p2). This is not a scaffold; it is the mode.

All lengths are in simulation units.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Tuple, Dict, Any

import numpy as np
import torch

try:
    from scipy.special import sph_harm, spherical_jn
    from scipy.optimize import brentq
except Exception as e:
    raise ImportError(
        "electron initialization requires SciPy (scipy.special, scipy.optimize). "
        "Install with: pip install scipy"
    ) from e

from branesim.core.state import BraneState, Dimensionality
from branesim.core.grid import BraneGrid


@dataclass(frozen=True)
class ElectronModeSpec:
    """Specification for an electron spherical-harmonic seed (3D brane only).

    Parameters
    ----------
    l, m : int
        Spherical harmonic indices (ℓ,m).
    n : int
        Radial index (nth positive zero of j_ℓ).
    radius : float
        Cavity radius a in simulation units.
    amplitude : float
        Target peak magnitude of |ψ| after normalization (simulation units).
    center : (float,float,float)
        Center of the mode in reference coordinates.
    wave_speed : float
        Wave speed in simulation units (typically c_sim=1). We set ω = wave_speed · k.
    polarization : str
        Polarization preset for the embedding-space plane. Common:
        - 'spatial'     : plane in X^1..X^3 only
        - 'spatial_x4'  : plane includes X^4 (recommended)
        - 'all'         : mixes all components (if embed_dim>=4)
        - 'xy','xz','yz': simple planes
    polarization_p1, polarization_p2 : Optional[Tuple[float,...]]
        Optional custom embedding-space vectors (length = embedding dimension).
    smooth_edge : float
        Smooth edge thickness for radial cutoff. If <=0, use hard mask r<=radius.
    """

    l: int = 1
    m: int = 1
    n: int = 1

    radius: float = 1.0
    amplitude: float = 0.25
    center: Tuple[float, float, float] = (0.0, 0.0, 0.0)
    wave_speed: float = 1.0

    polarization: str = "spatial_x4"
    polarization_p1: Optional[Tuple[float, ...]] = None
    polarization_p2: Optional[Tuple[float, ...]] = None

    smooth_edge: float = 2.0


def _orthonormalize_plane(p1: np.ndarray, p2: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    p1 = p1.astype(np.float64, copy=False)
    p2 = p2.astype(np.float64, copy=False)

    n1 = float(np.linalg.norm(p1))
    if n1 == 0.0:
        raise ValueError("polarization_p1 must not be the zero vector")
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
        p1[0:3] = 1.0
        p2[0] = 1.0
        p2[1] = -1.0
    elif mode == "spatial_x4":
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

    roots = []
    x = 1e-6
    step = float(np.pi)

    def f(t: float) -> float:
        return float(spherical_jn(l, t))

    for _ in range(100000):
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
    return float(roots[n - 1])


def _compute_spherical_coordinates(
    coords_xyz: np.ndarray,
    center: Tuple[float, float, float],
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Compute (r, theta, phi) from Euclidean xyz coordinates.

    theta is polar angle in [0, π], phi is azimuth in (-π, π].
    """
    x = coords_xyz[:, 0] - float(center[0])
    y = coords_xyz[:, 1] - float(center[1])
    z = coords_xyz[:, 2] - float(center[2])

    r = np.sqrt(x * x + y * y + z * z)
    r_safe = np.where(r > 0.0, r, 1.0)
    theta = np.arccos(np.clip(z / r_safe, -1.0, 1.0))
    phi = np.arctan2(y, x)
    return r, theta, phi


def _radial_window(r: np.ndarray, radius: float, smooth_edge: float) -> np.ndarray:
    """Smooth cutoff window ~1 for r<<radius, ~0 for r>>radius."""
    radius = max(float(radius), 1e-12)
    if smooth_edge <= 0.0:
        return (r <= radius).astype(np.float64)
    return 0.5 * (1.0 - np.tanh((r - radius) / (radius / float(smooth_edge))))


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
    """
    if state.dimension != Dimensionality.THREE_D:
        raise ValueError("initialize_electron_mode_3d requires a 3D brane state")
    if grid.dimension != Dimensionality.THREE_D:
        raise ValueError("initialize_electron_mode_3d requires a 3D brane grid")
    if state.grid_shape != grid.grid_shape:
        raise ValueError("state.grid_shape and grid.grid_shape must match")

    embed_dim = int(state.positions.shape[1])

    # Reference-space coordinates (simulation units)
    coords_t = grid.get_spatial_coordinates().detach().cpu().to(torch.float64)
    coords = coords_t.numpy()  # (N, 3)

    r, theta, phi = _compute_spherical_coordinates(coords, spec.center)

    # Radial wavenumber from cavity boundary j_l(k a)=0
    alpha_ln = _spherical_bessel_zero(int(spec.l), int(spec.n))
    k = float(alpha_ln) / max(float(spec.radius), 1e-12)
    omega = float(spec.wave_speed) * k

    window = _radial_window(r, spec.radius, spec.smooth_edge)
    R = spherical_jn(int(spec.l), k * r) * window

    # SciPy convention: sph_harm(m, n, theta, phi) where theta=azimuth, phi=polar
    Y = sph_harm(int(spec.m), int(spec.l), phi, theta)
    psi = R * Y
    Re = np.real(psi).astype(np.float64)
    Im = np.imag(psi).astype(np.float64)

    scale = 1.0
    if normalize:
        mag = np.sqrt(Re * Re + Im * Im)
        max_mag = float(np.max(mag))
        if max_mag < 1e-12:
            max_mag = 1.0
        scale = float(spec.amplitude) / max_mag
        Re *= scale
        Im *= scale

    # Polarization plane
    if spec.polarization_p1 is not None and spec.polarization_p2 is not None:
        p1 = np.asarray(spec.polarization_p1, dtype=np.float64)
        p2 = np.asarray(spec.polarization_p2, dtype=np.float64)
        if p1.shape[0] != embed_dim or p2.shape[0] != embed_dim:
            raise ValueError("polarization_p1/p2 must have length == embedding dimension")
        p1, p2 = _orthonormalize_plane(p1, p2)
    else:
        p1, p2 = _default_plane(embed_dim, spec.polarization)

    # u = Re*p1 + Im*p2
    # v = ω*(Im*p1 - Re*p2)
    u_np = (Re[:, None] * p1[None, :]) + (Im[:, None] * p2[None, :])
    v_np = omega * ((Im[:, None] * p1[None, :]) - (Re[:, None] * p2[None, :]))

    # Create torch tensors with state dtype/device.
    # Prefer generating numpy arrays in the matching precision to avoid extra copies.
    np_dtype = np.float32 if state.dtype == torch.float32 else np.float64

    u = torch.from_numpy(u_np.astype(np_dtype, copy=False)).to(device=state.device, dtype=state.dtype)
    if set_velocities:
        v = torch.from_numpy(v_np.astype(np_dtype, copy=False)).to(device=state.device, dtype=state.dtype)
    else:
        v = torch.zeros_like(u)

    state.set_kinematics(u=u, v=v)

    if not return_debug:
        return None

    return {
        "alpha_ln": float(alpha_ln),
        "k": float(k),
        "omega": float(omega),
        "scale": float(scale),
        "radius": float(spec.radius),
        "l": int(spec.l),
        "m": int(spec.m),
        "n": int(spec.n),
        "polarization": spec.polarization,
        "p1": p1,
        "p2": p2,
    }
