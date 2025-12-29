# electron/electron_initialization.py
from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Tuple, Dict, Any

import numpy as np

try:
    from scipy.special import sph_harm, spherical_jn
    from scipy.optimize import brentq
except Exception as e:
    raise ImportError(
        "electron initialization requires SciPy (scipy.special, scipy.optimize). "
        "Install with: pip install scipy"
    ) from e


@dataclass
class ElectronModeSpec:
    # Angular mode
    l: int = 1
    m: int = 1
    n: int = 1  # radial index (nth zero)

    # Spatial support / scale in reference coordinates
    radius: float = 1.0
    center: Tuple[float, float, float] = (0.0, 0.0, 0.0)

    # Overall amplitude (max vector displacement magnitude after normalization)
    amplitude: float = 0.25

    # Wave speed (in your simulation units), omega = c * k
    wave_speed: float = 1.0

    # --- Polarization plane for the rotating mode (embedding space) ---
    # If not provided, we choose a convenient default based on `polarization`.
    polarization: str = "xy"  # "xy", "xz", "yz", "spatial", "all"
    polarization_p1: Optional[Tuple[float, ...]] = None
    polarization_p2: Optional[Tuple[float, ...]] = None

    # --- Containment deformation (static) in a chosen embedding component ---
    # This is the "X^4 trap shape" part, kept separate from the oscillatory polarization.
    containment_component: int = 3  # default: X^4 in a 4D embedding
    containment_depth: float = 0.15
    containment_sigma: float = 0.5  # interpreted relative to `radius` if <= 0 -> disabled

    # Smooth radial taper at boundary (avoid sharp cutoffs)
    smooth_edge: float = 2.0  # larger => softer edge

    # Backward-compatibility alias (older patch used field_component)
    field_component: Optional[int] = None


def _orthonormalize_plane(p1: np.ndarray, p2: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    p1 = p1.astype(np.float64, copy=False)
    p2 = p2.astype(np.float64, copy=False)

    n1 = np.linalg.norm(p1)
    if n1 == 0:
        raise ValueError("polarization_p1 must not be the zero vector.")
    p1 = p1 / n1

    # Gram–Schmidt
    p2 = p2 - np.dot(p2, p1) * p1
    n2 = np.linalg.norm(p2)
    if n2 == 0:
        raise ValueError("polarization_p2 must not be collinear with polarization_p1.")
    p2 = p2 / n2
    return p1, p2


def _default_polarization_plane(embed_dim: int, mode: str) -> Tuple[np.ndarray, np.ndarray]:
    # Plane vectors in embedding space
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
        # spreads energy across all three spatial axes while still being a plane
        # p1 ~ (1,1,1,0), p2 ~ (1,-1,0,0) (then orthonormalize)
        if embed_dim < 3:
            raise ValueError("polarization='spatial' requires at least 3 embedding dims.")
        p1[0:3] = 1.0
        p2[0] = 1.0
        p2[1] = -1.0
    elif mode == "all":
        # include X^4 as well (if available): p1 ~ (1,1,1,1)
        if embed_dim < 4:
            # fall back to spatial
            return _default_polarization_plane(embed_dim, "spatial")
        p1[:] = 1.0
        p2[0] = 1.0
        p2[1] = -1.0
    else:
        raise ValueError(f"Unknown polarization mode: {mode!r}")
    return _orthonormalize_plane(p1, p2)


def _spherical_coords(x: np.ndarray, y: np.ndarray, z: np.ndarray):
    r = np.sqrt(x * x + y * y + z * z)
    # avoid division by zero
    r_safe = np.where(r > 0, r, 1.0)
    theta = np.arccos(np.clip(z / r_safe, -1.0, 1.0))  # polar angle [0, pi]
    phi = np.arctan2(y, x)  # azimuth [-pi, pi]
    return r, theta, phi


def _smooth_taper(r: np.ndarray, radius: float, smooth_edge: float) -> np.ndarray:
    # 1 inside, 0 outside, smooth around r ~ radius
    if smooth_edge <= 0:
        return (r <= radius).astype(np.float64)
    # tanh ramp: width controlled by smooth_edge
    return 0.5 * (1.0 - np.tanh((r - radius) / (radius / smooth_edge)))


def _spherical_jn_zero(l: int, n: int) -> float:
    """
    Find the nth positive zero of spherical_jn(l, x).
    Robust bracketing based on asymptotic ~ (n + l/2)pi, with fallback scan.
    """
    if n < 1:
        raise ValueError("n must be >= 1 for spherical_jn zeros.")
    # initial guess interval
    pi = np.pi
    a = (n + 0.5 * l - 0.5) * pi
    b = (n + 0.5 * l + 0.5) * pi
    a = max(a, 1e-6)

    f = lambda t: spherical_jn(l, t)

    # Try bracketing; if sign doesn't change, scan outward
    fa, fb = f(a), f(b)
    if fa == 0.0:
        return a
    if fb == 0.0:
        return b

    if fa * fb > 0:
        # scan intervals of length pi until we find a sign change
        left = a
        right = b
        for _ in range(80):
            left = right
            right = right + pi
            fl, fr = f(left), f(right)
            if fl * fr <= 0:
                a, b = left, right
                break
        else:
            raise RuntimeError(f"Could not bracket zero for spherical_jn(l={l}, n={n}).")
    return float(brentq(f, a, b))


def _get_XV(state):
    # Get positions and velocities from BraneState
    # Try both attribute names for compatibility
    X = getattr(state, "X", None)
    V = getattr(state, "V", None)

    # If X/V not found, try positions/velocities (BraneState convention)
    if X is None:
        X = getattr(state, "positions", None)
    if V is None:
        V = getattr(state, "velocities", None)

    if X is None or V is None:
        raise AttributeError("State must expose state.X and state.V (or state.positions and state.velocities) arrays.")
    return X, V


def _get_reference_xyz(grid, state) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Try to get reference coordinates for initialization.
    Fallback to current X[...,0:3] if no explicit grid coords exist.
    """
    # Common patterns in similar codebases
    for attr in ("positions", "coords", "X0", "X_ref", "X"):
        if hasattr(grid, attr):
            arr = getattr(grid, attr)
            arr = np.asarray(arr)
            if arr.shape[-1] >= 3:
                return arr[..., 0], arr[..., 1], arr[..., 2]

    # fallback: use current embedding (works at t=0)
    X, _ = _get_XV(state)
    return X[..., 0], X[..., 1], X[..., 2]


def initialize_electron_mode_3d(
    state,
    grid,
    spec: ElectronModeSpec,
    set_velocities: bool = True,
    normalize: bool = True,
    return_debug: bool = False,
) -> Dict[str, Any] | None:
    """
    Initialize an electron-like rotating standing wave on a 3D brane embedded in (>=3) dims.

    - Oscillatory part: circularly polarized vector field in a chosen embedding plane (p1, p2)
      derived from complex scalar mode psi = R(r) Y_lm(theta, phi).
    - Containment part: static Gaussian well in `containment_component` (default X^4).

    This ensures the electron energy is not confined to X^4 at t=0 and that a rotating phase
    field exists (for m != 0).
    """
    X, V = _get_XV(state)
    embed_dim = X.shape[-1]
    if embed_dim < 3:
        raise ValueError("Embedding dimension must be >= 3.")

    # Back-compat: allow old field_component to act as containment_component
    if spec.field_component is not None:
        spec.containment_component = int(spec.field_component)

    # Choose polarization plane
    if spec.polarization_p1 is not None and spec.polarization_p2 is not None:
        p1 = np.asarray(spec.polarization_p1, dtype=np.float64)
        p2 = np.asarray(spec.polarization_p2, dtype=np.float64)
        if p1.shape[0] != embed_dim or p2.shape[0] != embed_dim:
            raise ValueError("polarization_p1/p2 must have length == embedding dimension.")
        p1, p2 = _orthonormalize_plane(p1, p2)
    else:
        p1, p2 = _default_polarization_plane(embed_dim, spec.polarization)

    # Reference coordinates (in which we define r, theta, phi)
    x, y, z = _get_reference_xyz(grid, state)
    cx, cy, cz = spec.center
    dx = x - cx
    dy = y - cy
    dz = z - cz

    r, theta, phi = _spherical_coords(dx, dy, dz)

    # Radial profile: cavity-like via spherical Bessel with smooth taper
    alpha_ln = _spherical_jn_zero(spec.l, spec.n)
    k = alpha_ln / max(spec.radius, 1e-12)
    omega = spec.wave_speed * k

    taper = _smooth_taper(r, spec.radius, spec.smooth_edge)
    R = spherical_jn(spec.l, k * r) * taper

    # Spherical harmonic: SciPy uses sph_harm(m, l, phi, theta)
    Y = sph_harm(spec.m, spec.l, phi, theta)
    psi = R * Y

    Re = np.real(psi)
    Im = np.imag(psi)

    # Build *raw* vector displacement and velocity fields (no amplitude yet)
    raw_disp = Re[..., None] * p1[None, ...] + Im[..., None] * p2[None, ...]
    raw_vel = omega * (Im[..., None] * p1[None, ...] - Re[..., None] * p2[None, ...])

    scale = 1.0
    if normalize:
        max_norm = float(np.max(np.linalg.norm(raw_disp, axis=-1)))
        if max_norm > 0:
            scale = spec.amplitude / max_norm
        else:
            scale = spec.amplitude

    disp = scale * raw_disp
    X[...] = X[...] + disp

    if set_velocities:
        V[...] = V[...] + scale * raw_vel

    # Static containment deformation in X^4 (or chosen component)
    if 0 <= spec.containment_component < embed_dim and spec.containment_depth != 0.0:
        sigma = spec.containment_sigma
        if sigma <= 0:
            sigma = 0.5 * spec.radius
        if sigma > 0:
            well = -float(spec.containment_depth) * np.exp(-(r * r) / (2.0 * sigma * sigma))
            X[..., spec.containment_component] = X[..., spec.containment_component] + well

    if not return_debug:
        return None

    # Phase of the complex scalar mode (defined where |psi| > 0)
    amp = np.sqrt(Re * Re + Im * Im)
    phase = np.arctan2(Im, Re)

    return {
        "k": k,
        "omega": omega,
        "alpha_ln": alpha_ln,
        "p1": p1.copy(),
        "p2": p2.copy(),
        "amplitude_field": amp,
        "phase_field": phase,
        "taper": taper,
    }
