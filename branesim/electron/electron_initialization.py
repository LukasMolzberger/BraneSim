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
    """
    Electron initialization = rotating spherical-harmonic cavity seed.

    IMPORTANT: The previous 'containment scaffold' (static Gaussian well in X^4) has been removed.
    This initializer now seeds ONLY the rotating standing-wave mode. X^4 can participate via
    the polarization plane (e.g. polarization='spatial_x4' or 'all').
    """

    # Angular mode
    l: int = 1
    m: int = 1
    n: int = 1  # radial index (nth zero)

    # Spatial support / scale in reference coordinates
    radius: float = 1.0
    center: Tuple[float, float, float] = (0.0, 0.0, 0.0)

    # Overall amplitude (max vector displacement magnitude after normalization)
    amplitude: float = 0.25

    # Wave speed (simulation units), omega = c * k
    wave_speed: float = 1.0

    # Polarization plane for rotating mode in embedding space
    # Recommended defaults:
    # - "spatial"     : rotates within X^1..X^3 only (X^4 stays exactly zero in pure spring model)
    # - "spatial_x4"  : rotates in a plane that includes X^4 (X^4 is dynamically excited as part of the mode)
    # - "all"         : plane includes X^4 and mixes spatial components
    polarization: str = "spatial_x4"
    polarization_p1: Optional[Tuple[float, ...]] = None
    polarization_p2: Optional[Tuple[float, ...]] = None

    # Smooth radial taper at boundary (avoid sharp cutoffs)
    smooth_edge: float = 2.0

    # --- Deprecated (removed) containment parameters ---
    # Kept only to fail fast if some older experiment still passes them.
    containment_component: Optional[int] = None
    containment_depth: float = 0.0
    containment_sigma: float = 0.0
    field_component: Optional[int] = None  # legacy alias


def _orthonormalize_plane(p1: np.ndarray, p2: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    p1 = p1.astype(np.float64, copy=False)
    p2 = p2.astype(np.float64, copy=False)

    n1 = np.linalg.norm(p1)
    if n1 == 0:
        raise ValueError("polarization_p1 must not be the zero vector.")
    p1 = p1 / n1

    p2 = p2 - np.dot(p2, p1) * p1
    n2 = np.linalg.norm(p2)
    if n2 == 0:
        raise ValueError("polarization_p2 must not be collinear with polarization_p1.")
    p2 = p2 / n2
    return p1, p2


def _default_polarization_plane(embed_dim: int, mode: str) -> Tuple[np.ndarray, np.ndarray]:
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
        # plane entirely in X^1..X^3
        if embed_dim < 3:
            raise ValueError("polarization='spatial' requires at least 3 embedding dims.")
        p1[0:3] = 1.0
        p2[0] = 1.0
        p2[1] = -1.0
    elif mode == "spatial_x4":
        # plane includes X^4 so it is excited by the mode itself (no scaffold)
        if embed_dim < 4:
            # fallback to spatial if no X^4 exists
            return _default_polarization_plane(embed_dim, "spatial")
        p1[0:3] = 1.0
        p1[3] = 1.0
        p2[0] = 1.0
        p2[1] = -1.0
        p2[3] = 1.0
    elif mode == "all":
        # broad mixing, includes X^4 if available
        if embed_dim < 4:
            return _default_polarization_plane(embed_dim, "spatial")
        p1[:] = 1.0
        p2[0] = 1.0
        p2[1] = -1.0
        p2[2] = 0.5
        p2[3] = 1.0
    else:
        raise ValueError(f"Unknown polarization mode: {mode!r}")

    return _orthonormalize_plane(p1, p2)


def _spherical_coords(x: np.ndarray, y: np.ndarray, z: np.ndarray):
    r = np.sqrt(x * x + y * y + z * z)
    r_safe = np.where(r > 0, r, 1.0)
    theta = np.arccos(np.clip(z / r_safe, -1.0, 1.0))
    phi = np.arctan2(y, x)
    return r, theta, phi


def _smooth_taper(r: np.ndarray, radius: float, smooth_edge: float) -> np.ndarray:
    if smooth_edge <= 0:
        return (r <= radius).astype(np.float64)
    return 0.5 * (1.0 - np.tanh((r - radius) / (radius / smooth_edge)))


def _spherical_jn_zero(l: int, n: int) -> float:
    if n < 1:
        raise ValueError("n must be >= 1 for spherical_jn zeros.")
    pi = np.pi
    a = (n + 0.5 * l - 0.5) * pi
    b = (n + 0.5 * l + 0.5) * pi
    a = max(a, 1e-6)

    f = lambda t: spherical_jn(l, t)

    fa, fb = f(a), f(b)
    if fa == 0.0:
        return float(a)
    if fb == 0.0:
        return float(b)

    if fa * fb > 0:
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
    X = getattr(state, "X", None)
    V = getattr(state, "V", None)
    if X is None or V is None:
        raise AttributeError("State must expose state.X and state.V arrays.")
    return X, V


def _get_reference_xyz(grid, state) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    for attr in ("positions", "coords", "X0", "X_ref", "X"):
        if hasattr(grid, attr):
            arr = np.asarray(getattr(grid, attr))
            if arr.shape[-1] >= 3:
                return arr[..., 0], arr[..., 1], arr[..., 2]

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
    Initialize a rotating spherical-harmonic cavity seed on a 3D brane embedded in >=3 dims.

    This seeds ONLY the oscillatory mode. Any previous static "containment scaffold" has been removed.
    """

    # Fail fast if old scaffold params are still being used
    if spec.field_component is not None or spec.containment_component is not None:
        raise ValueError(
            "Containment scaffold parameters were removed. "
            "Do not pass field_component/containment_component anymore."
        )
    if abs(float(spec.containment_depth)) > 0 or abs(float(spec.containment_sigma)) > 0:
        raise ValueError(
            "Containment scaffold parameters were removed. "
            "Do not pass containment_depth/containment_sigma anymore."
        )

    X, V = _get_XV(state)
    embed_dim = X.shape[-1]
    if embed_dim < 3:
        raise ValueError("Embedding dimension must be >= 3.")

    # Polarization plane
    if spec.polarization_p1 is not None and spec.polarization_p2 is not None:
        p1 = np.asarray(spec.polarization_p1, dtype=np.float64)
        p2 = np.asarray(spec.polarization_p2, dtype=np.float64)
        if p1.shape[0] != embed_dim or p2.shape[0] != embed_dim:
            raise ValueError("polarization_p1/p2 must have length == embedding dimension.")
        p1, p2 = _orthonormalize_plane(p1, p2)
    else:
        p1, p2 = _default_polarization_plane(embed_dim, spec.polarization)

    # Reference coordinates
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

    # Vector displacement/velocity from quadrature embedding
    raw_disp = Re[..., None] * p1[None, ...] + Im[..., None] * p2[None, ...]
    raw_vel = omega * (Im[..., None] * p1[None, ...] - Re[..., None] * p2[None, ...])

    scale = 1.0
    if normalize:
        max_norm = float(np.max(np.linalg.norm(raw_disp, axis=-1)))
        scale = (spec.amplitude / max_norm) if max_norm > 0 else spec.amplitude

    X[...] = X[...] + scale * raw_disp
    if set_velocities:
        V[...] = V[...] + scale * raw_vel

    if not return_debug:
        return None

    amp = np.sqrt(Re * Re + Im * Im)
    phase = np.arctan2(Im, Re)

    return {
        "k": float(k),
        "omega": float(omega),
        "alpha_ln": float(alpha_ln),
        "p1": p1.copy(),
        "p2": p2.copy(),
        "scale": float(scale),
        "amplitude_field": amp,
        "phase_field": phase,
        "taper": taper,
    }
