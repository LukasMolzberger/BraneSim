"""Baryon initialization component (spherical coordinates + spherical harmonics)."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
import torch

try:
    from scipy.optimize import brentq
    from scipy.special import sph_harm, spherical_jn
except Exception as exc:  # pragma: no cover - explicit runtime error path
    raise ImportError(
        "baryon initialization requires SciPy (scipy.optimize, scipy.special). "
        "Install with: pip install scipy"
    ) from exc

from branesim.core.grid import BraneGrid
from branesim.core.state import BraneState, Dimensionality

from .io import InitialStatePackage
from .models import BaryonSeedConfig, DynamicsConfig, LatticeConfig


@dataclass(frozen=True)
class InitializationResult:
    """Return type for initialization component."""

    package: InitialStatePackage
    debug: dict[str, Any]


def _choose_device(device: str) -> torch.device:
    if device == "auto":
        if torch.cuda.is_available():
            return torch.device("cuda")
        if hasattr(torch.backends, "mps") and torch.backends.mps.is_available():
            return torch.device("mps")
        return torch.device("cpu")
    return torch.device(device)


def _spherical_bessel_zero(l: int, n: int) -> float:
    if l < 0:
        raise ValueError("l must be >= 0")
    if n <= 0:
        raise ValueError("n must be >= 1")

    roots: list[float] = []
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
        raise RuntimeError(f"Failed to find {n} roots for j_{l}. Found {len(roots)}")
    return roots[n - 1]


def _compute_spherical_coordinates(coords_xyz: np.ndarray, center: tuple[float, float, float]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    x = coords_xyz[:, 0] - float(center[0])
    y = coords_xyz[:, 1] - float(center[1])
    z = coords_xyz[:, 2] - float(center[2])

    r = np.sqrt(x * x + y * y + z * z)
    r_safe = np.where(r > 0.0, r, 1.0)
    theta = np.arccos(np.clip(z / r_safe, -1.0, 1.0))
    phi = np.arctan2(y, x)
    return r, theta, phi


def _radial_window(r: np.ndarray, radius: float, smooth_edge: float) -> np.ndarray:
    radius = max(float(radius), 1e-12)
    if smooth_edge <= 0.0:
        return (r <= radius).astype(np.float64)
    return 0.5 * (1.0 - np.tanh((r - radius) / (radius / float(smooth_edge))))


def _default_center(lattice: LatticeConfig) -> tuple[float, float, float]:
    return tuple(0.5 * (n - 1) * lattice.spacing for n in lattice.grid_shape)


def initialize_baryon_triplet_state(
    seed: BaryonSeedConfig,
    lattice: LatticeConfig,
    dynamics: DynamicsConfig,
    *,
    dtype: torch.dtype = torch.float64,
    device: str = "auto",
) -> InitializationResult:
    """Create a baryon-like triplet seed as a compressed-state package."""

    if len(lattice.grid_shape) != 3:
        raise ValueError("Baryon initializer expects a 3D lattice grid_shape=(nx,ny,nz)")

    torch_device = _choose_device(device)
    if torch_device.type == "mps" and dtype == torch.float64:
        dtype = torch.float32
    state = BraneState(lattice.grid_shape, Dimensionality.THREE_D, torch_device, dtype)
    state.initialize_flat_configuration(lattice.spacing)
    grid = BraneGrid(
        lattice.grid_shape,
        Dimensionality.THREE_D,
        lattice.spacing,
        torch_device,
        periodic_axes=lattice.periodic_axes,
    )

    if lattice.fixed_boundaries:
        state.set_fixed_boundaries()

    coords = grid.get_spatial_coordinates().detach().cpu().to(torch.float64).numpy()
    center = seed.center if seed.center is not None else _default_center(lattice)
    r, theta, phi = _compute_spherical_coordinates(coords, center)

    alpha_ln = _spherical_bessel_zero(seed.l, seed.n)
    k = float(alpha_ln) / max(float(seed.radius), 1e-12)
    omega = float(seed.wave_speed) * k

    window = _radial_window(r, seed.radius, seed.smooth_edge)
    radial = spherical_jn(seed.l, k * r) * window

    # SciPy convention: sph_harm(m, l, azimuth, polar)
    harmonic = sph_harm(seed.m, seed.l, phi, theta)
    base_mode = radial * harmonic

    axis_amplitudes = np.asarray(seed.axis_amplitudes, dtype=np.float64)
    axis_phases = np.asarray(seed.axis_phase_offsets, dtype=np.float64)
    if axis_amplitudes.shape != (3,) or axis_phases.shape != (3,):
        raise ValueError("axis_amplitudes and axis_phase_offsets must have exactly 3 entries")

    channels = np.stack(
        [
            axis_amplitudes[i] * base_mode * np.exp(1j * axis_phases[i])
            for i in range(3)
        ],
        axis=1,
    )

    mixing = float(np.clip(seed.mixing_strength, 0.0, 1.0))
    mixing_matrix = (1.0 - mixing) * np.eye(3, dtype=np.float64) + (mixing / 3.0) * np.ones((3, 3), dtype=np.float64)
    channels = channels @ mixing_matrix.T

    mode_mag = np.sqrt(np.sum(np.abs(channels) ** 2, axis=1))
    max_mag = float(np.max(mode_mag))
    if max_mag < 1e-12:
        max_mag = 1.0
    channels *= float(seed.amplitude) / max_mag

    u_np = np.zeros((coords.shape[0], 4), dtype=np.float64)
    v_np = np.zeros((coords.shape[0], 4), dtype=np.float64)

    # Triplet channels occupy the three in-brane embedding components.
    u_np[:, 0:3] = np.real(channels)
    v_np[:, 0:3] = omega * np.imag(channels)

    # U(1)-trace-like sector in X^4; this supports proton vs neutron hypotheses.
    trace_mode = np.sum(channels, axis=1) / np.sqrt(3.0)
    u_np[:, 3] = float(seed.x4_trace_weight) * np.real(trace_mode)
    v_np[:, 3] = float(seed.x4_trace_weight) * omega * np.imag(trace_mode)

    np_dtype = np.float32 if dtype == torch.float32 else np.float64
    u = torch.from_numpy(u_np.astype(np_dtype, copy=False)).to(device=torch_device, dtype=dtype)
    v = torch.from_numpy(v_np.astype(np_dtype, copy=False)).to(device=torch_device, dtype=dtype)

    state.set_kinematics(u=u, v=v)
    state.apply_fixed_boundaries()

    package = InitialStatePackage(
        positions=state.positions.detach().cpu().numpy(),
        velocities=state.velocities.detach().cpu().numpy(),
        rest_positions=state.rest_positions.detach().cpu().numpy(),
        lattice=lattice,
        dynamics=dynamics,
        metadata={
            "component": "baryon_initialization",
            "seed": seed.to_dict(),
            "debug": {
                "alpha_ln": alpha_ln,
                "k": k,
                "omega": omega,
                "mixing_strength": mixing,
                "center": list(center),
            },
        },
    )

    debug = {
        "alpha_ln": alpha_ln,
        "k": k,
        "omega": omega,
        "center": center,
        "max_channel_amplitude": float(np.max(np.abs(channels))),
        "trace_weight": seed.x4_trace_weight,
    }

    return InitializationResult(package=package, debug=debug)
