"""Component 1: baryon initialization.

Creates one initial-state file consumed by simulation.
"""

from __future__ import annotations

import argparse
import warnings
from dataclasses import asdict, dataclass

import numpy as np
import torch

try:
    from scipy.optimize import brentq
    from scipy.special import sph_harm, spherical_jn
except Exception as exc:  # pragma: no cover
    raise ImportError("Initialization requires SciPy (pip install scipy)") from exc

from components.shared import (
    DynamicsConfig,
    LatticeConfig,
    BraneState3D,
    choose_device,
    choose_dtype,
    parse_bool_triple,
    save_initial_state,
)


@dataclass(frozen=True)
class SeedConfig:
    l: int
    m: int
    n: int
    radius: float
    amplitude: float
    wave_speed: float
    smooth_edge: float
    axis_amplitudes: tuple[float, float, float]
    axis_phase_offsets: tuple[float, float, float]
    mixing_strength: float
    x4_trace_weight: float


def _triple_floats(value: str) -> tuple[float, float, float]:
    parts = [float(v.strip()) for v in value.split(",") if v.strip()]
    if len(parts) != 3:
        raise argparse.ArgumentTypeError("Expected 3 comma-separated floats")
    return parts[0], parts[1], parts[2]


def _spherical_bessel_zero(l: int, n: int) -> float:
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
        raise RuntimeError("Could not find required spherical Bessel root")
    return roots[n - 1]


def _compute_spherical(coords_xyz: np.ndarray, center: tuple[float, float, float]):
    x = coords_xyz[:, 0] - center[0]
    y = coords_xyz[:, 1] - center[1]
    z = coords_xyz[:, 2] - center[2]
    r = np.sqrt(x * x + y * y + z * z)
    r_safe = np.where(r > 0.0, r, 1.0)
    theta = np.arccos(np.clip(z / r_safe, -1.0, 1.0))
    phi = np.arctan2(y, x)
    return r, theta, phi


def _radial_window(r: np.ndarray, radius: float, smooth_edge: float) -> np.ndarray:
    if smooth_edge <= 0.0:
        return (r <= radius).astype(np.float64)
    radius = max(radius, 1e-12)
    return 0.5 * (1.0 - np.tanh((r - radius) / (radius / smooth_edge)))


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Component 1: baryon initialization")
    p.add_argument("--output", required=True)

    p.add_argument("--nx", type=int, default=64)
    p.add_argument("--ny", type=int, default=64)
    p.add_argument("--nz", type=int, default=64)
    p.add_argument("--spacing", type=float, default=1.0)
    p.add_argument("--periodic-axes", type=str, default="false,false,false")
    p.add_argument("--shell-weights", type=_triple_floats, default=(1.0, 1.0, 1.0),
                   help="Effective shell weights for axial, face-diagonal, and body-diagonal links.")
    p.add_argument("--free-boundaries", action="store_true")

    p.add_argument("--spring-constant", type=float, required=True)
    p.add_argument("--rest-length", type=float, required=True)
    p.add_argument("--mass-density", type=float, required=True)
    p.add_argument("--dt", type=float, required=True)
    p.add_argument("--num-steps", type=int, required=True)
    p.add_argument("--checkpoint-interval", type=int, default=1)

    p.add_argument("--l", type=int, default=1)
    p.add_argument("--m", type=int, default=1)
    p.add_argument("--n", type=int, default=1)
    p.add_argument("--radius", type=float, default=10.0)
    p.add_argument("--amplitude", type=float, default=0.25)
    p.add_argument("--wave-speed", type=float, default=1.0)
    p.add_argument("--smooth-edge", type=float, default=2.0)
    p.add_argument("--axis-amplitudes", type=_triple_floats, default=(1.0, 1.0, 1.0))
    p.add_argument("--axis-phase-offsets", type=_triple_floats, default=(0.0, 2.0943951023931953, 4.1887902047863905))
    p.add_argument("--mixing-strength", type=float, default=0.15)
    p.add_argument("--x4-trace-weight", type=float, default=0.35)

    p.add_argument("--device", type=str, default="auto")
    p.add_argument("--dtype", type=str, default="float64", choices=("float32", "float64"))
    return p.parse_args()


def main() -> None:
    args = parse_args()

    lattice = LatticeConfig(
        grid_shape=(args.nx, args.ny, args.nz),
        spacing=float(args.spacing),
        periodic_axes=parse_bool_triple(args.periodic_axes),
        fixed_boundaries=not args.free_boundaries,
        shell_weights=args.shell_weights,
    )
    dynamics = DynamicsConfig(
        spring_constant=float(args.spring_constant),
        rest_length=float(args.rest_length),
        mass_density=float(args.mass_density),
        dt=float(args.dt),
        num_steps=int(args.num_steps),
        checkpoint_interval=int(args.checkpoint_interval),
    )
    seed = SeedConfig(
        l=int(args.l),
        m=int(args.m),
        n=int(args.n),
        radius=float(args.radius),
        amplitude=float(args.amplitude),
        wave_speed=float(args.wave_speed),
        smooth_edge=float(args.smooth_edge),
        axis_amplitudes=args.axis_amplitudes,
        axis_phase_offsets=args.axis_phase_offsets,
        mixing_strength=float(args.mixing_strength),
        x4_trace_weight=float(args.x4_trace_weight),
    )

    device = choose_device(args.device)
    dtype = choose_dtype(args.dtype, device)

    state = BraneState3D(lattice.grid_shape, device=device, dtype=dtype)
    state.initialize_flat_configuration(lattice.spacing)
    if lattice.fixed_boundaries:
        state.set_fixed_boundaries()

    coords = state.positions[:, :3].detach().cpu().to(torch.float64).numpy()
    center = tuple(0.5 * (n - 1) * lattice.spacing for n in lattice.grid_shape)
    r, theta, phi = _compute_spherical(coords, center)

    alpha_ln = _spherical_bessel_zero(seed.l, seed.n)
    k = float(alpha_ln) / max(seed.radius, 1e-12)
    omega = seed.wave_speed * k

    radial = spherical_jn(seed.l, k * r) * _radial_window(r, seed.radius, seed.smooth_edge)
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message="`scipy.special.sph_harm` is deprecated",
            category=DeprecationWarning,
        )
        base = radial * sph_harm(seed.m, seed.l, phi, theta)

    amps = np.asarray(seed.axis_amplitudes, dtype=np.float64)
    phases = np.asarray(seed.axis_phase_offsets, dtype=np.float64)
    channels = np.stack([amps[i] * base * np.exp(1j * phases[i]) for i in range(3)], axis=1)

    mix = float(np.clip(seed.mixing_strength, 0.0, 1.0))
    mix_matrix = (1.0 - mix) * np.eye(3, dtype=np.float64) + (mix / 3.0) * np.ones((3, 3), dtype=np.float64)
    channels = channels @ mix_matrix.T

    mag = np.sqrt(np.sum(np.abs(channels) ** 2, axis=1))
    scale = seed.amplitude / max(float(np.max(mag)), 1e-12)
    channels *= scale

    u_np = np.zeros((coords.shape[0], 4), dtype=np.float64)
    v_np = np.zeros((coords.shape[0], 4), dtype=np.float64)
    u_np[:, 0:3] = np.real(channels)
    v_np[:, 0:3] = omega * np.imag(channels)

    trace_mode = np.sum(channels, axis=1) / np.sqrt(3.0)
    u_np[:, 3] = seed.x4_trace_weight * np.real(trace_mode)
    v_np[:, 3] = seed.x4_trace_weight * omega * np.imag(trace_mode)

    np_dtype = np.float32 if dtype == torch.float32 else np.float64
    state.set_kinematics(
        torch.from_numpy(u_np.astype(np_dtype, copy=False)).to(device=device, dtype=dtype),
        torch.from_numpy(v_np.astype(np_dtype, copy=False)).to(device=device, dtype=dtype),
    )
    state.apply_fixed_boundaries()

    save_initial_state(
        args.output,
        positions=state.positions.detach().cpu().numpy(),
        velocities=state.velocities.detach().cpu().numpy(),
        rest_positions=state.rest_positions.detach().cpu().numpy(),
        lattice=lattice,
        dynamics=dynamics,
        seed=asdict(seed),
        metadata={
            "component": "initialization",
            "debug": {
                "alpha_ln": float(alpha_ln),
                "k": float(k),
                "omega": float(omega),
                "scale": float(scale),
                "center": list(center),
            },
        },
    )

    print("Initialization complete")
    print(f"  output: {args.output}")
    print(f"  grid: {lattice.grid_shape}")
    print(f"  shell_weights: {tuple(round(v, 6) for v in lattice.shell_weights)}")
    print(f"  device: {device}")


if __name__ == "__main__":
    main()
