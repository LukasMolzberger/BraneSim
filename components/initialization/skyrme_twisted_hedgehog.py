"""Skyrme-twisted hedgehog initializer (Candidate 2, Phase 2).

Ansatz:
    xi^i(x) = u0 * x_hat^i * sin(F(r))
    dX4(x)  = u0 * cos(F(r))
    F(0) = pi,  F(inf) = 0,  F(w) = pi/2

Profile shapes
--------------
power2 : F(r) = pi / (1 + (r/w)^2)
    Algebraic 1/r^2 tail. F(0) = pi, F(w) = pi/2, F(inf) = 0 exactly.

tanh   : F(r) = pi * (1 - tanh(s*(r/w - 1))) / 2,  steepness s >= 1
    Exponentially decaying tail. For any s, F(w) = pi*(1-tanh(0))/2 = pi/2
    and F(inf) = 0 exactly. F(0) = pi*(1 - tanh(-s))/2 approaches pi as s
    increases; s=3 gives F(0) ~ 0.998*pi (sufficient for B=1 winding).

Topological winding
-------------------
The map (xi^1, xi^2, xi^3, dX4) = u0*(sin(F)*x_hat, cos(F)) traces S^3 as r
goes from 0 to infinity (Hopf-Skyrme style):
    r = 0   : sin(F) = 0, cos(F) = cos(pi) = -1   -> -X4 pole of S^3
    r = inf : sin(F) = 0, cos(F) = cos(0)  = +1   -> +X4 pole of S^3
This realizes B = 1 winding.

Derrick prediction
------------------
R_h / a ~ (ell0 / u0)^(2/3)  with  ell0 = alpha * a = 0.2 * a (at alpha=0.2).
For u0/ell0 = 0.03, R_h/a ~ 10.  See paper/derivations/vsh_channel_decomposition.md
section 2.5.

Usage
-----
    python -m components.initialization.skyrme_twisted_hedgehog \\
        --alpha 0.2 --u0 0.006 --w 10 --profile-shape tanh \\
        --grid-size 80 --spacing 1.0 --dt 0.05 --num-steps 10000 \\
        --checkpoint-interval 100 --output /path/to/initial_state.npz
"""

from __future__ import annotations

import argparse
import math
from dataclasses import asdict, dataclass

import numpy as np
import torch

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
class SkyrmeTwistedHedgehogSeed:
    alpha: float
    u0: float
    w: float
    profile_shape: str
    tanh_steepness: float
    J: int = 0
    L: int = 1
    B_winding: int = 1
    ansatz: str = "skyrme_twisted_hedgehog"


def _profile_power2(r: np.ndarray, w: float) -> np.ndarray:
    """F(r) = pi / (1 + (r/w)^2).  F(0)=pi, F(w)=pi/2, F(inf)=0."""
    return math.pi / (1.0 + (r / w) ** 2)


def _profile_tanh(r: np.ndarray, w: float, steepness: float = 3.0) -> np.ndarray:
    """F(r) = pi*(1 - tanh(steepness*(r/w - 1))) / 2.

    F(w) = pi/2 exactly for any steepness.
    F(0) ~ pi for large steepness (steepness=3 gives F(0)~0.998*pi).
    F(inf) = 0.
    """
    return math.pi * (1.0 - np.tanh(steepness * (r / w - 1.0))) / 2.0


def _derrick_prediction(u0: float, alpha: float, a: float = 1.0) -> float:
    """R_h/a ~ (ell0/u0)^(2/3), ell0 = alpha*a."""
    ell0 = alpha * a
    return (ell0 / max(u0, 1e-30)) ** (2.0 / 3.0)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Component 1 (Candidate 2): Skyrme-twisted hedgehog initializer"
    )
    p.add_argument("--output", required=True, help="Path to output .npz file")

    # Physics parameters
    p.add_argument("--alpha", type=float, default=0.20,
                   help="Prestress: rest_length / spacing. Default 0.20.")
    p.add_argument("--u0", type=float, required=True,
                   help="X4-amplitude (and lateral amplitude) scale. In lattice units.")
    p.add_argument("--w", type=float, required=True,
                   help="Profile half-width in lattice units; F(w) = pi/2.")
    p.add_argument("--profile-shape", choices=["power2", "tanh"], default="tanh",
                   help="Radial profile F(r). 'power2': pi/(1+(r/w)^2). 'tanh': exponential.")
    p.add_argument("--tanh-steepness", type=float, default=3.0,
                   help="Steepness factor s in tanh profile: tanh(s*(r/w - 1)). Default 3.")

    # Grid
    p.add_argument("--grid-size", type=int, default=80,
                   help="Cubic grid side length N (NxNxN). Overrides nx/ny/nz.")
    p.add_argument("--nx", type=int, default=None)
    p.add_argument("--ny", type=int, default=None)
    p.add_argument("--nz", type=int, default=None)
    p.add_argument("--spacing", type=float, default=1.0, help="Lattice spacing a.")

    # Dynamics (carried through to DynamicsConfig)
    p.add_argument("--spring-constant", type=float, default=1.0)
    p.add_argument("--mass-density", type=float, default=1.0)
    p.add_argument("--dt", type=float, default=0.05)
    p.add_argument("--num-steps", type=int, default=10000)
    p.add_argument("--checkpoint-interval", type=int, default=100)

    # Boundary conditions
    p.add_argument("--periodic-axes", type=str, default="true,true,true",
                   help="Periodic boundary conditions per axis. Default all-periodic.")
    p.add_argument("--fixed-boundaries", action="store_true",
                   help="If set, clamp boundary nodes to rest (only sensible with non-periodic axes). "
                        "Default is free boundaries, which is required for periodic soliton hunts.")
    p.add_argument("--axial-weight", type=float, default=1.0)

    # Device
    p.add_argument("--device", type=str, default="auto")
    p.add_argument("--dtype", type=str, default="float64", choices=("float32", "float64"))

    return p.parse_args()


def build_initial_state(
    alpha: float,
    u0: float,
    w: float,
    profile_shape: str,
    tanh_steepness: float,
    grid_shape: tuple[int, int, int],
    spacing: float,
    spring_constant: float,
    mass_density: float,
    dt: float,
    num_steps: int,
    checkpoint_interval: int,
    periodic_axes: tuple[bool, bool, bool],
    fixed_boundaries: bool,
    axial_weight: float,
    device: torch.device,
    dtype: torch.dtype,
) -> tuple[BraneState3D, LatticeConfig, DynamicsConfig, SkyrmeTwistedHedgehogSeed]:
    """Build initial positions/velocities for the Skyrme-twisted hedgehog."""
    rest_length = alpha * spacing

    lattice = LatticeConfig(
        grid_shape=grid_shape,
        spacing=spacing,
        periodic_axes=periodic_axes,
        fixed_boundaries=fixed_boundaries,
        axial_weight=axial_weight,
    )
    dynamics = DynamicsConfig(
        spring_constant=spring_constant,
        rest_length=rest_length,
        mass_density=mass_density,
        dt=dt,
        num_steps=num_steps,
        checkpoint_interval=checkpoint_interval,
    )
    seed = SkyrmeTwistedHedgehogSeed(
        alpha=alpha,
        u0=u0,
        w=w,
        profile_shape=profile_shape,
        tanh_steepness=tanh_steepness,
    )

    state = BraneState3D(grid_shape, device=device, dtype=dtype)
    state.initialize_flat_configuration(spacing)
    if fixed_boundaries:
        state.set_fixed_boundaries()

    # Node coordinates (Cartesian, in lattice units)
    coords = state.positions[:, :3].detach().cpu().to(torch.float64).numpy()
    nx, ny, nz = grid_shape
    center = np.array([0.5 * (nx - 1) * spacing,
                       0.5 * (ny - 1) * spacing,
                       0.5 * (nz - 1) * spacing], dtype=np.float64)
    dx = coords - center[np.newaxis, :]  # (N_nodes, 3)
    r = np.sqrt(np.sum(dx ** 2, axis=1))  # (N_nodes,)
    r_safe = np.where(r > 1e-12, r, 1e-12)

    # Unit radial vector x_hat^i = dx^i / r
    x_hat = dx / r_safe[:, np.newaxis]  # (N_nodes, 3)

    # Profile function F(r)
    if profile_shape == "power2":
        F = _profile_power2(r, w)
    elif profile_shape == "tanh":
        F = _profile_tanh(r, w, steepness=tanh_steepness)
    else:
        raise ValueError(f"Unknown profile_shape: {profile_shape!r}")

    sin_F = np.sin(F)
    cos_F = np.cos(F)

    # Lateral displacement: xi^i = u0 * x_hat^i * sin(F(r))
    # X4 displacement:      dX4  = u0 * cos(F(r))
    # Velocities: v=0 (static seed; soliton will relax dynamically)
    u_np = np.zeros((coords.shape[0], 4), dtype=np.float64)
    v_np = np.zeros((coords.shape[0], 4), dtype=np.float64)

    u_np[:, 0] = u0 * x_hat[:, 0] * sin_F
    u_np[:, 1] = u0 * x_hat[:, 1] * sin_F
    u_np[:, 2] = u0 * x_hat[:, 2] * sin_F
    u_np[:, 3] = u0 * cos_F

    np_dtype = np.float32 if dtype == torch.float32 else np.float64
    state.set_kinematics(
        torch.from_numpy(u_np.astype(np_dtype, copy=False)).to(device=device, dtype=dtype),
        torch.from_numpy(v_np.astype(np_dtype, copy=False)).to(device=device, dtype=dtype),
    )
    state.apply_fixed_boundaries()

    return state, lattice, dynamics, seed


def main() -> None:
    args = parse_args()

    # Resolve grid shape
    if args.nx is not None or args.ny is not None or args.nz is not None:
        nx = args.nx or args.grid_size
        ny = args.ny or args.grid_size
        nz = args.nz or args.grid_size
    else:
        nx = ny = nz = args.grid_size
    grid_shape = (nx, ny, nz)

    periodic = parse_bool_triple(args.periodic_axes)
    fixed_boundaries = bool(args.fixed_boundaries)
    if fixed_boundaries and any(periodic):
        raise ValueError(
            "--fixed-boundaries is incompatible with any periodic axis. "
            "Boundary clamping on a periodic domain imposes a hard positional "
            "constraint that violates backbone #5 (no clamps / cutoffs)."
        )

    device = choose_device(args.device)
    dtype = choose_dtype(args.dtype, device)

    state, lattice, dynamics, seed = build_initial_state(
        alpha=args.alpha,
        u0=args.u0,
        w=args.w,
        profile_shape=args.profile_shape,
        tanh_steepness=args.tanh_steepness,
        grid_shape=grid_shape,
        spacing=args.spacing,
        spring_constant=args.spring_constant,
        mass_density=args.mass_density,
        dt=args.dt,
        num_steps=args.num_steps,
        checkpoint_interval=args.checkpoint_interval,
        periodic_axes=periodic,
        fixed_boundaries=fixed_boundaries,
        axial_weight=args.axial_weight,
        device=device,
        dtype=dtype,
    )

    derrick_w = _derrick_prediction(args.u0, args.alpha, a=args.spacing)

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
            "initializer": "skyrme_twisted_hedgehog",
            "ansatz": "skyrme_twisted_hedgehog",
            "J": 0,
            "L": 1,
            "B_winding": 1,
            "profile_shape": args.profile_shape,
            "u0": args.u0,
            "w": args.w,
            "alpha": args.alpha,
            "derrick_predicted_w": derrick_w,
            "ell0": args.alpha * args.spacing,
            "u0_over_ell0": args.u0 / (args.alpha * args.spacing),
            "w_over_a": args.w / args.spacing,
            "debug": {
                "grid_shape": list(grid_shape),
                "spacing": args.spacing,
                "rest_length": args.alpha * args.spacing,
                "tanh_steepness": args.tanh_steepness,
            },
        },
    )

    print("Skyrme-twisted hedgehog initialization complete")
    print(f"  output:       {args.output}")
    print(f"  grid:         {grid_shape}")
    print(f"  alpha:        {args.alpha}")
    print(f"  u0:           {args.u0}  (u0/ell0 = {args.u0 / (args.alpha * args.spacing):.4f})")
    print(f"  w (seeded):   {args.w} a")
    print(f"  Derrick pred: w ~ {derrick_w:.2f} a  (= (ell0/u0)^(2/3))")
    print(f"  profile:      {args.profile_shape}", end="")
    if args.profile_shape == "tanh":
        print(f" (steepness={args.tanh_steepness})", end="")
    print()
    print(f"  F(0):         pi  (winding B=1 start)")
    print(f"  F(w):         pi/2  (profile half-width)")
    print(f"  F(inf):       0  (vacuum)")
    print(f"  device:       {device}")


if __name__ == "__main__":
    main()