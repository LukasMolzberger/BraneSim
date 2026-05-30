"""Axis-triplet initializer (Candidate 4, Phase 2) — NEGATIVE CONTROL.

Ansatz (J=1, L=0; three independent axis-polarized scalars):
    xi^i(x) = c_i * g(r)      i in {1,2,3},  c_i = u0 * weight_i
    dX4(x)  = 0

**This is deliberately NOT a baryon ansatz.** The color index ``i`` is *not*
locked to any spatial angular direction: each lateral component is an
independent isotropic bump with its own amplitude ``c_i``. There is no
color-spatial locking and no partial-wave structure, so the mode is expected
to radiate.

Why it exists
-------------
Per the baryon roadmap, this is the falsification baseline for the soliton
search. The confinement claim for the hedgehog / Skyrme-twisted candidates is
only interpretable *relative* to this control:

    If the axis-triplet confines but the locked candidates (1-3) do not, the
    color-spatial-locking framework is wrong.

The independent ``weight_i`` knobs let the control be run isotropically
(1,1,1) or with a deliberate axis asymmetry to probe how the substrate treats
unlocked anisotropic seeds.

Profile shapes (all give g(0)=1, g(inf)=0)
------------------------------------------
gaussian : g(r) = exp(-(r/w)^2)        Gaussian tail.
sech     : g(r) = 1/cosh(r/w)          Exponential tail.
power2   : g(r) = 1/(1+(r/w)^2)        Algebraic 1/r^2 tail.

Usage
-----
    python -m components.initialization.axis_triplet \\
        --alpha 0.2 --u0 0.006 --w 10 --weights 1,1,1 --profile-shape sech \\
        --grid-size 60 --spacing 1.0 --dt 0.05 --num-steps 10000 \\
        --checkpoint-interval 100 --output /path/to/initial_state.npz
"""

from __future__ import annotations

import argparse
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
class AxisTripletSeed:
    alpha: float
    u0: float
    w: float
    weights: tuple[float, float, float]
    profile_shape: str
    J: int = 1
    L: int = 0
    B_winding: int = 0
    ansatz: str = "axis_triplet"
    role: str = "negative_control"


def _profile_gaussian(r: np.ndarray, w: float) -> np.ndarray:
    """g(r) = exp(-(r/w)^2).  g(0)=1, g(inf)=0."""
    return np.exp(-((r / w) ** 2))


def _profile_sech(r: np.ndarray, w: float) -> np.ndarray:
    """g(r) = 1/cosh(r/w).  g(0)=1, g(inf)=0."""
    return 1.0 / np.cosh(r / w)


def _profile_power2(r: np.ndarray, w: float) -> np.ndarray:
    """g(r) = 1/(1+(r/w)^2).  g(0)=1, g(inf)=0."""
    return 1.0 / (1.0 + (r / w) ** 2)


_PROFILES = {
    "gaussian": _profile_gaussian,
    "sech": _profile_sech,
    "power2": _profile_power2,
}


def _parse_weights(s: str) -> tuple[float, float, float]:
    parts = [p.strip() for p in s.split(",") if p.strip() != ""]
    if len(parts) != 3:
        raise ValueError(f"--weights expects three comma-separated values, got {s!r}")
    return (float(parts[0]), float(parts[1]), float(parts[2]))


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Component 1 (Candidate 4): axis-triplet negative control (J=1, L=0)"
    )
    p.add_argument("--output", required=True, help="Path to output .npz file")

    # Physics parameters
    p.add_argument("--alpha", type=float, default=0.20,
                   help="Prestress: rest_length / spacing. Default 0.20.")
    p.add_argument("--u0", type=float, required=True,
                   help="Overall amplitude scale; c_i = u0 * weight_i. Lattice units.")
    p.add_argument("--w", type=float, required=True,
                   help="Profile width in lattice units.")
    p.add_argument("--weights", type=str, default="1,1,1",
                   help="Per-axis amplitude weights w1,w2,w3 (default isotropic 1,1,1).")
    p.add_argument("--profile-shape", choices=list(_PROFILES), default="sech",
                   help="Radial profile g(r). Default 'sech'.")

    # Grid
    p.add_argument("--grid-size", type=int, default=60,
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
    weights: tuple[float, float, float],
    profile_shape: str,
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
) -> tuple[BraneState3D, LatticeConfig, DynamicsConfig, AxisTripletSeed]:
    """Build initial positions/velocities for the axis-triplet control seed."""
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
    seed = AxisTripletSeed(
        alpha=alpha,
        u0=u0,
        w=w,
        weights=weights,
        profile_shape=profile_shape,
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

    # Profile function g(r) (dimensionless shape; amplitudes applied below)
    profile = _PROFILES[profile_shape]
    g = profile(r, w)

    # Independent axis-polarized scalars: xi^i = u0 * weight_i * g(r) ; X4 = 0
    # NO color-spatial locking (this is the negative control).
    u_np = np.zeros((coords.shape[0], 4), dtype=np.float64)
    v_np = np.zeros((coords.shape[0], 4), dtype=np.float64)

    u_np[:, 0] = u0 * weights[0] * g
    u_np[:, 1] = u0 * weights[1] * g
    u_np[:, 2] = u0 * weights[2] * g
    # u_np[:, 3] stays 0

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

    weights = _parse_weights(args.weights)
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
        weights=weights,
        profile_shape=args.profile_shape,
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
            "initializer": "axis_triplet",
            "ansatz": "axis_triplet",
            "role": "negative_control",
            "J": 1,
            "L": 0,
            "B_winding": 0,
            "profile_shape": args.profile_shape,
            "u0": args.u0,
            "w": args.w,
            "weights": list(weights),
            "alpha": args.alpha,
            "ell0": args.alpha * args.spacing,
            "u0_over_ell0": args.u0 / (args.alpha * args.spacing),
            "w_over_a": args.w / args.spacing,
            "color_spatial_locking": False,
            "topological_protection": False,
            "debug": {
                "grid_shape": list(grid_shape),
                "spacing": args.spacing,
                "rest_length": args.alpha * args.spacing,
            },
        },
    )

    print("Axis-triplet (negative control) initialization complete")
    print(f"  output:       {args.output}")
    print(f"  grid:         {grid_shape}")
    print(f"  alpha:        {args.alpha}")
    print(f"  u0:           {args.u0}  (u0/ell0 = {args.u0 / (args.alpha * args.spacing):.4f})")
    print(f"  weights:      {weights}")
    print(f"  w (seeded):   {args.w} a")
    print(f"  profile:      {args.profile_shape}")
    print(f"  (J,L):        (1,0)  NO color-spatial locking -> expected to radiate")
    print(f"  role:         negative control / falsification baseline")
    print(f"  device:       {device}")


if __name__ == "__main__":
    main()