"""Hedgehog initializer (Candidate 1, Phase 2).

Ansatz (J=0, L=1 vector spherical harmonic):
    xi^i(x) = f(r) * x_hat^i        f(0) = u0,  f(inf) = 0
    dX4(x)  = 0                     (no gravity-channel / no winding)

The color index ``i`` is locked to the spatial angular direction ``x_hat^i``;
the entire structure lives in the SU(3) traceless sector, so the U(1) trace
``sum_i xi^i`` averages to zero over angle and the far field is trace-neutral.
This is the canonical baryon seed (backbone #20, paper sec:soliton-labels) and
the working *neutron-class* hypothesis (trace-cancelled far field).

Topology / stability
---------------------
Winding number 0 (no S^3 wrap, unlike the Skyrme-twisted hedgehog). There is
therefore **no topological protection**: the mode is expected to radiate unless
the geometric quartic alone suffices to confine it. No Derrick-stable size is
predicted; the seeded width ``w`` is a search parameter, not a prediction.

Profile shapes (all give f(0)=u0, f(inf)=0)
-------------------------------------------
gaussian : f(r) = u0 * exp(-(r/w)^2)        Gaussian tail.
sech     : f(r) = u0 / cosh(r/w)            Exponential tail.
power2   : f(r) = u0 / (1 + (r/w)^2)        Algebraic 1/r^2 tail.

Note on the r=0 lock
--------------------
``xi^i = f(r) x_hat^i`` has an undefined direction at the exact center. For an
even grid no node lands on r=0 (the center falls between nodes), so the
singularity is never sampled; a lone center node (odd grid) is set to zero
displacement via the r_safe guard. The peak amplitude u0 is the radial-field
magnitude near the center.

Usage
-----
    python -m components.initialization.hedgehog \\
        --alpha 0.2 --u0 0.006 --w 10 --profile-shape sech \\
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
class HedgehogSeed:
    alpha: float
    u0: float
    w: float
    profile_shape: str
    J: int = 0
    L: int = 1
    B_winding: int = 0
    ansatz: str = "hedgehog"


def _profile_gaussian(r: np.ndarray, w: float) -> np.ndarray:
    """f(r) = exp(-(r/w)^2).  f(0)=1, f(inf)=0."""
    return np.exp(-((r / w) ** 2))


def _profile_sech(r: np.ndarray, w: float) -> np.ndarray:
    """f(r) = 1/cosh(r/w).  f(0)=1, f(inf)=0."""
    return 1.0 / np.cosh(r / w)


def _profile_power2(r: np.ndarray, w: float) -> np.ndarray:
    """f(r) = 1/(1+(r/w)^2).  f(0)=1, f(inf)=0."""
    return 1.0 / (1.0 + (r / w) ** 2)


_PROFILES = {
    "gaussian": _profile_gaussian,
    "sech": _profile_sech,
    "power2": _profile_power2,
}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Component 1 (Candidate 1): hedgehog initializer (J=0, L=1)"
    )
    p.add_argument("--output", required=True, help="Path to output .npz file")

    # Physics parameters
    p.add_argument("--alpha", type=float, default=0.20,
                   help="Prestress: rest_length / spacing. Default 0.20.")
    p.add_argument("--u0", type=float, required=True,
                   help="Peak lateral amplitude f(0), in lattice units.")
    p.add_argument("--w", type=float, required=True,
                   help="Profile width in lattice units.")
    p.add_argument("--profile-shape", choices=list(_PROFILES), default="sech",
                   help="Radial profile f(r). Default 'sech'.")

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
) -> tuple[BraneState3D, LatticeConfig, DynamicsConfig, HedgehogSeed]:
    """Build initial positions/velocities for the hedgehog seed."""
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
    seed = HedgehogSeed(
        alpha=alpha,
        u0=u0,
        w=w,
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

    # Unit radial vector x_hat^i = dx^i / r. An exact-center node (odd grid)
    # has dx=0, so dividing by max(r, eps) already yields 0 there; the explicit
    # zero-fill documents that the radial direction is undefined at r=0.
    x_hat = dx / np.maximum(r, 1e-12)[:, np.newaxis]  # (N_nodes, 3)
    x_hat[r <= 1e-12] = 0.0

    # Profile function f(r) (dimensionless shape; amplitude applied below)
    profile = _PROFILES[profile_shape]
    f = profile(r, w)

    # Lateral displacement: xi^i = u0 * f(r) * x_hat^i ; X4 = 0
    u_np = np.zeros((coords.shape[0], 4), dtype=np.float64)
    v_np = np.zeros((coords.shape[0], 4), dtype=np.float64)

    u_np[:, 0] = u0 * f * x_hat[:, 0]
    u_np[:, 1] = u0 * f * x_hat[:, 1]
    u_np[:, 2] = u0 * f * x_hat[:, 2]
    # u_np[:, 3] stays 0 (no gravity-channel / no winding)

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
            "initializer": "hedgehog",
            "ansatz": "hedgehog",
            "J": 0,
            "L": 1,
            "B_winding": 0,
            "profile_shape": args.profile_shape,
            "u0": args.u0,
            "w": args.w,
            "alpha": args.alpha,
            "ell0": args.alpha * args.spacing,
            "u0_over_ell0": args.u0 / (args.alpha * args.spacing),
            "w_over_a": args.w / args.spacing,
            "trace_sector": "traceless (SU(3)); trace-neutral far field",
            "topological_protection": False,
            "debug": {
                "grid_shape": list(grid_shape),
                "spacing": args.spacing,
                "rest_length": args.alpha * args.spacing,
            },
        },
    )

    print("Hedgehog initialization complete")
    print(f"  output:       {args.output}")
    print(f"  grid:         {grid_shape}")
    print(f"  alpha:        {args.alpha}")
    print(f"  u0:           {args.u0}  (u0/ell0 = {args.u0 / (args.alpha * args.spacing):.4f})")
    print(f"  w (seeded):   {args.w} a")
    print(f"  profile:      {args.profile_shape}")
    print(f"  (J,L):        (0,1)  winding B=0  -> no topological protection")
    print(f"  far field:    trace-neutral (neutron-class hypothesis)")
    print(f"  device:       {device}")


if __name__ == "__main__":
    main()