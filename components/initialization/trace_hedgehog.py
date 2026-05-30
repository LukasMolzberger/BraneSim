"""Trace-admixture hedgehog initializer (Candidate 3, Phase 2).

Ansatz (J=0, L=1 + L=0; hedgehog plus a U(1) trace scalar):
    xi^i(x) = u0 * [ f1(r) * x_hat^i  +  trace_frac * (1/sqrt(3)) * f0(r) ]
    dX4(x)  = 0

This is Candidate 1 (the pure hedgehog) with a scalar **trace admixture**
layered on top:
- The hedgehog part ``f1(r) x_hat^i`` is the L=1 traceless (SU(3)) sector; its
  U(1) trace ``sum_i x_hat^i`` averages to zero over angle -> trace-neutral.
- The trace part ``(1/sqrt(3))(1,1,1)^i f0(r)`` is an L=0 scalar pointing along
  the colour-diagonal (1,1,1) direction. Its trace is ``sqrt(3) f0(r)``, which
  does **not** angularly cancel -> a **nonzero U(1) far field**.

The nonzero trace far field is what makes this the working **proton-class**
hypothesis (vs the trace-neutral neutron-class hedgehog / Skyrme-twisted seeds).
Setting ``trace_frac = 0`` recovers the pure hedgehog (Candidate 1) exactly.

Topology / stability
---------------------
Winding number 0 (no S^3 wrap; X4 unused). No topological protection on its own
-- expected to radiate unless the geometric quartic confines it. The roadmap
notes this can be *combined* with the Skyrme twist of Candidate 2 to give a
topologically stable proton-class seed; that combined seed is a separate
follow-on (it would add the C2 winding + X4 = cos F to the hedgehog part here).

Profile shapes (all give f(0)=1, f(inf)=0)
------------------------------------------
gaussian : exp(-(r/w)^2)        Gaussian tail.
sech     : 1/cosh(r/w)          Exponential tail.
power2   : 1/(1+(r/w)^2)        Algebraic 1/r^2 tail.

The hedgehog part uses width ``w``; the trace part uses ``w_trace`` (defaults to
``w``). ``trace_frac`` is the amplitude ratio f0(0)/f1(0).

Usage
-----
    python -m components.initialization.trace_hedgehog \\
        --alpha 0.2 --u0 0.006 --w 10 --trace-frac 0.5 --profile-shape sech \\
        --grid-size 60 --spacing 1.0 --dt 0.05 --num-steps 10000 \\
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
class TraceHedgehogSeed:
    alpha: float
    u0: float
    w: float
    w_trace: float
    trace_frac: float
    profile_shape: str
    J: int = 0
    L: str = "1+0"
    B_winding: int = 0
    ansatz: str = "trace_hedgehog"


def _profile_gaussian(r: np.ndarray, w: float) -> np.ndarray:
    """exp(-(r/w)^2).  f(0)=1, f(inf)=0."""
    return np.exp(-((r / w) ** 2))


def _profile_sech(r: np.ndarray, w: float) -> np.ndarray:
    """1/cosh(r/w).  f(0)=1, f(inf)=0."""
    return 1.0 / np.cosh(r / w)


def _profile_power2(r: np.ndarray, w: float) -> np.ndarray:
    """1/(1+(r/w)^2).  f(0)=1, f(inf)=0."""
    return 1.0 / (1.0 + (r / w) ** 2)


_PROFILES = {
    "gaussian": _profile_gaussian,
    "sech": _profile_sech,
    "power2": _profile_power2,
}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Component 1 (Candidate 3): trace-admixture hedgehog (J=0, L=1+0)"
    )
    p.add_argument("--output", required=True, help="Path to output .npz file")

    # Physics parameters
    p.add_argument("--alpha", type=float, default=0.20,
                   help="Prestress: rest_length / spacing. Default 0.20.")
    p.add_argument("--u0", type=float, required=True,
                   help="Peak hedgehog amplitude f1(0), in lattice units.")
    p.add_argument("--w", type=float, required=True,
                   help="Hedgehog profile width in lattice units.")
    p.add_argument("--trace-frac", type=float, required=True,
                   help="Trace admixture ratio f0(0)/f1(0). 0 -> pure hedgehog (C1).")
    p.add_argument("--w-trace", type=float, default=None,
                   help="Trace profile width (defaults to --w if unset).")
    p.add_argument("--profile-shape", choices=list(_PROFILES), default="sech",
                   help="Radial profile for both f1 and f0. Default 'sech'.")

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
    w_trace: float,
    trace_frac: float,
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
) -> tuple[BraneState3D, LatticeConfig, DynamicsConfig, TraceHedgehogSeed]:
    """Build initial positions/velocities for the trace-admixture hedgehog."""
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
    seed = TraceHedgehogSeed(
        alpha=alpha,
        u0=u0,
        w=w,
        w_trace=w_trace,
        trace_frac=trace_frac,
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

    # Unit radial vector x_hat^i = dx^i / r (zero at an exact-center node)
    x_hat = dx / np.maximum(r, 1e-12)[:, np.newaxis]  # (N_nodes, 3)
    x_hat[r <= 1e-12] = 0.0

    profile = _PROFILES[profile_shape]
    f1 = profile(r, w)        # hedgehog (traceless L=1) radial profile
    f0 = profile(r, w_trace)  # trace (L=0) scalar profile

    # xi^i = u0 * [ f1 * x_hat^i + trace_frac * (1/sqrt(3)) * f0 ] ; X4 = 0
    trace_amp = u0 * trace_frac / math.sqrt(3.0)
    u_np = np.zeros((coords.shape[0], 4), dtype=np.float64)
    v_np = np.zeros((coords.shape[0], 4), dtype=np.float64)

    u_np[:, 0] = u0 * f1 * x_hat[:, 0] + trace_amp * f0
    u_np[:, 1] = u0 * f1 * x_hat[:, 1] + trace_amp * f0
    u_np[:, 2] = u0 * f1 * x_hat[:, 2] + trace_amp * f0
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

    w_trace = args.w_trace if args.w_trace is not None else args.w
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
        w_trace=w_trace,
        trace_frac=args.trace_frac,
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
            "initializer": "trace_hedgehog",
            "ansatz": "trace_hedgehog",
            "J": 0,
            "L": "1+0",
            "B_winding": 0,
            "profile_shape": args.profile_shape,
            "u0": args.u0,
            "w": args.w,
            "w_trace": w_trace,
            "trace_frac": args.trace_frac,
            "alpha": args.alpha,
            "ell0": args.alpha * args.spacing,
            "u0_over_ell0": args.u0 / (args.alpha * args.spacing),
            "w_over_a": args.w / args.spacing,
            "trace_sector": "nonzero U(1) far field (proton-class)",
            "topological_protection": False,
            "debug": {
                "grid_shape": list(grid_shape),
                "spacing": args.spacing,
                "rest_length": args.alpha * args.spacing,
            },
        },
    )

    print("Trace-admixture hedgehog initialization complete")
    print(f"  output:       {args.output}")
    print(f"  grid:         {grid_shape}")
    print(f"  alpha:        {args.alpha}")
    print(f"  u0:           {args.u0}  (u0/ell0 = {args.u0 / (args.alpha * args.spacing):.4f})")
    print(f"  w (hedgehog): {args.w} a")
    print(f"  w_trace:      {w_trace} a")
    print(f"  trace_frac:   {args.trace_frac}  (f0(0)/f1(0); 0 => pure hedgehog C1)")
    print(f"  profile:      {args.profile_shape}")
    print(f"  (J,L):        (0, 1+0)  winding B=0  -> no topological protection")
    print(f"  far field:    nonzero U(1) trace (proton-class hypothesis)")
    print(f"  device:       {device}")


if __name__ == "__main__":
    main()
