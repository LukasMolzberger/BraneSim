"""Component 1: baryon initialization (compressed file output)."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from branesim.baryon_pipeline import (
    BaryonSeedConfig,
    DynamicsConfig,
    LatticeConfig,
    initialize_baryon_triplet_state,
    save_initial_state_package,
)


def _triple_floats(value: str) -> tuple[float, float, float]:
    parts = [float(v.strip()) for v in value.split(",") if v.strip()]
    if len(parts) != 3:
        raise argparse.ArgumentTypeError("Expected exactly 3 comma-separated floats")
    return parts[0], parts[1], parts[2]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Initialize baryon-like spherical-harmonic state.")

    parser.add_argument("--output", required=True, help="Output .npz initial-state package")

    parser.add_argument("--nx", type=int, default=64)
    parser.add_argument("--ny", type=int, default=64)
    parser.add_argument("--nz", type=int, default=64)
    parser.add_argument("--spacing", type=float, default=1.0)
    parser.add_argument("--periodic-axes", type=str, default="false,false,false", help="Three booleans, e.g. true,true,true")
    parser.add_argument("--free-boundaries", action="store_true", help="Disable fixed boundary faces")

    parser.add_argument("--l", type=int, default=1)
    parser.add_argument("--m", type=int, default=1)
    parser.add_argument("--n", type=int, default=1)
    parser.add_argument("--radius", type=float, default=10.0)
    parser.add_argument("--amplitude", type=float, default=0.25)
    parser.add_argument("--wave-speed", type=float, default=1.0)
    parser.add_argument("--smooth-edge", type=float, default=2.0)
    parser.add_argument("--axis-amplitudes", type=_triple_floats, default=(1.0, 1.0, 1.0))
    parser.add_argument("--axis-phase-offsets", type=_triple_floats, default=(0.0, 2.0943951023931953, 4.1887902047863905))
    parser.add_argument("--mixing-strength", type=float, default=0.15)
    parser.add_argument("--x4-trace-weight", type=float, default=0.35)

    parser.add_argument("--spring-constant", type=float, default=1.0)
    parser.add_argument("--rest-length", type=float, default=0.2)
    parser.add_argument("--mass-density", type=float, default=1.0)
    parser.add_argument("--dt", type=float, default=0.01)
    parser.add_argument("--num-steps", type=int, default=200)
    parser.add_argument("--checkpoint-interval", type=int, default=1)

    parser.add_argument("--device", type=str, default="auto", help="auto/cpu/cuda/mps")
    parser.add_argument("--dtype", type=str, default="float64", choices=("float32", "float64"))

    return parser.parse_args()


def main() -> None:
    args = parse_args()

    periodic_axes = tuple(v.strip().lower() == "true" for v in args.periodic_axes.split(","))
    if len(periodic_axes) != 3:
        raise ValueError("--periodic-axes must contain 3 comma-separated boolean values")

    lattice = LatticeConfig(
        grid_shape=(args.nx, args.ny, args.nz),
        spacing=args.spacing,
        periodic_axes=periodic_axes,
        fixed_boundaries=not args.free_boundaries,
    )

    seed = BaryonSeedConfig(
        l=args.l,
        m=args.m,
        n=args.n,
        radius=args.radius,
        amplitude=args.amplitude,
        wave_speed=args.wave_speed,
        smooth_edge=args.smooth_edge,
        axis_amplitudes=args.axis_amplitudes,
        axis_phase_offsets=args.axis_phase_offsets,
        mixing_strength=args.mixing_strength,
        x4_trace_weight=args.x4_trace_weight,
    )

    dynamics = DynamicsConfig(
        spring_constant=args.spring_constant,
        rest_length=args.rest_length,
        mass_density=args.mass_density,
        dt=args.dt,
        num_steps=args.num_steps,
        checkpoint_interval=args.checkpoint_interval,
    )

    import torch

    dtype = torch.float32 if args.dtype == "float32" else torch.float64

    result = initialize_baryon_triplet_state(
        seed=seed,
        lattice=lattice,
        dynamics=dynamics,
        dtype=dtype,
        device=args.device,
    )

    output = save_initial_state_package(args.output, result.package)

    print("Baryon initialization complete")
    print(f"  output: {output}")
    print(f"  debug: {json.dumps(result.debug, indent=2)}")


if __name__ == "__main__":
    main()
