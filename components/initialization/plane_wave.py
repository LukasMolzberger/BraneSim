"""Plane-wave initializer for dispersion / isotropy experiments (sprint 1).

Sets up a low-amplitude standing-wave perturbation `u(x, 0) = ε p̂ cos(k·x)`
with `v(x, 0) = 0` on a periodic cubic lattice. The standing wave is a
superposition of forward and backward traveling waves with the same |k|, so
the amplitude at fixed wavevector oscillates at the dispersion frequency
`ω(k)`. Diagnostics extract `ω` from that oscillation.

Independently callable:
    python -m components.initialization.plane_wave --output ... [args]

Output is the same `initial_state.npz` format consumed by
`components.simulation.run`.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass

import numpy as np
import torch

from components.shared import (
    BraneState3D,
    DynamicsConfig,
    LatticeConfig,
    choose_device,
    choose_dtype,
    parse_bool_triple,
    save_initial_state,
)


@dataclass(frozen=True)
class PlaneWaveSeed:
    k_index: tuple[int, int, int]
    polarization: tuple[float, float, float]
    amplitude: float

    def to_dict(self) -> dict:
        return {
            "kind": "plane_wave",
            "k_index": list(self.k_index),
            "polarization": list(self.polarization),
            "amplitude": float(self.amplitude),
        }


def _triple_floats(value: str) -> tuple[float, float, float]:
    parts = [float(v.strip()) for v in value.split(",") if v.strip()]
    if len(parts) != 3:
        raise argparse.ArgumentTypeError("Expected 3 comma-separated floats")
    return parts[0], parts[1], parts[2]


def _triple_ints(value: str) -> tuple[int, int, int]:
    parts = [int(v.strip()) for v in value.split(",") if v.strip()]
    if len(parts) != 3:
        raise argparse.ArgumentTypeError("Expected 3 comma-separated ints")
    return parts[0], parts[1], parts[2]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Plane-wave initializer (dispersion / isotropy)")
    p.add_argument("--output", required=True)

    p.add_argument("--nx", type=int, default=32)
    p.add_argument("--ny", type=int, default=32)
    p.add_argument("--nz", type=int, default=32)
    p.add_argument("--spacing", type=float, default=1.0)
    p.add_argument("--periodic-axes", type=str, default="true,true,true")
    p.add_argument("--free-boundaries", action="store_true",
                   help="If set, do NOT clamp boundary nodes; needed for periodic dispersion measurements.")

    p.add_argument("--spring-constant", type=float, required=True)
    p.add_argument("--rest-length", type=float, required=True)
    p.add_argument("--mass-density", type=float, required=True)
    p.add_argument("--dt", type=float, required=True)
    p.add_argument("--num-steps", type=int, required=True)
    p.add_argument("--checkpoint-interval", type=int, default=1)

    p.add_argument("--k-index", type=_triple_ints, required=True,
                   help="Wavevector in lattice mode units: k = 2π * (kx/Lx, ky/Ly, kz/Lz). "
                        "Use comma-separated ints, e.g. '2,0,0' for [100] direction at mode 2.")
    p.add_argument("--polarization", type=_triple_floats, required=True,
                   help="Polarization unit vector p̂ for the lateral displacement. "
                        "Comma-separated floats; will be normalized.")
    p.add_argument("--amplitude", type=float, default=1e-3,
                   help="Displacement amplitude ε. Keep in linear regime: ε·|k|·spacing ≪ 1.")

    p.add_argument("--device", type=str, default="auto")
    p.add_argument("--dtype", type=str, default="float64", choices=("float32", "float64"))
    return p.parse_args()


def main() -> None:
    args = parse_args()

    periodic = parse_bool_triple(args.periodic_axes)
    if not all(periodic) and not args.free_boundaries:
        # Fixed boundaries on a non-periodic axis would clamp the wave at the edges.
        # Allow it (user might want a closed-box mode), but warn loudly.
        print("WARNING: non-periodic axis with fixed boundaries will distort plane-wave dispersion.")

    lattice = LatticeConfig(
        grid_shape=(args.nx, args.ny, args.nz),
        spacing=float(args.spacing),
        periodic_axes=periodic,
        fixed_boundaries=not args.free_boundaries,
    )
    dynamics = DynamicsConfig(
        spring_constant=float(args.spring_constant),
        rest_length=float(args.rest_length),
        mass_density=float(args.mass_density),
        dt=float(args.dt),
        num_steps=int(args.num_steps),
        checkpoint_interval=int(args.checkpoint_interval),
    )

    polarization = np.asarray(args.polarization, dtype=np.float64)
    pol_norm = float(np.linalg.norm(polarization))
    if pol_norm < 1e-12:
        raise ValueError("polarization vector has zero norm")
    polarization = polarization / pol_norm

    seed = PlaneWaveSeed(
        k_index=tuple(int(k) for k in args.k_index),
        polarization=(float(polarization[0]), float(polarization[1]), float(polarization[2])),
        amplitude=float(args.amplitude),
    )

    device = choose_device(args.device)
    dtype = choose_dtype(args.dtype, device)

    state = BraneState3D(lattice.grid_shape, device=device, dtype=dtype)
    state.initialize_flat_configuration(lattice.spacing)
    if lattice.fixed_boundaries:
        state.set_fixed_boundaries()

    coords = state.positions[:, :3].detach().cpu().to(torch.float64).numpy()
    nx, ny, nz = lattice.grid_shape
    Lx, Ly, Lz = nx * lattice.spacing, ny * lattice.spacing, nz * lattice.spacing
    kx = 2.0 * np.pi * seed.k_index[0] / Lx
    ky = 2.0 * np.pi * seed.k_index[1] / Ly
    kz = 2.0 * np.pi * seed.k_index[2] / Lz
    k_vec = np.array([kx, ky, kz], dtype=np.float64)

    phase = coords @ k_vec
    envelope = seed.amplitude * np.cos(phase)

    u_np = np.zeros((coords.shape[0], 4), dtype=np.float64)
    v_np = np.zeros((coords.shape[0], 4), dtype=np.float64)
    u_np[:, 0] = polarization[0] * envelope
    u_np[:, 1] = polarization[1] * envelope
    u_np[:, 2] = polarization[2] * envelope
    # X4 (amplitude direction) stays at zero — we are testing the in-brane lateral elasticity only.

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
        seed=seed.to_dict(),
        metadata={
            "component": "initialization",
            "seed_kind": "plane_wave",
            "debug": {
                "k_vec": [float(kx), float(ky), float(kz)],
                "k_magnitude": float(np.linalg.norm(k_vec)),
                "k_index": list(seed.k_index),
                "polarization": list(polarization),
                "linear_parameter_ka": float(np.linalg.norm(k_vec) * lattice.spacing),
                "amplitude_strain": float(seed.amplitude * np.linalg.norm(k_vec)),
            },
        },
    )

    print("Plane-wave initialization complete")
    print(f"  output: {args.output}")
    print(f"  grid: {lattice.grid_shape}, spacing={lattice.spacing}")
    print(f"  periodic: {periodic}, fixed_boundaries={lattice.fixed_boundaries}")
    print(f"  k_index: {seed.k_index}, |k|·a = {np.linalg.norm(k_vec) * lattice.spacing:.4f}")
    print(f"  polarization: {tuple(round(v, 4) for v in polarization)}")
    print(f"  amplitude: {seed.amplitude}, max strain ε|k| = {seed.amplitude * np.linalg.norm(k_vec):.2e}")
    print(f"  device: {device}")


if __name__ == "__main__":
    main()
