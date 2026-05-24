"""Component 2: simulation.

Reads initial-state package and writes trajectory package.
"""

from __future__ import annotations

import argparse

import torch

from components.shared import (
    BraneState3D,
    DynamicsConfig,
    choose_device,
    load_initial_state,
    TrajectoryWriter,
)
from components.simulation import BraneGrid3D, NodeMassModel, SpringForceComputer, VelocityVerletSolver


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Component 2: simulation")
    p.add_argument("--input", required=True)
    p.add_argument("--output", required=True)
    p.add_argument("--num-steps", type=int, default=None)
    p.add_argument("--checkpoint-interval", type=int, default=None)
    p.add_argument("--device", type=str, default="auto")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    data = load_initial_state(args.input)

    lattice = data["lattice"]
    dynamics: DynamicsConfig = data["dynamics"]

    num_steps = int(args.num_steps) if args.num_steps is not None else dynamics.num_steps
    checkpoint_interval = int(args.checkpoint_interval) if args.checkpoint_interval is not None else dynamics.checkpoint_interval
    if num_steps <= 0:
        raise ValueError("num_steps must be > 0")
    if checkpoint_interval <= 0:
        raise ValueError("checkpoint_interval must be > 0")

    device = choose_device(args.device)
    dtype = torch.float32 if device.type == "mps" else torch.float64

    state = BraneState3D(lattice.grid_shape, device=device, dtype=dtype)
    state.initialize_flat_configuration(lattice.spacing)

    state.rest_positions = torch.as_tensor(data["rest_positions"], device=device, dtype=dtype)
    state.positions = torch.as_tensor(data["positions"], device=device, dtype=dtype)
    state.velocities = torch.as_tensor(data["velocities"], device=device, dtype=dtype)

    if lattice.fixed_boundaries:
        state.set_fixed_boundaries()
        state.apply_fixed_boundaries()

    grid = BraneGrid3D(
        grid_shape=lattice.grid_shape,
        spacing=lattice.spacing,
        device=device,
        periodic_axes=lattice.periodic_axes,
        shell_weights=lattice.shell_weights,
    )
    physics = SpringForceComputer(
        spring_constant=dynamics.spring_constant,
        rest_length=dynamics.rest_length,
    )
    mass_model = NodeMassModel.from_density(density=dynamics.mass_density, spacing=lattice.spacing)
    solver = VelocityVerletSolver(dt=dynamics.dt, mass_model=mass_model, physics=physics, grid=grid)
    solver.initialize_accelerations(state)

    manifest = {
        "component": "simulation",
        "source_initial_state": args.input,
        "lattice": lattice.to_dict(),
        "dynamics": {
            **dynamics.to_dict(),
            "num_steps": num_steps,
            "checkpoint_interval": checkpoint_interval,
        },
        "metadata": data.get("metadata", {}),
    }

    with TrajectoryWriter(args.output, manifest=manifest) as writer:
        writer.write_npy("aux/rest_positions.npy", state.rest_positions.detach().cpu().numpy())
        writer.write_npy("aux/grid_coords.npy", state.grid_coords.detach().cpu().numpy())

        writer.write_frame(
            step=0,
            time=0.0,
            positions=state.positions.detach().cpu().numpy(),
            velocities=state.velocities.detach().cpu().numpy(),
        )

        for step in range(1, num_steps + 1):
            solver.step(state)
            if step % checkpoint_interval == 0 or step == num_steps:
                writer.write_frame(
                    step=step,
                    time=solver.time,
                    positions=state.positions.detach().cpu().numpy(),
                    velocities=state.velocities.detach().cpu().numpy(),
                )

    saved = 1 + (num_steps // checkpoint_interval) + (1 if (num_steps % checkpoint_interval) else 0)
    print("Simulation complete")
    print(f"  output: {args.output}")
    print(f"  frames: {saved}")
    print(f"  device: {device}")


if __name__ == "__main__":
    main()
