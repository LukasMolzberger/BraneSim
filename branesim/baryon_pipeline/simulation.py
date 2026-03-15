"""Pure simulation component: evolve lattice from file input to compressed trajectory."""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path
from typing import Any

import numpy as np
import torch

from branesim.core.dimensions import MassModel
from branesim.core.grid import BraneGrid
from branesim.core.solver import VelocityVerletSolver
from branesim.core.state import BraneState, Dimensionality
from branesim.physics.forces import SpringForceComputer

from .io import CompressedTrajectoryWriter, load_initial_state_package
from .models import DynamicsConfig


def _choose_device(device: str) -> torch.device:
    if device == "auto":
        if torch.cuda.is_available():
            return torch.device("cuda")
        if hasattr(torch.backends, "mps") and torch.backends.mps.is_available():
            return torch.device("mps")
        return torch.device("cpu")
    return torch.device(device)


def run_simulation_component(
    initial_state_path: str | Path,
    output_trajectory_path: str | Path,
    *,
    num_steps: int | None = None,
    checkpoint_interval: int | None = None,
    device: str = "auto",
) -> dict[str, Any]:
    """Run lattice simulation with no diagnostics or visualization side-effects."""

    package = load_initial_state_package(initial_state_path)
    lattice = package.lattice

    dynamics = package.dynamics
    if num_steps is not None:
        dynamics = replace(dynamics, num_steps=int(num_steps))
    if checkpoint_interval is not None:
        dynamics = replace(dynamics, checkpoint_interval=int(checkpoint_interval))

    if dynamics.num_steps <= 0:
        raise ValueError("num_steps must be positive")
    if dynamics.checkpoint_interval <= 0:
        raise ValueError("checkpoint_interval must be positive")

    torch_device = _choose_device(device)

    dtype = torch.float32 if torch_device.type == "mps" else torch.float64
    state = BraneState(lattice.grid_shape, Dimensionality.THREE_D, torch_device, dtype)
    state.initialize_flat_configuration(lattice.spacing)

    state.rest_positions = torch.as_tensor(package.rest_positions, device=torch_device, dtype=dtype)
    state.positions = torch.as_tensor(package.positions, device=torch_device, dtype=dtype)
    state.velocities = torch.as_tensor(package.velocities, device=torch_device, dtype=dtype)

    if lattice.fixed_boundaries:
        state.set_fixed_boundaries()
        state.apply_fixed_boundaries()

    grid = BraneGrid(
        lattice.grid_shape,
        Dimensionality.THREE_D,
        lattice.spacing,
        torch_device,
        periodic_axes=lattice.periodic_axes,
    )

    physics = SpringForceComputer(
        spring_constant=dynamics.spring_constant,
        rest_length=dynamics.rest_length,
    )
    mass_model = MassModel.from_density(
        density=dynamics.mass_density,
        intrinsic_dim=3,
        spacing=lattice.spacing,
    )
    solver = VelocityVerletSolver(dynamics.dt, mass_model, physics, grid)
    solver.initialize_accelerations(state)

    manifest = {
        "component": "baryon_simulation",
        "lattice": lattice.to_dict(),
        "dynamics": dynamics.to_dict(),
        "source_initial_state": str(initial_state_path),
        "metadata": package.metadata,
    }

    output_trajectory_path = Path(output_trajectory_path)
    with CompressedTrajectoryWriter(output_trajectory_path, manifest=manifest) as writer:
        writer.write_numpy("aux/rest_positions.npy", package.rest_positions)
        writer.write_numpy("aux/grid_coords.npy", state.grid_coords.detach().cpu().numpy())

        writer.write_frame(
            step=0,
            time=0.0,
            positions=state.positions.detach().cpu().numpy(),
            velocities=state.velocities.detach().cpu().numpy(),
        )

        for step in range(1, dynamics.num_steps + 1):
            solver.step(state)
            if step % dynamics.checkpoint_interval == 0 or step == dynamics.num_steps:
                writer.write_frame(
                    step=step,
                    time=solver.time,
                    positions=state.positions.detach().cpu().numpy(),
                    velocities=state.velocities.detach().cpu().numpy(),
                )

    summary = {
        "output_trajectory": str(output_trajectory_path),
        "num_steps": dynamics.num_steps,
        "checkpoint_interval": dynamics.checkpoint_interval,
        "num_saved_frames": (
            1
            + (dynamics.num_steps // dynamics.checkpoint_interval)
            + (1 if (dynamics.num_steps % dynamics.checkpoint_interval) != 0 else 0)
        ),
        "device": str(torch_device),
        "dtype": str(dtype),
    }
    return summary
