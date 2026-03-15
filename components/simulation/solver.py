"""Simulation-owned Velocity Verlet integrator."""

from __future__ import annotations

from dataclasses import dataclass

import torch

from components.simulation.forces import SpringForceComputer
from components.simulation.grid import BraneGrid3D
from components.shared.state import BraneState3D


@dataclass(frozen=True)
class NodeMassModel:
    mass_per_point: float

    @staticmethod
    def from_density(density: float, spacing: float) -> "NodeMassModel":
        return NodeMassModel(mass_per_point=float(density) * float(spacing) ** 3)


class VelocityVerletSolver:
    def __init__(
        self,
        dt: float,
        mass_model: NodeMassModel,
        physics: SpringForceComputer,
        grid: BraneGrid3D,
    ):
        self.dt = float(dt)
        self.mass_model = mass_model
        self.physics = physics
        self.grid = grid
        self.time = 0.0
        self.step_count = 0

    def initialize_accelerations(self, state: BraneState3D) -> None:
        forces = self.physics.compute_forces(state, self.grid)
        state.accelerations = forces / self.mass_model.mass_per_point
        state.new_accelerations = state.accelerations.clone()

    def step(self, state: BraneState3D) -> None:
        state.positions += state.velocities * self.dt + 0.5 * state.accelerations * (self.dt ** 2)
        state.apply_fixed_boundaries()

        forces = self.physics.compute_forces(state, self.grid)
        state.new_accelerations = forces / self.mass_model.mass_per_point
        state.apply_fixed_boundaries()

        state.velocities += 0.5 * (state.accelerations + state.new_accelerations) * self.dt
        state.apply_fixed_boundaries()

        state.accelerations = state.new_accelerations.clone()

        self.time += self.dt
        self.step_count += 1

    def compute_energy(self, state: BraneState3D) -> dict[str, float]:
        kinetic = 0.5 * self.mass_model.mass_per_point * torch.sum(state.velocities ** 2)
        potential = self.physics.compute_potential_energy(state, self.grid)
        return {
            "kinetic": float(kinetic.item()),
            "potential": float(potential.item()),
            "total": float((kinetic + potential).item()),
        }
