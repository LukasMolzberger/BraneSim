"""Simulation-owned spring force implementation."""

from __future__ import annotations

import torch

from components.simulation.grid import BraneGrid3D
from components.shared.state import BraneState3D


class SpringForceComputer:
    def __init__(self, spring_constant: float, rest_length: float):
        self.spring_constant = float(spring_constant)
        self.rest_length = float(rest_length)

    def _bond_delta(
        self,
        state: BraneState3D,
        grid: BraneGrid3D,
        ids: torch.Tensor,
        neighbor_idx: int,
        p_pos: torch.Tensor,
        q_pos: torch.Tensor,
    ) -> torch.Tensor:
        delta = q_pos - p_pos
        if not any(grid.periodic_axes):
            return delta

        offset = grid.neighbor_offsets[neighbor_idx]
        coords = state.grid_coords[ids]
        image_shift = torch.zeros_like(delta)

        for axis, is_periodic in enumerate(grid.periodic_axes):
            if not is_periodic:
                continue

            step = int(offset[axis].item())
            if step == -1:
                mask = coords[:, axis] == 0
                shift = -float(grid.grid_shape[axis]) * float(grid.spacing)
            elif step == 1:
                mask = coords[:, axis] == grid.grid_shape[axis] - 1
                shift = float(grid.grid_shape[axis]) * float(grid.spacing)
            else:
                continue

            image_shift[mask, axis] = shift

        return delta + image_shift

    def compute_forces(self, state: BraneState3D, grid: BraneGrid3D) -> torch.Tensor:
        forces = torch.zeros_like(state.positions)

        for neighbor_idx in range(grid.max_neighbors):
            neighbor_ids = grid.neighbors[:, neighbor_idx]
            valid_mask = neighbor_ids >= 0
            if not valid_mask.any():
                continue

            ids = torch.nonzero(valid_mask, as_tuple=False).squeeze(1)
            p_pos = state.positions[ids]
            q_pos = state.positions[neighbor_ids[ids]]

            delta = self._bond_delta(state, grid, ids, neighbor_idx, p_pos, q_pos)
            distance = torch.norm(delta, dim=1, keepdim=True)

            rest_length = self.rest_length * grid.neighbor_offset_norms[neighbor_idx].to(state.dtype)
            strain = distance - rest_length

            k_eff = self.spring_constant * grid.neighbor_weights[neighbor_idx].to(state.dtype)
            force_mag = k_eff * strain
            force_dir = delta / (distance + 1e-20)
            link_forces = force_mag * force_dir

            forces.index_add_(0, ids, link_forces)

        return forces

    def compute_potential_energy(self, state: BraneState3D, grid: BraneGrid3D) -> torch.Tensor:
        total_energy = torch.tensor(0.0, device=state.device, dtype=state.dtype)

        for neighbor_idx in range(grid.max_neighbors):
            neighbor_ids = grid.neighbors[:, neighbor_idx]
            valid_mask = neighbor_ids >= 0
            if not valid_mask.any():
                continue

            ids = torch.nonzero(valid_mask, as_tuple=False).squeeze(1)
            p_pos = state.positions[ids]
            q_pos = state.positions[neighbor_ids[ids]]

            delta = self._bond_delta(state, grid, ids, neighbor_idx, p_pos, q_pos)
            distance = torch.norm(delta, dim=1)
            rest_length = self.rest_length * grid.neighbor_offset_norms[neighbor_idx].to(state.dtype)
            strain = distance - rest_length

            k_eff = self.spring_constant * grid.neighbor_weights[neighbor_idx].to(state.dtype)
            link_energy = 0.5 * k_eff * strain ** 2
            total_energy += torch.sum(link_energy)

        return total_energy / 2.0
