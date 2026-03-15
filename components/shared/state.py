"""Shared core brane state datastructure (3D lattice in 4D embedding only)."""

from __future__ import annotations

from dataclasses import dataclass

import torch


@dataclass(frozen=True)
class BraneShape3D:
    nx: int
    ny: int
    nz: int

    @property
    def grid_shape(self) -> tuple[int, int, int]:
        return (self.nx, self.ny, self.nz)

    @property
    def num_points(self) -> int:
        return self.nx * self.ny * self.nz


class BraneState3D:
    """Simulation state container for a 3D lattice embedded in 4D."""

    def __init__(self, grid_shape: tuple[int, int, int], device: torch.device, dtype: torch.dtype):
        self.shape = BraneShape3D(*grid_shape)
        self.grid_shape = grid_shape
        self.num_points = self.shape.num_points
        self.device = device
        self.dtype = dtype

        self.positions = torch.zeros((self.num_points, 4), device=device, dtype=dtype)
        self.velocities = torch.zeros((self.num_points, 4), device=device, dtype=dtype)
        self.accelerations = torch.zeros((self.num_points, 4), device=device, dtype=dtype)
        self.new_accelerations = torch.zeros((self.num_points, 4), device=device, dtype=dtype)

        self.rest_positions: torch.Tensor | None = None
        self.grid_coords = self._build_grid_coords()

        self.fixed_mask = torch.zeros(self.num_points, device=device, dtype=torch.bool)
        self.fixed_positions = torch.zeros((self.num_points, 4), device=device, dtype=dtype)

    def _build_grid_coords(self) -> torch.Tensor:
        nx, ny, nz = self.grid_shape
        x_idx = torch.arange(nx, device=self.device)
        y_idx = torch.arange(ny, device=self.device)
        z_idx = torch.arange(nz, device=self.device)
        xx, yy, zz = torch.meshgrid(x_idx, y_idx, z_idx, indexing="ij")
        return torch.stack([xx.flatten(), yy.flatten(), zz.flatten()], dim=1).to(torch.int64)

    def initialize_flat_configuration(self, spacing: float) -> None:
        nx, ny, nz = self.grid_shape

        x_coords = torch.zeros(nx, device=self.device, dtype=self.dtype)
        y_coords = torch.zeros(ny, device=self.device, dtype=self.dtype)
        z_coords = torch.zeros(nz, device=self.device, dtype=self.dtype)

        for i in range(1, nx):
            x_coords[i] = x_coords[i - 1] + spacing
        for j in range(1, ny):
            y_coords[j] = y_coords[j - 1] + spacing
        for k in range(1, nz):
            z_coords[k] = z_coords[k - 1] + spacing

        xx, yy, zz = torch.meshgrid(x_coords, y_coords, z_coords, indexing="ij")
        self.positions[:, 0] = xx.flatten()
        self.positions[:, 1] = yy.flatten()
        self.positions[:, 2] = zz.flatten()
        self.positions[:, 3] = 0.0

        self.rest_positions = self.positions.clone()

    def set_kinematics(self, u: torch.Tensor, v: torch.Tensor) -> None:
        if self.rest_positions is None:
            raise ValueError("rest_positions is not initialized")
        self.positions = self.rest_positions + u
        self.velocities = v

    def set_fixed_boundaries(self) -> None:
        nx, ny, nz = self.grid_shape
        x = self.grid_coords[:, 0]
        y = self.grid_coords[:, 1]
        z = self.grid_coords[:, 2]
        face = (
            (x == 0)
            | (x == nx - 1)
            | (y == 0)
            | (y == ny - 1)
            | (z == 0)
            | (z == nz - 1)
        )
        self.fixed_mask[face] = True
        ids = torch.nonzero(self.fixed_mask, as_tuple=False).squeeze(1)
        self.fixed_positions.index_copy_(0, ids, torch.index_select(self.positions, 0, ids).clone())

    def apply_fixed_boundaries(self) -> None:
        if not self.fixed_mask.any():
            return
        ids = torch.nonzero(self.fixed_mask, as_tuple=False).squeeze(1)
        self.velocities.index_fill_(0, ids, 0.0)
        self.accelerations.index_fill_(0, ids, 0.0)
        self.new_accelerations.index_fill_(0, ids, 0.0)
        self.positions.index_copy_(0, ids, torch.index_select(self.fixed_positions, 0, ids))
