"""Simulation-owned 3D lattice topology."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import torch


@dataclass(frozen=True)
class BraneGrid3D:
    grid_shape: tuple[int, int, int]
    spacing: float
    device: torch.device
    periodic_axes: tuple[bool, bool, bool] = (False, False, False)

    def __post_init__(self):
        nx, ny, nz = self.grid_shape
        if nx <= 1 or ny <= 1 or nz <= 1:
            raise ValueError("grid_shape must have all dimensions > 1")

        offsets = []
        for di in (-1, 0, 1):
            for dj in (-1, 0, 1):
                for dk in (-1, 0, 1):
                    if not (di == 0 and dj == 0 and dk == 0):
                        offsets.append((di, dj, dk))
        object.__setattr__(self, "neighbor_offsets", torch.tensor(offsets, device=self.device, dtype=torch.int64))

        norms = torch.norm(self.neighbor_offsets.to(torch.float32), dim=1)
        object.__setattr__(self, "neighbor_offset_norms", norms)
        object.__setattr__(self, "neighbor_weights", 1.0 / (norms ** 2))

        neighbors = self._build_neighbors()
        object.__setattr__(self, "neighbors", neighbors)
        object.__setattr__(self, "max_neighbors", neighbors.shape[1])

    @property
    def num_points(self) -> int:
        nx, ny, nz = self.grid_shape
        return nx * ny * nz

    def _build_neighbors(self) -> torch.Tensor:
        nx, ny, nz = self.grid_shape
        periodic_x, periodic_y, periodic_z = self.periodic_axes

        neighbors = -np.ones((self.num_points, 26), dtype=np.int64)
        offsets = self.neighbor_offsets.detach().cpu().numpy()

        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    idx = (i * ny + j) * nz + k
                    for n_idx, (di, dj, dk) in enumerate(offsets):
                        ni, nj, nk = i + di, j + dj, k + dk

                        valid = True
                        if ni < 0 or ni >= nx:
                            if periodic_x:
                                ni = ni % nx
                            else:
                                valid = False
                        if nj < 0 or nj >= ny:
                            if periodic_y:
                                nj = nj % ny
                            else:
                                valid = False
                        if nk < 0 or nk >= nz:
                            if periodic_z:
                                nk = nk % nz
                            else:
                                valid = False

                        if valid:
                            neighbors[idx, n_idx] = (ni * ny + nj) * nz + nk

        return torch.from_numpy(neighbors).to(self.device)

    def get_spatial_coordinates(self) -> torch.Tensor:
        nx, ny, nz = self.grid_shape
        x_idx = torch.arange(nx, device=self.device, dtype=torch.float32)
        y_idx = torch.arange(ny, device=self.device, dtype=torch.float32)
        z_idx = torch.arange(nz, device=self.device, dtype=torch.float32)
        xx, yy, zz = torch.meshgrid(x_idx, y_idx, z_idx, indexing="ij")
        return torch.stack([xx.flatten(), yy.flatten(), zz.flatten()], dim=1) * self.spacing
