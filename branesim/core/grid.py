"""
BraneGrid: Grid topology and neighbor connectivity.

This module implements the grid structure and neighbor indexing for brane
simulations, supporting 1D, 2D, and 3D grids with dimension-agnostic operations.
"""

import torch
import numpy as np
from typing import Tuple
from branesim.core.state import Dimensionality


class BraneGrid:
    """
    Grid topology and neighbor connectivity manager.

    Precomputes neighbor indices for efficient force computation. Uses -1 to
    indicate invalid neighbors (boundary points).

    Attributes:
        grid_shape: Tuple specifying grid dimensions
        dimension: Dimensionality enum
        spacing: float, grid spacing in meters
        num_points: int, total number of points
        max_neighbors: int, maximum neighbors per point (2/8/26 for 1D/2D/3D)
        neighbors: [N, max_neighbors] tensor of neighbor indices (-1 for invalid)
        is_boundary: [N] boolean tensor marking boundary points
        periodic_axes: Tuple[bool, ...] indicating which axes are periodic
        device: torch.device
    """

    def __init__(
        self,
        grid_shape: Tuple[int, ...],
        dimension: Dimensionality,
        spacing: float,
        device: torch.device,
        periodic_axes: Tuple[bool, ...] = None
    ):
        """
        Initialize grid topology.

        Args:
            grid_shape: Tuple like (nx,) or (nx, ny) or (nx, ny, nz)
            dimension: Dimensionality enum
            spacing: Grid spacing h in meters
            device: torch.device for tensor storage
            periodic_axes: Tuple of booleans indicating which axes are periodic
                          (default: None = all non-periodic)
                          Example: (True, False, False) makes x-axis periodic
        """
        self.grid_shape = grid_shape
        self.dimension = dimension
        self.spacing = spacing
        self.device = device
        self.num_points = int(np.prod(grid_shape))

        # Set periodic axes (default: all non-periodic)
        ndim = len(grid_shape)
        if periodic_axes is None:
            self.periodic_axes = tuple([False] * ndim)
        else:
            if len(periodic_axes) != ndim:
                raise ValueError(f"periodic_axes length {len(periodic_axes)} must match grid dimensionality {ndim}")
            self.periodic_axes = tuple(periodic_axes)

        # Set max neighbors based on dimension
        if dimension == Dimensionality.ONE_D:
            self.max_neighbors = 2  # Left and right
        elif dimension == Dimensionality.TWO_D:
            self.max_neighbors = 8  # 3x3 - 1 (exclude center)
        else:  # THREE_D
            self.max_neighbors = 26  # 3x3x3 - 1

        # Build neighbor connectivity
        self.neighbors = self._build_neighbor_indices()

        # Compute boundary mask
        self.is_boundary = self._compute_boundary_mask()

    def _build_neighbor_indices(self) -> torch.Tensor:
        """
        Build neighbor index tensor [N, max_neighbors].

        Uses -1 for invalid neighbors (out of bounds).

        Returns:
            Tensor of neighbor indices
        """
        if self.dimension == Dimensionality.ONE_D:
            return self._build_neighbors_1d()
        elif self.dimension == Dimensionality.TWO_D:
            return self._build_neighbors_2d()
        else:  # THREE_D
            return self._build_neighbors_3d()

    def _build_neighbors_1d(self) -> torch.Tensor:
        """
        Build neighbors for 1D grid.

        Neighbors: [i-1, i+1]
        With periodic BCs, wraps around at boundaries.

        Returns:
            [N, 2] tensor of neighbor indices
        """
        nx = self.grid_shape[0]
        neighbors = -np.ones((nx, 2), dtype=np.int64)
        periodic_x = self.periodic_axes[0]

        for i in range(nx):
            # Left neighbor
            if i > 0:
                neighbors[i, 0] = i - 1
            elif periodic_x:
                neighbors[i, 0] = nx - 1  # Wrap to right end

            # Right neighbor
            if i < nx - 1:
                neighbors[i, 1] = i + 1
            elif periodic_x:
                neighbors[i, 1] = 0  # Wrap to left end

        return torch.from_numpy(neighbors).to(self.device)

    def _build_neighbors_2d(self) -> torch.Tensor:
        """
        Build neighbors for 2D grid.

        Neighbors: 8-connectivity (3x3 - 1).
        Order: [(-1,-1), (-1,0), (-1,1), (0,-1), (0,1), (1,-1), (1,0), (1,1)]

        Supports periodic boundary conditions via self.periodic_axes.

        Returns:
            [N, 8] tensor of neighbor indices
        """
        nx, ny = self.grid_shape
        N = nx * ny
        neighbors = -np.ones((N, 8), dtype=np.int64)

        periodic_x, periodic_y = self.periodic_axes

        # Neighbor offsets (di, dj)
        offsets = [
            (-1, -1), (-1, 0), (-1, 1),
            (0, -1),           (0, 1),
            (1, -1),  (1, 0),  (1, 1)
        ]

        for i in range(nx):
            for j in range(ny):
                idx = i * ny + j

                for k, (di, dj) in enumerate(offsets):
                    ni, nj = i + di, j + dj

                    # Handle periodic boundaries
                    valid = True
                    if ni < 0 or ni >= nx:
                        if periodic_x:
                            ni = ni % nx  # Wrap around
                        else:
                            valid = False
                    if nj < 0 or nj >= ny:
                        if periodic_y:
                            nj = nj % ny  # Wrap around
                        else:
                            valid = False

                    if valid:
                        neighbors[idx, k] = ni * ny + nj

        return torch.from_numpy(neighbors).to(self.device)

    def _build_neighbors_3d(self) -> torch.Tensor:
        """
        Build neighbors for 3D grid.

        Neighbors: 26-connectivity (3x3x3 - 1).
        All combinations of (di, dj, dk) ∈ {-1, 0, 1}³ except (0, 0, 0).

        Supports periodic boundary conditions via self.periodic_axes.

        Returns:
            [N, 26] tensor of neighbor indices
        """
        nx, ny, nz = self.grid_shape
        N = nx * ny * nz
        neighbors = -np.ones((N, 26), dtype=np.int64)

        periodic_x, periodic_y, periodic_z = self.periodic_axes

        # Generate all 26 offsets
        offsets = []
        for di in [-1, 0, 1]:
            for dj in [-1, 0, 1]:
                for dk in [-1, 0, 1]:
                    if not (di == 0 and dj == 0 and dk == 0):
                        offsets.append((di, dj, dk))

        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    idx = (i * ny + j) * nz + k

                    for n_idx, (di, dj, dk) in enumerate(offsets):
                        ni, nj, nk = i + di, j + dj, k + dk

                        # Handle periodic boundaries
                        valid = True
                        if ni < 0 or ni >= nx:
                            if periodic_x:
                                ni = ni % nx  # Wrap around
                            else:
                                valid = False
                        if nj < 0 or nj >= ny:
                            if periodic_y:
                                nj = nj % ny  # Wrap around
                            else:
                                valid = False
                        if nk < 0 or nk >= nz:
                            if periodic_z:
                                nk = nk % nz  # Wrap around
                            else:
                                valid = False

                        if valid:
                            neighbors[idx, n_idx] = (ni * ny + nj) * nz + nk

        return torch.from_numpy(neighbors).to(self.device)

    def _compute_boundary_mask(self) -> torch.Tensor:
        """
        Compute boolean mask [N] for boundary points.

        Returns:
            Boolean tensor marking boundary points
        """
        if self.dimension == Dimensionality.ONE_D:
            mask = np.zeros(self.grid_shape[0], dtype=bool)
            mask[0] = mask[-1] = True

        elif self.dimension == Dimensionality.TWO_D:
            nx, ny = self.grid_shape
            mask = np.zeros((nx, ny), dtype=bool)
            mask[0, :] = mask[-1, :] = True
            mask[:, 0] = mask[:, -1] = True
            mask = mask.flatten()

        else:  # THREE_D
            nx, ny, nz = self.grid_shape
            mask = np.zeros((nx, ny, nz), dtype=bool)
            mask[0, :, :] = mask[-1, :, :] = True
            mask[:, 0, :] = mask[:, -1, :] = True
            mask[:, :, 0] = mask[:, :, -1] = True
            mask = mask.flatten()

        return torch.from_numpy(mask).to(self.device)

    def get_spatial_coordinates(self) -> torch.Tensor:
        """
        Get spatial coordinates [N, D] in physical units (meters).

        Returns:
            Tensor of shape [N, D] with coordinates multiplied by spacing
        """
        if self.dimension == Dimensionality.ONE_D:
            coords = torch.arange(self.grid_shape[0], device=self.device, dtype=torch.float32)
            coords = coords.unsqueeze(1) * self.spacing

        elif self.dimension == Dimensionality.TWO_D:
            nx, ny = self.grid_shape
            x_idx = torch.arange(nx, device=self.device, dtype=torch.float32)
            y_idx = torch.arange(ny, device=self.device, dtype=torch.float32)
            xx, yy = torch.meshgrid(x_idx, y_idx, indexing='ij')
            coords = torch.stack([xx.flatten(), yy.flatten()], dim=1) * self.spacing

        else:  # THREE_D
            nx, ny, nz = self.grid_shape
            x_idx = torch.arange(nx, device=self.device, dtype=torch.float32)
            y_idx = torch.arange(ny, device=self.device, dtype=torch.float32)
            z_idx = torch.arange(nz, device=self.device, dtype=torch.float32)
            xx, yy, zz = torch.meshgrid(x_idx, y_idx, z_idx, indexing='ij')
            coords = torch.stack([xx.flatten(), yy.flatten(), zz.flatten()], dim=1) * self.spacing

        return coords

    def get_num_valid_neighbors(self, point_idx: int) -> int:
        """
        Count valid neighbors for a given point.

        Args:
            point_idx: Index of the point

        Returns:
            Number of valid (non -1) neighbors
        """
        return int((self.neighbors[point_idx] >= 0).sum().item())

    def __repr__(self) -> str:
        """String representation."""
        return (
            f"BraneGrid({self.dimension.name}, "
            f"grid_shape={self.grid_shape}, "
            f"spacing={self.spacing:.2e}m, "
            f"N={self.num_points}, "
            f"max_neighbors={self.max_neighbors})"
        )
