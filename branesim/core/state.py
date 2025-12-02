"""
BraneState: Tensor-based state representation for brane simulations.

This module implements the state container for brane points, storing positions,
velocities, and accelerations as PyTorch tensors for efficient GPU computation.
"""

import torch
from enum import Enum
from typing import Tuple, Optional


class Dimensionality(Enum):
    """Spatial dimensionality of the brane."""
    ONE_D = 1
    TWO_D = 2
    THREE_D = 3


class BraneState:
    """
    Container for brane simulation state.

    Stores the state of all brane points as PyTorch tensors. Each point has
    a 4D embedding position (X^0, X^1, X^2, X^3) where the fourth component
    represents the amplitude/normal displacement from the base 3D brane.

    Attributes:
        positions: [N, 4] tensor of 4D embedding coordinates
        velocities: [N, 4] tensor of velocity vectors
        accelerations: [N, 4] tensor of current accelerations
        new_accelerations: [N, 4] tensor for Velocity Verlet update
        grid_coords: [N, D] tensor of grid coordinates (D = 1, 2, or 3)
        dimension: Dimensionality enum
        num_points: int, total number of grid points
        device: torch.device (cpu or cuda)
        dtype: torch.dtype (float32 or float64)
    """

    def __init__(
        self,
        grid_shape: Tuple[int, ...],
        dimension: Dimensionality,
        device: Optional[torch.device] = None,
        dtype: torch.dtype = torch.float32
    ):
        """
        Initialize brane state.

        Args:
            grid_shape: Tuple specifying grid dimensions, e.g., (nx,) or (nx, ny, nz)
            dimension: Dimensionality enum (ONE_D, TWO_D, or THREE_D)
            device: torch.device, defaults to 'cuda' if available else 'cpu'
            dtype: torch.dtype, defaults to float32

        Raises:
            ValueError: If grid_shape doesn't match dimension
        """
        # Validate inputs
        if len(grid_shape) != dimension.value:
            raise ValueError(
                f"Grid shape {grid_shape} has {len(grid_shape)} dimensions, "
                f"but dimensionality is {dimension.value}D"
            )

        self.grid_shape = grid_shape
        self.dimension = dimension
        self.num_points = int(torch.prod(torch.tensor(grid_shape)).item())

        # Set device
        if device is None:
            self.device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
        else:
            self.device = device

        self.dtype = dtype

        # Initialize state tensors (all zeros initially)
        self.positions = torch.zeros((self.num_points, 4), device=self.device, dtype=dtype)
        self.velocities = torch.zeros((self.num_points, 4), device=self.device, dtype=dtype)
        self.accelerations = torch.zeros((self.num_points, 4), device=self.device, dtype=dtype)
        self.new_accelerations = torch.zeros((self.num_points, 4), device=self.device, dtype=dtype)

        # Grid coordinates [N, D] for topology lookups
        self.grid_coords = self._create_grid_coords()

        # Fixed boundary mask [N] - True for points that should not move
        self.fixed_mask = torch.zeros(self.num_points, device=self.device, dtype=torch.bool)

        # Store initial positions for fixed points
        self.fixed_positions = torch.zeros((self.num_points, 4), device=self.device, dtype=dtype)

    def _create_grid_coords(self) -> torch.Tensor:
        """
        Create grid coordinate tensor [N, D].

        For a 2D grid (nx, ny), creates coordinates:
        [[0, 0], [0, 1], ..., [0, ny-1], [1, 0], [1, 1], ..., [nx-1, ny-1]]

        Returns:
            Tensor of shape [N, D] with integer grid coordinates
        """
        if self.dimension == Dimensionality.ONE_D:
            # 1D: just indices 0 to nx-1
            coords = torch.arange(self.grid_shape[0], device=self.device).unsqueeze(1)

        elif self.dimension == Dimensionality.TWO_D:
            # 2D: meshgrid and flatten
            nx, ny = self.grid_shape
            x_idx = torch.arange(nx, device=self.device)
            y_idx = torch.arange(ny, device=self.device)
            xx, yy = torch.meshgrid(x_idx, y_idx, indexing='ij')
            coords = torch.stack([xx.flatten(), yy.flatten()], dim=1)

        else:  # THREE_D
            # 3D: meshgrid and flatten
            nx, ny, nz = self.grid_shape
            x_idx = torch.arange(nx, device=self.device)
            y_idx = torch.arange(ny, device=self.device)
            z_idx = torch.arange(nz, device=self.device)
            xx, yy, zz = torch.meshgrid(x_idx, y_idx, z_idx, indexing='ij')
            coords = torch.stack([xx.flatten(), yy.flatten(), zz.flatten()], dim=1)

        return coords.to(torch.int64)

    def to_device(self, device: torch.device) -> 'BraneState':
        """
        Move all tensors to specified device.

        Args:
            device: Target torch.device

        Returns:
            Self for method chaining
        """
        self.device = device
        self.positions = self.positions.to(device)
        self.velocities = self.velocities.to(device)
        self.accelerations = self.accelerations.to(device)
        self.new_accelerations = self.new_accelerations.to(device)
        self.grid_coords = self.grid_coords.to(device)
        self.fixed_mask = self.fixed_mask.to(device)
        self.fixed_positions = self.fixed_positions.to(device)
        return self

    def clone(self) -> 'BraneState':
        """
        Create a deep copy of the state.

        Returns:
            New BraneState with cloned tensors
        """
        new_state = BraneState(
            self.grid_shape,
            self.dimension,
            self.device,
            self.dtype
        )
        new_state.positions = self.positions.clone()
        new_state.velocities = self.velocities.clone()
        new_state.accelerations = self.accelerations.clone()
        new_state.new_accelerations = self.new_accelerations.clone()
        new_state.grid_coords = self.grid_coords.clone()
        new_state.fixed_mask = self.fixed_mask.clone()
        new_state.fixed_positions = self.fixed_positions.clone()
        return new_state

    def get_field_component(self, component_idx: int) -> torch.Tensor:
        """
        Extract single embedding component for visualization.

        Args:
            component_idx: Index 0-3 for X^0, X^1, X^2, X^3

        Returns:
            Tensor of shape [N] with component values
        """
        if component_idx < 0 or component_idx >= 4:
            raise ValueError(f"Component index must be 0-3, got {component_idx}")

        return self.positions[:, component_idx]

    def get_kinetic_energy(self) -> torch.Tensor:
        """
        Compute total kinetic energy (not per-unit-mass, just velocity magnitude).

        Returns:
            Scalar tensor with kinetic energy: 0.5 * Σ |v_p|²
        """
        return 0.5 * torch.sum(self.velocities ** 2)

    def get_velocity_magnitude(self) -> torch.Tensor:
        """
        Compute velocity magnitude for each point.

        Returns:
            Tensor of shape [N] with |v_p| for each point
        """
        return torch.norm(self.velocities, dim=1)

    def initialize_flat_configuration(self, spacing: float):
        """
        Initialize positions to flat grid configuration.

        Sets the first D components (X^0, ..., X^{D-1}) to spatial grid positions
        with given spacing, and leaves the amplitude (X^3) at zero.

        Uses explicit additive loop: position[i] = position[i-1] + spacing
        This ensures maximum uniformity - all distances between adjacent nodes
        are the same bit pattern, minimizing force imbalances at t=0.

        Args:
            spacing: Grid spacing in meters
        """
        if self.dimension == Dimensionality.ONE_D:
            # Build positions additively in loop to ensure uniform spacing
            nx = self.grid_shape[0]
            self.positions[0, 0] = 0.0
            for i in range(1, nx):
                self.positions[i, 0] = self.positions[i-1, 0] + spacing
            # X^1, X^2, X^3 remain zero

        elif self.dimension == Dimensionality.TWO_D:
            # Build coordinate arrays additively
            nx, ny = self.grid_shape

            # X-coordinates via loop
            x_coords = torch.zeros(nx, device=self.device, dtype=self.dtype)
            for i in range(1, nx):
                x_coords[i] = x_coords[i-1] + spacing

            # Y-coordinates via loop
            y_coords = torch.zeros(ny, device=self.device, dtype=self.dtype)
            for j in range(1, ny):
                y_coords[j] = y_coords[j-1] + spacing

            # Create meshgrid and flatten
            xx, yy = torch.meshgrid(x_coords, y_coords, indexing='ij')
            self.positions[:, 0] = xx.flatten()
            self.positions[:, 1] = yy.flatten()
            # X^2, X^3 remain zero

        else:  # THREE_D
            # Build coordinate arrays additively
            nx, ny, nz = self.grid_shape

            # X-coordinates via loop
            x_coords = torch.zeros(nx, device=self.device, dtype=self.dtype)
            for i in range(1, nx):
                x_coords[i] = x_coords[i-1] + spacing

            # Y-coordinates via loop
            y_coords = torch.zeros(ny, device=self.device, dtype=self.dtype)
            for j in range(1, ny):
                y_coords[j] = y_coords[j-1] + spacing

            # Z-coordinates via loop
            z_coords = torch.zeros(nz, device=self.device, dtype=self.dtype)
            for k in range(1, nz):
                z_coords[k] = z_coords[k-1] + spacing

            # Create meshgrid and flatten
            xx, yy, zz = torch.meshgrid(x_coords, y_coords, z_coords, indexing='ij')
            self.positions[:, 0] = xx.flatten()
            self.positions[:, 1] = yy.flatten()
            self.positions[:, 2] = zz.flatten()
            # X^3 remains zero

    def set_fixed_boundaries(self):
        """
        Set boundary points as fixed (immovable).

        For 1D: fixes the first and last points
        For 2D: fixes all edge points
        For 3D: fixes all face points

        Stores current positions as the fixed positions.
        """
        if self.dimension == Dimensionality.ONE_D:
            # Fix first and last points
            self.fixed_mask[0] = True
            self.fixed_mask[-1] = True

        elif self.dimension == Dimensionality.TWO_D:
            # Fix all edge points (vectorized)
            nx, ny = self.grid_shape
            x_coords = self.grid_coords[:, 0]
            y_coords = self.grid_coords[:, 1]
            edge_mask = (x_coords == 0) | (x_coords == nx - 1) | (y_coords == 0) | (y_coords == ny - 1)
            self.fixed_mask[edge_mask] = True

        else:  # THREE_D
            # Fix all face points (vectorized)
            nx, ny, nz = self.grid_shape
            x_coords = self.grid_coords[:, 0]
            y_coords = self.grid_coords[:, 1]
            z_coords = self.grid_coords[:, 2]
            face_mask = (
                (x_coords == 0) | (x_coords == nx - 1) |
                (y_coords == 0) | (y_coords == ny - 1) |
                (z_coords == 0) | (z_coords == nz - 1)
            )
            self.fixed_mask[face_mask] = True

        # Store initial positions for fixed points
        self.fixed_positions[self.fixed_mask] = self.positions[self.fixed_mask].clone()

    def apply_fixed_boundaries(self):
        """
        Enforce fixed boundary conditions by zeroing velocities/accelerations
        and restoring positions for fixed points.
        """
        if self.fixed_mask.any():
            # Zero out dynamics for fixed points
            self.velocities[self.fixed_mask] = 0.0
            self.accelerations[self.fixed_mask] = 0.0
            self.new_accelerations[self.fixed_mask] = 0.0
            # Restore fixed positions
            self.positions[self.fixed_mask] = self.fixed_positions[self.fixed_mask]

    def __repr__(self) -> str:
        """String representation."""
        return (
            f"BraneState({self.dimension.name}, "
            f"grid_shape={self.grid_shape}, "
            f"N={self.num_points}, "
            f"device={self.device}, "
            f"dtype={self.dtype})"
        )
