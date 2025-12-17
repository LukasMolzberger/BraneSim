"""
EMState: Tensor-based state representation for EM field simulations.

This module implements the state container for the four-potential A^μ,
storing potentials, velocities (time derivatives), and accelerations as
PyTorch tensors for efficient GPU computation.
"""

import torch
from typing import Tuple, Optional
from branesim.core.state import Dimensionality


class EMState:
    """
    Container for electromagnetic field state using four-potential formalism.

    Stores A^μ = (A^0, A^1, A^2, A^3) = (Φ/c, Ax, Ay, Az) on a spatial grid.
    Tensors are flattened as [N, 4] to match brane style, where N is the total
    number of grid points.

    Attributes:
        potential: [N, 4] tensor of A^μ components
        velocity: [N, 4] tensor of ∂_t A^μ (time derivatives)
        accel: [N, 4] tensor of ∂_tt A^μ (current accelerations)
        new_accel: [N, 4] tensor for Velocity Verlet update
        fixed_mask: [N] boolean tensor for boundary conditions
        fixed_potential: [N, 4] tensor storing fixed values
        grid_shape: Tuple specifying grid dimensions
        dimension: Dimensionality enum
        num_points: int, total number of grid points
        device: torch.device
        dtype: torch.dtype
    """

    def __init__(
        self,
        grid_shape: Tuple[int, ...],
        dimension: Dimensionality,
        device: Optional[torch.device] = None,
        dtype: torch.dtype = torch.float32,
    ):
        """
        Initialize EM state.

        Args:
            grid_shape: Tuple specifying grid dimensions, e.g., (nx,) or (nx, ny, nz)
            dimension: Dimensionality enum (ONE_D, TWO_D, or THREE_D)
            device: torch.device, defaults to 'cuda' if available else 'cpu'
            dtype: torch.dtype, defaults to float32

        Raises:
            ValueError: If grid_shape doesn't match dimension
        """
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

        # Initialize state tensors: A^μ, ∂_t A^μ, ∂_tt A^μ
        self.potential = torch.zeros((self.num_points, 4), device=self.device, dtype=dtype)
        self.velocity = torch.zeros((self.num_points, 4), device=self.device, dtype=dtype)
        self.accel = torch.zeros((self.num_points, 4), device=self.device, dtype=dtype)
        self.new_accel = torch.zeros((self.num_points, 4), device=self.device, dtype=dtype)

        # Optional fixed-mask boundaries (Dirichlet-like)
        self.fixed_mask = torch.zeros(self.num_points, device=self.device, dtype=torch.bool)
        self.fixed_potential = torch.zeros((self.num_points, 4), device=self.device, dtype=dtype)

    def view_grid(self, x: torch.Tensor) -> torch.Tensor:
        """
        Reshape flattened tensor to grid shape.

        Args:
            x: Tensor of shape [N, C] where C is number of components

        Returns:
            Tensor of shape [*grid_shape, C]
        """
        return x.view(*self.grid_shape, x.shape[-1])

    def flatten_grid(self, x: torch.Tensor) -> torch.Tensor:
        """
        Flatten grid tensor to [N, C] format.

        Args:
            x: Tensor of shape [*grid_shape, C]

        Returns:
            Tensor of shape [N, C]
        """
        return x.view(self.num_points, x.shape[-1])

    def apply_fixed_boundaries(self) -> None:
        """
        Enforce fixed boundary conditions by zeroing velocities/accelerations
        and restoring fixed potentials.
        """
        if self.fixed_mask is None:
            return
        if self.fixed_mask.any():
            self.potential[self.fixed_mask] = self.fixed_potential[self.fixed_mask]
            self.velocity[self.fixed_mask] = 0.0
            self.accel[self.fixed_mask] = 0.0
            self.new_accel[self.fixed_mask] = 0.0

    def set_fixed_from_mask(self, mask: torch.Tensor) -> None:
        """
        Set fixed boundary points from a boolean mask.

        Args:
            mask: Boolean tensor of shape [N] marking fixed points
        """
        mask = mask.to(device=self.device)
        self.fixed_mask = mask.bool()
        self.fixed_potential = self.potential.clone()

    def enforce_lorenz_gauge(self, grid, c: float, bc: str = "periodic") -> torch.Tensor:
        """
        Enforce Lorenz gauge condition: ∇·A + (1/c) ∂_t A^0 = 0.

        Simple gauge cleaning: adjust A^0_dot so that div(A) + (1/c) A^0_dot = 0.

        Args:
            grid: BraneGrid instance with spacing and dimension info
            c: Speed of light in m/s
            bc: Boundary condition type ("periodic" or "dirichlet")

        Returns:
            Gauge residual field (flattened [N]) after correction
        """
        from branesim.em.derivatives import divergence

        Agrid = self.view_grid(self.potential)     # [..., 4]
        Vgrid = self.view_grid(self.velocity)      # [..., 4]

        Avec = Agrid[..., 1:4]                     # [..., 3]
        divA = divergence(Avec, grid, bc=bc)        # [...]
        gauge = divA + (1.0 / c) * Vgrid[..., 0]    # [...]

        # Correct A^0_dot: A^0_dot <- A^0_dot - c * gauge
        Vgrid[..., 0] = Vgrid[..., 0] - c * gauge

        self.velocity = self.flatten_grid(Vgrid)
        self.apply_fixed_boundaries()

        # Recompute residual for monitoring
        Vgrid2 = self.view_grid(self.velocity)
        gauge2 = divA + (1.0 / c) * Vgrid2[..., 0]
        return gauge2.reshape(self.num_points)

    def to_device(self, device: torch.device) -> 'EMState':
        """
        Move all tensors to specified device.

        Args:
            device: Target torch.device

        Returns:
            Self for method chaining
        """
        self.device = device
        self.potential = self.potential.to(device)
        self.velocity = self.velocity.to(device)
        self.accel = self.accel.to(device)
        self.new_accel = self.new_accel.to(device)
        self.fixed_mask = self.fixed_mask.to(device)
        self.fixed_potential = self.fixed_potential.to(device)
        return self

    def clone(self) -> 'EMState':
        """
        Create a deep copy of the state.

        Returns:
            New EMState with cloned tensors
        """
        new_state = EMState(
            self.grid_shape,
            self.dimension,
            self.device,
            self.dtype
        )
        new_state.potential = self.potential.clone()
        new_state.velocity = self.velocity.clone()
        new_state.accel = self.accel.clone()
        new_state.new_accel = self.new_accel.clone()
        new_state.fixed_mask = self.fixed_mask.clone()
        new_state.fixed_potential = self.fixed_potential.clone()
        return new_state

    def __repr__(self) -> str:
        """String representation."""
        return (
            f"EMState({self.dimension.name}, "
            f"grid_shape={self.grid_shape}, "
            f"N={self.num_points}, "
            f"device={self.device}, "
            f"dtype={self.dtype})"
        )