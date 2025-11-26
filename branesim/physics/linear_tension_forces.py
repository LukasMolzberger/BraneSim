"""
LinearTensionForceComputer: Simplified force model for wave propagation.

Instead of computing full geometric spring forces with pre-tension,
this implements a linearized tension model that directly gives the wave equation:

    m · ∂²ξ/∂t² = T · ∇²ξ

For a discrete system:
    F_i = (T/h) · Σⱼ (ξⱼ - ξᵢ)

where the sum is over all neighbors j of point i.

This avoids numerical stiffness from geometric nonlinearity while maintaining
correct wave propagation physics in 1D, 2D, and 3D.
"""

import torch
import numpy as np
from branesim.core.state import BraneState, Dimensionality
from branesim.core.grid import BraneGrid


class LinearTensionForceComputer:
    """
    Linearized tension force for wave propagation in 1D/2D/3D.

    Implements F = T · ∇²ξ in discrete form, giving clean wave equation
    without geometric nonlinearity or stiffness issues.

    For each point i:
        F_i = (T/h) · Σⱼ (ξⱼ - ξᵢ)

    where j ranges over all neighbors of i.

    Attributes:
        tension: Tension T [N]
        grid_spacing: Grid spacing h [m]
    """

    def __init__(self, tension: float, grid_spacing: float):
        """
        Initialize linear tension force computer.

        Args:
            tension: Tension T [N]
            grid_spacing: Grid spacing h [m]
        """
        self.tension = tension
        self.grid_spacing = grid_spacing

        # Effective spring constant for compatibility
        self.spring_constant = tension / grid_spacing
        self.rest_length = grid_spacing  # For potential energy calculation

    def compute_forces(self, state: BraneState, grid: BraneGrid) -> torch.Tensor:
        """
        Compute linear tension forces with proper neighbor weighting.

        Forces are computed in all (D+1) dimensions where D is the brane dimensionality:
        - 1D brane: uses (x, ξ) - 2 dimensions
        - 2D brane: uses (x, y, ξ) - 3 dimensions
        - 3D brane: uses (x, y, z, ξ) - 4 dimensions

        Distances are calculated in the full (D+1)-dimensional space.

        Args:
            state: BraneState with current positions
            grid: BraneGrid with neighbor connectivity

        Returns:
            Forces tensor [N, 4] in 4D embedding space
        """
        forces = torch.zeros_like(state.positions)

        # Determine which dimensions to use based on brane dimensionality
        if grid.dimension == Dimensionality.ONE_D:
            # 1D brane: use (x, ξ) = dimensions 0 and 3
            active_dims = [0, 3]
            # Force scaling: F = T·h^(d-2) = T/h
            force_scale = self.tension / self.grid_spacing
            # 1D: 2 neighbors (left, right)
            weights = torch.ones(2, device=state.device, dtype=state.dtype)
        elif grid.dimension == Dimensionality.TWO_D:
            # 2D brane: use (x, y, ξ) = dimensions 0, 1, and 3
            active_dims = [0, 1, 3]
            # Force scaling: F = T·h^(d-2) = T
            force_scale = self.tension
            # 2D: 8 neighbors
            # Diagonal neighbors in grid space have different 3D distances
            inv_sqrt2 = 1.0 / np.sqrt(2.0)
            weights = torch.tensor([inv_sqrt2, 1.0, inv_sqrt2, 1.0, 1.0, inv_sqrt2, 1.0, inv_sqrt2],
                                   device=state.device, dtype=state.dtype)
        else:  # THREE_D
            # 3D brane: use (x, y, z, ξ) = all 4 dimensions
            active_dims = [0, 1, 2, 3]
            # Force scaling: F = T·h^(d-2) = T·h
            force_scale = self.tension * self.grid_spacing
            # 3D: 26 neighbors with different 4D distances
            inv_sqrt2 = 1.0 / np.sqrt(2.0)
            inv_sqrt3 = 1.0 / np.sqrt(3.0)
            weights_3d = []
            for di in [-1, 0, 1]:
                for dj in [-1, 0, 1]:
                    for dk in [-1, 0, 1]:
                        if not (di == 0 and dj == 0 and dk == 0):
                            num_nonzero = (di != 0) + (dj != 0) + (dk != 0)
                            if num_nonzero == 1:  # Face neighbor
                                weights_3d.append(1.0)
                            elif num_nonzero == 2:  # Edge neighbor
                                weights_3d.append(inv_sqrt2)
                            else:  # Corner neighbor
                                weights_3d.append(inv_sqrt3)
            weights = torch.tensor(weights_3d, device=state.device, dtype=state.dtype)

        # Vectorized computation over all neighbor slots
        for neighbor_idx in range(grid.max_neighbors):
            # Get neighbor indices [N]
            neighbor_ids = grid.neighbors[:, neighbor_idx]

            # Mask valid neighbors (not -1)
            valid_mask = neighbor_ids >= 0

            if not valid_mask.any():
                continue

            # Get positions in the active dimensions [N_valid, D+1]
            pos_i = state.positions[valid_mask][:, active_dims]
            pos_j = state.positions[neighbor_ids[valid_mask]][:, active_dims]

            # Compute displacement vector in (D+1) dimensions [N_valid, D+1]
            delta = pos_j - pos_i

            # Compute distance in (D+1)-dimensional space [N_valid]
            distance = torch.norm(delta, dim=1, keepdim=True)

            # Force direction (normalized displacement) [N_valid, D+1]
            force_dir = delta / (distance + 1e-10)

            # Force magnitude with proper weight and scaling [N_valid, 1]
            weight = weights[neighbor_idx]
            force_mag = weight * force_scale * distance

            # Total force in (D+1) dimensions [N_valid, D+1]
            link_forces = force_mag * force_dir

            # Accumulate forces in the active dimensions
            for dim_idx, dim in enumerate(active_dims):
                forces[valid_mask, dim] += link_forces[:, dim_idx]

        return forces

    def compute_potential_energy(self, state: BraneState, grid: BraneGrid) -> torch.Tensor:
        """
        Compute potential energy using (D+1)-dimensional distances.

        For linear tension model with proper dimensionality:
            PE = (T/2) · Σ_links |R_j - R_i|²

        where |R_j - R_i| is computed in the (D+1)-dimensional space.

        Args:
            state: BraneState
            grid: BraneGrid

        Returns:
            Total potential energy [J]
        """
        total_energy = torch.tensor(0.0, device=state.device, dtype=state.dtype)

        # Determine active dimensions based on brane dimensionality
        if grid.dimension == Dimensionality.ONE_D:
            active_dims = [0, 3]  # (x, ξ)
        elif grid.dimension == Dimensionality.TWO_D:
            active_dims = [0, 1, 3]  # (x, y, ξ)
        else:  # THREE_D
            active_dims = [0, 1, 2, 3]  # (x, y, z, ξ)

        # Sum over all neighbor connections
        for neighbor_idx in range(grid.max_neighbors):
            neighbor_ids = grid.neighbors[:, neighbor_idx]
            valid = neighbor_ids >= 0

            if not valid.any():
                continue

            # Get positions in active dimensions [N_valid, D+1]
            pos_i = state.positions[valid][:, active_dims]
            pos_j = state.positions[neighbor_ids[valid]][:, active_dims]

            # Compute displacement in (D+1) dimensions [N_valid, D+1]
            delta = pos_j - pos_i

            # Distance squared in (D+1)-dimensional space [N_valid]
            distance_sq = torch.sum(delta ** 2, dim=1)

            # Strain energy: (T/2) · |ΔR|²
            total_energy += 0.5 * self.tension / self.grid_spacing * torch.sum(distance_sq)

        # Divide by 2 to account for double-counting (each link counted from both ends)
        return total_energy / 2.0

    def __repr__(self) -> str:
        """String representation."""
        return (
            f"LinearTensionForceComputer("
            f"T={self.tension:.2e} N, "
            f"h={self.grid_spacing:.2e} m)"
        )
