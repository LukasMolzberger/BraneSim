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

        For discrete Laplacian:
        - 1D: 2 cardinal neighbors, weight = 1.0
        - 2D: 4 cardinal neighbors (weight = 1.0), 4 diagonal (weight = 0.5)
        - 3D: 6 cardinal neighbors (weight = 1.0), 12 edge/face diagonal (weight varies)

        This approximates: ∂²ξ/∂t² = (T/ρ) · ∇²ξ

        Args:
            state: BraneState with current positions
            grid: BraneGrid with neighbor connectivity

        Returns:
            Forces tensor [N, 4] in 4D embedding space
        """
        forces = torch.zeros_like(state.positions)

        # Only compute for transverse direction (ξ³)
        xi = state.positions[:, 3]

        # Force scaling: F = T·h^(d-2) · Σⱼ (ξⱼ - ξᵢ) where d = dimensionality
        # This accounts for mass scaling as h^d and Laplacian as 1/h²
        if grid.dimension == Dimensionality.ONE_D:
            # 1D: force_scale = T·h^(1-2) = T/h
            force_scale = self.tension / self.grid_spacing
            # 1D: 2 neighbors (left, right) at distance h
            weights = torch.ones(2, device=state.device, dtype=state.dtype)
        elif grid.dimension == Dimensionality.TWO_D:
            # 2D: force_scale = T·h^(2-2) = T
            force_scale = self.tension
            # 2D: 8 neighbors in order [(-1,-1), (-1,0), (-1,1), (0,-1), (0,1), (1,-1), (1,0), (1,1)]
            # Diagonal neighbors at distance √2·h: indices 0, 2, 5, 7 (weight 1/√2)
            # Cardinal neighbors at distance h: indices 1, 3, 4, 6 (weight 1.0)
            inv_sqrt2 = 1.0 / np.sqrt(2.0)
            weights = torch.tensor([inv_sqrt2, 1.0, inv_sqrt2, 1.0, 1.0, inv_sqrt2, 1.0, inv_sqrt2],
                                   device=state.device, dtype=state.dtype)
        else:  # THREE_D
            # 3D: force_scale = T·h^(3-2) = T·h
            force_scale = self.tension * self.grid_spacing
            # 3D: 26 neighbors at different distances
            # Face neighbors (distance h): weight 1.0
            # Edge neighbors (distance √2·h): weight 1/√2
            # Corner neighbors (distance √3·h): weight 1/√3
            # Order from grid: all combos of (di,dj,dk) ∈ {-1,0,1}³ except (0,0,0)
            inv_sqrt2 = 1.0 / np.sqrt(2.0)
            inv_sqrt3 = 1.0 / np.sqrt(3.0)
            weights_3d = []
            for di in [-1, 0, 1]:
                for dj in [-1, 0, 1]:
                    for dk in [-1, 0, 1]:
                        if not (di == 0 and dj == 0 and dk == 0):
                            num_nonzero = (di != 0) + (dj != 0) + (dk != 0)
                            if num_nonzero == 1:  # Face neighbor (distance h)
                                weights_3d.append(1.0)
                            elif num_nonzero == 2:  # Edge neighbor (distance √2·h)
                                weights_3d.append(inv_sqrt2)
                            else:  # Corner neighbor (distance √3·h, num_nonzero == 3)
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

            # For valid neighbors: F += weight · force_scale · (ξⱼ - ξᵢ)
            xi_i = xi[valid_mask]
            xi_j = xi[neighbor_ids[valid_mask]]

            # Force contribution from this neighbor with proper weight and scaling
            weight = weights[neighbor_idx]
            force_contribution = weight * force_scale * (xi_j - xi_i)

            # Accumulate forces
            forces[valid_mask, 3] += force_contribution

        return forces

    def compute_potential_energy(self, state: BraneState, grid: BraneGrid) -> torch.Tensor:
        """
        Compute potential energy.

        For linear tension model:
            PE = (T/4h) · Σ_links (ξⱼ - ξᵢ)²

        Factor of 1/4 accounts for double-counting of links in neighbor iteration.

        Args:
            state: BraneState
            grid: BraneGrid

        Returns:
            Total potential energy [J]
        """
        xi = state.positions[:, 3]
        total_energy = torch.tensor(0.0, device=state.device, dtype=state.dtype)

        # Sum over all neighbor connections
        for neighbor_idx in range(grid.max_neighbors):
            neighbor_ids = grid.neighbors[:, neighbor_idx]
            valid = neighbor_ids >= 0

            if not valid.any():
                continue

            xi_i = xi[valid]
            xi_j = xi[neighbor_ids[valid]]

            # Strain energy: (T/2h) · (Δξ)²
            # Factor of 1/2 from each spring
            delta_xi = xi_j - xi_i
            total_energy += 0.5 * self.tension / self.grid_spacing * torch.sum(delta_xi ** 2)

        # Divide by 2 to account for double-counting (each link counted from both ends)
        return total_energy / 2.0

    def __repr__(self) -> str:
        """String representation."""
        return (
            f"LinearTensionForceComputer("
            f"T={self.tension:.2e} N, "
            f"h={self.grid_spacing:.2e} m)"
        )
