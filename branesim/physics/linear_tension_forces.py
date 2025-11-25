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
from branesim.core.state import BraneState
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
        Compute linear tension forces: F_i = (T/h) · Σⱼ (ξⱼ - ξᵢ).

        This approximates the continuous wave equation:
            ∂²ξ/∂t² = (T/ρ) · ∇²ξ

        Works for 1D, 2D, and 3D by summing over all neighbors.

        Args:
            state: BraneState with current positions
            grid: BraneGrid with neighbor connectivity

        Returns:
            Forces tensor [N, 4] in 4D embedding space
        """
        forces = torch.zeros_like(state.positions)

        # Only compute for transverse direction (ξ³)
        xi = state.positions[:, 3]

        # Vectorized computation over all neighbor slots
        for neighbor_idx in range(grid.max_neighbors):
            # Get neighbor indices [N]
            neighbor_ids = grid.neighbors[:, neighbor_idx]

            # Mask valid neighbors (not -1)
            valid_mask = neighbor_ids >= 0

            if not valid_mask.any():
                continue

            # For valid neighbors: F += (T/h) · (ξⱼ - ξᵢ)
            xi_i = xi[valid_mask]
            xi_j = xi[neighbor_ids[valid_mask]]

            # Force contribution from this neighbor
            force_contribution = (self.tension / self.grid_spacing) * (xi_j - xi_i)

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
