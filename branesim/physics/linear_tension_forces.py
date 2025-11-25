"""
LinearTensionForceComputer: Simplified force model for 1D wave propagation.

Instead of computing full geometric spring forces with pre-tension,
this implements a linearized tension model that directly gives the wave equation:

    m · ∂²ξ/∂t² = T · ∂²ξ/∂x²

For a discrete 1D system:
    F_i = T/h² · (ξᵢ₊₁ + ξᵢ₋₁ - 2ξᵢ)

This avoids numerical stiffness from geometric nonlinearity while maintaining
correct wave propagation physics.
"""

import torch
from branesim.core.state import BraneState
from branesim.core.grid import BraneGrid


class LinearTensionForceComputer:
    """
    Linearized tension force for 1D wave propagation.

    Implements F = T · ∇²ξ in discrete form, giving clean wave equation
    without geometric nonlinearity or stiffness issues.

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
        Compute linear tension forces: F = T/h² · (ξᵢ₊₁ + ξᵢ₋₁ - 2ξᵢ).

        This gives the discrete wave equation:
            m · ∂²ξ/∂t² = T/h² · (ξᵢ₊₁ + ξᵢ₋₁ - 2ξᵢ)

        which approximates the continuous wave equation:
            ∂²ξ/∂t² = (T/μ) · ∂²ξ/∂x²

        Args:
            state: BraneState with current positions
            grid: BraneGrid with neighbor connectivity

        Returns:
            Forces tensor [N, 4] in 4D embedding space
        """
        forces = torch.zeros_like(state.positions)

        # Only compute for transverse direction (ξ³)
        xi = state.positions[:, 3]

        # For 1D, each interior point has exactly 2 neighbors (left and right)
        # neighbors shape: [N, 2]
        left_ids = grid.neighbors[:, 0]   # Left neighbor
        right_ids = grid.neighbors[:, 1]  # Right neighbor

        left_valid = left_ids >= 0
        right_valid = right_ids >= 0

        # Interior points: have both neighbors
        interior = left_valid & right_valid

        if interior.any():
            xi_left = xi[left_ids[interior]]
            xi_center = xi[interior]
            xi_right = xi[right_ids[interior]]

            # Discrete second difference: Δ²ξ = (ξᵢ₊₁ + ξᵢ₋₁ - 2ξᵢ)
            delta_2_xi = xi_right + xi_left - 2 * xi_center

            # Force: F = (T/h) · Δ²ξ
            # Spring model: k_eff = T/h connecting adjacent masses
            # This gives wave speed c = √(T/μ) in the long-wavelength limit
            force_magnitude = (self.tension / self.grid_spacing) * delta_2_xi

            # Apply only in transverse direction (ξ³)
            forces[interior, 3] = force_magnitude

        return forces

    def compute_potential_energy(self, state: BraneState, grid: BraneGrid) -> torch.Tensor:
        """
        Compute potential energy.

        For linear tension model:
            PE = (T/2h) · Σ (ξᵢ₊₁ - ξᵢ)²

        This approximates the continuous string energy:
            PE = (T/2) ∫ (∂ξ/∂x)² dx

        Args:
            state: BraneState
            grid: BraneGrid

        Returns:
            Total potential energy [J]
        """
        xi = state.positions[:, 3]

        # Sum over all links
        right_ids = grid.neighbors[:, 1]
        valid = right_ids >= 0

        if not valid.any():
            return torch.tensor(0.0, device=state.device, dtype=state.dtype)

        xi_center = xi[valid]
        xi_right = xi[right_ids[valid]]

        # Strain energy: (T/2h) · (Δξ)²
        delta_xi = xi_right - xi_center
        potential = 0.5 * self.tension / self.grid_spacing * torch.sum(delta_xi ** 2)

        return potential

    def __repr__(self) -> str:
        """String representation."""
        return (
            f"LinearTensionForceComputer("
            f"T={self.tension:.2e} N, "
            f"h={self.grid_spacing:.2e} m)"
        )
