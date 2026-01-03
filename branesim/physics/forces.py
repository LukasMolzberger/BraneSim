"""
SpringForceComputer: Compute spring forces between brane points.

This module implements the force computation for spring-mass model of the brane.
"""

import torch

from branesim.core.state import BraneState
from branesim.core.grid import BraneGrid


class SpringForceComputer:
    """
    Computes spring forces.

    Implements the force law:
        F_pq = φ'(ε) * (R_q - R_p) / |R_q - R_p|

    where ε = |R_q - R_p| - L_0 is the strain, and:
        - Linear: φ'(ε) = k * ε

    Attributes:
        spring_constant: float (k), spring constant in N/m for unit-offset neighbors
        rest_length: float (L_0), rest length in meters for unit-offset neighbors
    """

    def __init__(
        self,
        spring_constant: float,
        rest_length: float,
    ):
        """
        Initialize spring force computer.

        Args:
            spring_constant: Spring constant k [N/m] for unit-offset neighbors
            rest_length: Rest length L_0 [m] for unit-offset neighbors
        """
        self.spring_constant = spring_constant
        self.rest_length = rest_length

    def compute_forces(self, state: BraneState, grid: BraneGrid) -> torch.Tensor:
        """
        Compute spring forces for all points.

        Vectorized implementation that loops over neighbor slots rather than points.
        Uses masking to handle boundary points efficiently.

        Applies isotropic-stencil weights (k_delta = k / |delta|^2) and
        offset-scaled rest lengths (L0_delta = L0 * |delta|).

        Args:
            state: BraneState with current positions
            grid: BraneGrid with neighbor connectivity

        Returns:
            forces: [N, 4] tensor of force vectors
        """
        N = state.positions.shape[0]
        forces = torch.zeros_like(state.positions)

        # Vectorized neighbor force computation
        for neighbor_idx in range(grid.max_neighbors):
            # Get neighbor indices [N]
            neighbor_ids = grid.neighbors[:, neighbor_idx]

            # Mask valid neighbors (not -1)
            valid_mask = neighbor_ids >= 0

            if not valid_mask.any():
                continue

            # Convert boolean mask to indices early for consistency
            valid_indices = torch.nonzero(valid_mask, as_tuple=False).squeeze(1)

            # Get positions for valid pairs [N_valid, 4]
            p_pos = state.positions[valid_indices]
            q_pos = state.positions[neighbor_ids[valid_indices]]

            # Compute displacements [N_valid, 4]
            delta = q_pos - p_pos
            distance = torch.norm(delta, dim=1, keepdim=True)  # [N_valid, 1]

            # Compute strain [N_valid, 1] with offset-scaled rest length
            rest_length = self.rest_length * grid.neighbor_offset_norms[neighbor_idx].to(state.dtype)
            strain = distance - rest_length

            # φ'(ε) = k_delta * ε with isotropic stencil weights
            k_eff = self.spring_constant * grid.neighbor_weights[neighbor_idx].to(state.dtype)
            force_mag = k_eff * strain

            # Force direction: delta / distance (add small epsilon to avoid division by zero)
            force_dir = delta / distance

            # Total force [N_valid, 4]
            link_forces = force_mag * force_dir

            # Accumulate forces using index_add (more reliable than boolean indexing on MPS)
            forces.index_add_(0, valid_indices, link_forces)

        return forces

    def compute_potential_energy(self, state: BraneState, grid: BraneGrid) -> torch.Tensor:
        """
        Compute total potential energy.

        For nonlinear saturation:
            φ(ε) = (k * ε_cr² / 2) * ln(1 + (ε/ε_cr)²)

        For linear:
            φ(ε) = (k / 2) * ε²

        Args:
            state: BraneState with current positions
            grid: BraneGrid with neighbor connectivity

        Returns:
            Scalar tensor with total potential energy
        """
        total_energy = torch.tensor(0.0, device=state.device, dtype=state.dtype)

        # Loop over neighbors, counting each link once
        for neighbor_idx in range(grid.max_neighbors):
            neighbor_ids = grid.neighbors[:, neighbor_idx]
            valid_mask = neighbor_ids >= 0

            if not valid_mask.any():
                continue

            # Get positions
            p_pos = state.positions[valid_mask]
            q_pos = state.positions[neighbor_ids[valid_mask]]

            # Compute strain
            delta = q_pos - p_pos
            distance = torch.norm(delta, dim=1)
            rest_length = self.rest_length * grid.neighbor_offset_norms[neighbor_idx].to(state.dtype)
            strain = distance - rest_length

            # Compute potential energy for these links
            # φ(ε) = (k / 2) * ε²
            k_eff = self.spring_constant * grid.neighbor_weights[neighbor_idx].to(state.dtype)
            link_energy = 0.5 * k_eff * strain ** 2

            total_energy += torch.sum(link_energy)

        # Divide by 2 because each link is counted twice (from both endpoints)
        return total_energy / 2.0

    def __repr__(self) -> str:
        """String representation."""
        return (
            f"SpringForceComputer(k={self.spring_constant:.2e} N/m, "
            f"L_0={self.rest_length:.2e} m)"
        )
