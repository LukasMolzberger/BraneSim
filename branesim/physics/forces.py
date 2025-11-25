"""
SpringForceComputer: Compute spring forces between brane points.

This module implements the force computation for spring-mass model of the brane,
with optional nonlinear saturation.
"""

import torch
from typing import Optional

from branesim.core.state import BraneState
from branesim.core.grid import BraneGrid


class SpringForceComputer:
    """
    Computes spring forces with optional nonlinear saturation.

    Implements the force law:
        F_pq = φ'(ε) * (R_q - R_p) / |R_q - R_p|

    where ε = |R_q - R_p| - L_0 is the strain, and:
        - Linear: φ'(ε) = k * ε
        - Nonlinear saturation: φ'(ε) = k * ε / (1 + (ε/ε_cr)²)

    Attributes:
        spring_constant: float (k), spring constant in N/m
        rest_length: float (L_0), rest length in meters
        critical_strain: float or None (ε_cr), for nonlinear saturation
        use_saturation: bool, whether to use nonlinear saturation
        apply_boundary_tension: bool, whether to apply boundary forces
    """

    def __init__(
        self,
        spring_constant: float,
        rest_length: float,
        critical_strain: Optional[float] = None,
        apply_boundary_tension: bool = False
    ):
        """
        Initialize spring force computer.

        Args:
            spring_constant: Spring constant k [N/m]
            rest_length: Rest length L_0 [m]
            critical_strain: Critical strain ε_cr for saturation (None for linear)
            apply_boundary_tension: Whether to apply boundary tension forces
        """
        self.spring_constant = spring_constant
        self.rest_length = rest_length
        self.critical_strain = critical_strain
        self.use_saturation = (critical_strain is not None)
        self.apply_boundary_tension = apply_boundary_tension

    def compute_forces(self, state: BraneState, grid: BraneGrid) -> torch.Tensor:
        """
        Compute spring forces for all points.

        Vectorized implementation that loops over neighbor slots rather than points.
        Uses masking to handle boundary points efficiently.

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

            # Get positions for valid pairs [N_valid, 4]
            p_pos = state.positions[valid_mask]
            q_pos = state.positions[neighbor_ids[valid_mask]]

            # Compute displacements [N_valid, 4]
            delta = q_pos - p_pos
            distance = torch.norm(delta, dim=1, keepdim=True)  # [N_valid, 1]

            # Compute strain [N_valid, 1]
            strain = distance - self.rest_length

            # Compute force magnitude with optional saturation
            if self.use_saturation:
                # φ'(ε) = k*ε/(1 + (ε/ε_cr)²)
                force_mag = self.spring_constant * strain / (
                    1.0 + (strain / self.critical_strain) ** 2
                )
            else:
                # φ'(ε) = k*ε
                force_mag = self.spring_constant * strain

            # Force direction: delta / distance (add small epsilon to avoid division by zero)
            force_dir = delta / (distance + 1e-10)

            # Total force [N_valid, 4]
            link_forces = force_mag * force_dir

            # Accumulate forces
            forces[valid_mask] += link_forces

        # Apply boundary tension if enabled
        if self.apply_boundary_tension:
            forces = self._apply_boundary_tension(forces, state, grid)

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
            strain = distance - self.rest_length

            # Compute potential energy for these links
            if self.use_saturation:
                # φ(ε) = (k * ε_cr² / 2) * ln(1 + (ε/ε_cr)²)
                link_energy = (
                    0.5 * self.spring_constant * self.critical_strain ** 2 *
                    torch.log(1.0 + (strain / self.critical_strain) ** 2)
                )
            else:
                # φ(ε) = (k / 2) * ε²
                link_energy = 0.5 * self.spring_constant * strain ** 2

            total_energy += torch.sum(link_energy)

        # Divide by 2 because each link is counted twice (from both endpoints)
        return total_energy / 2.0

    def _apply_boundary_tension(
        self,
        forces: torch.Tensor,
        state: BraneState,
        grid: BraneGrid
    ) -> torch.Tensor:
        """
        Apply tension forces to boundary points.

        Boundary points are pulled back toward their rest configuration
        to prevent the brane from collapsing or expanding at the edges.

        Args:
            forces: Current force tensor [N, 4]
            state: BraneState with current positions
            grid: BraneGrid with boundary mask

        Returns:
            Updated forces with boundary tension
        """
        # Get boundary points
        boundary_mask = grid.is_boundary

        if not boundary_mask.any():
            return forces

        # Get rest positions for boundary points (flat configuration)
        rest_positions = torch.zeros_like(state.positions)
        spatial_coords = grid.get_spatial_coordinates()

        # Set spatial components to grid positions
        if state.dimension.value == 1:
            rest_positions[:, 0] = spatial_coords[:, 0]
        elif state.dimension.value == 2:
            rest_positions[:, :2] = spatial_coords
        else:  # 3D
            rest_positions[:, :3] = spatial_coords

        # Compute restoring force for boundary points
        displacement = state.positions[boundary_mask] - rest_positions[boundary_mask]
        restoring_force = -self.spring_constant * displacement

        # Add to total forces
        forces[boundary_mask] += restoring_force

        return forces

    def __repr__(self) -> str:
        """String representation."""
        sat_str = f"ε_cr={self.critical_strain:.2e}" if self.use_saturation else "linear"
        return (
            f"SpringForceComputer(k={self.spring_constant:.2e} N/m, "
            f"L_0={self.rest_length:.2e} m, {sat_str})"
        )
