"""
VelocityVerletSolver: Velocity Verlet time integration for brane simulations.

This module implements the 2nd-order symplectic Velocity Verlet integrator,
which provides excellent long-term energy conservation for Hamiltonian systems.
"""

import torch
from typing import Dict

from branesim.core.state import BraneState
from branesim.core.grid import BraneGrid


class VelocityVerletSolver:
    """
    Velocity Verlet integrator for brane wave equation.

    Implements the three-phase update:
        Phase 1: R_new = R + v * dt + 0.5 * a * dt²
        Phase 2: a_new = F(R_new) / m
        Phase 3: v_new = v + 0.5 * (a + a_new) * dt

    This scheme is:
    - 2nd order accurate in both position and velocity
    - Symplectic (preserves phase space volume)
    - Time-reversible
    - Requires only 1 force evaluation per step

    Attributes:
        dt: float, time step in seconds
        mass_density: float, mass per unit volume (ρ_m) [kg/m³]
        physics: Force computer instance
        grid: BraneGrid instance
        time: float, current simulation time
        step_count: int, number of steps taken
    """

    def __init__(
        self,
        dt: float,
        mass_density: float,
        physics,  # SpringForceComputer
        grid: BraneGrid
    ):
        """
        Initialize Velocity Verlet solver.

        Args:
            dt: Time step [s]
            mass_density: Mass per unit volume ρ_m [kg/m³] (or per unit length in 1D, per unit area in 2D)
            physics: Force computer (e.g., SpringForceComputer)
            grid: BraneGrid with topology
        """
        self.dt = dt
        self.mass_density = mass_density
        self.physics = physics
        self.grid = grid
        self.time = 0.0
        self.step_count = 0

        # Compute actual mass per point based on dimensionality
        # 1D: m = ρ * h
        # 2D: m = ρ * h²
        # 3D: m = ρ * h³
        if grid.dimension.value == 1:
            self.mass_per_point = mass_density * grid.spacing
        elif grid.dimension.value == 2:
            self.mass_per_point = mass_density * grid.spacing ** 2
        else:  # 3D
            self.mass_per_point = mass_density * grid.spacing ** 3

    def step(self, state: BraneState) -> BraneState:
        """
        Perform single Velocity Verlet time step.

        Updates state in-place and returns it for method chaining.

        Args:
            state: BraneState with current positions, velocities, accelerations

        Returns:
            Updated state (same object)
        """
        # Phase 1: Update positions using current velocity and acceleration
        # R_new = R + v * dt + 0.5 * a * dt²
        state.positions += (
            state.velocities * self.dt +
            0.5 * state.accelerations * self.dt ** 2
        )

        # Apply fixed boundary conditions (restore fixed positions)
        state.apply_fixed_boundaries()

        # Phase 2: Compute new accelerations at new positions
        # F_new = compute_forces(R_new)
        # a_new = F_new / m_point
        forces = self.physics.compute_forces(state, self.grid)
        state.new_accelerations = forces / self.mass_per_point

        # Apply fixed boundary conditions (zero accelerations for fixed points)
        state.apply_fixed_boundaries()

        # Phase 3: Update velocities using average of old and new accelerations
        # v_new = v + 0.5 * (a + a_new) * dt
        state.velocities += 0.5 * (state.accelerations + state.new_accelerations) * self.dt

        # Apply fixed boundary conditions (zero velocities for fixed points)
        state.apply_fixed_boundaries()

        # Swap accelerations: a = a_new for next step
        state.accelerations = state.new_accelerations.clone()

        # Update time and step count
        self.time += self.dt
        self.step_count += 1

        return state

    def initialize_accelerations(self, state: BraneState):
        """
        Compute initial accelerations before first step.

        This ensures that acceleration is properly initialized for the
        first Verlet step.

        Args:
            state: BraneState to initialize
        """
        forces = self.physics.compute_forces(state, self.grid)
        state.accelerations = forces / self.mass_per_point
        state.new_accelerations = state.accelerations.clone()

    def compute_energy(self, state: BraneState) -> Dict[str, float]:
        """
        Compute kinetic, potential, and total energy.

        Energy:
            KE = 0.5 * Σ m_point * |v_p|²
            PE = Σ φ(ε_link)
            E_total = KE + PE

        Args:
            state: BraneState with current positions and velocities

        Returns:
            Dictionary with keys: 'kinetic', 'potential', 'total'
        """
        # Kinetic energy: KE = 0.5 * m_point * Σ |v_p|²
        kinetic = 0.5 * self.mass_per_point * torch.sum(state.velocities ** 2)

        # Potential energy: PE = Σ φ(ε)
        potential = self.physics.compute_potential_energy(state, self.grid)

        return {
            'kinetic': kinetic.item(),
            'potential': potential.item(),
            'total': (kinetic + potential).item()
        }

    def verify_wave_speed(self) -> tuple[float, float, float]:
        """
        Verify that the configured parameters reproduce the physical speed of light.

        CRITICAL: c = 3×10⁸ m/s is a PHYSICAL CONSTANT, not something we compute.
        This method checks if the discrete model parameters reproduce this value.

        Returns:
            Tuple of (expected_c, computed_c, relative_error)
            where:
                expected_c = 3×10⁸ m/s (physical constant)
                computed_c = wave speed from current parameters
                relative_error = |computed_c - expected_c| / expected_c
        """
        from branesim.config.simulation_config import PhysicalConstants

        constants = PhysicalConstants()
        expected_c = constants.c

        # Compute effective wave speed from SpringForceComputer parameters
        k = self.physics.spring_constant
        L_0 = self.physics.rest_length
        h = self.grid.spacing

        if self.grid.dimension.value == 1:
            T_0 = k * (h - L_0)
            mu = self.mass_density
            computed_c = (T_0 / mu) ** 0.5
        else:
            tension = k * L_0
            computed_c = (tension / self.mass_density) ** 0.5

        relative_error = abs(computed_c - expected_c) / expected_c

        return expected_c, computed_c, relative_error

    def _get_effective_wave_speed(self) -> float:
        """
        Internal helper to get effective wave speed for CFL checking.

        Returns:
            Computed wave speed in m/s
        """
        _, computed_c, _ = self.verify_wave_speed()
        return computed_c

    def reset_time(self):
        """Reset simulation time and step count to zero."""
        self.time = 0.0
        self.step_count = 0

    def __repr__(self) -> str:
        """String representation."""
        expected_c, computed_c, error = self.verify_wave_speed()

        # Show warning if wave speed deviates significantly from physical constant
        if error > 0.01:  # > 1% error
            warning = f" [WARNING: {error*100:.1f}% from c_physical]"
        else:
            warning = ""

        return (
            f"VelocityVerletSolver(dt={self.dt:.2e}s, "
            f"ρ_m={self.mass_density:.2e} kg/m³, "
            f"c_eff={computed_c:.2e} m/s{warning}, "
            f"t={self.time:.2e}s, steps={self.step_count})"
        )
