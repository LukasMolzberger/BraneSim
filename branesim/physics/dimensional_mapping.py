"""
Dimensional mapping between physical and simulation units.

This module provides a clean interface for converting between physical (SI) units
and dimensionless simulation units using three fundamental scaling factors.
"""

import torch
import numpy as np
from typing import Union

# Type alias for values that can be scalars or arrays
ScalarOrArray = Union[float, int, torch.Tensor, np.ndarray]


class DimensionalMapper:
    """
    Transparent scaling layer for converting between physical and simulation units.

    The dimensional mapping is purely mechanical - it provides conversion factors
    to scale physical quantities into numerically tractable dimensionless units.

    Scaling factors:
        length_scale: Chosen so h_sim = 1.0
        time_scale: Chosen so c_light_sim = 1.0 (time_scale = length_scale / c_light)
        mass_scale: Fixed reference mass (independent of actual system mass)

    All other quantities are derived from these three fundamental scales using
    dimensional analysis.
    """

    def __init__(self, h_phys: float, c_light: float, mass_reference: float):
        """
        Initialize dimensional mapper with three fundamental scales.

        Parameters
        ----------
        h_phys : float
            Physical grid spacing [m]. Defines length_scale.
        c_light : float
            Speed of light [m/s]. Used with length_scale to define time_scale.
        mass_reference : float
            Fixed reference mass [kg]. Defines mass_scale.
        """
        # Three fundamental scaling factors
        self.length_scale = h_phys
        self.time_scale = self.length_scale / c_light
        self.mass_scale = mass_reference

    # ========================================================================
    # LENGTH conversions
    # ========================================================================

    def to_sim_length(self, length_phys: ScalarOrArray) -> ScalarOrArray:
        """Convert physical length [m] to simulation units."""
        return length_phys / self.length_scale

    def to_phys_length(self, length_sim: ScalarOrArray) -> ScalarOrArray:
        """Convert simulation length to physical units [m]."""
        return length_sim * self.length_scale

    # ========================================================================
    # TIME conversions
    # ========================================================================

    def to_sim_time(self, time_phys: ScalarOrArray) -> ScalarOrArray:
        """Convert physical time [s] to simulation units."""
        return time_phys / self.time_scale

    def to_phys_time(self, time_sim: ScalarOrArray) -> ScalarOrArray:
        """Convert simulation time to physical units [s]."""
        return time_sim * self.time_scale

    # ========================================================================
    # MASS conversions
    # ========================================================================

    def to_sim_mass(self, mass_phys: ScalarOrArray) -> ScalarOrArray:
        """Convert physical mass [kg] to simulation units."""
        return mass_phys / self.mass_scale

    def to_phys_mass(self, mass_sim: ScalarOrArray) -> ScalarOrArray:
        """Convert simulation mass to physical units [kg]."""
        return mass_sim * self.mass_scale

    # ========================================================================
    # FREQUENCY conversions (1/time)
    # ========================================================================

    def to_sim_frequency(self, frequency_phys: ScalarOrArray) -> ScalarOrArray:
        """Convert physical frequency [rad/s or Hz] to simulation units."""
        return frequency_phys * self.time_scale

    def to_phys_frequency(self, frequency_sim: ScalarOrArray) -> ScalarOrArray:
        """Convert simulation frequency to physical units [rad/s or Hz]."""
        return frequency_sim / self.time_scale

    # ========================================================================
    # VELOCITY conversions (length/time)
    # ========================================================================

    def to_sim_velocity(self, velocity_phys: ScalarOrArray) -> ScalarOrArray:
        """Convert physical velocity [m/s] to simulation units."""
        velocity_scale = self.length_scale / self.time_scale
        return velocity_phys / velocity_scale

    def to_phys_velocity(self, velocity_sim: ScalarOrArray) -> ScalarOrArray:
        """Convert simulation velocity to physical units [m/s]."""
        velocity_scale = self.length_scale / self.time_scale
        return velocity_sim * velocity_scale

    # ========================================================================
    # ACCELERATION conversions (length/time²)
    # ========================================================================

    def to_sim_acceleration(self, accel_phys: ScalarOrArray) -> ScalarOrArray:
        """Convert physical acceleration [m/s²] to simulation units."""
        accel_scale = self.length_scale / (self.time_scale ** 2)
        return accel_phys / accel_scale

    def to_phys_acceleration(self, accel_sim: ScalarOrArray) -> ScalarOrArray:
        """Convert simulation acceleration to physical units [m/s²]."""
        accel_scale = self.length_scale / (self.time_scale ** 2)
        return accel_sim * accel_scale

    # ========================================================================
    # FORCE conversions (mass × length/time²)
    # ========================================================================

    def to_sim_force(self, force_phys: ScalarOrArray) -> ScalarOrArray:
        """Convert physical force [N] to simulation units."""
        force_scale = self.mass_scale * self.length_scale / (self.time_scale ** 2)
        return force_phys / force_scale

    def to_phys_force(self, force_sim: ScalarOrArray) -> ScalarOrArray:
        """Convert simulation force to physical units [N]."""
        force_scale = self.mass_scale * self.length_scale / (self.time_scale ** 2)
        return force_sim * force_scale

    # ========================================================================
    # ENERGY conversions (mass × length²/time²)
    # ========================================================================

    def to_sim_energy(self, energy_phys: ScalarOrArray) -> ScalarOrArray:
        """Convert physical energy [J] to simulation units."""
        energy_scale = self.mass_scale * (self.length_scale ** 2) / (self.time_scale ** 2)
        return energy_phys / energy_scale

    def to_phys_energy(self, energy_sim: ScalarOrArray) -> ScalarOrArray:
        """Convert simulation energy to physical units [J]."""
        energy_scale = self.mass_scale * (self.length_scale ** 2) / (self.time_scale ** 2)
        return energy_sim * energy_scale

    # ========================================================================
    # SPRING CONSTANT conversions (force/length = mass/time²)
    # ========================================================================

    def to_sim_spring_constant(self, k_phys: ScalarOrArray) -> ScalarOrArray:
        """Convert physical spring constant [N/m] to simulation units."""
        k_scale = self.mass_scale / (self.time_scale ** 2)
        return k_phys / k_scale

    def to_phys_spring_constant(self, k_sim: ScalarOrArray) -> ScalarOrArray:
        """Convert simulation spring constant to physical units [N/m]."""
        k_scale = self.mass_scale / (self.time_scale ** 2)
        return k_sim * k_scale

    # ========================================================================
    # BENDING STIFFNESS conversions (energy × length = mass × length³/time²)
    # ========================================================================

    def to_sim_bending_stiffness(self, kappa_phys: ScalarOrArray) -> ScalarOrArray:
        """
        Convert physical bending stiffness [J·m] to simulation units.

        Bending stiffness has dimensions [energy × length] = [kg·m³/s²].
        """
        kappa_scale = self.mass_scale * (self.length_scale ** 3) / (self.time_scale ** 2)
        return kappa_phys / kappa_scale

    def to_phys_bending_stiffness(self, kappa_sim: ScalarOrArray) -> ScalarOrArray:
        """
        Convert simulation bending stiffness to physical units [J·m].

        Bending stiffness has dimensions [energy × length] = [kg·m³/s²].
        """
        kappa_scale = self.mass_scale * (self.length_scale ** 3) / (self.time_scale ** 2)
        return kappa_sim * kappa_scale

    # ========================================================================
    # TENSION conversions (force/length = mass/time²)
    # ========================================================================

    def to_sim_tension(self, T_phys: ScalarOrArray) -> ScalarOrArray:
        """
        Convert physical tension [N/m] to simulation units.

        Tension has dimensions [force/length] = [mass/time²].
        Note: This has the same dimensions as spring constant k.
        """
        T_scale = self.mass_scale / (self.time_scale ** 2)
        return T_phys / T_scale

    def to_phys_tension(self, T_sim: ScalarOrArray) -> ScalarOrArray:
        """
        Convert simulation tension to physical units [N/m].

        Tension has dimensions [force/length] = [mass/time²].
        """
        T_scale = self.mass_scale / (self.time_scale ** 2)
        return T_sim * T_scale

    # ========================================================================
    # Convenience methods for common conversions
    # ========================================================================

    def to_nanometers(self, length_sim: ScalarOrArray) -> ScalarOrArray:
        """Convert simulation length to nanometers."""
        return self.to_phys_length(length_sim) * 1e9

    def to_femtoseconds(self, time_sim: ScalarOrArray) -> ScalarOrArray:
        """Convert simulation time to femtoseconds."""
        return self.to_phys_time(time_sim) * 1e15

    def to_picometers(self, length_sim: ScalarOrArray) -> ScalarOrArray:
        """Convert simulation length to picometers."""
        return self.to_phys_length(length_sim) * 1e12

    # ========================================================================
    # Position/state conversions (handle [N, 4] tensors)
    # ========================================================================

    def to_sim_positions(self, positions_phys: torch.Tensor) -> torch.Tensor:
        """
        Convert physical positions [N, 4] to simulation units.

        Parameters
        ----------
        positions_phys : torch.Tensor
            Physical positions with shape [N, 4] in meters.

        Returns
        -------
        torch.Tensor
            Positions in simulation units [N, 4].
        """
        return positions_phys / self.length_scale

    def to_phys_positions(self, positions_sim: torch.Tensor) -> torch.Tensor:
        """
        Convert simulation positions [N, 4] to physical units [m].

        Parameters
        ----------
        positions_sim : torch.Tensor
            Simulation positions with shape [N, 4].

        Returns
        -------
        torch.Tensor
            Positions in physical units [N, 4] in meters.
        """
        return positions_sim * self.length_scale

    def to_sim_velocities(self, velocities_phys: torch.Tensor) -> torch.Tensor:
        """Convert physical velocities [N, 4] to simulation units."""
        velocity_scale = self.length_scale / self.time_scale
        return velocities_phys / velocity_scale

    def to_phys_velocities(self, velocities_sim: torch.Tensor) -> torch.Tensor:
        """Convert simulation velocities [N, 4] to physical units [m/s]."""
        velocity_scale = self.length_scale / self.time_scale
        return velocities_sim * velocity_scale

    def to_sim_accelerations(self, accelerations_phys: torch.Tensor) -> torch.Tensor:
        """Convert physical accelerations [N, 4] to simulation units."""
        accel_scale = self.length_scale / (self.time_scale ** 2)
        return accelerations_phys / accel_scale

    def to_phys_accelerations(self, accelerations_sim: torch.Tensor) -> torch.Tensor:
        """Convert simulation accelerations [N, 4] to physical units [m/s²]."""
        accel_scale = self.length_scale / (self.time_scale ** 2)
        return accelerations_sim * accel_scale

    def to_sim_forces(self, forces_phys: torch.Tensor) -> torch.Tensor:
        """Convert physical forces [N, 4] to simulation units."""
        force_scale = self.mass_scale * self.length_scale / (self.time_scale ** 2)
        return forces_phys / force_scale

    def to_phys_forces(self, forces_sim: torch.Tensor) -> torch.Tensor:
        """Convert simulation forces [N, 4] to physical units [N]."""
        force_scale = self.mass_scale * self.length_scale / (self.time_scale ** 2)
        return forces_sim * force_scale

    def __repr__(self) -> str:
        """String representation showing scaling factors."""
        return (
            f"DimensionalMapper(\n"
            f"  length_scale = {self.length_scale:.6e} m,\n"
            f"  time_scale   = {self.time_scale:.6e} s,\n"
            f"  mass_scale   = {self.mass_scale:.6e} kg\n"
            f")"
        )