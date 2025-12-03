"""
Dimensional mapping between physical and simulation units.

This module provides a clean interface for converting between physical (SI) units
and dimensionless simulation units. It encapsulates all scaling factors and
conversion logic in a single class.
"""

import torch
import numpy as np
from typing import Union

# Type alias for values that can be scalars or arrays
ScalarOrArray = Union[float, int, torch.Tensor, np.ndarray]


class DimensionalMapper:
    """
    Bidirectional converter between physical and simulation units.

    This class encapsulates the scaling factors and provides methods for
    converting all physical quantities between SI units and dimensionless
    simulation units.

    The simulation uses clean O(1) dimensionless units where:
        h_sim = 1.0 (grid spacing)
        m_sim = 1.0 (point mass)
        k_sim = 1.0 (spring constant)
        c_sim = 1.0 (wave speed)

    Physical quantities are scaled by:
        L0: length scale [m]
        T0: time scale [s]
        M0: mass scale [kg]
        E0: energy scale [J]

    Attributes:
        L0: Length scale [m]
        T0: Time scale [s]
        M0: Mass scale [kg]
        E0: Energy scale [J]
        c_phys: Physical wave speed [m/s]
        h_phys: Physical grid spacing [m]
    """

    def __init__(self, phys_params: dict, cfl_factor: float = 0.1):
        """
        Initialize dimensional mapper from physical parameters.

        Parameters
        ----------
        phys_params : dict
            Physical parameter dictionary containing:
            - "h_phys": grid spacing [m]
            - "m_point": mass per lattice point [kg]
            - "c_phys": wave speed [m/s]
        cfl_factor : float, optional
            CFL factor for time step calculation, default is 0.1.
        """
        self.h_phys = phys_params["h_phys"]
        self.m_point_phys = phys_params["m_point"]
        self.c_phys = phys_params["c_phys"]

        # Compute scaling factors
        self.L0 = self.h_phys
        self.T0 = self.L0 / self.c_phys
        self.M0 = self.m_point_phys
        self.E0 = self.M0 * (self.L0 / self.T0) ** 2

        # Store CFL factor
        self.cfl_factor = cfl_factor

    # ========================================================================
    # LENGTH conversions
    # ========================================================================

    def to_sim_length(self, length_phys: ScalarOrArray) -> ScalarOrArray:
        """Convert physical length [m] to simulation units."""
        return length_phys / self.L0

    def to_phys_length(self, length_sim: ScalarOrArray) -> ScalarOrArray:
        """Convert simulation length to physical units [m]."""
        return length_sim * self.L0

    # ========================================================================
    # TIME conversions
    # ========================================================================

    def to_sim_time(self, time_phys: ScalarOrArray) -> ScalarOrArray:
        """Convert physical time [s] to simulation units."""
        return time_phys / self.T0

    def to_phys_time(self, time_sim: ScalarOrArray) -> ScalarOrArray:
        """Convert simulation time to physical units [s]."""
        return time_sim * self.T0

    # ========================================================================
    # MASS conversions
    # ========================================================================

    def to_sim_mass(self, mass_phys: ScalarOrArray) -> ScalarOrArray:
        """Convert physical mass [kg] to simulation units."""
        return mass_phys / self.M0

    def to_phys_mass(self, mass_sim: ScalarOrArray) -> ScalarOrArray:
        """Convert simulation mass to physical units [kg]."""
        return mass_sim * self.M0

    # ========================================================================
    # VELOCITY conversions (length/time)
    # ========================================================================

    def to_sim_velocity(self, velocity_phys: ScalarOrArray) -> ScalarOrArray:
        """Convert physical velocity [m/s] to simulation units."""
        return velocity_phys / (self.L0 / self.T0)

    def to_phys_velocity(self, velocity_sim: ScalarOrArray) -> ScalarOrArray:
        """Convert simulation velocity to physical units [m/s]."""
        return velocity_sim * (self.L0 / self.T0)

    # ========================================================================
    # ACCELERATION conversions (length/time²)
    # ========================================================================

    def to_sim_acceleration(self, accel_phys: ScalarOrArray) -> ScalarOrArray:
        """Convert physical acceleration [m/s²] to simulation units."""
        return accel_phys / (self.L0 / self.T0**2)

    def to_phys_acceleration(self, accel_sim: ScalarOrArray) -> ScalarOrArray:
        """Convert simulation acceleration to physical units [m/s²]."""
        return accel_sim * (self.L0 / self.T0**2)

    # ========================================================================
    # FORCE conversions (mass × length/time²)
    # ========================================================================

    def to_sim_force(self, force_phys: ScalarOrArray) -> ScalarOrArray:
        """Convert physical force [N] to simulation units."""
        return force_phys / (self.M0 * self.L0 / self.T0**2)

    def to_phys_force(self, force_sim: ScalarOrArray) -> ScalarOrArray:
        """Convert simulation force to physical units [N]."""
        return force_sim * (self.M0 * self.L0 / self.T0**2)

    # ========================================================================
    # ENERGY conversions (mass × length²/time²)
    # ========================================================================

    def to_sim_energy(self, energy_phys: ScalarOrArray) -> ScalarOrArray:
        """Convert physical energy [J] to simulation units."""
        return energy_phys / self.E0

    def to_phys_energy(self, energy_sim: ScalarOrArray) -> ScalarOrArray:
        """Convert simulation energy to physical units [J]."""
        return energy_sim * self.E0

    # ========================================================================
    # SPRING CONSTANT conversions (force/length = N/m)
    # ========================================================================

    def to_sim_spring_constant(self, k_phys: ScalarOrArray) -> ScalarOrArray:
        """Convert physical spring constant [N/m] to simulation units."""
        return k_phys / (self.M0 / self.T0**2)

    def to_phys_spring_constant(self, k_sim: ScalarOrArray) -> ScalarOrArray:
        """Convert simulation spring constant to physical units [N/m]."""
        return k_sim * (self.M0 / self.T0**2)

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
        return positions_phys / self.L0

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
        return positions_sim * self.L0

    def to_sim_velocities(self, velocities_phys: torch.Tensor) -> torch.Tensor:
        """Convert physical velocities [N, 4] to simulation units."""
        return velocities_phys / (self.L0 / self.T0)

    def to_phys_velocities(self, velocities_sim: torch.Tensor) -> torch.Tensor:
        """Convert simulation velocities [N, 4] to physical units [m/s]."""
        return velocities_sim * (self.L0 / self.T0)

    def to_sim_accelerations(self, accelerations_phys: torch.Tensor) -> torch.Tensor:
        """Convert physical accelerations [N, 4] to simulation units."""
        return accelerations_phys / (self.L0 / self.T0**2)

    def to_phys_accelerations(self, accelerations_sim: torch.Tensor) -> torch.Tensor:
        """Convert simulation accelerations [N, 4] to physical units [m/s²]."""
        return accelerations_sim * (self.L0 / self.T0**2)

    def to_sim_forces(self, forces_phys: torch.Tensor) -> torch.Tensor:
        """Convert physical forces [N, 4] to simulation units."""
        return forces_phys / (self.M0 * self.L0 / self.T0**2)

    def to_phys_forces(self, forces_sim: torch.Tensor) -> torch.Tensor:
        """Convert simulation forces [N, 4] to physical units [N]."""
        return forces_sim * (self.M0 * self.L0 / self.T0**2)

    # ========================================================================
    # Utility methods
    # ========================================================================

    def get_sim_time_step(self) -> float:
        """Get dimensionless time step based on CFL condition."""
        dt_phys = self.cfl_factor * self.h_phys / self.c_phys
        return dt_phys / self.T0

    def get_phys_time_step(self) -> float:
        """Get physical time step [s] based on CFL condition."""
        return self.cfl_factor * self.h_phys / self.c_phys

    def __repr__(self) -> str:
        """String representation showing scaling factors."""
        return (
            f"DimensionalMapper(\n"
            f"  L0 = {self.L0:.6e} m,\n"
            f"  T0 = {self.T0:.6e} s,\n"
            f"  M0 = {self.M0:.6e} kg,\n"
            f"  E0 = {self.E0:.6e} J,\n"
            f"  c_phys = {self.c_phys:.6e} m/s,\n"
            f"  h_phys = {self.h_phys:.6e} m\n"
            f")"
        )


def create_mapper_from_params(params: dict) -> DimensionalMapper:
    """
    Create a DimensionalMapper from parameter dictionary.

    This is a convenience function that extracts the necessary fields
    from a parameter dict (typically from get_dimensionless_params)
    and creates a mapper.

    Parameters
    ----------
    params : dict
        Parameter dictionary containing h_phys, m_point/m_point_phys,
        c_phys, and optionally cfl_factor.

    Returns
    -------
    DimensionalMapper
        Initialized dimensional mapper.
    """
    # Build minimal phys_params dict for mapper
    phys_params = {
        "h_phys": params["h_phys"],
        "m_point": params.get("m_point_phys", params.get("m_point")),
        "c_phys": params["c_phys"],
    }

    cfl_factor = params.get("cfl_factor", 0.1)

    return DimensionalMapper(phys_params, cfl_factor)


# Example usage
if __name__ == "__main__":
    from branesim.config.simulation_config import PhysicalConstants
    from branesim.physics.parameters import compton_calibrated_brane_lattice_params

    # Setup
    constants = PhysicalConstants()
    h_phys = 10.0 * constants.lambda_C

    # Get physical parameters
    phys_params = compton_calibrated_brane_lattice_params(
        grid_spacing_m=h_phys,
        dimensionality=1,
        c=constants.c
    )
    phys_params["h_phys"] = h_phys
    phys_params["c_phys"] = constants.c

    # Create mapper
    mapper = DimensionalMapper(phys_params, cfl_factor=0.1)
    print(mapper)

    # Test conversions
    print("\n" + "=" * 60)
    print("Testing bidirectional conversions:")
    print("=" * 60)

    wavelength_phys = 40 * h_phys
    wavelength_sim = mapper.to_sim_length(wavelength_phys)
    wavelength_back = mapper.to_phys_length(wavelength_sim)

    print(f"Length: {wavelength_phys:.6e} m → {wavelength_sim:.1f} sim → {wavelength_back:.6e} m")
    print(f"Match: {abs(wavelength_phys - wavelength_back) < 1e-20}")

    time_phys = 1e-15  # 1 fs
    time_sim = mapper.to_sim_time(time_phys)
    time_back = mapper.to_phys_time(time_sim)

    print(f"Time: {time_phys:.6e} s → {time_sim:.6f} sim → {time_back:.6e} s")
    print(f"Match: {abs(time_phys - time_back) < 1e-25}")

    # Test convenience methods
    print(f"\nConvenience: {wavelength_sim:.1f} sim = {mapper.to_nanometers(wavelength_sim):.3f} nm")
    print(f"Convenience: {time_sim:.6f} sim = {mapper.to_femtoseconds(time_sim):.3f} fs")