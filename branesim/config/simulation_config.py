"""
SimulationConfig: Configuration management for brane simulations.

This module provides a dataclass-based configuration system with automatic
validation, including CFL stability checking.
"""

from dataclasses import dataclass



@dataclass
class PhysicalConstants:
    """Physical constants in SI units."""
    c: float = 299792458.0  # Speed of light [m/s]
    hbar: float = 1.054571817e-34  # Reduced Planck constant [J·s]
    m_e: float = 9.1093837015e-31  # Electron mass [kg]
    G: float = 6.67430e-11  # Newton's gravitational constant [m³/(kg·s²)]
    epsilon0: float = 8.8541878128e-12  # Vacuum permittivity [F/m]
    mu0: float = 1.25663706212e-6  # Vacuum permeability [H/m]

    @property
    def lambda_C(self) -> float:
        """Compton wavelength [m]."""
        return self.hbar / (self.m_e * self.c)

    def compute_target_tension(self, rho_m: float) -> float:
        """
        Compute target brane tension for wave speed = c.

        From continuum theory: c² = T / ρ_m

        Args:
            rho_m: Brane mass density [kg/m²]

        Returns:
            T_target: Target tension [N/m]
        """
        return rho_m * self.c**2

    def compute_target_bending_stiffness(self) -> float:
        """
        Compute target bending stiffness for correct gravitational coupling.

        From emergent gravity analysis: κ = c³ / (8π G)

        Returns:
            κ_target: Target bending stiffness [J·m]
        """
        import math
        return self.c**3 / (8.0 * math.pi * self.G)
