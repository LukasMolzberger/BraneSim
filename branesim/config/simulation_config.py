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

    @property
    def lambda_C(self) -> float:
        """Compton wavelength [m]."""
        return self.hbar / (self.m_e * self.c)
