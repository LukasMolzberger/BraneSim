"""
Physical Constants for BraneSim

This module defines all fundamental physical constants and calibrated parameters
used throughout the simulation. Constants are organized by their source and role
in the theoretical framework.

All values are in SI units unless otherwise specified.
"""

import math
from dataclasses import dataclass


@dataclass
class PhysicalConstants:
    """
    Fundamental physical constants and calibrated simulation parameters.

    This class contains:
    1. Universal physical constants (c, ℏ, G, etc.)
    2. Calibrated microscopic parameters (rest_length_frac)
    3. Derived quantities (λ_C, target tension, bending stiffness)
    """

    # ========================================================================
    # Universal Physical Constants
    # ========================================================================
    # Source: CODATA 2018 recommended values
    # https://physics.nist.gov/cuu/Constants/

    c: float = 299792458.0
    """Speed of light in vacuum [m/s] - EXACT (defining constant since 2019)"""

    hbar: float = 1.054571817e-34
    """Reduced Planck constant ℏ [J·s] - CODATA 2018"""

    m_e: float = 9.1093837015e-31
    """Electron mass [kg] - CODATA 2018"""

    G: float = 6.67430e-11
    """Newton's gravitational constant [m³/(kg·s²)] - CODATA 2018"""

    epsilon0: float = 8.8541878128e-12
    """Vacuum permittivity ε₀ [F/m] - CODATA 2018"""

    mu0: float = 1.25663706212e-6
    """Vacuum permeability μ₀ [H/m] - CODATA 2018"""

    # ========================================================================
    # Calibrated Microscopic Parameters
    # ========================================================================
    # These parameters are determined by the physical calibration procedure
    # described in branesim/tools/calibrate_physical_rest_length.py
    #
    # The calibration ensures that:
    # 1. Wave speed equals c (from T = ρ_m c²)
    # 2. Gravitational coupling matches G (from κ = c³/(8πG))
    # 3. Rest length is physically consistent with continuum theory

    rest_length_frac: float = 2.321269e-07
    """
    Dimensionless rest length ratio L₀/a.

    This is the ratio of spring rest length L₀ to lattice spacing a.

    Source: Physical calibration with parameters:
      - ρ_m = 3.975878e+14 kg/m³ (from m_point = 2.861821e-27 kg, h = λ_C/20)
      - k = 6.899394e+17 N/m (spring constant)
      - a = 1.930796e-14 m (lattice spacing = h)

    Derivation: From 1D chain model
      T = (k/a)(1 - L₀/a)
      T_target = ρ_m c²

      Therefore: L₀/a = 1 - (T_target × a)/k

    This value ensures that the microscopic lattice model reproduces the
    continuum wave speed c = √(T/ρ_m) exactly.

    NOTE: This is a PHYSICAL parameter, not a tuning knob. The dimensional
    mapping layer must INHERIT this value, not redefine it.
    """

    # ========================================================================
    # Derived Quantities
    # ========================================================================

    @property
    def lambda_C(self) -> float:
        """
        Compton wavelength of the electron λ_C = ℏ/(m_e c) [m].

        This is the characteristic length scale at which quantum effects
        become important for the electron. In the brane model, this sets
        the length scale for toroidal electron structures.

        Value: ≈ 2.426e-12 m
        """
        return self.hbar / (self.m_e * self.c)

    def compute_target_tension(self, rho_m: float) -> float:
        """
        Compute target brane tension for wave speed = c.

        From continuum theory: c² = T / ρ_m
        Therefore: T = ρ_m c²

        This relation ensures that transverse waves on the brane propagate
        at the speed of light, which is a fundamental requirement for
        emergent Lorentz invariance.

        Args:
            rho_m: Brane mass density [kg/m³] for 3D
                                      [kg/m²] for 2D
                                      [kg/m]  for 1D

        Returns:
            T_target: Target tension [J/m³] for 3D (elastic modulus)
                                    [N/m]   for 2D (tension)
                                    [N]     for 1D (force)

        Reference: PROJECT_PRINCIPLES.md, Section 2.1
        """
        return rho_m * self.c**2

    def compute_target_bending_stiffness(self) -> float:
        """
        Compute target bending stiffness for correct gravitational coupling.

        From emergent gravity analysis: κ = c³ / (8π G)

        This relation ensures that the effective gravitational constant
        emerging from brane curvature matches Newton's G.

        Returns:
            κ_target: Target bending stiffness [J·m]

        Value: ≈ 1.606e34 J·m

        Reference: PROJECT_PRINCIPLES.md, Section 2.1
        Paper: "A 3D-Brane Based Model of a Non-Classical Aether"
        """
        return self.c**3 / (8.0 * math.pi * self.G)

    def compute_rest_length(self, lattice_spacing: float) -> float:
        """
        Compute physical rest length L₀ from lattice spacing.

        Args:
            lattice_spacing: Lattice spacing a [m]

        Returns:
            L₀: Rest length [m]

        Usage:
            L₀ = rest_length_frac × a
        """
        return self.rest_length_frac * lattice_spacing

    def __repr__(self) -> str:
        """String representation showing key constants."""
        return (
            f"PhysicalConstants(\n"
            f"  c = {self.c:.6e} m/s\n"
            f"  ℏ = {self.hbar:.6e} J·s\n"
            f"  m_e = {self.m_e:.6e} kg\n"
            f"  G = {self.G:.6e} m³/(kg·s²)\n"
            f"  λ_C = {self.lambda_C:.6e} m\n"
            f"  rest_length_frac = {self.rest_length_frac:.6e}\n"
            f")"
        )


# ============================================================================
# Module-level convenience instance
# ============================================================================

# Default instance for convenient access
CONSTANTS = PhysicalConstants()

# Expose common constants at module level for convenience
c = CONSTANTS.c
hbar = CONSTANTS.hbar
m_e = CONSTANTS.m_e
G = CONSTANTS.G
epsilon0 = CONSTANTS.epsilon0
mu0 = CONSTANTS.mu0
lambda_C = CONSTANTS.lambda_C
rest_length_frac = CONSTANTS.rest_length_frac


# ============================================================================
# Notes on Usage
# ============================================================================
"""
Usage Examples:

1. Import the dataclass (recommended for most code):
   ```python
   from branesim.config.physical_constants import PhysicalConstants

   constants = PhysicalConstants()
   wave_speed = math.sqrt(T / rho_m)
   assert abs(wave_speed - constants.c) < 1e-6
   ```

2. Import module-level constants (for quick access):
   ```python
   from branesim.config.physical_constants import c, lambda_C, rest_length_frac

   h = lambda_C / 20  # Grid spacing
   L0 = rest_length_frac * h  # Rest length
   ```

3. Compute derived quantities:
   ```python
   constants = PhysicalConstants()
   T_target = constants.compute_target_tension(rho_m)
   kappa_target = constants.compute_target_bending_stiffness()
   ```

Calibration Procedure:
  To recalibrate rest_length_frac for different simulation parameters:

  ```bash
  python -m branesim.tools.calibrate_physical_rest_length \
    --rho-m <your_rho_m> \
    --spring-k <your_k> \
    --lattice-spacing <your_a> \
    --output-json config/rest_length_physical.json
  ```

  Then update the rest_length_frac value in this file with the calibrated result.

Architecture Notes:
  - This module contains ONLY physical constants and calibrated parameters
  - These values are INPUTS to the simulation, not outputs
  - The dimensional mapping layer (dimensional_mapping.py) inherits these
    values and provides unit conversion, but does NOT redefine them
  - All physics flows from these fundamental constants

References:
  - CODATA 2018: https://physics.nist.gov/cuu/Constants/
  - PROJECT_PRINCIPLES.md: Theoretical foundation
  - branesim/tools/calibrate_physical_rest_length.py: Calibration procedure
"""