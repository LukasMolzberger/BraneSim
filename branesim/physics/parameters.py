"""
Physical parameter calibration for 3D brane simulations.

This module implements the Compton-cell-based calibration procedure described
in the paper (Section "Amplitude scale calibration"), connecting continuum
brane parameters to discrete lattice parameters.
"""

import math
from typing import Dict, Optional


def brane_lattice_params_3d(
    grid_spacing_m: float,
    use_compton_default: bool = True,
    rho_mass_density: Optional[float] = None,
    c: float = 2.99792458e8
) -> Dict[str, float]:
    """
    Compute brane point mass and spring constant for a 3D cubic lattice,
    using the Compton-cell-based physical estimate if desired.

    This implements Route (i) from the paper: Compton-cell mass calibration.
    The key assumptions are:
        - Reduced Compton wavelength: λ_C = ℏ/(m_e c)
        - Brane volume mass density: ρ = m_e / λ_C³  (Compton-cell assumption)
        - Brane bulk modulus: K = ρ c²

    For a discrete cubic lattice with spacing h:
        - Point mass: m_point = ρ h³
        - Axial spring constant: k_spring = K h

    This mapping ensures the lattice reproduces the continuum wave speed
    c² = K/ρ in the long-wavelength limit.

    Parameters
    ----------
    grid_spacing_m : float
        Lattice spacing h in meters.
    use_compton_default : bool, optional
        If True, use the Compton-cell assumption: ρ = m_e / λ_C³, K = ρ c².
        If False, you must supply rho_mass_density explicitly.
    rho_mass_density : float, optional
        Brane volume mass density ρ in kg/m³.
        Only used if use_compton_default is False.
    c : float, optional
        Wave speed in the brane (m/s), default is speed of light.

    Returns
    -------
    params : dict
        {
            "lambda_C": reduced Compton wavelength (m),
            "rho": mass density ρ (kg/m³),
            "K": bulk modulus K (Pa = N/m²),
            "m_point": mass per lattice point (kg),
            "k_spring": spring constant per axial spring (N/m)
        }

    Examples
    --------
    >>> # Use Compton-cell calibration with grid spacing 1e-13 m
    >>> params = brane_lattice_params_3d(1e-13)
    >>> print(f"Point mass: {params['m_point']:.2e} kg")
    >>> print(f"Spring constant: {params['k_spring']:.2e} N/m")

    >>> # Use custom mass density
    >>> params = brane_lattice_params_3d(
    ...     1e-13,
    ...     use_compton_default=False,
    ...     rho_mass_density=1e10
    ... )

    References
    ----------
    See paper Section "Amplitude scale calibration" (experimental-setting.tex)
    for the theoretical derivation of these relations.
    """
    # Physical constants (CODATA-like)
    hbar = 1.054571817e-34     # J·s
    m_e = 9.1093837015e-31     # kg

    # Reduced Compton wavelength of the electron
    lambda_C = hbar / (m_e * c)

    if use_compton_default:
        # Compton-cell assumption: ρ · λ_C³ ≈ m_e
        rho = m_e / (lambda_C ** 3)
    else:
        if rho_mass_density is None:
            raise ValueError(
                "rho_mass_density must be provided when use_compton_default=False."
            )
        rho = rho_mass_density

    # Effective bulk modulus / stiffness of the brane
    K = rho * c ** 2

    # Discrete mapping for a 3D cubic lattice:
    # - each grid cell ~ volume h³
    # - point mass m = ρ h³
    # - axial spring constant k = K h
    h = grid_spacing_m
    m_point = rho * h ** 3
    k_spring = K * h

    return {
        "lambda_C": lambda_C,
        "rho": rho,
        "K": K,
        "m_point": m_point,
        "k_spring": k_spring,
    }


def print_calibration_summary(params: Dict[str, float], grid_spacing: float) -> None:
    """
    Print a human-readable summary of brane calibration parameters.

    Parameters
    ----------
    params : dict
        Output from brane_lattice_params_3d()
    grid_spacing : float
        Lattice spacing h in meters
    """
    print("\n" + "=" * 60)
    print("3D Brane Lattice Calibration (Compton-cell based)")
    print("=" * 60)
    print(f"Reduced Compton wavelength λ_C = {params['lambda_C']:.4e} m")
    print(f"Grid spacing h                 = {grid_spacing:.4e} m")
    print(f"Grid spacing / λ_C             = {grid_spacing / params['lambda_C']:.2f}")
    print()
    print("Continuum parameters:")
    print(f"  Mass density ρ                = {params['rho']:.4e} kg/m³")
    print(f"  Bulk modulus K                = {params['K']:.4e} Pa")
    print(f"  Wave speed c = √(K/ρ)         = {math.sqrt(params['K'] / params['rho']):.4e} m/s")
    print()
    print("Discrete lattice parameters:")
    print(f"  Point mass m = ρ h³           = {params['m_point']:.4e} kg")
    print(f"  Spring constant k = K h       = {params['k_spring']:.4e} N/m")
    print("=" * 60 + "\n")


# Example usage (can be removed or kept for testing)
if __name__ == "__main__":
    # Example: lattice spacing at 10× Compton wavelength
    from branesim.config.simulation_config import PhysicalConstants

    constants = PhysicalConstants()
    h = 10.0 * constants.lambda_C

    params = brane_lattice_params_3d(h)
    print_calibration_summary(params, h)

    # Show how to use in simulation setup
    print("Usage in simulation:")
    print(f"  mass_per_point   = {params['m_point']:.4e}  # kg")
    print(f"  spring_constant  = {params['k_spring']:.4e}  # N/m")