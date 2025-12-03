"""
Physical parameter calibration for brane simulations using Compton-cell calibration.

This module implements the Compton-cell-based calibration procedure described
in the paper (Section "Amplitude scale calibration"), connecting continuum
brane parameters to discrete lattice parameters.

Note: For unit conversion between physical and simulation units, see
branesim.physics.dimensional_mapping.DimensionalMapper and map_to_dimensionless_params
"""

import math
from branesim.physics.dimensional_mapping import map_to_dimensionless_params


def compton_calibrated_brane_lattice_params(
    grid_spacing_m: float,
    dimensionality: int,
    c: float = 2.99792458e8,
) -> dict:
    """
    Compton-calibrated brane lattice parameters for a D-dimensional lattice.

    This uses the simple Compton-cell assumption:

        (1) Reduced Compton wavelength:
            lambda_C = hbar / (m_e * c)

        (2) Effective brane mass density in D dimensions:
            rho_D = m_e / lambda_C^D

        (3) Effective D-dimensional "tension" / stiffness:
            T_D = rho_D * c^2

    For a hypercubic lattice with spacing h:

        mass per lattice point:
            m_point = rho_D * h^D

        spring constant of an axial spring:
            k_spring = T_D * h^(D - 2)

    which reproduces c^2 = T_D / rho_D in the long-wavelength limit for
    D = 1, 2, 3:

        D = 1:  rho_1 = m_e / lambda_C      (kg/m)
                T_1   = rho_1 * c^2         (N)
                m     = rho_1 * h           (kg)
                k     = T_1 / h             (N/m)

        D = 2:  rho_2 = m_e / lambda_C^2    (kg/m^2)
                T_2   = rho_2 * c^2         (N/m)
                m     = rho_2 * h^2         (kg)
                k     = T_2                 (N/m)

        D = 3:  rho_3 = m_e / lambda_C^3    (kg/m^3)
                T_3   = rho_3 * c^2         (N/m^2)
                m     = rho_3 * h^3         (kg)
                k     = T_3 * h             (N/m)

    Parameters
    ----------
    grid_spacing_m : float
        Lattice spacing h in meters.
    dimensionality : int
        Spatial dimensionality D of the brane model (1, 2, or 3).
    c : float, optional
        Wave speed in the brane (m/s), default is the speed of light.

    Returns
    -------
    params : dict
        {
            "dim":        D,
            "lambda_C":   reduced Compton wavelength (m),
            "rho_D":      D-dim mass density (kg/m^D),
            "T_D":        D-dim tension / stiffness (SI: N/m^(D-1)),
            "m_point":    mass per lattice point (kg),
            "k_spring":   axial spring constant (N/m),
            "A_estimate": amplitude scale estimate ~ lambda_C / sqrt(pi) (m),
            "rest_length": spring rest length for pretension implementation (m)
        }
    """
    if dimensionality not in (1, 2, 3):
        raise ValueError("dimensionality must be 1, 2, or 3.")

    # Physical constants (CODATA-like)
    hbar = 1.054571817e-34      # J*s
    m_e  = 9.1093837015e-31     # kg

    # Reduced Compton wavelength of the electron
    lambda_C = hbar / (m_e * c)

    D = dimensionality
    h = grid_spacing_m

    # D-dimensional Compton-based mass density:
    #   rho_D * lambda_C^D = m_e  =>  rho_D = m_e / lambda_C^D
    rho_D = m_e / (lambda_C ** D)

    # Effective D-dimensional "tension" / stiffness:
    #   T_D = rho_D * c^2
    T_D = rho_D * c**2

    # Discrete mass per lattice point:
    m_point = rho_D * (h ** D)

    # Axial spring constant following k = T_D * h^(D - 2)
    k_spring = T_D * (h ** (D - 2))

    # Compton-scale amplitude estimate for one-quantum excitation
    # A ~ lambda_C / sqrt(pi)  (from earlier derivation)
    A_estimate = lambda_C / math.sqrt(math.pi)

    # Rest length for pretension implementation
    # CRITICAL: To implement the continuum pretension κ = ρc² in the discrete model,
    # springs must be pre-stretched. For Compton calibration, this gives L_0 = 0.
    #
    # Physical reasoning:
    #   - Continuum has pretension/tension T_D = ρ_D * c²
    #   - Each discrete spring must carry background force to realize this pretension
    #   - For spacing h and spring constant k_spring:
    #     * D=1: F_0 = T_1,      k = T_1/h    → L_0 = h - T_1/k = h - h = 0
    #     * D=2: F_0 = T_2*h,    k = T_2      → L_0 = h - (T_2*h)/T_2 = 0
    #     * D=3: F_0 = T_3*h²,   k = T_3*h    → L_0 = h - (T_3*h²)/(T_3*h) = 0
    #   - Without this pretension, the discrete model has zero background tension,
    #     violating the κ = ρc² assumption used in the continuum derivation.
    rest_length = 0.0

    return {
        "dim": D,
        "lambda_C": lambda_C,
        "rho_D": rho_D,
        "T_D": T_D,
        "m_point": m_point,
        "k_spring": k_spring,
        "A_estimate": A_estimate,
        "rest_length": rest_length,
    }


def print_calibration_summary(params: dict, grid_spacing: float) -> None:
    """
    Print a human-readable summary of brane calibration parameters.

    Parameters
    ----------
    params : dict
        Output from compton_calibrated_brane_lattice_params()
    grid_spacing : float
        Lattice spacing h in meters
    """
    D = params['dim']
    dim_label = {1: "1D", 2: "2D", 3: "3D"}[D]

    print("\n" + "=" * 60)
    print(f"{dim_label} Brane Lattice Calibration (Compton-cell based)")
    print("=" * 60)
    print(f"Reduced Compton wavelength λ_C = {params['lambda_C']:.4e} m")
    print(f"Grid spacing h                 = {grid_spacing:.4e} m")
    print(f"Grid spacing / λ_C             = {grid_spacing / params['lambda_C']:.2f}")
    print()
    print("Continuum parameters:")
    print(f"  Mass density ρ_{D}             = {params['rho_D']:.4e} kg/m^{D}")
    print(f"  Tension/stiffness T_{D}        = {params['T_D']:.4e} N/m^{D-1}")
    print(f"  Wave speed c = √(T_{D}/ρ_{D})  = {math.sqrt(params['T_D'] / params['rho_D']):.4e} m/s")
    print()
    print("Discrete lattice parameters:")
    print(f"  Point mass m = ρ_{D} h^{D}      = {params['m_point']:.4e} kg")
    print(f"  Spring constant k              = {params['k_spring']:.4e} N/m")
    print(f"  Rest length L_0 (pretension)   = {params['rest_length']:.4e} m")
    print(f"  Amplitude estimate A           = {params['A_estimate']:.4e} m")
    print("=" * 60 + "\n")
