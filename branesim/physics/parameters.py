"""
Physical parameter calibration for brane simulations using Compton-cell calibration.

This module implements the Compton-cell-based calibration procedure described
in the paper (Section "Amplitude scale calibration"), connecting continuum
brane parameters to discrete lattice parameters.

Note: For unit conversion between physical and simulation units, see
branesim.physics.dimensional_mapping.DimensionalMapper
"""

import math


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


def manual_brane_lattice_params(
    grid_spacing_m: float,
    dimensionality: int,
    mass_scale_multiplier: float = 1.0,
    c: float = 2.99792458e8,
) -> dict:
    """
    Manually specify brane lattice parameters with custom mass scale.

    This function allows you to specify a custom mass scale to explore
    how mass/energy density affects wave propagation and lateral coupling.

    The stiffness T_D is FIXED to the Compton-calibrated value:
        T_D = (m_e / lambda_C^D) × c²  (reference stiffness)

    The mass density ρ_D is VARIED:
        ρ_D = mass_scale_multiplier × (m_e / lambda_C^D)

    This makes the wave speed VARY:
        c_wave = √(T_D/ρ_D) = c / √(mass_scale_multiplier)

    Derived parameters:
        m_point = ρ_D × h^D  (point mass)
        k_spring = T_D × h^(D-2)  (spring constant, FIXED)

    Parameters
    ----------
    grid_spacing_m : float
        Lattice spacing h in meters.
    dimensionality : int
        Spatial dimensionality D of the brane model (1, 2, or 3).
    mass_scale_multiplier : float, optional
        Multiplier for the mass density relative to Compton calibration.
        - mass_scale_multiplier = 1.0: Compton calibration (default)
        - mass_scale_multiplier < 1.0: Lower mass density (weaker lateral coupling)
        - mass_scale_multiplier > 1.0: Higher mass density (stronger lateral coupling)
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
            "rest_length": spring rest length for pretension implementation (m),
            "mass_scale_multiplier": the multiplier used
        }

    Examples
    --------
    # Compton calibration (baseline)
    params_baseline = manual_brane_lattice_params(h, D=1, mass_scale_multiplier=1.0)

    # 100x lighter (explore weak lateral coupling)
    params_light = manual_brane_lattice_params(h, D=1, mass_scale_multiplier=0.01)

    # 100x heavier (explore strong lateral coupling)
    params_heavy = manual_brane_lattice_params(h, D=1, mass_scale_multiplier=100.0)
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

    # D-dimensional Compton-based mass density (reference)
    rho_D_compton = m_e / (lambda_C ** D)

    # FIXED stiffness (based on Compton calibration, independent of mass_scale_multiplier)
    T_D_reference = rho_D_compton * c**2
    T_D = T_D_reference  # FIXED

    # VARIED mass density (depends on mass_scale_multiplier)
    rho_D = mass_scale_multiplier * rho_D_compton

    # Actual wave speed (varies with mass_scale_multiplier)
    # c_wave = √(T_D/ρ_D) = c / √(mass_scale_multiplier)
    c_wave = math.sqrt(T_D / rho_D)

    # Discrete mass per lattice point (varies with mass_scale_multiplier)
    m_point = rho_D * (h ** D)

    # Axial spring constant (FIXED, since T_D is fixed)
    k_spring = T_D * (h ** (D - 2))

    # Compton-scale amplitude estimate for one-quantum excitation
    # Note: This is based on Compton wavelength, independent of mass scale multiplier
    A_estimate = lambda_C / math.sqrt(math.pi)

    # Rest length for pretension implementation (always 0 for this model)
    rest_length = 0.0

    # Reference mass for dimensional mapping (Compton-calibrated value)
    # This ensures m_sim and k_sim vary with mass_scale_multiplier
    m_point_reference = rho_D_compton * (h ** D)

    return {
        "dim": D,
        "lambda_C": lambda_C,
        "rho_D": rho_D,
        "T_D": T_D,
        "m_point": m_point,
        "k_spring": k_spring,
        "A_estimate": A_estimate,
        "rest_length": rest_length,
        "mass_scale_multiplier": mass_scale_multiplier,
        "m_point_reference": m_point_reference,
        "c_wave": c_wave,  # Actual wave speed in the brane
        "c_reference": c,  # Speed of light (reference)
    }


def print_calibration_summary(params: dict, grid_spacing: float) -> None:
    """
    Print a human-readable summary of brane calibration parameters.

    Parameters
    ----------
    params : dict
        Output from compton_calibrated_brane_lattice_params() or manual_brane_lattice_params()
    grid_spacing : float
        Lattice spacing h in meters
    """
    D = params['dim']
    dim_label = {1: "1D", 2: "2D", 3: "3D"}[D]

    # Check if this is manual or Compton calibration
    is_manual = "mass_scale_multiplier" in params
    calibration_type = "Manual mass scale" if is_manual else "Compton-cell based"

    print("\n" + "=" * 60)
    print(f"{dim_label} Brane Lattice Calibration ({calibration_type})")
    print("=" * 60)
    print(f"Reduced Compton wavelength λ_C = {params['lambda_C']:.4e} m")
    print(f"Grid spacing h                 = {grid_spacing:.4e} m")
    print(f"Grid spacing / λ_C             = {grid_spacing / params['lambda_C']:.2f}")

    if is_manual:
        multiplier = params['mass_scale_multiplier']
        print(f"Mass scale multiplier          = {multiplier:.4e}")
        if multiplier < 1.0:
            print(f"  → {1.0/multiplier:.2f}× lighter than Compton calibration")
        elif multiplier > 1.0:
            print(f"  → {multiplier:.2f}× heavier than Compton calibration")
        else:
            print(f"  → Same as Compton calibration")

    print()
    print("Continuum parameters:")
    print(f"  Mass density ρ_{D}             = {params['rho_D']:.4e} kg/m^{D}")
    print(f"  Tension/stiffness T_{D}        = {params['T_D']:.4e} N/m^{D-1}")

    c_wave = math.sqrt(params['T_D'] / params['rho_D'])
    print(f"  Wave speed c_wave              = {c_wave:.4e} m/s")

    # Show wave speed ratio for manual calibration
    if is_manual and 'c_reference' in params:
        c_ref = params['c_reference']
        ratio = c_wave / c_ref
        print(f"  → c_wave / c_light              = {ratio:.6f}")
        if ratio < 1.0:
            print(f"  → Wave propagates at {ratio*100:.2f}% speed of light")
        elif ratio > 1.0:
            print(f"  → Wave propagates at {ratio*100:.2f}% speed of light")
        else:
            print(f"  → Wave propagates at speed of light")

    print()
    print("Discrete lattice parameters:")
    print(f"  Point mass m = ρ_{D} h^{D}      = {params['m_point']:.4e} kg")
    print(f"  Spring constant k              = {params['k_spring']:.4e} N/m")
    print(f"  Rest length L_0 (pretension)   = {params['rest_length']:.4e} m")
    print(f"  Amplitude estimate A           = {params['A_estimate']:.4e} m")
    print("=" * 60 + "\n")
