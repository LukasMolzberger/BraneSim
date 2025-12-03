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


def map_to_dimensionless_params(
    phys_params: dict,
    cfl_factor: float = 0.1,
) -> dict:
    """
    Map physical parameters to dimensionless simulation parameters.

    This is the core scaling layer that converts any physical parameter set
    to clean dimensionless units for numerical simulation:
        h_sim = 1.0
        m_point_sim = 1.0
        k_spring_sim = 1.0
        c_sim = 1.0

    The solver runs in these dimensionless units, and all physical quantities
    are mapped via scaling factors:
        - Length scale L0 = h_phys
        - Time scale T0 = L0 / c_phys
        - Mass scale M0 = m_point_phys

    This ensures clean O(1) numerics without changing any physical ratios,
    regardless of which calibration method produced the physical parameters.

    Parameters
    ----------
    phys_params : dict
        Physical parameter dictionary containing at minimum:
        {
            "h_phys": grid spacing [m],
            "m_point": mass per lattice point [kg],
            "k_spring": spring constant [N/m],
            "c_phys": wave speed [m/s],
            "rest_length": spring rest length [m],
            "dim": dimensionality (1, 2, or 3),
        }
        May also contain rho_D, T_D, lambda_C, etc. depending on calibration method.
    cfl_factor : float, optional
        CFL factor for time step calculation, default is 0.1.

    Returns
    -------
    params : dict
        Dictionary containing physical, dimensionless, and scaling parameters:
        {
            # Physical parameters (passed through)
            "h_phys": grid spacing [m],
            "m_point_phys": point mass [kg],
            "k_spring_phys": spring constant [N/m],
            "c_phys": wave speed [m/s],
            "rest_length": spring rest length [m],
            "dt_phys": physical time step [s],
            ... (other params from phys_params)

            # Scaling factors
            "L0": length scale [m],
            "T0": time scale [s],
            "M0": mass scale [kg],
            "E0": energy scale [J],

            # Dimensionless simulation parameters
            "h_sim": grid spacing in sim units (always 1.0),
            "m_point_sim": point mass in sim units (always 1.0),
            "k_spring_sim": spring constant in sim units (always 1.0),
            "c_sim": wave speed in sim units (always 1.0),
            "dt_sim": time step in sim units,

            # Metadata
            "cfl_factor": CFL factor,
        }
    """
    # Extract required physical parameters
    h_phys = phys_params["h_phys"]
    m_point_phys = phys_params["m_point"]
    k_spring_phys = phys_params["k_spring"]
    c_phys = phys_params["c_phys"]
    rest_length = phys_params["rest_length"]

    # Physical time step (from CFL condition)
    dt_phys = cfl_factor * h_phys / c_phys

    # Scaling choices
    # L0: Use the physical grid spacing
    L0 = h_phys
    # T0: Choose so that c_sim = 1
    T0 = L0 / c_phys
    # M0: Choose so that m_point_sim = 1
    M0 = m_point_phys
    # E0: Natural energy scale
    E0 = M0 * (L0 / T0) ** 2

    # Dimensionless simulation parameters
    h_sim = 1.0                  # Grid spacing in units of L0
    m_point_sim = 1.0            # Point mass in units of M0
    k_spring_sim = 1.0           # Spring constant chosen so c_sim = 1
    c_sim = 1.0                  # Wave speed in sim units (by construction)
    dt_sim = dt_phys / T0        # Dimensionless time step

    # Build result by copying all physical params and adding simulation params
    result = phys_params.copy()

    # Standardize naming (ensure both m_point and m_point_phys exist)
    result["m_point_phys"] = m_point_phys
    result["k_spring_phys"] = k_spring_phys
    result["dt_phys"] = dt_phys

    # Add scaling factors
    result.update({
        "L0": L0,
        "T0": T0,
        "M0": M0,
        "E0": E0,
    })

    # Add dimensionless simulation parameters
    result.update({
        "h_sim": h_sim,
        "m_point_sim": m_point_sim,
        "k_spring_sim": k_spring_sim,
        "c_sim": c_sim,
        "dt_sim": dt_sim,
        "cfl_factor": cfl_factor,
    })

    return result


def get_dimensionless_params(
    grid_spacing_m: float,
    dimensionality: int,
    c: float = 2.99792458e8,
    lambda_C_multiplier: float = 10.0,
    cfl_factor: float = 0.1,
) -> dict:
    """
    Get dimensionless simulation parameters using Compton calibration.

    This is a convenience function that combines Compton-based physical parameter
    configuration with dimensionless mapping. It's equivalent to:
        phys_params = compton_calibrated_brane_lattice_params(...)
        sim_params = map_to_dimensionless_params(phys_params, cfl_factor)

    For other calibration methods, use map_to_dimensionless_params() directly.

    For unit conversions between physical and simulation units, create a
    DimensionalMapper:
        from branesim.physics.dimensional_mapping import create_mapper_from_params
        params = get_dimensionless_params(...)
        mapper = create_mapper_from_params(params)
        length_sim = mapper.to_sim_length(length_phys)
        length_phys = mapper.to_phys_length(length_sim)

    Parameters
    ----------
    grid_spacing_m : float
        Physical lattice spacing h in meters.
    dimensionality : int
        Spatial dimensionality D of the brane model (1, 2, or 3).
    c : float, optional
        Speed of light in m/s, default is 2.99792458e8.
    lambda_C_multiplier : float, optional
        Grid spacing as multiple of Compton wavelength, default is 10.0.
    cfl_factor : float, optional
        CFL factor for time step calculation, default is 0.1.

    Returns
    -------
    params : dict
        Dictionary containing both physical and dimensionless parameters.
        See map_to_dimensionless_params() for details.
    """
    # Step 1: Configure physical parameters using Compton calibration
    phys_params = compton_calibrated_brane_lattice_params(
        grid_spacing_m=grid_spacing_m,
        dimensionality=dimensionality,
        c=c
    )

    # Add required fields for mapping
    phys_params["h_phys"] = grid_spacing_m
    phys_params["c_phys"] = c

    # Step 2: Map to dimensionless parameters
    sim_params = map_to_dimensionless_params(phys_params, cfl_factor)

    # Add Compton-specific metadata
    sim_params["lambda_C_multiplier"] = lambda_C_multiplier

    # Add physical constants if not already present
    if "hbar" not in sim_params:
        sim_params["hbar"] = 1.054571817e-34  # J*s
    if "m_e" not in sim_params:
        sim_params["m_e"] = 9.1093837015e-31  # kg

    return sim_params


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
