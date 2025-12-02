"""
Physical parameter calibration for brane simulations using Compton-cell calibration.

This module implements the Compton-cell-based calibration procedure described
in the paper (Section "Amplitude scale calibration"), connecting continuum
brane parameters to discrete lattice parameters.
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


def get_dimensionless_params(
    grid_spacing_m: float,
    dimensionality: int,
    c: float = 2.99792458e8,
    lambda_C_multiplier: float = 10.0,
    cfl_factor: float = 0.1,
) -> dict:
    """
    Get dimensionless simulation parameters with scaling layer.

    This function wraps the physical parameter computation and adds a scaling
    layer that converts everything to clean dimensionless units for the solver:
        h_sim = 1.0
        m_point_sim = 1.0
        k_spring_sim = 1.0
        c_sim = 1.0

    The solver runs in these dimensionless units, and all physical quantities
    are mapped via:
        - Length scale L0 (typically lambda_C * lambda_C_multiplier)
        - Time scale T0 = L0 / c_phys
        - Mass scale M0 = m_point_phys

    This ensures clean O(1) numerics without changing any physical ratios.

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
        Dictionary containing both physical and dimensionless parameters:
        {
            # Physical constants
            "c_phys": speed of light [m/s],
            "lambda_C": Compton wavelength [m],
            "m_e": electron mass [kg],
            "hbar": reduced Planck constant [J·s],

            # Physical parameters
            "h_phys": grid spacing [m],
            "m_point_phys": point mass [kg],
            "k_spring_phys": spring constant [N/m],
            "rho_D": D-dim mass density [kg/m^D],
            "T_D": D-dim tension [N/m^(D-1)],
            "dt_phys": physical time step [s],
            "rest_length": spring rest length [m],

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

            # Other
            "dim": dimensionality,
            "A_estimate": amplitude scale estimate [m],
            "lambda_C_multiplier": grid spacing / lambda_C,
            "cfl_factor": CFL factor,
        }
    """
    # Get physical parameters
    phys_params = compton_calibrated_brane_lattice_params(
        grid_spacing_m=grid_spacing_m,
        dimensionality=dimensionality,
        c=c
    )

    # Physical constants
    hbar = 1.054571817e-34      # J*s
    m_e  = 9.1093837015e-31     # kg
    lambda_C = phys_params["lambda_C"]

    # Physical parameters
    h_phys = grid_spacing_m
    m_point_phys = phys_params["m_point"]
    k_spring_phys = phys_params["k_spring"]
    rho_D = phys_params["rho_D"]
    T_D = phys_params["T_D"]

    # Physical time step (from CFL condition)
    dt_phys = cfl_factor * h_phys / c

    # Scaling choices
    # L0: Use the physical grid spacing (which is lambda_C * lambda_C_multiplier)
    L0 = h_phys
    # T0: Choose so that c_sim = 1
    T0 = L0 / c
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

    return {
        # Physical constants
        "c_phys": c,
        "lambda_C": lambda_C,
        "m_e": m_e,
        "hbar": hbar,

        # Physical parameters
        "h_phys": h_phys,
        "m_point_phys": m_point_phys,
        "k_spring_phys": k_spring_phys,
        "rho_D": rho_D,
        "T_D": T_D,
        "dt_phys": dt_phys,
        "rest_length": phys_params["rest_length"],

        # Scaling factors
        "L0": L0,
        "T0": T0,
        "M0": M0,
        "E0": E0,

        # Dimensionless simulation parameters
        "h_sim": h_sim,
        "m_point_sim": m_point_sim,
        "k_spring_sim": k_spring_sim,
        "c_sim": c_sim,
        "dt_sim": dt_sim,

        # Other
        "dim": dimensionality,
        "A_estimate": phys_params["A_estimate"],
        "lambda_C_multiplier": lambda_C_multiplier,
        "cfl_factor": cfl_factor,
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


# Example usage (can be removed or kept for testing)
if __name__ == "__main__":
    # Example: lattice spacing at 10× Compton wavelength
    from branesim.config.simulation_config import PhysicalConstants

    constants = PhysicalConstants()
    h = 10.0 * constants.lambda_C

    # Test all dimensions
    for dim in [1, 2, 3]:
        params = compton_calibrated_brane_lattice_params(h, dimensionality=dim)
        print_calibration_summary(params, h)

    # Show how to use in simulation setup
    params_3d = compton_calibrated_brane_lattice_params(h, dimensionality=3)
    print("Usage in 3D simulation:")
    print(f"  mass_per_point   = {params_3d['m_point']:.4e}  # kg")
    print(f"  spring_constant  = {params_3d['k_spring']:.4e}  # N/m")