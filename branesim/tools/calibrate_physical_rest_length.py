"""
Physical calibration of spring rest_length from continuum theory.

This script works entirely in SI units to determine a physically plausible
rest_length = L0/a for the microscopic lattice model. It:

1. Takes physical inputs: ρ_m [kg/m³], k [N/m], a [m]
2. Computes target continuum constants: T = ρ_m c², κ = c³/(8π G)
3. Inverts the 1D chain model to get rest_length that yields T_target
4. Provides optional gravitational diagnostics

This calibration is INDEPENDENT of the dimensionless simulation mapping.
The mapping layer should INHERIT this rest_length as input, not redefine it.
"""

import argparse
import json
import math
from dataclasses import dataclass

from branesim.config.physical_constants import PhysicalConstants


@dataclass
class MicroParams:
    """Microscopic physical parameters in SI units."""
    rho_m: float           # brane mass density [kg/m³]
    spring_k: float        # microscopic spring constant k [N/m]
    lattice_spacing: float # lattice spacing a [m]

    # optional: for the gravitational "diagnostic" only
    mu_G: float = 1.0          # dimensionless O(1) factor for bending
    ell_star_factor: float = 1.0  # ℓ_* = ell_star_factor * a


def compute_target_continuum_constants(rho_m: float,
                                      constants: PhysicalConstants = None):
    """
    Compute physical target tension T and bending stiffness κ
    from ρ_m, c, G.

    Uses continuum relations from the paper:
    - T_target = ρ_m c²
    - κ_target = c³ / (8π G)

    Parameters
    ----------
    rho_m : float
        Brane mass density [kg/m³]
    constants : PhysicalConstants, optional
        Physical constants. If None, uses default values.

    Returns
    -------
    T_target : float
        Target tension [J/m³]
    kappa_target : float
        Target bending stiffness [J/m]
    """
    if constants is None:
        constants = PhysicalConstants()

    T_target = constants.compute_target_tension(rho_m)
    kappa_target = constants.compute_target_bending_stiffness()

    return T_target, kappa_target


def tension_from_rest_length(rest_length: float,
                             spring_k: float,
                             lattice_spacing: float) -> float:
    """
    Physical continuum tension T [J/m³] from given dimensionless
    rest_length = L0/a, with a, k in SI units.

    From 1D chain model:
    T = (k/a) * (1 - rest_length)

    Parameters
    ----------
    rest_length : float
        Dimensionless rest length L0/a, must be in (0, 1)
    spring_k : float
        Spring constant [N/m]
    lattice_spacing : float
        Lattice spacing a [m]

    Returns
    -------
    T : float
        Continuum tension [J/m³]
    """
    if not (0.0 < rest_length < 1.0):
        raise ValueError(f"rest_length must be in (0,1), got {rest_length}")

    a = lattice_spacing
    return (spring_k / a) * (1.0 - rest_length)


def rest_length_from_tension(T_target: float,
                             spring_k: float,
                             lattice_spacing: float) -> float:
    """
    Invert T = (k/a) * (1 - rest_length) to obtain a physical
    rest_length = L0/a that yields T_target.

    Parameters
    ----------
    T_target : float
        Target continuum tension [J/m³]
    spring_k : float
        Spring constant [N/m]
    lattice_spacing : float
        Lattice spacing a [m]

    Returns
    -------
    rest_length : float
        Dimensionless rest length L0/a

    Raises
    ------
    ValueError
        If computed rest_length is outside (0, 1)
    """
    a = lattice_spacing
    rest_length = 1.0 - (T_target * a) / spring_k

    if rest_length <= 0.0 or rest_length >= 1.0:
        raise ValueError(
            f"Computed rest_length={rest_length} outside (0,1). "
            "Adjust rho_m, spring_k, or lattice_spacing."
        )
    return rest_length


def estimate_micro_kappa(T: float, params: MicroParams) -> float:
    """
    Estimate microscopic bending stiffness from tension and length scale.

    Simple ansatz: κ_micro ≈ μ_G * T * ℓ_*²
    where ℓ_* = ell_star_factor * a

    Parameters
    ----------
    T : float
        Tension [J/m³]
    params : MicroParams
        Microscopic parameters including mu_G and ell_star_factor

    Returns
    -------
    kappa_micro : float
        Microscopic bending stiffness [J/m]
    """
    a = params.lattice_spacing
    ell_star = params.ell_star_factor * a
    return params.mu_G * T * (ell_star ** 2)


def effective_zeta_G(kappa_micro: float, ell_star: float,
                     constants: PhysicalConstants = None) -> float:
    """
    Compute effective ζ_G that would match Einstein's G.

    From: G = c³ / (16π ζ_G κ_micro ℓ_*²)
    Therefore: ζ_G^eff = c³ / (16π G κ_micro ℓ_*²)

    Parameters
    ----------
    kappa_micro : float
        Microscopic bending stiffness [J/m]
    ell_star : float
        Length scale [m]
    constants : PhysicalConstants, optional
        Physical constants. If None, uses default values.

    Returns
    -------
    zeta_G_eff : float
        Effective geometric factor (dimensionless)
    """
    if constants is None:
        constants = PhysicalConstants()

    numerator = constants.c ** 3
    denominator = 16.0 * math.pi * constants.G * kappa_micro * (ell_star ** 2)
    return numerator / denominator


def calibrate_rest_length(params: MicroParams,
                          constants: PhysicalConstants = None) -> dict:
    """
    Perform full physical calibration of rest_length.

    Parameters
    ----------
    params : MicroParams
        Microscopic physical parameters
    constants : PhysicalConstants, optional
        Physical constants. If None, uses default values.

    Returns
    -------
    result : dict
        Dictionary containing all calibration results
    """
    if constants is None:
        constants = PhysicalConstants()

    # 1) Target continuum physics
    T_target, kappa_target = compute_target_continuum_constants(
        params.rho_m, constants
    )

    # 2) Physical rest_length (L0/a)
    rest_length = rest_length_from_tension(
        T_target=T_target,
        spring_k=params.spring_k,
        lattice_spacing=params.lattice_spacing,
    )

    # 3) Diagnostics
    T_micro = tension_from_rest_length(
        rest_length=rest_length,
        spring_k=params.spring_k,
        lattice_spacing=params.lattice_spacing,
    )
    kappa_micro = estimate_micro_kappa(T_micro, params)
    ell_star = params.ell_star_factor * params.lattice_spacing
    zeta_G = effective_zeta_G(kappa_micro, ell_star, constants)

    L0 = rest_length * params.lattice_spacing

    return {
        "rho_m": params.rho_m,
        "spring_k": params.spring_k,
        "lattice_spacing": params.lattice_spacing,
        "rest_length": rest_length,  # physically calibrated L0/a
        "L0": L0,
        "T_target": T_target,
        "kappa_target": kappa_target,
        "T_micro": T_micro,
        "kappa_micro": kappa_micro,
        "mu_G": params.mu_G,
        "ell_star_factor": params.ell_star_factor,
        "ell_star": ell_star,
        "zeta_G_effective": zeta_G,
    }


def print_results(result: dict) -> None:
    """Print calibration results in a readable format."""
    print("\n=== Physical rest_length calibration (SI) ===")
    print(f"ρ_m             : {result['rho_m']:.6e}  kg/m³")
    print(f"k (spring)      : {result['spring_k']:.6e}  N/m")
    print(f"a (lattice)     : {result['lattice_spacing']:.6e}  m")
    print(f"T_target        : {result['T_target']:.6e}  J/m³")
    print(f"κ_target        : {result['kappa_target']:.6e}  J/m")
    print(f"rest_length     : {result['rest_length']:.6e}  (L0/a, physical)")
    print(f"L0              : {result['L0']:.6e}  m")
    print(f"T_micro (check) : {result['T_micro']:.6e}  J/m³")
    print(f"κ_micro (ansatz): {result['kappa_micro']:.6e}  J/m")
    print(f"ℓ_*             : {result['ell_star']:.6e}  m")
    print(f"ζ_G (effective) : {result['zeta_G_effective']:.6e}")
    print("============================================\n")


def main(argv=None):
    """
    Command-line interface for physical rest_length calibration.

    Example usage:
        python -m branesim.tools.calibrate_physical_rest_length \
          --rho-m 1.0e18 \
          --spring-k 1.0e5 \
          --lattice-spacing 1.0e-9 \
          --output-json config/rest_length_physical.json
    """
    parser = argparse.ArgumentParser(
        description="Calibrate a physically plausible spring rest_length "
                    "for the brane model (pure SI units)."
    )
    parser.add_argument("--rho-m", type=float, required=True,
                        help="Brane mass density ρ_m [kg/m³].")
    parser.add_argument("--spring-k", type=float, required=True,
                        help="Microscopic spring constant k [N/m].")
    parser.add_argument("--lattice-spacing", type=float, required=True,
                        help="Lattice spacing a [m].")
    parser.add_argument("--mu-G", type=float, default=1.0,
                        help="Dimensionless bending factor μ_G (diagnostic only).")
    parser.add_argument("--ell-star-factor", type=float, default=1.0,
                        help="ℓ_* = ell_star_factor * a (diagnostic only).")
    parser.add_argument("--output-json", type=str, default=None,
                        help="Optional path to write calibration result as JSON.")
    args = parser.parse_args(argv)

    params = MicroParams(
        rho_m=args.rho_m,
        spring_k=args.spring_k,
        lattice_spacing=args.lattice_spacing,
        mu_G=args.mu_G,
        ell_star_factor=args.ell_star_factor,
    )

    try:
        result = calibrate_rest_length(params)
        print_results(result)

        if args.output_json:
            with open(args.output_json, "w") as f:
                json.dump(result, f, indent=2)
            print(f"Wrote calibration result to {args.output_json}")

        return 0

    except ValueError as e:
        print(f"\nError: {e}")
        print("\nSuggestions:")
        print("  - Increase spring_k to support higher tension")
        print("  - Decrease rho_m to lower the required tension")
        print("  - Adjust lattice_spacing")
        return 1


if __name__ == "__main__":
    raise SystemExit(main())