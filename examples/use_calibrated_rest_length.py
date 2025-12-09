"""
Example: Using physically calibrated rest_length in a simulation.

This demonstrates the proper workflow:
1. Physical calibration determines rest_length_frac = L0/a from continuum theory
2. Simulation inherits this value as input (doesn't redefine it)
3. Dimensional mapping provides numerical convenience without changing physics
"""

import json
from branesim.tools.calibrate_physical_rest_length import (
    MicroParams,
    calibrate_rest_length,
    print_results
)
from branesim.config.physical_constants import PhysicalConstants


def example_calibration_workflow():
    """
    Demonstrate the two-stage workflow: physical calibration → simulation setup.
    """

    print("=" * 70)
    print("STAGE 1: Physical Calibration (SI Units)")
    print("=" * 70)

    # Define physical micro-parameters
    params = MicroParams(
        rho_m=1.0e18,           # kg/m³ - brane mass density
        spring_k=1.0e26,        # N/m - microscopic spring constant
        lattice_spacing=1.0e-9, # m - lattice spacing
        mu_G=1.0,               # dimensionless bending factor
        ell_star_factor=1.0,    # ℓ_* = a
    )

    # Perform physical calibration
    result = calibrate_rest_length(params)
    print_results(result)

    # Extract the physically calibrated rest_length
    rest_length_frac_physical = result['rest_length']  # This is L0/a
    L0_physical = result['L0']  # This is L0 in meters

    print("=" * 70)
    print("STAGE 2: Simulation Setup")
    print("=" * 70)

    # Simulation parameters (these can be chosen for numerical convenience)
    h_sim = 1.0e-9  # Grid spacing in simulation [m] (same as lattice_spacing)

    # IMPORTANT: The simulation INHERITS the physical rest_length_frac
    # It does NOT redefine or recompute it
    rest_length_frac_sim = rest_length_frac_physical  # Simply copy it
    rest_length_sim = rest_length_frac_sim * h_sim    # L0 = (L0/a) * a

    print(f"\nSimulation inherits physical calibration:")
    print(f"  h_sim               : {h_sim:.6e} m")
    print(f"  rest_length_frac_sim: {rest_length_frac_sim:.6f} (inherited from physics)")
    print(f"  rest_length_sim     : {rest_length_sim:.6e} m")

    # Verify consistency
    assert abs(rest_length_sim - L0_physical) < 1e-20, "Inconsistency detected!"
    print(f"\n✓ Consistency check passed: rest_length_sim = L0_physical")

    print("\n" + "=" * 70)
    print("STAGE 3: Dimensional Mapping (Optional)")
    print("=" * 70)

    # If using dimensionless units, the mapping layer just rescales
    # BUT it does NOT change rest_length_frac (already dimensionless)
    print(f"\nDimensionless mapping:")
    print(f"  length_scale: {h_sim:.6e} m (makes h_sim → 1.0)")
    print(f"  rest_length_frac: {rest_length_frac_sim:.6f} (UNCHANGED - already dimensionless)")
    print(f"  rest_length_dimensionless: 1.0 * {rest_length_frac_sim:.6f} = {rest_length_frac_sim:.6f}")

    print("\n" + "=" * 70)
    print("Summary")
    print("=" * 70)
    print(f"\nThe key principle:")
    print(f"  • Physical calibration DETERMINES rest_length_frac from continuum theory")
    print(f"  • Simulation INHERITS this value as input")
    print(f"  • Dimensional mapping PRESERVES dimensionless ratios like rest_length_frac")
    print(f"\nResult: rest_length_frac = {rest_length_frac_physical:.6f}")
    print(f"        (This would replace the hardcoded 0.1 in experiments)")

    return result


def save_calibration_for_experiment(output_path: str = "config/rest_length_calibrated.json"):
    """
    Save calibration result for use in experiments.

    This creates a JSON file that experiments can load to get the
    physically correct rest_length_frac value.
    """
    params = MicroParams(
        rho_m=1.0e18,
        spring_k=1.0e26,
        lattice_spacing=1.0e-9,
    )

    result = calibrate_rest_length(params)

    # Create experiment-ready config
    experiment_config = {
        "rest_length_frac": result['rest_length'],  # The key output
        "rest_length_meters": result['L0'],
        "calibration_metadata": {
            "rho_m": result['rho_m'],
            "spring_k": result['spring_k'],
            "lattice_spacing": result['lattice_spacing'],
            "T_target": result['T_target'],
            "kappa_target": result['kappa_target'],
        }
    }

    with open(output_path, 'w') as f:
        json.dump(experiment_config, f, indent=2)

    print(f"\nSaved experiment config to {output_path}")
    print(f"Experiments can now use: rest_length_frac = {result['rest_length']:.6f}")

    return experiment_config


if __name__ == "__main__":
    print("\n" + "=" * 70)
    print("Example: Physical Calibration → Simulation Workflow")
    print("=" * 70 + "\n")

    # Run the complete workflow example
    result = example_calibration_workflow()

    print("\n\nTo use in actual experiments, run:")
    print("  python -m branesim.tools.calibrate_physical_rest_length \\")
    print("    --rho-m 1.0e18 \\")
    print("    --spring-k 1.0e26 \\")
    print("    --lattice-spacing 1.0e-9 \\")
    print("    --output-json config/rest_length_physical.json")
    print("\nThen in your experiment:")
    print("  import json")
    print("  with open('config/rest_length_physical.json') as f:")
    print("      calib = json.load(f)")
    print("  rest_length_frac = calib['rest_length']")
    print("  rest_length = rest_length_frac * h  # h is your grid spacing")