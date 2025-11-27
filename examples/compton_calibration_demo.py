"""
Example: Using Compton-cell calibration for 3D brane simulations.

This demonstrates how to use the amplitude scale calibration from the paper
(Section "Amplitude scale calibration") to set up a physically-motivated
3D brane simulation.
"""

from branesim.config.simulation_config import SimulationConfig
from branesim.core.state import Dimensionality
from branesim.physics.parameters import brane_lattice_params_3d, print_calibration_summary


def example_direct_calibration():
    """Example 1: Using the calibration function directly."""
    print("=" * 70)
    print("EXAMPLE 1: Direct use of brane_lattice_params_3d()")
    print("=" * 70)

    # Set grid spacing at 10× Compton wavelength
    lambda_C = 2.4263102367e-12  # meters
    h = 10.0 * lambda_C

    # Compute calibrated parameters
    params = brane_lattice_params_3d(
        grid_spacing_m=h,
        use_compton_default=True  # Use Compton-cell assumption
    )

    # Print summary
    print_calibration_summary(params, h)

    # Show how to use in simulation
    print("Usage in simulation setup:")
    print(f"  mass_per_point   = {params['m_point']:.4e} kg")
    print(f"  spring_constant  = {params['k_spring']:.4e} N/m")
    print()


def example_simulation_config():
    """Example 2: Using SimulationConfig.from_compton_calibration()."""
    print("=" * 70)
    print("EXAMPLE 2: Creating SimulationConfig from Compton calibration")
    print("=" * 70)

    # Create configuration for a 32³ grid
    config = SimulationConfig.from_compton_calibration(
        grid_shape=(32, 32, 32),
        dimension=Dimensionality.THREE_D,
        lambda_C_multiplier=10.0,  # h = 10 λ_C
        cfl_factor=0.4,             # CFL number
        critical_strain=0.1,        # Optional: ε_cr for saturation
        device='cpu',
        dtype='float64'
    )

    # The config is ready to use!
    print(config)
    print()

    # Verify wave speed
    expected_c, computed_c, error = config.verify_wave_speed()
    print(f"Wave speed verification:")
    print(f"  Expected (physical): {expected_c:.6e} m/s")
    print(f"  Computed (discrete): {computed_c:.6e} m/s")
    print(f"  Relative error:      {error*100:.4f}%")
    print()


def example_parameter_scan():
    """Example 3: Scan over different grid spacings."""
    print("=" * 70)
    print("EXAMPLE 3: Parameter scan at different grid spacings")
    print("=" * 70)

    lambda_C = 2.4263102367e-12  # meters

    print(f"{'Multiplier':<12} {'h (m)':<15} {'m_point (kg)':<18} {'k_spring (N/m)':<18}")
    print("-" * 70)

    for multiplier in [5.0, 10.0, 20.0, 50.0, 100.0]:
        h = multiplier * lambda_C
        params = brane_lattice_params_3d(h)

        print(f"{multiplier:<12.1f} {h:<15.4e} {params['m_point']:<18.4e} {params['k_spring']:<18.4e}")

    print()


def example_custom_density():
    """Example 4: Using custom mass density instead of Compton-cell assumption."""
    print("=" * 70)
    print("EXAMPLE 4: Custom mass density (not using Compton-cell assumption)")
    print("=" * 70)

    h = 1e-13  # meters
    custom_rho = 1e10  # kg/m³ (custom choice)

    params = brane_lattice_params_3d(
        grid_spacing_m=h,
        use_compton_default=False,  # Use custom density
        rho_mass_density=custom_rho
    )

    print(f"Grid spacing h = {h:.4e} m")
    print(f"Custom mass density ρ = {custom_rho:.4e} kg/m³")
    print(f"Bulk modulus K = ρ c² = {params['K']:.4e} Pa")
    print(f"Point mass m = ρ h³ = {params['m_point']:.4e} kg")
    print(f"Spring constant k = K h = {params['k_spring']:.4e} N/m")
    print()


if __name__ == "__main__":
    # Run all examples
    example_direct_calibration()
    print("\n\n")

    example_simulation_config()
    print("\n\n")

    example_parameter_scan()
    print("\n\n")

    example_custom_density()