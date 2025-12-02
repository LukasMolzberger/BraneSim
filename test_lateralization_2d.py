"""
Quick test to demonstrate lateralization measurement works in 2D.

In 2D, there ARE lateral directions (y-direction perpendicular to x-propagation),
so we should see non-zero R_lat unlike the 1D case.
"""

import sys
import os
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

import torch
import numpy as np

from branesim.core.state import BraneState, Dimensionality
from branesim.core.grid import BraneGrid
from branesim.diagnostics.lateralization import (
    LateralizationMeasurement,
    LateralizationConfig,
)

def test_lateralization_2d():
    """Test lateralization measurement in 2D with realistic photon-like setup."""

    # Small 2D grid
    nx, ny = 20, 10
    h = 1.0  # Spacing

    grid = BraneGrid((nx, ny), Dimensionality.TWO_D, h, torch.device('cpu'))
    state = BraneState((nx, ny), Dimensionality.TWO_D, torch.device('cpu'), torch.float64)
    state.initialize_flat_configuration(h)

    # Get coordinates
    coords = grid.get_spatial_coordinates()
    x = coords[:, 0]
    y = coords[:, 1]

    # Create a Gaussian wave packet centered at x=10, y=5
    center_x = 10.0
    center_y = 5.0
    sigma_x = 3.0
    sigma_y = 2.0
    wavelength = 4.0
    k = 2 * np.pi / wavelength

    # Position field (amplitude in X⁴)
    envelope = torch.exp(-((x - center_x) ** 2) / (2 * sigma_x ** 2)
                        -((y - center_y) ** 2) / (2 * sigma_y ** 2))
    state.positions[:, 3] = 0.1 * envelope * torch.cos(k * (x - center_x))

    # Set velocities:
    # - x-velocities (parallel to propagation)
    # - y-velocities (PERPENDICULAR to propagation - this is lateral!)
    # - amplitude velocities

    # Give some x-velocity (longitudinal)
    state.velocities[:, 0] = 0.5 * envelope

    # Give some y-velocity (LATERAL) - stronger in certain regions
    state.velocities[:, 1] = 1.0 * envelope * torch.sin(k * (x - center_x))

    # Give amplitude velocity
    state.velocities[:, 3] = 0.3 * envelope

    # Set up lateralization measurement
    photon_tube_radius = 2.0 * sigma_x  # Smaller than extent

    lat_config = LateralizationConfig(
        propagation_axis=0,  # x-axis
        amplitude_component=3,  # X⁴
        photon_tube_radius=photon_tube_radius,
        count_amplitude_outside_tube_as_lateral=True
    )
    lateralization = LateralizationMeasurement(lat_config, grid, m_point=1.0)

    # Create dummy physics for potential energy (won't use it for this test)
    class DummyPhysics:
        spring_constant = 1.0
        rest_length = 1.0
    physics = DummyPhysics()

    # Measure
    R_lat_local, R_lat_global, diagnostics = lateralization.measure(state, physics)

    # Print results
    print("=" * 60)
    print("2D Lateralization Measurement Test")
    print("=" * 60)
    print(f"\nConfiguration:")
    print(f"  Grid: {nx} × {ny} points")
    print(f"  Propagation axis: x (index 0)")
    print(f"  Photon tube radius: {photon_tube_radius:.2f}")
    print(f"  σ_x: {sigma_x:.2f}, σ_y: {sigma_y:.2f}")

    print(f"\nVelocity setup:")
    print(f"  v_x (longitudinal): 0.5 × envelope")
    print(f"  v_y (LATERAL): 1.0 × envelope × sin(kx)")
    print(f"  v_amplitude: 0.3 × envelope")

    print(f"\nMeasurement Results:")
    print(f"  Global R_lat: {R_lat_global:.4f}")
    print(f"  Points in tube: {diagnostics['in_tube'].sum().item()} / {nx * ny}")

    # Find point with maximum R_lat
    max_idx = torch.argmax(R_lat_local)
    max_R_lat = R_lat_local[max_idx].item()

    print(f"  Max local R_lat: {max_R_lat:.4f}")
    print(f"  Min local R_lat: {R_lat_local.min().item():.4f}")
    print(f"  Mean local R_lat: {R_lat_local.mean().item():.4f}")

    # Show energy breakdown at the maximum point
    print(f"\nEnergy breakdown at max R_lat point (index {max_idx}):")
    print(f"  E_long_kin: {diagnostics['E_long_kin'][max_idx].item():.6e} J")
    print(f"  E_lat_kin: {diagnostics['E_lat_kin'][max_idx].item():.6e} J")
    print(f"  Ratio: {diagnostics['E_lat_kin'][max_idx].item() / (diagnostics['E_long_kin'][max_idx].item() + 1e-16):.4f}")

    # Check a few specific points
    print(f"\nSample of local R_lat values (first 10 points):")
    for i in range(min(10, nx * ny)):
        r_lat = R_lat_local[i].item()
        in_tube = diagnostics['in_tube'][i].item()
        print(f"  Point {i:2d}: R_lat = {r_lat:.4f}, in_tube = {in_tube}")

    print(f"\n{'=' * 60}")

    # Check that we have non-zero lateral energy (unlike 1D!)
    if R_lat_global > 0.01:
        print("✓ SUCCESS: Lateralization measurement is NON-ZERO in 2D!")
        print("  This is because we have genuine transverse (y) motion.")
    else:
        print("⚠ WARNING: R_lat is very low, expected higher values")

    print(f"{'=' * 60}\n")

    return R_lat_global

if __name__ == '__main__':
    R_lat = test_lateralization_2d()