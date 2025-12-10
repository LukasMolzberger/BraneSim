"""
Visualize Circularly Polarized Photon on Straight Path

This test visualizes the initial brane state created from EM→brane mapping
for a circularly polarized photon propagating along a straight line.

Creates 24 2D plots showing:
- 3 slices (XY, XZ, YZ)
- 4 components per slice (x, y, z, amplitude)
- 2 quantities (position, velocity)

Goal: Verify the EM-to-brane mapping produces physically reasonable initial conditions.
"""

import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import torch
import numpy as np

from branesim.core.state import BraneState, Dimensionality
from branesim.core.grid import BraneGrid
from branesim.config.physical_constants import PhysicalConstants
from branesim.physics.dimensional_mapping import DimensionalMapper
from branesim.physics.em_to_brane_mapping import initialize_brane_from_em_fields
from branesim.visualization import visualize_brane_state
from branesim.visualization.em_field_viz import (
    visualize_em_fields_along_centerline,
    visualize_em_field_components_2d
)
from branesim.utils import TestRunManager


def compute_straight_waveguide_em_fields(
    coords_phys: torch.Tensor,
    wavelength: float,
    sigma_transverse: float,
    peak_E_field: float,
    propagation_axis: int = 0,
    constants: PhysicalConstants = None,
) -> tuple[torch.Tensor, torch.Tensor]:
    """
    Compute EM fields for a circularly polarized photon in a straight waveguide.
    """
    if constants is None:
        constants = PhysicalConstants()

    coords_np = coords_phys.cpu().numpy()

    # Define transverse directions based on propagation axis
    if propagation_axis == 0:  # Propagating along x
        tangent_dir = np.array([1.0, 0.0, 0.0])
        normal_dir = np.array([0.0, 1.0, 0.0])
        binormal_dir = np.array([0.0, 0.0, 1.0])
        transverse_coords = coords_np[:, 1:]  # (y, z)
    elif propagation_axis == 1:  # Propagating along y
        tangent_dir = np.array([0.0, 1.0, 0.0])
        normal_dir = np.array([0.0, 0.0, 1.0])
        binormal_dir = np.array([1.0, 0.0, 0.0])
        transverse_coords = coords_np[:, [0, 2]]  # (x, z)
    else:  # Propagating along z
        tangent_dir = np.array([0.0, 0.0, 1.0])
        normal_dir = np.array([1.0, 0.0, 0.0])
        binormal_dir = np.array([0.0, 1.0, 0.0])
        transverse_coords = coords_np[:, :2]  # (x, y)

    k = 2.0 * np.pi / wavelength
    z_long = coords_np[:, propagation_axis]
    r_transverse = np.linalg.norm(transverse_coords, axis=1)

    envelope_transverse = np.exp(-0.5 * (r_transverse / sigma_transverse) ** 2)
    phase = k * z_long

    # For circularly polarized wave, amplitude magnitude is constant (only direction rotates)
    amplitude = peak_E_field * envelope_transverse

    cos_pol = np.cos(phase)
    sin_pol = np.sin(phase)

    E_field_np = amplitude[:, np.newaxis] * (
        cos_pol[:, np.newaxis] * normal_dir[np.newaxis, :] +
        sin_pol[:, np.newaxis] * binormal_dir[np.newaxis, :]
    )

    B_magnitude = amplitude / constants.c
    B_dir = (
        -sin_pol[:, np.newaxis] * normal_dir[np.newaxis, :] +
        cos_pol[:, np.newaxis] * binormal_dir[np.newaxis, :]
    )
    B_field_np = B_magnitude[:, np.newaxis] * B_dir

    device = coords_phys.device
    dtype = coords_phys.dtype
    E_field = torch.from_numpy(E_field_np).to(device=device, dtype=dtype)
    B_field = torch.from_numpy(B_field_np).to(device=device, dtype=dtype)

    return E_field, B_field


def main():
    """Visualize circularly polarized photon on straight path."""
    print("=" * 70)
    print("Circularly Polarized Photon - Straight Path Visualization")
    print("=" * 70)

    # Initialize test run manager
    run_manager = TestRunManager(experiment_name="visualize_circular_photon_straight")
    print(run_manager.get_summary())

    constants = PhysicalConstants()

    # Physical setup
    wavelength_phys = constants.lambda_C
    points_per_wavelength = 20
    h_phys = wavelength_phys / points_per_wavelength
    D = 3
    m_point = 2.861821e-27  # kg

    rho_D = m_point / (h_phys ** D)
    rest_length_phys = constants.rest_length_frac * h_phys

    mapper = DimensionalMapper(
        h_phys=h_phys,
        c_light=constants.c,
        mass_reference=m_point
    )

    h_sim = mapper.to_sim_length(h_phys)

    # Domain size
    nx = 100
    ny = 100
    nz = 100

    print(f"\nGrid Configuration:")
    print(f"  Domain: {nx} × {ny} × {nz} points")
    print(f"  Lattice spacing: {h_phys*1e12:.3f} pm")

    # Auto-select device
    if torch.cuda.is_available():
        device = torch.device('cuda')
        print(f"\n✓ Using NVIDIA GPU")
        dtype = torch.float64
    elif hasattr(torch.backends, 'mps') and torch.backends.mps.is_available():
        device = torch.device('mps')
        print(f"\n✓ Using Apple Silicon GPU (MPS)")
        dtype = torch.float32
    else:
        device = torch.device('cpu')
        print(f"\n✓ Using CPU")
        dtype = torch.float64

    # Create brane state
    state = BraneState((nx, ny, nz), Dimensionality.THREE_D, device, dtype)
    state.initialize_flat_configuration(h_sim)
    initial_positions = state.positions.clone()

    grid = BraneGrid((nx, ny, nz), Dimensionality.THREE_D, h_sim, device)

    # Generate EM fields
    print(f"\nGenerating EM fields...")
    coords_sim = grid.get_spatial_coordinates()
    coords_phys = mapper.to_phys_length(coords_sim)

    omega_phys = 2.0 * np.pi * constants.c / wavelength_phys
    sigma_transverse = 3.0 * wavelength_phys

    # Target amplitude: Use smaller value to avoid float32 overflow
    A_target = 0.005 * h_phys
    # From energy matching: u_EM = (1/2)ρω²A² gives E₀ = ωA√(ρ/(2ε₀))
    peak_E_field = omega_phys * A_target * np.sqrt(rho_D / (2.0 * constants.epsilon0))

    print(f"  Target amplitude: {A_target*1e12:.3f} pm")
    print(f"  Peak E-field: {peak_E_field:.6e} V/m")

    E_field_phys, B_field_phys = compute_straight_waveguide_em_fields(
        coords_phys=coords_phys,
        wavelength=wavelength_phys,
        sigma_transverse=sigma_transverse,
        peak_E_field=peak_E_field,
        propagation_axis=0,
        constants=constants,
    )

    # ========================================================================
    # VISUALIZE EM FIELDS ALONG CENTERLINE (BEFORE MAPPING)
    # ========================================================================
    print(f"\nGenerating EM field visualizations...")

    # Extract centerline: points along x-axis at center of domain (y=0, z=0)
    domain_length_x = nx * h_phys
    domain_length_y = ny * h_phys
    domain_length_z = nz * h_phys

    # Sample centerline at center of transverse directions
    num_centerline_points = 200
    x_centerline = np.linspace(0, domain_length_x, num_centerline_points)
    y_centerline = np.full(num_centerline_points, domain_length_y / 2.0)
    z_centerline = np.full(num_centerline_points, domain_length_z / 2.0)

    centerline_coords = np.stack([x_centerline, y_centerline, z_centerline], axis=1)
    centerline_coords_torch = torch.from_numpy(centerline_coords).to(
        device=device, dtype=dtype
    )

    # Compute E and B fields along centerline
    E_field_centerline, B_field_centerline = compute_straight_waveguide_em_fields(
        coords_phys=centerline_coords_torch,
        wavelength=wavelength_phys,
        sigma_transverse=sigma_transverse,
        peak_E_field=peak_E_field,
        propagation_axis=0,
        constants=constants,
    )

    # Convert to numpy for visualization
    E_field_np = E_field_centerline.cpu().numpy()
    B_field_np = B_field_centerline.cpu().numpy()

    # Generate 3D visualization of EM fields
    em_3d_path = run_manager.get_plot_path("em_fields_3d_views.png")
    visualize_em_fields_along_centerline(
        centerline=centerline_coords,
        E_field=E_field_np,
        B_field=B_field_np,
        output_path=em_3d_path,
        title="EM Fields: Circularly Polarized Photon (Straight Waveguide)",
        arrow_scale=1.5,
        subsample_step=5,
        figsize=(15, 5),
        dpi=150
    )
    print(f"  ✓ Saved EM 3D views: em_fields_3d_views.png")

    # Generate 2D component plots
    em_2d_path = run_manager.get_plot_path("em_field_components.png")
    visualize_em_field_components_2d(
        centerline=centerline_coords,
        E_field=E_field_np,
        B_field=B_field_np,
        output_path=em_2d_path,
        propagation_axis=0,
        title="EM Field Components Along Propagation (X-axis)",
        figsize=(12, 8),
        dpi=150
    )
    print(f"  ✓ Saved EM components: em_field_components.png")

    # Apply EM → Brane mapping
    print(f"\nApplying EM → Brane mapping...")
    initialize_brane_from_em_fields(
        state=state,
        grid=grid,
        mapper=mapper,
        m_point_phys=m_point,
        h_phys=h_phys,
        omega_phys=omega_phys,
        E_field_phys=E_field_phys,
        B_field_phys=B_field_phys,
        epsilon_eff=constants.epsilon0,
        mu_eff=constants.mu0,
        c_light=constants.c,
        field_component=3,
        max_amplitude_fraction_of_h=0.1,
        velocity_clip_to_c=True,
    )

    print(f"  ✓ Brane state initialized")

    # Save configuration
    config = {
        "experiment": "Circularly Polarized Photon - Straight Path",
        "grid_size": f"{nx}×{ny}×{nz}",
        "lattice_spacing": f"{h_phys*1e12:.3f} pm",
        "wavelength": f"{wavelength_phys*1e9:.3f} nm",
        "target_amplitude": f"{A_target*1e12:.3f} pm",
        "peak_E_field": f"{peak_E_field:.6e} V/m",
        "device": str(device),
    }
    run_manager.save_config(config)

    # ========================================================================
    # VISUALIZE BRANE STATE
    # ========================================================================
    print(f"\nGenerating brane state visualizations...")
    saved_files = visualize_brane_state(
        state=state,
        grid=grid,
        mapper=mapper,
        output_dir=run_manager.plots_dir,  # Use plots directory
        filename_prefix="polarized_photon",
        initial_positions=initial_positions,
        print_stats=True,
        dpi=150,
        csv_output_dir=run_manager.data_dir  # Use data directory for CSV files
    )

    print(f"\nExpectations:")
    print(f"  • Amplitude should show Gaussian transverse profile with λ_C oscillations along x")
    print(f"  • Lateral displacements (x,y,z) should encode Poynting vector (flow pattern)")
    print(f"  • Circular polarization: rotating pattern in y-z plane")
    print(f"  • Velocities encode time derivative of positions")

    print(f"\n{'=' * 70}")
    print(f"All outputs saved to: {run_manager.run_dir}")
    print(f"{'=' * 70}")


if __name__ == '__main__':
    main()