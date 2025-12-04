"""
3D Photon in Waveguide Geometry with Realistic Physical Scales

Simulates a photon wave packet propagating through a 3D waveguide.
Uses actual speed of light c = 299,792,458 m/s and physical length scales
based on the Compton wavelength.

Visualization options:
1. 2D slice through middle of volume (default)
2. Volumetric rendering (optional, requires appropriate backend)
3. Multiple orthogonal slices
"""

import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import math
import torch
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, FFMpegWriter
from matplotlib.colors import hsv_to_rgb
from mpl_toolkits.axes_grid1 import make_axes_locatable

from branesim.core.state import BraneState, Dimensionality
from branesim.core.grid import BraneGrid
from branesim.core.solver import VelocityVerletSolver
from branesim.physics.forces import SpringForceComputer
from branesim.config.simulation_config import PhysicalConstants
from branesim.physics.dimensional_mapping import DimensionalMapper


def displacement_to_rgb_3d(disp_x, disp_y, max_magnitude=None):
    """
    Convert 2D displacement vectors (in a slice) to RGB image using HSV color coding.

    Args:
        disp_x: x-component of displacement (2D array)
        disp_y: y-component of displacement (2D array)
        max_magnitude: Maximum magnitude for normalization (if None, uses max of data)

    Returns:
        RGB image (nx, ny, 3) where:
        - Hue encodes direction (0-360 degrees)
        - Saturation = 1.0 (full color)
        - Value encodes magnitude (0-1)
    """
    # Compute magnitude and angle
    magnitude = np.sqrt(disp_x**2 + disp_y**2)
    angle = np.arctan2(disp_y, disp_x)  # Returns angle in [-pi, pi]

    # Normalize angle to [0, 1] for hue
    hue = (angle + np.pi) / (2 * np.pi)

    # Normalize magnitude to [0, 1] for value
    if max_magnitude is None:
        max_magnitude = magnitude.max()

    if max_magnitude > 0:
        value = magnitude / max_magnitude
    else:
        value = np.zeros_like(magnitude)

    # Full saturation
    saturation = np.ones_like(magnitude)

    # Stack into HSV image
    hsv = np.stack([hue, saturation, value], axis=-1)

    # Convert to RGB
    rgb = hsv_to_rgb(hsv)

    return rgb, magnitude, angle


def initialize_waveguide_wave_shape_3d(state, grid, amplitude, center_x, width_x, wavelength, width_y=None, width_z=None):
    """
    Initialize ONLY the shape of a 3D localized photon wave packet.

    Uses Gaussian envelopes in x, y, and z directions to create a single
    localized wave packet centered in the domain. This represents a free-space
    photon propagating in 3D.

    Velocities are initialized separately.

    Args:
        state: BraneState
        grid: BraneGrid
        amplitude: Peak amplitude [m]
        center_x: Center x coordinate [m]
        width_x: Gaussian width σ_x [m]
        wavelength: Carrier wavelength [m]
        width_y: Gaussian width σ_y [m] (if None, uses same as width_x)
        width_z: Gaussian width σ_z [m] (if None, uses same as width_x)
    """
    # Get spatial coordinates
    coords = grid.get_spatial_coordinates()
    x = coords[:, 0]
    y = coords[:, 1]
    z = coords[:, 2]

    # Domain dimensions
    nx, ny, nz = grid.grid_shape
    domain_length_x = (nx - 1) * grid.spacing
    domain_length_y = (ny - 1) * grid.spacing
    domain_length_z = (nz - 1) * grid.spacing

    # Center the wave packet in y and z directions
    center_y = domain_length_y / 2.0
    center_z = domain_length_z / 2.0

    # Use same width in y and z as x if not specified
    if width_y is None:
        width_y = width_x
    if width_z is None:
        width_z = width_x

    # Gaussian envelope in x (propagation direction)
    envelope_x = torch.exp(-((x - center_x) ** 2) / (2 * width_x ** 2))

    # Hard-truncate Gaussian at 4σ to ensure compact support
    # This prevents far-field regions from showing lateral motion before photon arrival
    cutoff_x = 4.0 * width_x
    mask_x = torch.abs(x - center_x) <= cutoff_x
    envelope_x = envelope_x * mask_x

    # Gaussian envelope in y (transverse direction)
    envelope_y = torch.exp(-((y - center_y) ** 2) / (2 * width_y ** 2))

    # Hard-truncate in y as well
    cutoff_y = 4.0 * width_y
    mask_y = torch.abs(y - center_y) <= cutoff_y
    envelope_y = envelope_y * mask_y

    # Gaussian envelope in z (transverse direction)
    envelope_z = torch.exp(-((z - center_z) ** 2) / (2 * width_z ** 2))

    # Hard-truncate in z as well
    cutoff_z = 4.0 * width_z
    mask_z = torch.abs(z - center_z) <= cutoff_z
    envelope_z = envelope_z * mask_z

    # Full 3D Gaussian envelope
    envelope = amplitude * envelope_x * envelope_y * envelope_z

    # Wave number
    k = 2 * np.pi / wavelength

    # Position field only - velocities set separately
    # Carrier wave propagating in x-direction
    state.positions[:, 3] = envelope * torch.cos(k * (x - center_x))

    print(f"  Wavelength λ = {wavelength:.6e} m ({wavelength/grid.spacing:.1f} × h)")
    print(f"  Width σ_x = {width_x:.6e} m ({width_x/wavelength:.2f} × λ)")
    print(f"  Width σ_y = {width_y:.6e} m ({width_y/wavelength:.2f} × λ)")
    print(f"  Width σ_z = {width_z:.6e} m ({width_z/wavelength:.2f} × λ)")
    print(f"  Center (x, y, z) = ({center_x:.6e}, {center_y:.6e}, {center_z:.6e}) m")
    print(f"  Amplitude = {amplitude:.6e} m")
    print(f"  Wave number k = {k:.6e} rad/m")
    print(f"  Gaussian envelope: single localized wave packet")
    print(f"  Free-space propagation (no waveguide confinement)")


def extract_slice_xy(field_3d, grid_shape, z_index=None):
    """
    Extract a 2D slice from 3D field at constant z.

    Args:
        field_3d: 1D array of length nx*ny*nz (flattened 3D field)
        grid_shape: (nx, ny, nz)
        z_index: z-index for slice (default: middle)

    Returns:
        2D array of shape (nx, ny)
    """
    nx, ny, nz = grid_shape
    field = field_3d.reshape(nx, ny, nz)

    if z_index is None:
        z_index = nz // 2

    return field[:, :, z_index]


def extract_slice_xz(field_3d, grid_shape, y_index=None):
    """
    Extract a 2D slice from 3D field at constant y.

    Args:
        field_3d: 1D array of length nx*ny*nz (flattened 3D field)
        grid_shape: (nx, ny, nz)
        y_index: y-index for slice (default: middle)

    Returns:
        2D array of shape (nx, nz)
    """
    nx, ny, nz = grid_shape
    field = field_3d.reshape(nx, ny, nz)

    if y_index is None:
        y_index = ny // 2

    return field[:, y_index, :]


def extract_slice_yz(field_3d, grid_shape, x_index=None):
    """
    Extract a 2D slice from 3D field at constant x.

    Args:
        field_3d: 1D array of length nx*ny*nz (flattened 3D field)
        grid_shape: (nx, ny, nz)
        x_index: x-index for slice (default: middle)

    Returns:
        2D array of shape (ny, nz)
    """
    nx, ny, nz = grid_shape
    field = field_3d.reshape(nx, ny, nz)

    if x_index is None:
        x_index = nx // 2

    return field[x_index, :, :]


def main():
    """Run 3D photon simulation in waveguide geometry with realistic scales."""
    print("=" * 70)
    print("3D Photon in Waveguide - Realistic Physical Scales")
    print("=" * 70)

    # Physical constants
    constants = PhysicalConstants()

    print(f"\nPhysical Constants:")
    print(f"  Speed of light c = {constants.c:.6e} m/s")
    print(f"  Compton wavelength λ_C = {constants.lambda_C:.6e} m")
    print(f"  ℏ = {constants.hbar:.6e} J·s")
    print(f"  m_e = {constants.m_e:.6e} kg")

    # Configuration
    lambda_C_multiplier = 5.0  # Grid spacing = 5 × λ_C
    h_phys = constants.lambda_C * lambda_C_multiplier
    cfl_factor = 0.1

    D = 3

    # 3D brane parameters constrained to give wave speed = c
    m_e = constants.m_e
    rho_D = m_e / (constants.lambda_C ** 3)  # kg/m³ (volume mass density - derived from Compton scale)
    T_D = rho_D * constants.c**2  # Pa (elastic modulus - computed from c² = T_D/rho_D)
    rest_length_phys = 0.0

    # Wave speed (exactly equals c by construction)
    c_wave = constants.c

    # Discrete mass per lattice point
    m_point = rho_D * (h_phys ** D)

    # Axial spring constant (for 3D: k = T_D * h^(D-2) = T_D * h)
    k_spring = T_D * (h_phys ** (D - 2))

    # Create dimensional mapper for unit conversions
    mapper = DimensionalMapper(
        h_phys=h_phys,
        c_light=constants.c,
        mass_reference=m_point
    )

    # Simulation uses dimensionless units
    h_sim = mapper.to_sim_length(h_phys)  # = 1.0 always
    m_sim = mapper.to_sim_mass(m_point)
    k_sim = mapper.to_sim_spring_constant(k_spring)
    c_sim = mapper.to_sim_velocity(c_wave)  # = 1.0 always (c_wave = c)
    rest_length_sim = mapper.to_sim_length(rest_length_phys)

    # Time step calculation (CFL condition based on wave speed)
    dt_phys = cfl_factor * h_phys / c_wave
    dt_sim = mapper.to_sim_time(dt_phys)

    # Domain size - waveguide geometry (long in x, narrow in y and z)
    nx = 200  # Long waveguide
    ny = 40   # Narrow transverse
    nz = 40   # Narrow transverse
    domain_length_phys_x = nx * h_phys
    domain_length_phys_y = ny * h_phys
    domain_length_phys_z = nz * h_phys
    domain_length_sim = nx * h_sim  # = nx * 1.0 = nx

    # Verify wave speed
    expected_wave_speed = math.sqrt(T_D / rho_D)

    print(f"\nPhysical Parameters:")
    print(f"  3D volume mass density ρ_3 = {rho_D:.6e} kg/m³")
    print(f"  3D elastic modulus T_3 = {T_D:.6e} Pa")
    print(f"  Spring constant k = {k_spring:.6e} N/m")
    print(f"  Point mass m = {m_point:.6e} kg")
    print(f"  Time step dt = {dt_phys:.6e} s")
    print(f"  Expected wave speed = {expected_wave_speed:.6e} m/s")
    print(f"  Speed of light c = {constants.c:.6e} m/s")
    print(f"  Wave speed error = {abs(expected_wave_speed - constants.c)/constants.c:.6e}")

    print(f"\nScaling Factors:")
    print(mapper)

    print(f"\nDimensionless Simulation Parameters:")
    print(f"  h_sim = {h_sim:.6e}  (FIXED to 1.0, defines length scale L0)")
    print(f"  c_light_sim = 1.000000e+00  (FIXED to 1.0, defines time scale T0 = L0/c)")
    print(f"  m_point_sim = {m_sim:.6e}")
    print(f"  k_spring_sim = {k_sim:.6e}")
    print(f"  dt_sim = {dt_sim:.6e}  (time step)")
    print(f"  Wave propagation:")
    print(f"    c_wave_sim = √(k_sim/m_sim) = {(k_sim/m_sim)**0.5:.6e}  (= 1.0 always)")

    print(f"\nSimulation Configuration:")
    print(f"  Domain: {nx} × {ny} × {nz} points")
    print(f"  Domain (physical): {domain_length_phys_x:.6e} × {domain_length_phys_y:.6e} × {domain_length_phys_z:.6e} m")
    print(f"  Domain (sim units): {domain_length_sim:.1f} × {domain_length_sim * ny / nx:.1f} × {domain_length_sim * nz / nx:.1f}")
    print(f"  Total points: {nx*ny*nz:,}")
    print(f"  Aspect ratio: {nx}:{ny}:{nz} (waveguide geometry)")
    print(f"  CFL number = {cfl_factor:.3f}")

    # Create components using SIMULATION UNITS
    device = torch.device('cpu')
    dtype = torch.float64

    state = BraneState((nx, ny, nz), Dimensionality.THREE_D, device, dtype)
    state.initialize_flat_configuration(h_sim)  # Use sim spacing = 1.0

    # Store initial positions for lateral distortion tracking
    initial_positions = state.positions.clone()

    # Set fixed boundaries (all 6 walls)
    state.set_fixed_boundaries()
    print(f"\nBoundary Conditions:")
    print(f"  Fixed points: {state.fixed_mask.sum().item()} / {nx * ny * nz}")
    print(f"  All 6 walls: FIXED (waveguide)")

    grid = BraneGrid((nx, ny, nz), Dimensionality.THREE_D, h_sim, device)  # Sim spacing = 1.0

    print(f"\nPretension Implementation (Sim Units):")
    print(f"  Rest length L_0 (sim) = {rest_length_sim:.6e}")
    print(f"  Rest length L_0 (phys) = {rest_length_phys:.6e} m")
    print(f"  Actual spacing (sim) = {h_sim:.6e}")
    print(f"  Spring constant (sim) = {k_sim:.6e}")
    print(f"  Background tension F_0 (sim) = k×(h-L_0) = {k_sim * (h_sim - rest_length_sim):.6e}")

    physics = SpringForceComputer(k_sim, rest_length_sim)
    solver = VelocityVerletSolver(dt_sim, m_sim, physics, grid)

    # Initialize wave packet IN SIMULATION UNITS
    print(f"\nInitializing photon wave packet...")

    # Physical values (what we want in real units)
    wavelength_phys = 40 * h_phys  # 40 points per wavelength (matches 1D/2D)
    amplitude_phys = 0.1 * h_phys
    width_x_phys = 3 * wavelength_phys / (2 * np.pi)
    center_x_phys = domain_length_phys_x / 3.0

    # Convert to sim units using mapper
    wavelength_sim = mapper.to_sim_length(wavelength_phys)  # = 40.0
    amplitude_sim = mapper.to_sim_length(amplitude_phys)    # = 0.1
    width_x_sim = mapper.to_sim_length(width_x_phys)
    center_x_sim = mapper.to_sim_length(center_x_phys)

    print(f"  Physical wavelength: {wavelength_phys:.6e} m")
    print(f"  Sim wavelength: {wavelength_sim:.1f} grid units")
    print(f"  Physical amplitude: {amplitude_phys:.6e} m")
    print(f"  Sim amplitude: {amplitude_sim:.3f} grid units")

    # Initialize shape only (in sim units)
    print(f"\n[1] Initializing wave shape (sim units)...")
    initialize_waveguide_wave_shape_3d(state, grid, amplitude_sim, center_x_sim, width_x_sim, wavelength_sim)

    # Note: Velocities can be initialized here if needed using
    # initialize_right_moving_velocities() from initial_conditions.py

    # Compute initial accelerations
    solver.initialize_accelerations(state)
    state.apply_fixed_boundaries()

    # Initial measurements
    initial_energy = solver.compute_energy(state)

    print(f"\nInitial State:")
    print(f"  Energy = {initial_energy['total']:.6e} J")

    # Run simulation
    # Time for light to cross domain: t = L/c (in physical units)
    crossing_time_phys = domain_length_phys_x / constants.c
    simulation_time_phys = 0.5 * crossing_time_phys  # 0.5 crossings (reduced for faster demo)
    simulation_time_sim = mapper.to_sim_time(simulation_time_phys)

    num_steps = int(simulation_time_sim / dt_sim)

    print(f"\nRunning simulation...")
    print(f"  Light crossing time (phys) = {crossing_time_phys:.6e} s")
    print(f"  Simulation time (phys) = {simulation_time_phys:.6e} s")
    print(f"  Simulation time (sim) = {simulation_time_sim:.3f} time units")
    print(f"  Number of steps = {num_steps:,}")

    # Tracking
    times_phys = []  # Physical times for plotting
    energies = []

    # Take snapshots at regular intervals (in physical time)
    num_snapshots = 7
    snapshot_times_phys = np.linspace(0, simulation_time_phys, num_snapshots)
    snapshots = {}
    snapshots_lateral_x = {}  # x-component of lateral displacement
    snapshots_lateral_y = {}  # y-component of lateral displacement
    snapshots_lateral_z = {}  # z-component of lateral displacement
    snapshot_steps = {int(t / dt_phys): t for t in snapshot_times_phys}

    # For animation - save every N steps
    animation_frames = []
    animation_frames_lateral_x = []
    animation_frames_lateral_y = []
    animation_times = []
    frame_interval = max(1, num_steps // 200)  # ~200 frames total

    print_interval = max(1, num_steps // 20)

    for step in range(num_steps + 1):
        if step in snapshot_steps:
            field = state.get_field_component(3).cpu().numpy()
            snapshots[snapshot_steps[step]] = field.copy()

            # Store lateral displacements (x, y, z components)
            lateral_disp_x = (state.positions[:, 0] - initial_positions[:, 0]).cpu().numpy()
            lateral_disp_y = (state.positions[:, 1] - initial_positions[:, 1]).cpu().numpy()
            lateral_disp_z = (state.positions[:, 2] - initial_positions[:, 2]).cpu().numpy()
            snapshots_lateral_x[snapshot_steps[step]] = lateral_disp_x.copy()
            snapshots_lateral_y[snapshot_steps[step]] = lateral_disp_y.copy()
            snapshots_lateral_z[snapshot_steps[step]] = lateral_disp_z.copy()

        # Save frames for animation
        if step % frame_interval == 0:
            field = state.get_field_component(3).cpu().numpy()
            animation_frames.append(field.copy())

            # Save lateral displacement frames
            lateral_disp_x = (state.positions[:, 0] - initial_positions[:, 0]).cpu().numpy()
            lateral_disp_y = (state.positions[:, 1] - initial_positions[:, 1]).cpu().numpy()
            animation_frames_lateral_x.append(lateral_disp_x.copy())
            animation_frames_lateral_y.append(lateral_disp_y.copy())

            animation_times.append(solver.time)

        if step % max(1, num_steps // 100) == 0:  # Track 100 points
            energy = solver.compute_energy(state)
            time_phys = mapper.to_phys_time(solver.time)  # Convert sim time to physical
            times_phys.append(time_phys)
            energies.append(energy['total'])

        if step % print_interval == 0:
            time_phys = mapper.to_phys_time(solver.time)  # Convert for printing
            print(f"  Step {step:8d}/{num_steps}: t={time_phys:.6e}s, "
                  f"E={energy['total']:.6e}J")

        if step < num_steps:
            solver.step(state)

    # Final analysis
    final_energy = solver.compute_energy(state)

    energy_drift = abs(final_energy['total'] - initial_energy['total']) / initial_energy['total']

    print(f"\n{'=' * 70}")
    print("Results:")
    print(f"{'=' * 70}")
    print(f"\nWave Propagation:")
    print(f"  Photon propagates at c = {constants.c:.6e} m/s")
    print(f"  Waveguide dimensions: {domain_length_phys_x:.6e} × {domain_length_phys_y:.6e} × {domain_length_phys_z:.6e} m")

    print(f"\nEnergy Conservation:")
    print(f"  Initial: {initial_energy['total']:.6e} J")
    print(f"  Final:   {final_energy['total']:.6e} J")
    print(f"  Drift:   {energy_drift:.6e} ({energy_drift*100:.6f}%)")

    if energy_drift < 1e-2:
        print(f"  ✓ Good energy conservation")
    else:
        print(f"  ⚠ Energy drift: {energy_drift*100:.2f}%")

    # ========================================================================
    # VISUALIZATION
    # ========================================================================
    print(f"\nCreating visualizations...")
    print(f"  Using 2D slice through middle of volume (z = {nz//2})")

    # Coordinate arrays (in nanometers) - convert from physical units
    x_coords_phys = np.arange(nx) * h_phys
    y_coords_phys = np.arange(ny) * h_phys
    z_coords_phys = np.arange(nz) * h_phys
    x_coords = x_coords_phys * 1e9
    y_coords = y_coords_phys * 1e9
    z_coords = z_coords_phys * 1e9
    amplitude_nm = amplitude_phys * 1e9

    # ========================================================================
    # 1. AMPLITUDE FIELD VISUALIZATION (XY slice at middle z)
    # ========================================================================
    fig, axes = plt.subplots(num_snapshots, 1, figsize=(14, 12))
    fig.suptitle(f'3D Photon in Waveguide - Amplitude Field (XY slice, z={nz//2})',
                 fontsize=16, fontweight='bold')

    for idx, (step, t) in enumerate(snapshot_steps.items()):
        if t in snapshots:
            # Extract XY slice at middle z (sim units)
            field_slice_sim = extract_slice_xy(snapshots[t], (nx, ny, nz))
            field_slice_nm = mapper.to_phys_length(field_slice_sim) * 1e9  # Convert sim → phys → nm

            # Plot heatmap
            im = axes[idx].imshow(field_slice_nm.T, origin='lower',
                                 extent=[x_coords[0], x_coords[-1],
                                        y_coords[0], y_coords[-1]],
                                 cmap='RdBu_r',
                                 vmin=-amplitude_nm*1.2, vmax=amplitude_nm*1.2,
                                 aspect='auto')

            axes[idx].set_ylabel('y [nm]', fontsize=11)
            axes[idx].set_xlim(x_coords[0], x_coords[-1])
            axes[idx].set_ylim(y_coords[0], y_coords[-1])

            # Time in femtoseconds
            t_fs = t * 1e15
            axes[idx].text(0.02, 0.95, f't = {t_fs:.3f} fs',
                          transform=axes[idx].transAxes,
                          fontsize=12, verticalalignment='top',
                          bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

            # Add colorbar
            divider = make_axes_locatable(axes[idx])
            cax = divider.append_axes("right", size="3%", pad=0.05)
            plt.colorbar(im, cax=cax, label='ξ [nm]')

            if idx == num_snapshots - 1:
                axes[idx].set_xlabel('x [nm]', fontsize=12)

    plt.tight_layout()
    plt.savefig('photon_3d_example_propagation_xy.png', dpi=150, bbox_inches='tight')
    print(f"  ✓ Saved: photon_3d_example_propagation_xy.png")

    # ========================================================================
    # 2. ORTHOGONAL SLICES (XZ at middle y, YZ at middle x)
    # ========================================================================
    fig2, axes2 = plt.subplots(2, num_snapshots, figsize=(18, 6))
    fig2.suptitle(f'3D Photon - Orthogonal Slices', fontsize=16, fontweight='bold')

    for idx, (step, t) in enumerate(snapshot_steps.items()):
        if t in snapshots:
            # XZ slice at middle y (sim units)
            field_xz_sim = extract_slice_xz(snapshots[t], (nx, ny, nz))
            field_xz_nm = mapper.to_phys_length(field_xz_sim) * 1e9  # Convert sim → phys → nm

            # YZ slice at middle x (sim units)
            field_yz_sim = extract_slice_yz(snapshots[t], (nx, ny, nz))
            field_yz_nm = mapper.to_phys_length(field_yz_sim) * 1e9  # Convert sim → phys → nm

            # Plot XZ slice
            im_xz = axes2[0, idx].imshow(field_xz_nm.T, origin='lower',
                                         extent=[x_coords[0], x_coords[-1],
                                                z_coords[0], z_coords[-1]],
                                         cmap='RdBu_r',
                                         vmin=-amplitude_nm*1.2, vmax=amplitude_nm*1.2,
                                         aspect='auto')
            axes2[0, idx].set_ylabel('z [nm]', fontsize=9)
            axes2[0, idx].set_xlabel('x [nm]', fontsize=9)
            t_fs = t * 1e15
            axes2[0, idx].set_title(f't = {t_fs:.2f} fs', fontsize=10)

            # Plot YZ slice
            im_yz = axes2[1, idx].imshow(field_yz_nm.T, origin='lower',
                                         extent=[y_coords[0], y_coords[-1],
                                                z_coords[0], z_coords[-1]],
                                         cmap='RdBu_r',
                                         vmin=-amplitude_nm*1.2, vmax=amplitude_nm*1.2,
                                         aspect='equal')
            axes2[1, idx].set_ylabel('z [nm]', fontsize=9)
            axes2[1, idx].set_xlabel('y [nm]', fontsize=9)

    plt.tight_layout()
    plt.savefig('photon_3d_example_orthogonal_slices.png', dpi=150, bbox_inches='tight')
    print(f"  ✓ Saved: photon_3d_example_orthogonal_slices.png")

    # ========================================================================
    # 3. ENERGY CONSERVATION PLOT
    # ========================================================================
    fig3, ax = plt.subplots(figsize=(10, 6))

    times_fs = np.array(times_phys) * 1e15  # Already in physical units
    energy_array = np.array(energies)
    ax.plot(times_fs, energy_array / initial_energy['total'], 'g-', linewidth=2)
    ax.axhline(y=1.0, color='r', linestyle='--', linewidth=1, alpha=0.5)
    ax.set_xlabel('Time [fs]', fontsize=12)
    ax.set_ylabel('E(t) / E(0)', fontsize=12)
    ax.set_title('Energy Conservation', fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig('photon_3d_example_energy.png', dpi=150, bbox_inches='tight')
    print(f"  ✓ Saved: photon_3d_example_energy.png")

    # ========================================================================
    # 4. LATERAL DISTORTION (XY slice with color-coded direction)
    # ========================================================================
    print(f"\nCreating lateral distortion plots...")

    # Find max lateral displacement magnitude for consistent scaling (convert sim → phys)
    max_disp_mag_phys = 0
    for t in snapshots_lateral_x.keys():
        # Extract XY slice at middle z (sim units)
        disp_x_slice_sim = extract_slice_xy(snapshots_lateral_x[t], (nx, ny, nz))
        disp_y_slice_sim = extract_slice_xy(snapshots_lateral_y[t], (nx, ny, nz))
        # Convert sim → phys
        disp_x_slice_phys = mapper.to_phys_length(disp_x_slice_sim)
        disp_y_slice_phys = mapper.to_phys_length(disp_y_slice_sim)
        mag = np.sqrt(disp_x_slice_phys**2 + disp_y_slice_phys**2).max()
        max_disp_mag_phys = max(max_disp_mag_phys, mag)

    fig_lat, axes_lat = plt.subplots(num_snapshots, 1, figsize=(14, 12))
    fig_lat.suptitle(f'3D Photon - Lateral Distortion (XY slice, Color=Direction, Brightness=Magnitude)',
                     fontsize=16, fontweight='bold')

    for idx, (step, t) in enumerate(snapshot_steps.items()):
        if t in snapshots_lateral_x:
            # Extract XY slice at middle z (sim units)
            disp_x_slice_sim = extract_slice_xy(snapshots_lateral_x[t], (nx, ny, nz))
            disp_y_slice_sim = extract_slice_xy(snapshots_lateral_y[t], (nx, ny, nz))

            # Convert sim → phys
            disp_x_slice_phys = mapper.to_phys_length(disp_x_slice_sim)
            disp_y_slice_phys = mapper.to_phys_length(disp_y_slice_sim)

            # Convert to RGB image (using physical units)
            rgb_image, magnitude, angle = displacement_to_rgb_3d(disp_x_slice_phys, disp_y_slice_phys, max_disp_mag_phys)

            # Plot RGB image (transpose spatial dimensions only, keep color channel last)
            axes_lat[idx].imshow(np.transpose(rgb_image, (1, 0, 2)), origin='lower',
                                extent=[x_coords[0], x_coords[-1],
                                       y_coords[0], y_coords[-1]],
                                aspect='auto')

            axes_lat[idx].set_ylabel('y [nm]', fontsize=11)
            axes_lat[idx].set_xlim(x_coords[0], x_coords[-1])
            axes_lat[idx].set_ylim(y_coords[0], y_coords[-1])

            # Time in femtoseconds
            t_fs = t * 1e15
            axes_lat[idx].text(0.02, 0.95, f't = {t_fs:.3f} fs',
                              transform=axes_lat[idx].transAxes,
                              fontsize=12, verticalalignment='top',
                              bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

            # Add magnitude info (already in phys units)
            max_mag_pm = max_disp_mag_phys * 1e12
            axes_lat[idx].text(0.98, 0.95, f'max: {max_mag_pm:.2f} pm',
                              transform=axes_lat[idx].transAxes,
                              fontsize=10, verticalalignment='top',
                              horizontalalignment='right',
                              bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.8))

            if idx == num_snapshots - 1:
                axes_lat[idx].set_xlabel('x [nm]', fontsize=12)

    plt.tight_layout()
    plt.savefig('photon_3d_example_lateral_distortion.png', dpi=150, bbox_inches='tight')
    print(f"  ✓ Saved: photon_3d_example_lateral_distortion.png")

    # ========================================================================
    # 5. ANIMATION (XY slice)
    # ========================================================================
    print(f"\nCreating animation...")
    print(f"  Total frames: {len(animation_frames)}")

    fig_anim, ax_anim = plt.subplots(figsize=(12, 4))

    # Initial frame (convert sim → nm)
    field_init_sim = extract_slice_xy(animation_frames[0], (nx, ny, nz))
    field_init_nm = mapper.to_phys_length(field_init_sim) * 1e9  # Convert sim → phys → nm
    im_anim = ax_anim.imshow(field_init_nm.T, origin='lower',
                             extent=[x_coords[0], x_coords[-1],
                                    y_coords[0], y_coords[-1]],
                             cmap='RdBu_r', vmin=-amplitude_nm*1.2, vmax=amplitude_nm*1.2,
                             aspect='auto', animated=True)

    ax_anim.set_xlabel('x [nm]', fontsize=12)
    ax_anim.set_ylabel('y [nm]', fontsize=12)
    time_text = ax_anim.text(0.02, 0.95, '', transform=ax_anim.transAxes,
                            fontsize=12, verticalalignment='top',
                            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

    divider = make_axes_locatable(ax_anim)
    cax = divider.append_axes("right", size="3%", pad=0.05)
    plt.colorbar(im_anim, cax=cax, label='ξ [nm]')
    ax_anim.set_title('3D Photon in Waveguide (XY slice, z=middle) - Dimensionless Units',
                     fontsize=14, fontweight='bold')

    def animate(frame_idx):
        """Update function for animation."""
        field_sim = extract_slice_xy(animation_frames[frame_idx], (nx, ny, nz))
        field_nm = mapper.to_phys_length(field_sim) * 1e9  # Convert sim → phys → nm
        im_anim.set_array(field_nm.T)
        t_sim = animation_times[frame_idx]
        t_fs = mapper.to_phys_time(t_sim) * 1e15  # Convert sim → phys → fs
        time_text.set_text(f't = {t_fs:.3f} fs')
        return [im_anim, time_text]

    anim = FuncAnimation(fig_anim, animate, frames=len(animation_frames),
                        interval=50, blit=True, repeat=True)

    # Save animation
    writer = FFMpegWriter(fps=20, bitrate=2000)
    anim.save('photon_3d_example.mp4', writer=writer, dpi=100)
    print(f"  ✓ Saved: photon_3d_example.mp4")

    plt.close(fig_anim)

    print(f"\n{'=' * 70}")
    print("Simulation complete!")
    print(f"{'=' * 70}")
    print(f"\nPhysical Interpretation:")
    print(f"  Domain size: {domain_length_phys_x*1e9:.3f} × {domain_length_phys_y*1e9:.3f} × {domain_length_phys_z*1e9:.3f} nm")
    print(f"  Domain size: {domain_length_phys_x/constants.lambda_C:.0f} × {domain_length_phys_y/constants.lambda_C:.0f} × {domain_length_phys_z/constants.lambda_C:.0f} λ_C")
    print(f"  Wavelength: {wavelength_phys*1e9:.3f} nm")
    print(f"  Simulation time: {simulation_time_phys*1e15:.3f} femtoseconds")

    print(f"\nVisualization Notes:")
    print(f"  • Primary plots show XY slice at z = {nz//2} (middle of waveguide)")
    print(f"  • Orthogonal slices show XZ and YZ planes")
    print(f"  • Lateral distortion uses color to encode direction (hue) and brightness for magnitude")
    print(f"  • Wave propagates along x-axis in waveguide mode")


if __name__ == '__main__':
    main()