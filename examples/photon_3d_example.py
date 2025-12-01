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
from branesim.physics.parameters import compton_calibrated_brane_lattice_params


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

    # Configuration with Compton-cell calibration
    # Grid spacing as multiple of Compton wavelength
    lambda_C_multiplier = 10.0  # Grid spacing = 10 × λ_C
    h = constants.lambda_C * lambda_C_multiplier

    # Get 3D Compton-calibrated parameters
    params = compton_calibrated_brane_lattice_params(
        grid_spacing_m=h,
        dimensionality=3,
        c=constants.c
    )

    # Extract 3D parameters
    rho = params["rho_D"]  # 3D volume mass density [kg/m³]
    tension = params["T_D"]  # 3D "tension" [Pa = N/m²]

    # Domain size - waveguide geometry (long in x, narrow in y and z)
    # Note: Reduced size for faster computation (3D is computationally intensive)
    nx = 100  # Long waveguide
    ny = 20   # Narrow transverse
    nz = 20   # Narrow transverse
    domain_length_x = nx * h
    domain_length_y = ny * h
    domain_length_z = nz * h

    # Wave speed (speed of light)
    c = constants.c

    # CFL condition for stability
    cfl_factor = 0.1
    dt = cfl_factor * h / c

    # Verify that wave speed will be c
    # For 3D: c = √(T/ρ) where T is the 3D "tension" (actually pressure/stress)
    expected_wave_speed = np.sqrt(tension / rho)

    print(f"\nCompton-Cell Calibration (3D):")
    print(f"  Reduced Compton wavelength λ_C = {params['lambda_C']:.4e} m")
    print(f"  Grid spacing h = {lambda_C_multiplier:.0f} × λ_C = {h:.6e} m")
    print(f"  3D volume mass density ρ_3 = m_e/λ_C³ = {rho:.6e} kg/m³")
    print(f"  3D elastic modulus T_3 = ρ_3×c² = {tension:.6e} Pa")
    print(f"  Spring constant k = T_3×h = {params['k_spring']:.6e} N/m")
    print(f"  Point mass m = ρ_3×h³ = {params['m_point']:.6e} kg")
    print(f"  Expected wave speed = √(T_3/ρ_3) = {expected_wave_speed:.6e} m/s")
    print(f"  Speed of light c = {c:.6e} m/s")
    print(f"  Wave speed error = {abs(expected_wave_speed - c)/c:.6e}")

    print(f"\nSimulation Configuration:")
    print(f"  Domain: {nx} × {ny} × {nz} points")
    print(f"  Domain size: {domain_length_x:.6e} × {domain_length_y:.6e} × {domain_length_z:.6e} m")
    print(f"  Total points: {nx*ny*nz:,}")
    print(f"  Aspect ratio: {nx}:{ny}:{nz} (waveguide geometry)")
    print(f"  Time step dt = {dt:.6e} s")
    print(f"  CFL number = {cfl_factor:.3f}")

    # Create components
    device = torch.device('cpu')
    dtype = torch.float64

    state = BraneState((nx, ny, nz), Dimensionality.THREE_D, device, dtype)
    state.initialize_flat_configuration(h)

    # Store initial positions for lateral distortion tracking
    initial_positions = state.positions.clone()

    # Set fixed boundaries (all 6 walls)
    state.set_fixed_boundaries()
    print(f"\nBoundary Conditions:")
    print(f"  Fixed points: {state.fixed_mask.sum().item()} / {nx * ny * nz}")
    print(f"  All 6 walls: FIXED (waveguide)")

    grid = BraneGrid((nx, ny, nz), Dimensionality.THREE_D, h, device)

    # CRITICAL: Implement pretension κ = ρc²
    # In continuum: elastic modulus T_3 = ρ_3 * c² [Pa]
    # In discrete: each spring carries force F_0 = T_3 * h² [N]
    # With k_spring = T_3 * h: F_0 = k(h - L_0) must equal T_3 * h²
    # Solving: T_3*h(h - L_0) = T_3*h² → L_0 = 0
    rest_length = 0.0  # Springs pre-stretched to carry background tension
    spring_constant = params["k_spring"]  # k = T_3 * h [N/m]
    background_force = spring_constant * (h - rest_length)  # Force per spring [N]

    print(f"\nPretension Implementation:")
    print(f"  Rest length L_0 = {rest_length:.6e} m")
    print(f"  Actual spacing a = {h:.6e} m")
    print(f"  Background strain = (a - L_0)/L_0 = {'infinite (L_0=0)' if rest_length == 0 else f'{(h-rest_length)/rest_length:.6e}'}")
    print(f"  Background force per spring F_0 = k(a - L_0) = {background_force:.6e} N")
    print(f"  Expected force per spring = T_3*h² = {tension*h**2:.6e} N")
    print(f"  Match: {abs(background_force - tension*h**2)/(tension*h**2):.6e}")

    physics = SpringForceComputer(spring_constant, rest_length)
    solver = VelocityVerletSolver(dt, rho, physics, grid)

    # Initialize wave packet
    print(f"\nInitializing photon wave packet...")

    # Wavelength: multiple of grid spacing for good resolution
    wavelength = 20 * h  # 20 points per wavelength

    # Amplitude: small compared to domain
    amplitude = 0.1 * h

    # Width: several wavelengths (only in x)
    width_x = 3 * wavelength / (2 * np.pi)

    # Center: in the left third of domain
    center_x = domain_length_x / 3.0

    # Initialize shape only
    print(f"\n[1] Initializing wave shape...")
    initialize_waveguide_wave_shape_3d(state, grid, amplitude, center_x, width_x, wavelength)

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
    # Time for light to cross domain: t = L/c
    crossing_time = domain_length_x / c
    simulation_time = 0.5 * crossing_time  # 0.5 crossings (reduced for faster demo)

    num_steps = int(simulation_time / dt)

    print(f"\nRunning simulation...")
    print(f"  Light crossing time = {crossing_time:.6e} s")
    print(f"  Simulation time = {simulation_time:.6e} s")
    print(f"  Number of steps = {num_steps:,}")

    # Tracking
    times = []
    energies = []

    # Take snapshots at regular intervals
    num_snapshots = 7
    snapshot_times = np.linspace(0, simulation_time, num_snapshots)
    snapshots = {}
    snapshots_lateral_x = {}  # x-component of lateral displacement
    snapshots_lateral_y = {}  # y-component of lateral displacement
    snapshots_lateral_z = {}  # z-component of lateral displacement
    snapshot_steps = {int(t / dt): t for t in snapshot_times}

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
            times.append(solver.time)
            energies.append(energy['total'])

        if step % print_interval == 0:
            print(f"  Step {step:8d}/{num_steps}: t={solver.time:.6e}s, "
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
    print(f"  Photon propagates at c = {c:.6e} m/s")
    print(f"  Waveguide dimensions: {domain_length_x:.6e} × {domain_length_y:.6e} × {domain_length_z:.6e} m")

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

    # Coordinate arrays (in nanometers)
    x_coords = np.arange(nx) * h * 1e9
    y_coords = np.arange(ny) * h * 1e9
    z_coords = np.arange(nz) * h * 1e9
    amplitude_nm = amplitude * 1e9

    # ========================================================================
    # 1. AMPLITUDE FIELD VISUALIZATION (XY slice at middle z)
    # ========================================================================
    fig, axes = plt.subplots(num_snapshots, 1, figsize=(14, 12))
    fig.suptitle(f'3D Photon in Waveguide - Amplitude Field (XY slice, z={nz//2})',
                 fontsize=16, fontweight='bold')

    for idx, (step, t) in enumerate(snapshot_steps.items()):
        if t in snapshots:
            # Extract XY slice at middle z
            field_slice = extract_slice_xy(snapshots[t], (nx, ny, nz))
            field_slice_nm = field_slice * 1e9

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
            # XZ slice at middle y
            field_xz = extract_slice_xz(snapshots[t], (nx, ny, nz))
            field_xz_nm = field_xz * 1e9

            # YZ slice at middle x
            field_yz = extract_slice_yz(snapshots[t], (nx, ny, nz))
            field_yz_nm = field_yz * 1e9

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

    times_fs = np.array(times) * 1e15
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

    # Find max lateral displacement magnitude for consistent scaling
    max_disp_mag = 0
    for t in snapshots_lateral_x.keys():
        # Extract XY slice at middle z
        disp_x_slice = extract_slice_xy(snapshots_lateral_x[t], (nx, ny, nz))
        disp_y_slice = extract_slice_xy(snapshots_lateral_y[t], (nx, ny, nz))
        mag = np.sqrt(disp_x_slice**2 + disp_y_slice**2).max()
        max_disp_mag = max(max_disp_mag, mag)

    fig_lat, axes_lat = plt.subplots(num_snapshots, 1, figsize=(14, 12))
    fig_lat.suptitle(f'3D Photon - Lateral Distortion (XY slice, Color=Direction, Brightness=Magnitude)',
                     fontsize=16, fontweight='bold')

    for idx, (step, t) in enumerate(snapshot_steps.items()):
        if t in snapshots_lateral_x:
            # Extract XY slice at middle z
            disp_x_slice = extract_slice_xy(snapshots_lateral_x[t], (nx, ny, nz))
            disp_y_slice = extract_slice_xy(snapshots_lateral_y[t], (nx, ny, nz))

            # Convert to RGB image
            rgb_image, magnitude, angle = displacement_to_rgb_3d(disp_x_slice, disp_y_slice, max_disp_mag)

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

            # Add magnitude info
            max_mag_pm = max_disp_mag * 1e12
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

    # Initial frame
    field_init = extract_slice_xy(animation_frames[0], (nx, ny, nz)) * 1e9
    im_anim = ax_anim.imshow(field_init.T, origin='lower',
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
    ax_anim.set_title('3D Photon in Waveguide (XY slice, z=middle)',
                     fontsize=14, fontweight='bold')

    def animate(frame_idx):
        """Update function for animation."""
        field = extract_slice_xy(animation_frames[frame_idx], (nx, ny, nz)) * 1e9
        im_anim.set_array(field.T)
        t_fs = animation_times[frame_idx] * 1e15
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
    print(f"  Domain size: {domain_length_x*1e9:.3f} × {domain_length_y*1e9:.3f} × {domain_length_z*1e9:.3f} nm")
    print(f"  Domain size: {domain_length_x/constants.lambda_C:.0f} × {domain_length_y/constants.lambda_C:.0f} × {domain_length_z/constants.lambda_C:.0f} λ_C")
    print(f"  Wavelength: {wavelength*1e9:.3f} nm")
    print(f"  Simulation time: {simulation_time*1e15:.3f} femtoseconds")

    print(f"\nVisualization Notes:")
    print(f"  • Primary plots show XY slice at z = {nz//2} (middle of waveguide)")
    print(f"  • Orthogonal slices show XZ and YZ planes")
    print(f"  • Lateral distortion uses color to encode direction (hue) and brightness for magnitude")
    print(f"  • Wave propagates along x-axis in waveguide mode")


if __name__ == '__main__':
    main()