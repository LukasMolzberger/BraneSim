"""
2D Photon in Tunnel Geometry with Realistic Physical Scales

Uses actual speed of light c = 299,792,458 m/s and physical length scales
based on the Compton wavelength.
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

from branesim.core.state import BraneState, Dimensionality
from branesim.core.grid import BraneGrid
from branesim.core.solver import VelocityVerletSolver
from branesim.physics.forces import SpringForceComputer
from branesim.config.simulation_config import PhysicalConstants
from branesim.physics.dimensional_mapping import DimensionalMapper
from branesim.core.initial_conditions import (
    initialize_right_moving_velocities_time_reversed,
    verify_wave_propagation,
)


def displacement_to_rgb(disp_x, disp_y, max_magnitude=None):
    """
    Convert 2D displacement vectors to RGB image using HSV color coding.

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


def initialize_tunnel_wave_shape_2d(state, grid, amplitude, center_x, width_x, wavelength, width_y=None):
    """
    Initialize ONLY the shape of a 2D localized photon wave packet.

    Uses Gaussian envelopes in both x and y directions to create a single
    localized wave packet centered in the domain. This represents a free-space
    photon propagating in 2D.

    Velocities are initialized separately using initialize_right_moving_velocities().
    This allows dimension-independent velocity initialization.

    Args:
        state: BraneState
        grid: BraneGrid
        amplitude: Peak amplitude [m]
        center_x: Center x coordinate [m]
        width_x: Gaussian width σ_x [m]
        wavelength: Carrier wavelength [m]
        width_y: Gaussian width σ_y [m] (if None, uses same as width_x)
    """
    # Get spatial coordinates
    coords = grid.get_spatial_coordinates()
    x = coords[:, 0]
    y = coords[:, 1]

    # Domain dimensions
    nx, ny = grid.grid_shape
    domain_length_x = (nx - 1) * grid.spacing
    domain_length_y = (ny - 1) * grid.spacing

    # Center the wave packet in y-direction
    center_y = domain_length_y / 2.0

    # Use same width in y as x if not specified
    if width_y is None:
        width_y = width_x

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

    # Full 2D Gaussian envelope
    envelope = amplitude * envelope_x * envelope_y

    # Wave number
    k = 2 * np.pi / wavelength

    # Position field only - velocities set separately
    # Carrier wave propagating in x-direction
    state.positions[:, 3] = envelope * torch.cos(k * (x - center_x))

    print(f"  Wavelength λ = {wavelength:.6e} m ({wavelength/grid.spacing:.1f} × h)")
    print(f"  Width σ_x = {width_x:.6e} m ({width_x/wavelength:.2f} × λ)")
    print(f"  Width σ_y = {width_y:.6e} m ({width_y/wavelength:.2f} × λ)")
    print(f"  Center (x, y) = ({center_x:.6e}, {center_y:.6e}) m")
    print(f"  Amplitude = {amplitude:.6e} m")
    print(f"  Wave number k = {k:.6e} rad/m")
    print(f"  Gaussian envelope: single localized wave packet")
    print(f"  Free-space propagation (no waveguide confinement)")


def main():
    """Run 2D photon simulation in tunnel geometry with realistic scales."""
    print("=" * 70)
    print("2D Photon in Tunnel - Realistic Physical Scales")
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

    D = 2

    # 2D brane parameters constrained to give wave speed = c
    # Use same value as 1D for consistency (arbitrary choice)
    m_e = constants.m_e
    rho_D = m_e / (constants.lambda_C ** 2)  # kg/m² (surface mass density - derived from Compton scale)
    T_D = rho_D * constants.c**2  # N/m (tension - computed from c² = T_D/rho_D)
    rest_length_phys = 0.0 * h_phys

    # Wave speed (exactly equals c by construction)
    c_wave = constants.c

    # Discrete mass per lattice point
    m_point = rho_D * (h_phys ** D)

    # Axial spring constant (for 2D: k = T_D * h^(D-2) = T_D * h^0 = T_D)
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

    # Domain size - tunnel geometry (long in x, narrow in y)
    nx = 800  # Long tunnel
    ny = 200  # Narrow tunnel
    domain_length_phys = nx * h_phys
    domain_length_sim = nx * h_sim  # = nx * 1.0 = nx

    # Verify wave speed
    expected_wave_speed = math.sqrt(T_D / rho_D)

    print(f"\nPhysical Parameters:")
    print(f"  2D surface mass density ρ_2 = {rho_D:.6e} kg/m²")
    print(f"  2D tension T_2 = {T_D:.6e} N/m")
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
    print(f"  Domain: {nx} × {ny} points")
    print(f"  Domain (physical): {domain_length_phys:.6e} × {domain_length_phys * ny / nx:.6e} m")
    print(f"  Domain (sim units): {domain_length_sim:.1f} × {domain_length_sim * ny / nx:.1f}")
    print(f"  Aspect ratio: {nx/ny:.1f}:1 (tunnel geometry)")
    print(f"  CFL number = {cfl_factor:.3f}")

    # Create components using SIMULATION UNITS
    device = torch.device('cpu')
    dtype = torch.float64

    state = BraneState((nx, ny), Dimensionality.TWO_D, device, dtype)
    state.initialize_flat_configuration(h_sim)  # Use sim spacing = 1.0

    # Store initial positions for lateral distortion tracking
    initial_positions = state.positions.clone()

    # Set fixed boundaries (all edges)
    state.set_fixed_boundaries()
    print(f"\nBoundary Conditions:")
    print(f"  Fixed points: {state.fixed_mask.sum().item()} / {nx * ny}")
    print(f"  Top/bottom walls: FIXED (tunnel)")
    print(f"  Left/right walls: FIXED (end caps)")

    grid = BraneGrid((nx, ny), Dimensionality.TWO_D, h_sim, device)  # Sim spacing = 1.0

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
    wavelength_phys = 40 * h_phys  # 40 points per wavelength
    amplitude_phys = 10 * h_phys
    width_x_phys = 3 * wavelength_phys / (2 * np.pi)
    center_x_phys = domain_length_phys / 3.0

    # Convert to sim units using mapper
    wavelength_sim = mapper.to_sim_length(wavelength_phys)  # = 40.0
    amplitude_sim = mapper.to_sim_length(amplitude_phys)    # = 0.2
    width_x_sim = mapper.to_sim_length(width_x_phys)
    center_x_sim = mapper.to_sim_length(center_x_phys)

    print(f"  Physical wavelength: {wavelength_phys:.6e} m")
    print(f"  Sim wavelength: {wavelength_sim:.1f} grid units")
    print(f"  Physical amplitude: {amplitude_phys:.6e} m")
    print(f"  Sim amplitude: {amplitude_sim:.3f} grid units")

    # Step 1: Initialize shape only (in sim units)
    print(f"\n[1] Initializing wave shape (sim units)...")
    initialize_tunnel_wave_shape_2d(state, grid, amplitude_sim, center_x_sim, width_x_sim, wavelength_sim)

    # Step 2: Initialize velocities for right-moving wave at speed c
    print(f"\n[2] Initializing velocities for right-moving wave...")
    initialize_right_moving_velocities_time_reversed(
        state=state,
        grid=grid,
        physics=physics,
        m_point=m_sim,
        wave_speed=c_sim,  # Actual wave speed in sim units (= √(k_sim/m_sim) = 1.0)
        field_component=3,
        shift_cells=1,
    )

    # Step 3: Compute initial accelerations
    solver.initialize_accelerations(state)
    state.apply_fixed_boundaries()

    # Initial measurements
    initial_energy = solver.compute_energy(state)

    print(f"\nInitial State:")
    print(f"  Energy = {initial_energy['total']:.6e} J")

    # Run simulation
    # Time for light to cross domain: t = L/c (in physical units)
    crossing_time_phys = domain_length_phys / constants.c
    simulation_time_phys = 3.0 * crossing_time_phys  # 3 crossings
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
    snapshots_lateral_x = {}  # Store x-component of lateral displacement
    snapshots_lateral_y = {}  # Store y-component of lateral displacement
    snapshot_steps = {int(t / dt_phys): t for t in snapshot_times_phys}

    # For animation - save every N steps
    animation_frames = []
    animation_frames_lateral_x = []  # For lateral distortion animation
    animation_frames_lateral_y = []  # For lateral distortion animation
    animation_times = []
    frame_interval = max(1, num_steps // 300)  # ~300 frames total

    # Physical coordinates for plotting (convert from sim to physical)
    x_coords_phys = np.arange(nx) * h_phys
    y_coords_phys = np.arange(ny) * h_phys

    print_interval = max(1, num_steps // 20)

    for step in range(num_steps + 1):
        if step in snapshot_steps:
            field = state.get_field_component(3).cpu().numpy()
            snapshots[snapshot_steps[step]] = field.copy()

            # Store lateral displacements (x and y components)
            lateral_disp_x = (state.positions[:, 0] - initial_positions[:, 0]).cpu().numpy()
            lateral_disp_y = (state.positions[:, 1] - initial_positions[:, 1]).cpu().numpy()
            snapshots_lateral_x[snapshot_steps[step]] = lateral_disp_x.copy()
            snapshots_lateral_y[snapshot_steps[step]] = lateral_disp_y.copy()

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
    print(f"  Tunnel length: {domain_length_phys:.6e} m")
    print(f"  Tunnel height: {domain_length_phys * ny / nx:.6e} m")

    print(f"\nEnergy Conservation:")
    print(f"  Initial: {initial_energy['total']:.6e} J")
    print(f"  Final:   {final_energy['total']:.6e} J")
    print(f"  Drift:   {energy_drift:.6e} ({energy_drift*100:.6f}%)")

    if energy_drift < 1e-2:
        print(f"  ✓ Good energy conservation")
    else:
        print(f"  ⚠ Energy drift: {energy_drift*100:.2f}%")

    # Visualization
    print(f"\nCreating plots...")

    fig, axes = plt.subplots(num_snapshots, 1, figsize=(16, 12))
    fig.suptitle(f'2D Photon in Tunnel (c = {constants.c:.3e} m/s) - Dimensionless Units',
                 fontsize=16, fontweight='bold')

    # Convert to nanometers for better readability
    x_nm = x_coords_phys * 1e9
    y_nm = y_coords_phys * 1e9
    amplitude_nm = amplitude_phys * 1e9

    for idx, (step, t) in enumerate(snapshot_steps.items()):
        if t in snapshots:
            field_sim = snapshots[t].reshape(nx, ny)
            field_nm = mapper.to_phys_length(field_sim) * 1e9  # Convert sim → phys → nm

            # Plot heatmap
            im = axes[idx].imshow(field_nm.T, origin='lower',
                                 extent=[x_nm[0], x_nm[-1], y_nm[0], y_nm[-1]],
                                 cmap='RdBu_r', vmin=-amplitude_nm*1.2, vmax=amplitude_nm*1.2,
                                 aspect='auto')

            axes[idx].set_ylabel('y [nm]', fontsize=11)
            axes[idx].set_xlim(x_nm[0], x_nm[-1])
            axes[idx].set_ylim(y_nm[0], y_nm[-1])

            # Time in femtoseconds
            t_fs = t * 1e15
            axes[idx].text(0.02, 0.95, f't = {t_fs:.3f} fs',
                          transform=axes[idx].transAxes,
                          fontsize=12, verticalalignment='top',
                          bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

            # Add colorbar
            plt.colorbar(im, ax=axes[idx], label='ξ [nm]', fraction=0.046, pad=0.04)

            if idx == num_snapshots - 1:
                axes[idx].set_xlabel('x [nm]', fontsize=12)

    plt.tight_layout()
    plt.savefig('photon_2d_example_propagation.png', dpi=150, bbox_inches='tight')
    print(f"  ✓ Saved: photon_2d_example_propagation.png")

    # Energy conservation plot
    fig2, ax = plt.subplots(figsize=(10, 6))

    times_fs = np.array(times_phys) * 1e15
    energy_array = np.array(energies)
    ax.plot(times_fs, energy_array / initial_energy['total'], 'g-', linewidth=2)
    ax.axhline(y=1.0, color='r', linestyle='--', linewidth=1, alpha=0.5)
    ax.set_xlabel('Time [fs]', fontsize=12)
    ax.set_ylabel('E(t) / E(0)', fontsize=12)
    ax.set_title('Energy Conservation', fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig('photon_2d_example_energy.png', dpi=150, bbox_inches='tight')
    print(f"  ✓ Saved: photon_2d_example_energy.png")

    # Create animation
    print(f"\nCreating animation...")
    print(f"  Total frames: {len(animation_frames)}")

    fig_anim, ax_anim = plt.subplots(figsize=(12, 3))

    # Initial frame (convert sim → nm)
    field_init = mapper.to_phys_length(animation_frames[0].reshape(nx, ny)) * 1e9
    im_anim = ax_anim.imshow(field_init.T, origin='lower',
                             extent=[x_nm[0], x_nm[-1], y_nm[0], y_nm[-1]],
                             cmap='RdBu_r', vmin=-amplitude_nm*1.2, vmax=amplitude_nm*1.2,
                             aspect='auto', animated=True)

    ax_anim.set_xlabel('x [nm]', fontsize=12)
    ax_anim.set_ylabel('y [nm]', fontsize=12)
    time_text = ax_anim.text(0.02, 0.95, '', transform=ax_anim.transAxes,
                            fontsize=12, verticalalignment='top',
                            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

    plt.colorbar(im_anim, ax=ax_anim, label='ξ [nm]', fraction=0.046, pad=0.04)
    ax_anim.set_title('2D Photon in Tunnel (c = 299,792,458 m/s) - Dimensionless Units', fontsize=14, fontweight='bold')

    def animate(frame_idx):
        """Update function for animation."""
        field_sim = animation_frames[frame_idx].reshape(nx, ny)
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
    anim.save('photon_2d_example.mp4', writer=writer, dpi=100)
    print(f"  ✓ Saved: photon_2d_example.mp4")

    plt.close(fig_anim)

    # Lateral distortion visualization
    print(f"\nCreating lateral distortion plots...")

    # Find max lateral displacement magnitude for consistent scaling (convert sim → phys)
    max_disp_mag_phys = 0
    for t in snapshots_lateral_x.keys():
        disp_x_sim = snapshots_lateral_x[t].reshape(nx, ny)
        disp_y_sim = snapshots_lateral_y[t].reshape(nx, ny)
        disp_x_phys = mapper.to_phys_length(disp_x_sim)
        disp_y_phys = mapper.to_phys_length(disp_y_sim)
        mag = np.sqrt(disp_x_phys**2 + disp_y_phys**2).max()
        max_disp_mag_phys = max(max_disp_mag_phys, mag)

    # Create snapshot plots for lateral distortion
    fig_lat, axes_lat = plt.subplots(num_snapshots, 1, figsize=(16, 12))
    fig_lat.suptitle(f'2D Photon - Lateral Distortion (Color = Direction, Brightness = Magnitude)',
                     fontsize=16, fontweight='bold')

    for idx, (step, t) in enumerate(snapshot_steps.items()):
        if t in snapshots_lateral_x:
            disp_x_sim = snapshots_lateral_x[t].reshape(nx, ny)
            disp_y_sim = snapshots_lateral_y[t].reshape(nx, ny)

            # Convert sim → phys
            disp_x_phys = mapper.to_phys_length(disp_x_sim)
            disp_y_phys = mapper.to_phys_length(disp_y_sim)

            # Convert to RGB image (using physical units)
            rgb_image, magnitude, angle = displacement_to_rgb(disp_x_phys, disp_y_phys, max_disp_mag_phys)

            # Plot RGB image (transpose spatial dimensions only, keep color channel last)
            axes_lat[idx].imshow(np.transpose(rgb_image, (1, 0, 2)), origin='lower',
                                extent=[x_nm[0], x_nm[-1], y_nm[0], y_nm[-1]],
                                aspect='auto')

            axes_lat[idx].set_ylabel('y [nm]', fontsize=11)
            axes_lat[idx].set_xlim(x_nm[0], x_nm[-1])
            axes_lat[idx].set_ylim(y_nm[0], y_nm[-1])

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
    plt.savefig('photon_2d_example_lateral_distortion.png', dpi=150, bbox_inches='tight')
    print(f"  ✓ Saved: photon_2d_example_lateral_distortion.png")

    # Create lateral distortion animation
    print(f"\nCreating lateral distortion animation...")

    # Find max magnitude across all frames for consistent scaling (convert sim → phys)
    max_disp_mag_anim_phys = 0
    for idx in range(len(animation_frames_lateral_x)):
        disp_x_sim = animation_frames_lateral_x[idx].reshape(nx, ny)
        disp_y_sim = animation_frames_lateral_y[idx].reshape(nx, ny)
        disp_x_phys = mapper.to_phys_length(disp_x_sim)
        disp_y_phys = mapper.to_phys_length(disp_y_sim)
        mag = np.sqrt(disp_x_phys**2 + disp_y_phys**2).max()
        max_disp_mag_anim_phys = max(max_disp_mag_anim_phys, mag)

    fig_anim_lat, ax_anim_lat = plt.subplots(figsize=(12, 3))

    # Initial frame (convert sim → phys)
    disp_x_init_sim = animation_frames_lateral_x[0].reshape(nx, ny)
    disp_y_init_sim = animation_frames_lateral_y[0].reshape(nx, ny)
    disp_x_init_phys = mapper.to_phys_length(disp_x_init_sim)
    disp_y_init_phys = mapper.to_phys_length(disp_y_init_sim)
    rgb_init, _, _ = displacement_to_rgb(disp_x_init_phys, disp_y_init_phys, max_disp_mag_anim_phys)

    # Transpose spatial dimensions only, keep color channel last
    im_anim_lat = ax_anim_lat.imshow(np.transpose(rgb_init, (1, 0, 2)), origin='lower',
                                    extent=[x_nm[0], x_nm[-1], y_nm[0], y_nm[-1]],
                                    aspect='auto', animated=True)

    ax_anim_lat.set_xlabel('x [nm]', fontsize=12)
    ax_anim_lat.set_ylabel('y [nm]', fontsize=12)
    time_text_lat = ax_anim_lat.text(0.02, 0.95, '', transform=ax_anim_lat.transAxes,
                                    fontsize=12, verticalalignment='top',
                                    bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

    ax_anim_lat.set_title('2D Photon - Lateral Distortion (Color = Direction, Brightness = Magnitude)',
                         fontsize=14, fontweight='bold')

    def animate_lateral(frame_idx):
        """Update function for lateral distortion animation."""
        disp_x_sim = animation_frames_lateral_x[frame_idx].reshape(nx, ny)
        disp_y_sim = animation_frames_lateral_y[frame_idx].reshape(nx, ny)
        disp_x_phys = mapper.to_phys_length(disp_x_sim)
        disp_y_phys = mapper.to_phys_length(disp_y_sim)
        rgb_image, _, _ = displacement_to_rgb(disp_x_phys, disp_y_phys, max_disp_mag_anim_phys)

        # Transpose spatial dimensions only, keep color channel last
        im_anim_lat.set_array(np.transpose(rgb_image, (1, 0, 2)))
        t_sim = animation_times[frame_idx]
        t_fs = mapper.to_phys_time(t_sim) * 1e15  # Convert sim → phys → fs
        time_text_lat.set_text(f't = {t_fs:.3f} fs')
        return [im_anim_lat, time_text_lat]

    anim_lat = FuncAnimation(fig_anim_lat, animate_lateral, frames=len(animation_frames_lateral_x),
                            interval=50, blit=True, repeat=True)

    # Save animation
    writer_lat = FFMpegWriter(fps=20, bitrate=2000)
    anim_lat.save('photon_2d_example_lateral_distortion.mp4', writer=writer_lat, dpi=100)
    print(f"  ✓ Saved: photon_2d_example_lateral_distortion.mp4")

    plt.close(fig_anim_lat)

    print(f"\n{'=' * 70}")
    print("Simulation complete!")
    print(f"{'=' * 70}")
    print(f"\nPhysical Interpretation:")
    print(f"  Domain size: {domain_length_phys*1e9:.3f} × {domain_length_phys * ny / nx *1e9:.3f} nm")
    print(f"  Domain size: {domain_length_phys/constants.lambda_C:.0f} × {domain_length_phys * ny / nx / constants.lambda_C:.0f} λ_C")
    print(f"  Wavelength: {wavelength_phys*1e9:.3f} nm")
    print(f"  Simulation time: {simulation_time_phys*1e15:.3f} femtoseconds")


if __name__ == '__main__':
    main()