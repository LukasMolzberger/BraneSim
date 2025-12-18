"""
Polarized Photon Experiment - EM Four-Potential Solver

This experiment simulates a circularly polarized photon in a straight waveguide
using the four-potential electromagnetic field solver.

Key features:
- Uses A^μ = (Φ/c, Ax, Ay, Az) formulation
- Solves d'Alembertian wave equation in vacuum
- Circularly polarized tubular mode with transverse Gaussian profile
- Computes E and B fields from A^μ for analysis and visualization

The goal is to demonstrate:
1. Pure EM wave propagation in vacuum (no brane mapping)
2. Circular polarization encoded in rotating vector potential
3. Energy conservation in the EM solver
4. Proper |E| = c|B| relationship for electromagnetic waves
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

from branesim.core.state import Dimensionality
from branesim.core.grid import BraneGrid
from branesim.em.em_state import EMState
from branesim.em.potential_solver import FourPotentialVerletSolver
from branesim.config.physical_constants import PhysicalConstants
from branesim.utils import TestRunManager


def initialize_circular_polarized_waveguide(
    state: EMState,
    grid: BraneGrid,
    wavelength: float,
    sigma_transverse: float,
    amplitude: float,
    c: float,
    propagation_axis: int = 0,
) -> None:
    """
    Initialize four-potential for a circularly polarized photon in a waveguide.

    For a wave propagating along x-axis (propagation_axis=0):
        A_y(x,t=0) = A_0 * exp(-r_⊥²/(2σ²)) * sin(kx)
        A_z(x,t=0) = A_0 * exp(-r_⊥²/(2σ²)) * cos(kx)
        A^0 = 0 (no scalar potential)

        ∂A_y/∂t(t=0) = -ω * A_0 * exp(-r_⊥²/(2σ²)) * cos(kx)
        ∂A_z/∂t(t=0) = ω * A_0 * exp(-r_⊥²/(2σ²)) * sin(kx)

    This gives circularly polarized E and B fields.

    Args:
        state: EMState to initialize
        grid: BraneGrid with spatial coordinates
        wavelength: Wavelength λ in meters
        sigma_transverse: Gaussian width in transverse directions [m]
        amplitude: Peak vector potential magnitude A_0 [V⋅s/m]
        c: Speed of light [m/s]
        propagation_axis: Propagation direction (0=x, 1=y, 2=z)
    """
    k = 2.0 * math.pi / wavelength
    omega = c * k

    # Get spatial coordinates [N, D]
    coords = grid.get_spatial_coordinates().cpu().numpy()
    N = coords.shape[0]

    # Determine transverse coordinates based on propagation axis
    if propagation_axis == 0:  # Propagate along x
        z_long = coords[:, 0]  # longitudinal coordinate
        transverse_coords = coords[:, 1:]  # (y, z) if 3D, or (y,) if 2D
        pol_indices = [2, 3]  # A_y = index 2, A_z = index 3
    elif propagation_axis == 1:  # Propagate along y
        z_long = coords[:, 1]
        transverse_coords = np.column_stack([coords[:, 0], coords[:, 2]]) if coords.shape[1] > 2 else coords[:, [0]]
        pol_indices = [3, 1]  # A_z, A_x
    else:  # Propagate along z
        z_long = coords[:, 2]
        transverse_coords = coords[:, :2]
        pol_indices = [1, 2]  # A_x, A_y

    # Center Gaussian at middle of transverse domain
    if transverse_coords.ndim == 1:
        transverse_coords = transverse_coords.reshape(-1, 1)

    transverse_center = transverse_coords.mean(axis=0)
    r_transverse = np.linalg.norm(transverse_coords - transverse_center, axis=1)

    # Gaussian transverse envelope
    envelope = np.exp(-0.5 * (r_transverse / sigma_transverse) ** 2)

    # Phase along propagation direction
    phase = k * z_long

    # Vector potential components (circular polarization)
    A_pol1 = amplitude * envelope * np.sin(phase)  # First polarization component
    A_pol2 = amplitude * envelope * np.cos(phase)  # Second (90° shifted)

    # Time derivatives
    A_pol1_dot = -omega * amplitude * envelope * np.cos(phase)
    A_pol2_dot = omega * amplitude * envelope * np.sin(phase)

    # Set state (zero out first, then set components)
    state.potential.zero_()
    state.velocity.zero_()

    # Convert to torch and assign
    device = state.device
    dtype = state.dtype

    state.potential[:, pol_indices[0]] = torch.from_numpy(A_pol1).to(device=device, dtype=dtype)
    state.potential[:, pol_indices[1]] = torch.from_numpy(A_pol2).to(device=device, dtype=dtype)

    state.velocity[:, pol_indices[0]] = torch.from_numpy(A_pol1_dot).to(device=device, dtype=dtype)
    state.velocity[:, pol_indices[1]] = torch.from_numpy(A_pol2_dot).to(device=device, dtype=dtype)

    state.apply_fixed_boundaries()


def vector_field_to_rgb(vx, vy, max_magnitude=None):
    """
    Convert 2D vector field to RGB using HSV color coding.

    Args:
        vx: x-component (2D array)
        vy: y-component (2D array)
        max_magnitude: Maximum magnitude for normalization

    Returns:
        RGB image where hue=direction, value=magnitude
    """
    magnitude = np.sqrt(vx**2 + vy**2)
    angle = np.arctan2(vy, vx)

    hue = (angle + np.pi) / (2 * np.pi)

    if max_magnitude is None:
        max_magnitude = magnitude.max()

    if max_magnitude > 0:
        value = magnitude / max_magnitude
    else:
        value = np.zeros_like(magnitude)

    saturation = np.ones_like(magnitude)
    hsv = np.stack([hue, saturation, value], axis=-1)
    rgb = hsv_to_rgb(hsv)
    rgb = np.clip(rgb, 0.0, 1.0)

    return rgb, magnitude, angle


def extract_slice_xy(field_3d, grid_shape, z_index=None):
    """Extract 2D XY slice from 3D field at constant z."""
    nx, ny, nz = grid_shape
    field = field_3d.reshape(nx, ny, nz)
    if z_index is None:
        z_index = nz // 2
    return field[:, :, z_index]


def extract_slice_xz(field_3d, grid_shape, y_index=None):
    """Extract 2D XZ slice from 3D field at constant y."""
    nx, ny, nz = grid_shape
    field = field_3d.reshape(nx, ny, nz)
    if y_index is None:
        y_index = ny // 2
    return field[:, y_index, :]


def extract_slice_yz(field_3d, grid_shape, x_index=None):
    """Extract 2D YZ slice from 3D field at constant x."""
    nx, ny, nz = grid_shape
    field = field_3d.reshape(nx, ny, nz)
    if x_index is None:
        x_index = nx // 2
    return field[x_index, :, :]


def compute_em_energy(E_field, B_field, epsilon0, mu0, grid_spacing, ndim):
    """
    Compute total electromagnetic energy.

    U_EM = ∫ (ε₀/2 |E|² + 1/(2μ₀) |B|²) dV

    Args:
        E_field: [N, 3] electric field
        B_field: [N, 3] magnetic field
        epsilon0: Permittivity
        mu0: Permeability
        grid_spacing: Spatial grid spacing
        ndim: Spatial dimension (for volume element)

    Returns:
        Total EM energy in Joules
    """
    E_energy_density = 0.5 * epsilon0 * torch.sum(E_field ** 2, dim=1)
    B_energy_density = 0.5 / mu0 * torch.sum(B_field ** 2, dim=1)
    total_energy_density = E_energy_density + B_energy_density

    # Volume element
    dV = grid_spacing ** ndim

    return torch.sum(total_energy_density) * dV


def main():
    """Run polarized photon experiment with EM four-potential solver."""
    print("=" * 70)
    print("Polarized Photon Experiment - EM Four-Potential Solver")
    print("=" * 70)

    # Initialize test run manager
    run_manager = TestRunManager(experiment_name="polarized_photon_em_solver")
    print(run_manager.get_summary())

    constants = PhysicalConstants()
    print(f"\nPhysical Constants:")
    print(f"  c = {constants.c:.6e} m/s")
    print(f"  λ_C = {constants.lambda_C:.6e} m")
    print(f"  ε₀ = {constants.epsilon0:.6e} F/m")
    print(f"  μ₀ = {constants.mu0:.6e} H/m")

    # --- SIMULATION SETUP ---
    # Use a practical wavelength (visible light range) instead of Compton wavelength
    # to avoid numerical underflow issues
    wavelength_phys = 500e-9  # 500 nm (green light)
    points_per_wavelength = 20
    h_phys = wavelength_phys / points_per_wavelength
    cfl_factor = 0.3  # Conservative CFL for stability

    # Domain size
    nx = 80
    ny = 80
    nz = 80
    domain_length_x = nx * h_phys
    domain_length_y = ny * h_phys
    domain_length_z = nz * h_phys

    print(f"\nSimulation Configuration:")
    print(f"  Domain: {nx} × {ny} × {nz} points")
    print(f"  Domain (physical): {domain_length_x*1e9:.3f} × {domain_length_y*1e9:.3f} × {domain_length_z*1e9:.3f} nm")
    print(f"  Grid spacing: {h_phys*1e12:.3f} pm")
    print(f"  Wavelength: {wavelength_phys*1e9:.3f} nm")
    print(f"  Points per wavelength: {points_per_wavelength}")
    print(f"  Total points: {nx*ny*nz:,}")

    # Time step from CFL condition
    dt_phys = cfl_factor * h_phys / constants.c

    # Auto-select device
    if torch.cuda.is_available():
        device = torch.device('cuda')
        print(f"\n✓ Using NVIDIA GPU: {torch.cuda.get_device_name(0)}")
        dtype = torch.float64
    elif hasattr(torch.backends, 'mps') and torch.backends.mps.is_available():
        device = torch.device('mps')
        print(f"\n✓ Using Apple Silicon GPU (MPS)")
        dtype = torch.float32
    else:
        device = torch.device('cpu')
        print(f"\n⚠ Using CPU (no GPU detected)")
        dtype = torch.float64

    # Create EM state and grid
    grid = BraneGrid((nx, ny, nz), Dimensionality.THREE_D, spacing=h_phys, device=device)
    state = EMState((nx, ny, nz), Dimensionality.THREE_D, device=device, dtype=dtype)

    # Create solver
    solver = FourPotentialVerletSolver(
        dt=dt_phys,
        c=constants.c,
        mu0=constants.mu0,
        grid=grid,
        bc="periodic",  # Periodic boundaries for wave propagation
        enforce_lorenz=False,  # Initial conditions satisfy Lorenz gauge
        current_source=None,  # Vacuum (no sources)
    )

    print(f"\nSolver Configuration:")
    print(f"  Time step: {dt_phys:.6e} s ({dt_phys*1e15:.3f} fs)")
    print(f"  CFL condition: {solver.cfl_ok()}")
    print(f"  Boundary conditions: periodic")
    print(f"  Lorenz gauge enforcement: {solver.enforce_lorenz}")

    # --- INITIALIZE CIRCULARLY POLARIZED WAVE ---
    print(f"\n{'=' * 70}")
    print("Initializing Circularly Polarized Wave")
    print(f"{'=' * 70}")

    sigma_transverse = 3.0 * wavelength_phys

    # Choose amplitude for vector potential
    # For a plane wave: E₀ = ωA₀, so A₀ = E₀/ω
    # Choose E₀ such that energy density is reasonable
    E_target = 1e6  # 1 MV/m (typical for photon simulations)
    omega = 2.0 * math.pi * constants.c / wavelength_phys
    A_amplitude = E_target / omega

    print(f"  Wavelength: {wavelength_phys*1e9:.3f} nm")
    print(f"  Transverse width σ: {sigma_transverse*1e9:.3f} nm ({sigma_transverse/wavelength_phys:.2f} × λ)")
    print(f"  Vector potential amplitude: {A_amplitude:.6e} V⋅s/m")
    print(f"  Expected E-field magnitude: {E_target:.6e} V/m")
    print(f"  Expected B-field magnitude: {E_target/constants.c:.6e} T")
    print(f"  Propagation: along +x axis")
    print(f"  Polarization: circular (y-z plane)")

    initialize_circular_polarized_waveguide(
        state=state,
        grid=grid,
        wavelength=wavelength_phys,
        sigma_transverse=sigma_transverse,
        amplitude=A_amplitude,
        c=constants.c,
        propagation_axis=0,
    )

    # Initialize accelerations
    solver.initialize_accelerations(state)

    # Debug: Check initial potentials
    print(f"\nInitial Four-Potential Statistics:")
    print(f"  |A^μ| max: {torch.linalg.norm(state.potential, dim=1).max().item():.6e} V⋅s/m")
    print(f"  |∂A^μ/∂t| max: {torch.linalg.norm(state.velocity, dim=1).max().item():.6e} V⋅s/m/s")
    print(f"  A components: A0={state.potential[:, 0].abs().max().item():.2e}, "
          f"Ax={state.potential[:, 1].abs().max().item():.2e}, "
          f"Ay={state.potential[:, 2].abs().max().item():.2e}, "
          f"Az={state.potential[:, 3].abs().max().item():.2e}")

    # Compute initial fields
    E0, B0 = solver.compute_fields(state)
    E0_mag = torch.linalg.norm(E0, dim=1)
    B0_mag = torch.linalg.norm(B0, dim=1)

    print(f"\nInitial Field Statistics:")
    print(f"  |E| max: {E0_mag.max().item():.6e} V/m")
    print(f"  |E| mean: {E0_mag.mean().item():.6e} V/m")
    print(f"  |B| max: {B0_mag.max().item():.6e} T")
    print(f"  |B| mean: {B0_mag.mean().item():.6e} T")

    # Check if fields are non-zero
    if E0_mag.max().item() > 0 and B0_mag.max().item() > 0:
        print(f"  |E|/(c|B|) max: {(E0_mag / (constants.c * B0_mag + 1e-30)).max().item():.4f}")
    else:
        print(f"  ⚠ Warning: Fields are zero or near-zero!")

    # Initial energy
    initial_energy = compute_em_energy(E0, B0, constants.epsilon0, constants.mu0, h_phys, 3)
    print(f"  Initial EM energy: {initial_energy.item():.6e} J")

    if initial_energy.item() < 1e-30:
        print(f"\n⚠ ERROR: Initial energy is effectively zero!")
        print(f"  This suggests the fields were not initialized correctly.")
        print(f"  Check initialization function and field computation.")
        return

    # --- RUN SIMULATION ---
    num_steps = 500
    simulation_time_phys = num_steps * dt_phys
    crossing_time_phys = domain_length_x / constants.c

    print(f"\n{'=' * 70}")
    print("Running Simulation")
    print(f"{'=' * 70}")
    print(f"  Light crossing time: {crossing_time_phys*1e15:.3f} fs")
    print(f"  Simulation time: {simulation_time_phys*1e15:.3f} fs")
    print(f"  Number of steps: {num_steps:,}")
    print(f"  Crossing time / simulation time: {crossing_time_phys / simulation_time_phys:.2f}")

    # Tracking
    times_phys = []
    energies = []

    # Snapshots
    num_snapshots = 5
    snapshot_times = np.linspace(0, simulation_time_phys, num_snapshots)
    snapshots_E = {}
    snapshots_B = {}
    snapshot_steps = {int(t / dt_phys): t for t in snapshot_times}

    # Animation frames
    animation_frames_E = []
    animation_frames_B = []
    animation_times = []
    frame_interval = max(1, num_steps // 200)

    print_interval = max(1, num_steps // 20)

    for step in range(num_steps + 1):
        # Take snapshots
        if step in snapshot_steps:
            E_snap, B_snap = solver.compute_fields(state)
            snapshots_E[snapshot_steps[step]] = E_snap.cpu().numpy()
            snapshots_B[snapshot_steps[step]] = B_snap.cpu().numpy()

        # Animation frames
        if step % frame_interval == 0:
            E_frame, B_frame = solver.compute_fields(state)
            animation_frames_E.append(E_frame.cpu().numpy())
            animation_frames_B.append(B_frame.cpu().numpy())
            animation_times.append(solver.time)

        # Energy tracking
        if step % max(1, num_steps // 100) == 0:
            E_curr, B_curr = solver.compute_fields(state)
            energy = compute_em_energy(E_curr, B_curr, constants.epsilon0, constants.mu0, h_phys, 3)
            times_phys.append(solver.time)
            energies.append(energy.item())

        # Progress printing
        if step % print_interval == 0:
            print(f"  Step {step:8d}/{num_steps}: t={solver.time*1e15:.3f} fs, "
                  f"E={energies[-1] if energies else 0:.6e} J")

        # Advance
        if step < num_steps:
            solver.step(state)

    # Final analysis
    E_final, B_final = solver.compute_fields(state)
    final_energy = compute_em_energy(E_final, B_final, constants.epsilon0, constants.mu0, h_phys, 3).item()
    energy_drift = abs(final_energy - initial_energy.item()) / initial_energy.item()

    print(f"\n{'=' * 70}")
    print("Results")
    print(f"{'=' * 70}")
    print(f"\nEnergy Conservation:")
    print(f"  Initial: {initial_energy.item():.6e} J")
    print(f"  Final:   {final_energy:.6e} J")
    print(f"  Drift:   {energy_drift:.6e} ({energy_drift*100:.6f}%)")

    if energy_drift < 1e-2:
        print(f"  ✓ Good energy conservation")
    else:
        print(f"  ⚠ Energy drift: {energy_drift*100:.2f}%")

    # --- VISUALIZATION ---
    print(f"\n{'=' * 70}")
    print("Creating Visualizations")
    print(f"{'=' * 70}")

    x_coords_nm = np.arange(nx) * h_phys * 1e9
    y_coords_nm = np.arange(ny) * h_phys * 1e9
    z_coords_nm = np.arange(nz) * h_phys * 1e9

    # ========================================================================
    # E-FIELD MAGNITUDE SNAPSHOTS - XY SLICE
    # ========================================================================
    fig_E_xy, axes_E_xy = plt.subplots(num_snapshots, 1, figsize=(14, 12))
    if num_snapshots == 1:
        axes_E_xy = [axes_E_xy]
    fig_E_xy.suptitle(f'E-Field Magnitude (XY slice, z={nz//2})',
                      fontsize=16, fontweight='bold')

    E_max = max([np.linalg.norm(E, axis=1).max() for E in snapshots_E.values()])

    for idx, (step, t) in enumerate(snapshot_steps.items()):
        if t in snapshots_E:
            E_field = snapshots_E[t]
            E_mag = np.linalg.norm(E_field, axis=1)
            E_mag_slice = extract_slice_xy(E_mag, (nx, ny, nz))

            im = axes_E_xy[idx].imshow(E_mag_slice.T, origin='lower',
                                       extent=[x_coords_nm[0], x_coords_nm[-1],
                                              y_coords_nm[0], y_coords_nm[-1]],
                                       cmap='hot',
                                       vmin=0, vmax=E_max,
                                       aspect='equal')

            axes_E_xy[idx].set_ylabel('y [nm]', fontsize=11)
            t_fs = t * 1e15
            axes_E_xy[idx].text(0.02, 0.95, f't = {t_fs:.3f} fs',
                               transform=axes_E_xy[idx].transAxes,
                               fontsize=12, verticalalignment='top',
                               bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

            divider = make_axes_locatable(axes_E_xy[idx])
            cax = divider.append_axes("right", size="3%", pad=0.05)
            plt.colorbar(im, cax=cax, label='|E| [V/m]')

            if idx == num_snapshots - 1:
                axes_E_xy[idx].set_xlabel('x [nm]', fontsize=12)

    plt.tight_layout()
    plt.savefig(run_manager.get_plot_path('E_field_magnitude_xy.png'), dpi=150, bbox_inches='tight')
    print(f"  ✓ Saved: E_field_magnitude_xy.png")
    plt.close(fig_E_xy)

    # ========================================================================
    # B-FIELD MAGNITUDE SNAPSHOTS - XY SLICE
    # ========================================================================
    fig_B_xy, axes_B_xy = plt.subplots(num_snapshots, 1, figsize=(14, 12))
    if num_snapshots == 1:
        axes_B_xy = [axes_B_xy]
    fig_B_xy.suptitle(f'B-Field Magnitude (XY slice, z={nz//2})',
                      fontsize=16, fontweight='bold')

    B_max = max([np.linalg.norm(B, axis=1).max() for B in snapshots_B.values()])

    for idx, (step, t) in enumerate(snapshot_steps.items()):
        if t in snapshots_B:
            B_field = snapshots_B[t]
            B_mag = np.linalg.norm(B_field, axis=1)
            B_mag_slice = extract_slice_xy(B_mag, (nx, ny, nz))

            im = axes_B_xy[idx].imshow(B_mag_slice.T, origin='lower',
                                       extent=[x_coords_nm[0], x_coords_nm[-1],
                                              y_coords_nm[0], y_coords_nm[-1]],
                                       cmap='viridis',
                                       vmin=0, vmax=B_max,
                                       aspect='equal')

            axes_B_xy[idx].set_ylabel('y [nm]', fontsize=11)
            t_fs = t * 1e15
            axes_B_xy[idx].text(0.02, 0.95, f't = {t_fs:.3f} fs',
                               transform=axes_B_xy[idx].transAxes,
                               fontsize=12, verticalalignment='top',
                               bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

            divider = make_axes_locatable(axes_B_xy[idx])
            cax = divider.append_axes("right", size="3%", pad=0.05)
            plt.colorbar(im, cax=cax, label='|B| [T]')

            if idx == num_snapshots - 1:
                axes_B_xy[idx].set_xlabel('x [nm]', fontsize=12)

    plt.tight_layout()
    plt.savefig(run_manager.get_plot_path('B_field_magnitude_xy.png'), dpi=150, bbox_inches='tight')
    print(f"  ✓ Saved: B_field_magnitude_xy.png")
    plt.close(fig_B_xy)

    # ========================================================================
    # E-FIELD VECTOR (Y-Z) SNAPSHOTS - XY SLICE
    # ========================================================================
    fig_E_vec, axes_E_vec = plt.subplots(num_snapshots, 1, figsize=(14, 12))
    if num_snapshots == 1:
        axes_E_vec = [axes_E_vec]
    fig_E_vec.suptitle(f'E-Field Transverse Components (XY slice, z={nz//2})',
                       fontsize=16, fontweight='bold')

    # Find max magnitude for consistent coloring
    max_E_transverse = 0
    for t in snapshots_E.keys():
        E_field = snapshots_E[t]
        Ey_slice = extract_slice_xy(E_field[:, 1], (nx, ny, nz))
        Ez_slice = extract_slice_xy(E_field[:, 2], (nx, ny, nz))
        mag = np.sqrt(Ey_slice**2 + Ez_slice**2).max()
        max_E_transverse = max(max_E_transverse, mag)

    for idx, (step, t) in enumerate(snapshot_steps.items()):
        if t in snapshots_E:
            E_field = snapshots_E[t]
            Ey_slice = extract_slice_xy(E_field[:, 1], (nx, ny, nz))
            Ez_slice = extract_slice_xy(E_field[:, 2], (nx, ny, nz))

            rgb, magnitude, angle = vector_field_to_rgb(Ey_slice, Ez_slice, max_E_transverse)

            axes_E_vec[idx].imshow(np.transpose(rgb, (1, 0, 2)), origin='lower',
                                   extent=[x_coords_nm[0], x_coords_nm[-1],
                                          y_coords_nm[0], y_coords_nm[-1]],
                                   aspect='equal')

            axes_E_vec[idx].set_ylabel('y [nm]', fontsize=11)
            t_fs = t * 1e15
            axes_E_vec[idx].text(0.02, 0.95, f't = {t_fs:.3f} fs',
                                transform=axes_E_vec[idx].transAxes,
                                fontsize=12, verticalalignment='top',
                                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

            if idx == num_snapshots - 1:
                axes_E_vec[idx].set_xlabel('x [nm]', fontsize=12)

    plt.tight_layout()
    plt.savefig(run_manager.get_plot_path('E_field_transverse_xy.png'), dpi=150, bbox_inches='tight')
    print(f"  ✓ Saved: E_field_transverse_xy.png")
    plt.close(fig_E_vec)

    # ========================================================================
    # ENERGY CONSERVATION PLOT
    # ========================================================================
    fig_energy, ax_energy = plt.subplots(figsize=(10, 6))
    times_fs = np.array(times_phys) * 1e15
    energy_array = np.array(energies)
    ax_energy.plot(times_fs, energy_array / initial_energy.item(), 'g-', linewidth=2)
    ax_energy.axhline(y=1.0, color='r', linestyle='--', linewidth=1, alpha=0.5)
    ax_energy.set_xlabel('Time [fs]', fontsize=12)
    ax_energy.set_ylabel('E(t) / E(0)', fontsize=12)
    ax_energy.set_title('EM Energy Conservation', fontsize=14, fontweight='bold')
    ax_energy.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(run_manager.get_plot_path('em_energy_conservation.png'), dpi=150, bbox_inches='tight')
    print(f"  ✓ Saved: em_energy_conservation.png")
    plt.close(fig_energy)

    # ========================================================================
    # ANIMATION - E-FIELD MAGNITUDE (XY slice)
    # ========================================================================
    print(f"\nCreating E-field animation...")
    print(f"  Total frames: {len(animation_frames_E)}")

    fig_anim_E, ax_anim_E = plt.subplots(figsize=(12, 4))

    E_init_mag = np.linalg.norm(animation_frames_E[0], axis=1)
    E_init_slice = extract_slice_xy(E_init_mag, (nx, ny, nz))

    im_anim_E = ax_anim_E.imshow(E_init_slice.T, origin='lower',
                                 extent=[x_coords_nm[0], x_coords_nm[-1],
                                        y_coords_nm[0], y_coords_nm[-1]],
                                 cmap='hot', vmin=0, vmax=E_max,
                                 aspect='equal', animated=True)

    ax_anim_E.set_xlabel('x [nm]', fontsize=12)
    ax_anim_E.set_ylabel('y [nm]', fontsize=12)
    time_text_E = ax_anim_E.text(0.02, 0.95, '', transform=ax_anim_E.transAxes,
                                 fontsize=12, verticalalignment='top',
                                 bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

    divider_E = make_axes_locatable(ax_anim_E)
    cax_E = divider_E.append_axes("right", size="3%", pad=0.05)
    plt.colorbar(im_anim_E, cax=cax_E, label='|E| [V/m]')
    ax_anim_E.set_title('Circularly Polarized Photon - E-Field (XY slice)',
                        fontsize=14, fontweight='bold')

    def animate_E(frame_idx):
        E_mag = np.linalg.norm(animation_frames_E[frame_idx], axis=1)
        E_slice = extract_slice_xy(E_mag, (nx, ny, nz))
        im_anim_E.set_array(E_slice.T)
        t_fs = animation_times[frame_idx] * 1e15
        time_text_E.set_text(f't = {t_fs:.3f} fs')
        return [im_anim_E, time_text_E]

    anim_E = FuncAnimation(fig_anim_E, animate_E, frames=len(animation_frames_E),
                          interval=50, blit=True, repeat=True)

    writer_E = FFMpegWriter(fps=20, bitrate=2000)
    anim_E.save(run_manager.get_plot_path('E_field_magnitude_xy.mp4'), writer=writer_E, dpi=100)
    print(f"  ✓ Saved: E_field_magnitude_xy.mp4")
    plt.close(fig_anim_E)

    # Save configuration
    config = {
        "experiment": "Polarized Photon EM Solver",
        "grid_size": f"{nx}×{ny}×{nz}",
        "num_steps": num_steps,
        "dt_phys": f"{dt_phys:.6e} s",
        "wavelength": f"{wavelength_phys*1e9:.3f} nm",
        "E_target": f"{E_target:.6e} V/m",
        "energy_drift": f"{energy_drift*100:.6f}%",
        "device": str(device),
    }
    run_manager.save_config(config)

    print(f"\n{'=' * 70}")
    print("Experiment Complete!")
    print(f"{'=' * 70}")
    print(f"\nKey Results:")
    print(f"  Domain size: {domain_length_x*1e9:.3f} × {domain_length_y*1e9:.3f} × {domain_length_z*1e9:.3f} nm")
    print(f"  Wavelength: {wavelength_phys*1e9:.3f} nm")
    print(f"  Simulation time: {simulation_time_phys*1e15:.3f} fs")
    print(f"  Energy drift: {energy_drift*100:.6f}%")
    print(f"  |E| max: {E_max:.6e} V/m")
    print(f"  |B| max: {B_max:.6e} T")
    print(f"  |E|/(c|B|): ~{E_max / (constants.c * B_max):.4f} (should be ~1)")

    print(f"\nInterpretation:")
    print(f"  • Pure EM wave in vacuum (no brane mapping)")
    print(f"  • Circular polarization encoded in A_y and A_z rotating 90° out of phase")
    print(f"  • E and B fields computed from four-potential derivatives")
    print(f"  • Energy conservation validates solver accuracy")

    print(f"\nAll outputs saved to: {run_manager.run_dir}")


if __name__ == '__main__':
    main()