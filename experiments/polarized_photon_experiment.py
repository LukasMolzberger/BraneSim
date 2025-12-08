"""
Polarized Photon Experiment - Straight Tubular Waveguide

This experiment tests the inverted EM→brane mapping in a simplified setting:
- Straight waveguide (no torus, no electron)
- Tubular photon mode with transverse Gaussian profile
- Circularly polarized EM fields
- Direct mapping from (E, B) to brane state via energy matching

The goal is to verify that:
1. Amplitude encodes EM energy density correctly
2. Lateral velocities encode Poynting vector direction/magnitude
3. The mapping is physically reasonable before applying to torus geometry
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
from branesim.physics.em_to_brane_mapping import initialize_brane_from_em_fields
from branesim.geometry.tubular_photon_mode import (
    PhotonModeParameters,
    compute_circular_polarization_EB,
)


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

    This creates a tubular mode that:
    - Propagates along the specified axis (default: x-axis, index 0)
    - Has a Gaussian transverse profile in the perpendicular directions
    - Is circularly polarized
    - Has NO longitudinal Gaussian envelope (constant amplitude along propagation)

    Args:
        coords_phys: (N, 3) physical coordinates [m]
        wavelength: Photon wavelength [m]
        sigma_transverse: Gaussian width in transverse directions [m]
        peak_E_field: Peak electric field magnitude [V/m]
        propagation_axis: Axis along which photon propagates (0=x, 1=y, 2=z)
        constants: PhysicalConstants (for c, epsilon0, mu0)

    Returns:
        E_field: (N, 3) electric field [V/m]
        B_field: (N, 3) magnetic field [T]
    """
    if constants is None:
        constants = PhysicalConstants()

    # Convert to numpy for easier manipulation
    coords_np = coords_phys.cpu().numpy()
    N = coords_np.shape[0]

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

    # Wave number
    k = 2.0 * np.pi / wavelength
    omega = 2.0 * np.pi * constants.c / wavelength

    # Longitudinal coordinate (for phase)
    z_long = coords_np[:, propagation_axis]

    # Transverse distance from axis
    r_transverse = np.linalg.norm(transverse_coords, axis=1)

    # Gaussian transverse envelope
    envelope_transverse = np.exp(-0.5 * (r_transverse / sigma_transverse) ** 2)

    # Longitudinal phase (cosine wave)
    phase = k * z_long
    amplitude_long = np.cos(phase)

    # Total amplitude at each point
    amplitude = peak_E_field * envelope_transverse * amplitude_long

    # Circular polarization: E rotates in (normal, binormal) plane
    # Phase of rotation follows the longitudinal phase
    cos_pol = np.cos(phase)
    sin_pol = np.sin(phase)

    # Electric field vector at each point
    E_field_np = amplitude[:, np.newaxis] * (
        cos_pol[:, np.newaxis] * normal_dir[np.newaxis, :] +
        sin_pol[:, np.newaxis] * binormal_dir[np.newaxis, :]
    )

    # Magnetic field: B = (1/c) (t × E) for a plane wave
    # For circular polarization: B is orthogonal to both t and E
    B_magnitude = amplitude / constants.c
    B_dir = (
        -sin_pol[:, np.newaxis] * normal_dir[np.newaxis, :] +
        cos_pol[:, np.newaxis] * binormal_dir[np.newaxis, :]
    )
    B_field_np = B_magnitude[:, np.newaxis] * B_dir

    # Convert back to torch tensors
    device = coords_phys.device
    dtype = coords_phys.dtype
    E_field = torch.from_numpy(E_field_np).to(device=device, dtype=dtype)
    B_field = torch.from_numpy(B_field_np).to(device=device, dtype=dtype)

    return E_field, B_field


# Import visualization functions from photon_3d_experiment
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
    magnitude = np.sqrt(disp_x**2 + disp_y**2)
    angle = np.arctan2(disp_y, disp_x)

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


def main():
    """Run polarized photon experiment with EM→brane mapping."""
    print("=" * 70)
    print("Polarized Photon Experiment - Tubular Mode, Straight Waveguide")
    print("=" * 70)

    constants = PhysicalConstants()
    print(f"\nPhysical Constants:")
    print(f"  c = {constants.c:.6e} m/s")
    print(f"  λ_C = {constants.lambda_C:.6e} m")
    print(f"  ε₀ = {constants.epsilon0:.6e} F/m")
    print(f"  μ₀ = {constants.mu0:.6e} H/m")

    # --- PHYSICAL AND SIMULATION SETUP ---
    # Match photon_3d_experiment.py setup

    wavelength_phys = constants.lambda_C
    points_per_wavelength = 20
    h_phys = wavelength_phys / points_per_wavelength
    cfl_factor = 0.1
    D = 3

    # Universal point mass
    m_point = 2.861821e-27  # kg

    # 3D brane parameters
    rho_D = m_point / (h_phys ** D)
    T_D = rho_D * constants.c**2
    rest_length_phys = 0.0 * h_phys
    c_wave = constants.c
    k_spring = T_D * (h_phys ** (D - 2))

    # Dimensional mapper
    mapper = DimensionalMapper(
        h_phys=h_phys,
        c_light=constants.c,
        mass_reference=m_point
    )

    h_sim = mapper.to_sim_length(h_phys)
    m_sim = mapper.to_sim_mass(m_point)
    k_sim = mapper.to_sim_spring_constant(k_spring)
    c_sim = mapper.to_sim_velocity(c_wave)
    rest_length_sim = mapper.to_sim_length(rest_length_phys)

    # Time step
    dt_phys = cfl_factor * h_phys / c_wave
    dt_sim = mapper.to_sim_time(dt_phys)

    # Domain size
    nx = 100
    ny = 100
    nz = 100
    domain_length_phys_x = nx * h_phys
    domain_length_phys_y = ny * h_phys
    domain_length_phys_z = nz * h_phys

    print(f"\nPhysical Parameters:")
    print(f"  3D volume mass density ρ_3 = {rho_D:.6e} kg/m³")
    print(f"  3D elastic modulus T_3 = {T_D:.6e} Pa")
    print(f"  Spring constant k = {k_spring:.6e} N/m")
    print(f"  Point mass m = {m_point:.6e} kg")
    print(f"  Time step dt = {dt_phys:.6e} s")

    print(f"\nSimulation Configuration:")
    print(f"  Domain: {nx} × {ny} × {nz} points")
    print(f"  Domain (physical): {domain_length_phys_x:.6e} × {domain_length_phys_y:.6e} × {domain_length_phys_z:.6e} m")
    print(f"  Total points: {nx*ny*nz:,}")

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

    # Create brane components
    state = BraneState((nx, ny, nz), Dimensionality.THREE_D, device, dtype)
    state.initialize_flat_configuration(h_sim)
    initial_positions = state.positions.clone()

    state.set_fixed_boundaries()
    print(f"\nBoundary Conditions:")
    print(f"  Fixed points: {state.fixed_mask.sum().item()} / {nx * ny * nz}")

    grid = BraneGrid((nx, ny, nz), Dimensionality.THREE_D, h_sim, device)

    physics = SpringForceComputer(k_sim, rest_length_sim)
    solver = VelocityVerletSolver(dt_sim, m_sim, physics, grid)

    # --- NEW: GENERATE TUBULAR PHOTON MODE ---
    print(f"\n{'=' * 70}")
    print("Generating Tubular Photon Mode")
    print(f"{'=' * 70}")

    # Get grid coordinates in physical units
    coords_sim = grid.get_spatial_coordinates()  # (N, 3) in sim units
    coords_phys = mapper.to_phys_length(coords_sim)  # (N, 3) in meters

    # Photon parameters
    omega_phys = 2.0 * np.pi * constants.c / wavelength_phys
    sigma_transverse = 3.0 * wavelength_phys  # Transverse Gaussian width

    # Peak E field magnitude (arbitrary, will be scaled by energy matching)
    # Choose something reasonable: corresponds to ~1 photon worth of energy
    photon_energy = constants.hbar * omega_phys  # Energy of 1 photon [J]
    waveguide_volume = domain_length_phys_x * (np.pi * sigma_transverse**2)  # Approximate volume
    energy_density_target = photon_energy / waveguide_volume  # [J/m³]
    # For EM: u ≈ ε₀ E²  →  E ≈ √(u/ε₀)
    peak_E_field = np.sqrt(2.0 * energy_density_target / constants.epsilon0)  # [V/m]

    print(f"\nPhoton Mode Parameters:")
    print(f"  Wavelength: {wavelength_phys:.6e} m (= λ_C)")
    print(f"  Angular frequency ω: {omega_phys:.6e} rad/s")
    print(f"  Transverse width σ: {sigma_transverse:.6e} m ({sigma_transverse/wavelength_phys:.2f} × λ)")
    print(f"  Peak E field: {peak_E_field:.6e} V/m")
    print(f"  Propagation: along +x axis")
    print(f"  Polarization: circular")

    # Compute EM fields on grid points
    print(f"\nComputing EM fields on brane grid...")
    E_field_phys, B_field_phys = compute_straight_waveguide_em_fields(
        coords_phys=coords_phys,
        wavelength=wavelength_phys,
        sigma_transverse=sigma_transverse,
        peak_E_field=peak_E_field,
        propagation_axis=0,  # x-axis
        constants=constants,
    )

    # --- NEW: MAP EM FIELDS → BRANE INITIAL STATE ---
    print(f"\n{'=' * 70}")
    print("Applying EM → Brane Mapping")
    print(f"{'=' * 70}")

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

    # Initialize accelerations
    solver.initialize_accelerations(state)
    state.apply_fixed_boundaries()

    # Initial measurements
    initial_energy = solver.compute_energy(state)
    print(f"\nInitial State:")
    print(f"  Energy = {initial_energy['total']:.6e} J")

    # --- RUN SIMULATION ---
    num_steps = 1000
    simulation_time_sim = num_steps * dt_sim
    simulation_time_phys = mapper.to_phys_time(simulation_time_sim)
    crossing_time_phys = domain_length_phys_x / constants.c

    print(f"\n{'=' * 70}")
    print("Running Simulation")
    print(f"{'=' * 70}")
    print(f"  Light crossing time (phys) = {crossing_time_phys:.6e} s")
    print(f"  Simulation time (phys) = {simulation_time_phys:.6e} s")
    print(f"  Number of steps = {num_steps:,}")

    # Tracking
    times_phys = []
    energies = []

    # Snapshots
    num_snapshots = 7
    snapshot_times_phys = np.linspace(0, simulation_time_phys, num_snapshots)
    snapshots = {}
    snapshots_lateral_x = {}
    snapshots_lateral_y = {}
    snapshots_lateral_z = {}
    snapshot_steps = {int(t / dt_phys): t for t in snapshot_times_phys}

    # Animation frames
    animation_frames = []
    animation_frames_lateral_x = []
    animation_frames_lateral_y = []
    animation_times = []
    frame_interval = max(1, num_steps // 200)

    print_interval = max(1, num_steps // 20)

    for step in range(num_steps + 1):
        if step in snapshot_steps:
            field = state.get_field_component(3).cpu().numpy()
            snapshots[snapshot_steps[step]] = field.copy()

            lateral_disp_x = (state.positions[:, 0] - initial_positions[:, 0]).cpu().numpy()
            lateral_disp_y = (state.positions[:, 1] - initial_positions[:, 1]).cpu().numpy()
            lateral_disp_z = (state.positions[:, 2] - initial_positions[:, 2]).cpu().numpy()
            snapshots_lateral_x[snapshot_steps[step]] = lateral_disp_x.copy()
            snapshots_lateral_y[snapshot_steps[step]] = lateral_disp_y.copy()
            snapshots_lateral_z[snapshot_steps[step]] = lateral_disp_z.copy()

        if step % frame_interval == 0:
            field = state.get_field_component(3).cpu().numpy()
            animation_frames.append(field.copy())

            lateral_disp_x = (state.positions[:, 0] - initial_positions[:, 0]).cpu().numpy()
            lateral_disp_y = (state.positions[:, 1] - initial_positions[:, 1]).cpu().numpy()
            animation_frames_lateral_x.append(lateral_disp_x.copy())
            animation_frames_lateral_y.append(lateral_disp_y.copy())
            animation_times.append(solver.time)

        if step % max(1, num_steps // 100) == 0:
            energy = solver.compute_energy(state)
            time_phys = mapper.to_phys_time(solver.time)
            times_phys.append(time_phys)
            energies.append(energy['total'])

        if step % print_interval == 0:
            time_phys = mapper.to_phys_time(solver.time)
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
    print(f"\nEnergy Conservation:")
    print(f"  Initial: {initial_energy['total']:.6e} J")
    print(f"  Final:   {final_energy['total']:.6e} J")
    print(f"  Drift:   {energy_drift:.6e} ({energy_drift*100:.6f}%)")

    if energy_drift < 1e-2:
        print(f"  ✓ Good energy conservation")
    else:
        print(f"  ⚠ Energy drift: {energy_drift*100:.2f}%")

    # --- VISUALIZATION ---
    print(f"\n{'=' * 70}")
    print("Creating Visualizations")
    print(f"{'=' * 70}")

    x_coords_phys = np.arange(nx) * h_phys
    y_coords_phys = np.arange(ny) * h_phys
    z_coords_phys = np.arange(nz) * h_phys
    x_coords = x_coords_phys * 1e9  # nm
    y_coords = y_coords_phys * 1e9  # nm
    z_coords = z_coords_phys * 1e9  # nm

    # Get amplitude scale from initial state
    amplitude_phys = mapper.to_phys_length(state.positions[:, 3].max().item())
    amplitude_nm = amplitude_phys * 1e9

    # 1. Amplitude field snapshots (XY slice)
    fig, axes = plt.subplots(num_snapshots, 1, figsize=(14, 12))
    fig.suptitle(f'Polarized Photon - Amplitude Field (XY slice, z={nz//2})',
                 fontsize=16, fontweight='bold')

    for idx, (step, t) in enumerate(snapshot_steps.items()):
        if t in snapshots:
            field_slice_sim = extract_slice_xy(snapshots[t], (nx, ny, nz))
            field_slice_nm = mapper.to_phys_length(field_slice_sim) * 1e9

            im = axes[idx].imshow(field_slice_nm.T, origin='lower',
                                 extent=[x_coords[0], x_coords[-1],
                                        y_coords[0], y_coords[-1]],
                                 cmap='RdBu_r',
                                 vmin=-amplitude_nm*1.2, vmax=amplitude_nm*1.2,
                                 aspect='auto')

            axes[idx].set_ylabel('y [nm]', fontsize=11)
            t_fs = t * 1e15
            axes[idx].text(0.02, 0.95, f't = {t_fs:.3f} fs',
                          transform=axes[idx].transAxes,
                          fontsize=12, verticalalignment='top',
                          bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

            divider = make_axes_locatable(axes[idx])
            cax = divider.append_axes("right", size="3%", pad=0.05)
            plt.colorbar(im, cax=cax, label='ξ [nm]')

            if idx == num_snapshots - 1:
                axes[idx].set_xlabel('x [nm]', fontsize=12)

    plt.tight_layout()
    plt.savefig('polarized_photon_amplitude_xy.png', dpi=150, bbox_inches='tight')
    print(f"  ✓ Saved: polarized_photon_amplitude_xy.png")

    # 2. Lateral distortion (XY slice)
    max_disp_mag_phys = 0
    for t in snapshots_lateral_x.keys():
        disp_x_slice_sim = extract_slice_xy(snapshots_lateral_x[t], (nx, ny, nz))
        disp_y_slice_sim = extract_slice_xy(snapshots_lateral_y[t], (nx, ny, nz))
        disp_x_slice_phys = mapper.to_phys_length(disp_x_slice_sim)
        disp_y_slice_phys = mapper.to_phys_length(disp_y_slice_sim)
        mag = np.sqrt(disp_x_slice_phys**2 + disp_y_slice_phys**2).max()
        max_disp_mag_phys = max(max_disp_mag_phys, mag)

    fig_lat, axes_lat = plt.subplots(num_snapshots, 1, figsize=(14, 12))
    fig_lat.suptitle(f'Polarized Photon - Lateral Distortion (Color=Direction, Brightness=Magnitude)',
                     fontsize=16, fontweight='bold')

    for idx, (step, t) in enumerate(snapshot_steps.items()):
        if t in snapshots_lateral_x:
            disp_x_slice_sim = extract_slice_xy(snapshots_lateral_x[t], (nx, ny, nz))
            disp_y_slice_sim = extract_slice_xy(snapshots_lateral_y[t], (nx, ny, nz))
            disp_x_slice_phys = mapper.to_phys_length(disp_x_slice_sim)
            disp_y_slice_phys = mapper.to_phys_length(disp_y_slice_sim)

            rgb_image, magnitude, angle = displacement_to_rgb_3d(
                disp_x_slice_phys, disp_y_slice_phys, max_disp_mag_phys
            )

            axes_lat[idx].imshow(np.transpose(rgb_image, (1, 0, 2)), origin='lower',
                                extent=[x_coords[0], x_coords[-1],
                                       y_coords[0], y_coords[-1]],
                                aspect='auto')

            axes_lat[idx].set_ylabel('y [nm]', fontsize=11)
            t_fs = t * 1e15
            axes_lat[idx].text(0.02, 0.95, f't = {t_fs:.3f} fs',
                              transform=axes_lat[idx].transAxes,
                              fontsize=12, verticalalignment='top',
                              bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

            max_mag_pm = max_disp_mag_phys * 1e12
            axes_lat[idx].text(0.98, 0.95, f'max: {max_mag_pm:.2f} pm',
                              transform=axes_lat[idx].transAxes,
                              fontsize=10, verticalalignment='top',
                              horizontalalignment='right',
                              bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.8))

            if idx == num_snapshots - 1:
                axes_lat[idx].set_xlabel('x [nm]', fontsize=12)

    plt.tight_layout()
    plt.savefig('polarized_photon_lateral_distortion.png', dpi=150, bbox_inches='tight')
    print(f"  ✓ Saved: polarized_photon_lateral_distortion.png")

    # 3. Energy conservation
    fig3, ax = plt.subplots(figsize=(10, 6))
    times_fs = np.array(times_phys) * 1e15
    energy_array = np.array(energies)
    ax.plot(times_fs, energy_array / initial_energy['total'], 'g-', linewidth=2)
    ax.axhline(y=1.0, color='r', linestyle='--', linewidth=1, alpha=0.5)
    ax.set_xlabel('Time [fs]', fontsize=12)
    ax.set_ylabel('E(t) / E(0)', fontsize=12)
    ax.set_title('Energy Conservation', fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig('polarized_photon_energy.png', dpi=150, bbox_inches='tight')
    print(f"  ✓ Saved: polarized_photon_energy.png")

    print(f"\n{'=' * 70}")
    print("Experiment Complete!")
    print(f"{'=' * 70}")
    print(f"\nKey Results:")
    print(f"  Maximum amplitude: {amplitude_nm:.3f} nm")
    print(f"  Maximum lateral displacement: {max_disp_mag_phys*1e12:.3f} pm")
    print(f"  Energy drift: {energy_drift*100:.6f}%")
    print(f"\nInterpretation:")
    print(f"  • Amplitude encodes EM energy density via u = (1/2)ρω²A²")
    print(f"  • Lateral velocities encode Poynting vector S/u")
    print(f"  • Circular polarization should appear as rotating displacement pattern")


if __name__ == '__main__':
    main()