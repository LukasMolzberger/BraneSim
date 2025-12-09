"""
1D Photon Amplitude Sweep - Finding the Lateralization Threshold

This experiment sweeps through different photon amplitudes to find the threshold
where geometric nonlinearity causes excessive energy coupling into lateral modes,
preventing proper wave propagation at the speed of light.

Key metrics tracked:
- Wave propagation speed (should be c for small amplitudes)
- Global lateralization ratio R_lat (energy in lateral modes)
- Energy conservation
- Wave packet integrity
"""

import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import math
import torch
import numpy as np
import matplotlib.pyplot as plt
import csv

from branesim.core.state import BraneState, Dimensionality
from branesim.core.grid import BraneGrid
from branesim.core.solver import VelocityVerletSolver
from branesim.physics.forces import SpringForceComputer
from branesim.config.physical_constants import PhysicalConstants
from branesim.physics.dimensional_mapping import DimensionalMapper
from branesim.core.initial_conditions import initialize_right_moving_velocities_time_reversed
from branesim.diagnostics.lateralization import LateralizationMeasurement, LateralizationConfig


def initialize_wave_shape_1d(state, grid, wavelength, amplitude, center):
    """Initialize ONLY the shape of a 1D wave packet."""
    x = torch.arange(grid.grid_shape[0], device=state.device, dtype=state.dtype) * grid.spacing
    k = 2 * np.pi / wavelength
    sigma = 3 * wavelength / (2 * np.pi)
    envelope = amplitude * torch.exp(-((x - center) ** 2) / (2 * sigma ** 2))
    cutoff = 4.0 * sigma
    mask = torch.abs(x - center) <= cutoff
    envelope = envelope * mask
    state.positions[:, 3] = envelope * torch.cos(k * (x - center))


def track_wave_center(state, grid):
    """Track center of wave energy."""
    energy_density = state.velocities[:, 3] ** 2 + state.positions[:, 3] ** 2
    total = energy_density.sum()
    if total > 1e-10:
        x_coords = torch.arange(len(energy_density), device=energy_density.device,
                               dtype=energy_density.dtype)
        center = (x_coords * energy_density).sum() / total
        return center.item() * grid.spacing
    return 0.0


def export_csv_snapshot(filename, state, grid, initial_positions, spring_constant,
                       lateralization, physics, h, mapper, rest_length=0.0):
    """
    Export detailed CSV snapshot with all brane point data in SI units.

    All quantities are converted from simulation units to physical SI units before export.

    Columns: point_idx, x_position, xi_position, x_velocity, xi_velocity,
             x_acceleration, xi_acceleration, delta_x, delta_xi, h_spacing,
             L_to_next, delta_L_to_next,
             F_left_x, F_left_xi, F_right_x, F_right_xi,
             E_amp_kin, E_amp_pot, E_lat_kin, E_lat_pot, R_lat

    Units:
        - Positions, displacements, lengths, h_spacing: meters [m]
        - Velocities: meters per second [m/s]
        - Accelerations: meters per second squared [m/s²]
        - Forces: Newtons [N]
        - Energies: Joules [J] (already in physical units)
        - R_lat: dimensionless
    """
    nx = state.positions.shape[0]
    neighbors = grid.neighbors

    # Get lateralization measurement
    R_lat_local, R_lat_global, diagnostics = lateralization.measure(state, physics)

    # Convert to numpy (still in sim units)
    positions_sim = state.positions.cpu().numpy()
    velocities_sim = state.velocities.cpu().numpy()
    accelerations_sim = state.accelerations.cpu().numpy()
    initial_pos_sim = initial_positions.cpu().numpy()
    R_lat = R_lat_local.cpu().numpy()
    E_amp_kin = diagnostics['E_amp_kin'].cpu().numpy()  # Already in J
    E_amp_pot = diagnostics['E_amp_pot'].cpu().numpy()  # Already in J
    E_lat_kin = diagnostics['E_lat_kin'].cpu().numpy()  # Already in J
    E_lat_pot = diagnostics['E_lat_pot'].cpu().numpy()  # Already in J

    # Compute spring forces for each point (in sim units)
    k_sim = spring_constant
    rest_length_sim = rest_length

    # Store forces from left and right neighbors (in sim units)
    F_left_sim = np.zeros((nx, 4))
    F_right_sim = np.zeros((nx, 4))

    for i in range(nx):
        for j_idx in neighbors[i]:
            j = j_idx.item() if isinstance(j_idx, torch.Tensor) else int(j_idx)
            if j < 0:
                continue

            # Vector from i to j (sim units)
            d = positions_sim[j] - positions_sim[i]
            length = np.sqrt(np.sum(d**2))

            if length < 1e-20:
                continue

            # Spring force magnitude: F = k * (L - L0) (sim units)
            extension = length - rest_length_sim
            force_mag = k_sim * extension

            # Force direction (unit vector)
            direction = d / length

            # Force vector on point i from point j (sim units)
            force = force_mag * direction

            # Determine if this is left or right neighbor
            if j == i - 1:
                F_left_sim[i] = force
            elif j == i + 1:
                F_right_sim[i] = force

    with open(filename, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)

        # Write header with units
        writer.writerow([
            'point_idx',
            'x_position [m]', 'xi_position [m]',
            'x_velocity [m/s]', 'xi_velocity [m/s]',
            'x_acceleration [m/s^2]', 'xi_acceleration [m/s^2]',
            'delta_x [m]', 'delta_xi [m]',
            'h_spacing [m]',
            'L_to_next [m]', 'delta_L_to_next [m]',
            'F_left_x [N]', 'F_left_xi [N]', 'F_right_x [N]', 'F_right_xi [N]',
            'E_amp_kin [J]', 'E_amp_pot [J]', 'E_lat_kin [J]', 'E_lat_pot [J]', 'R_lat'
        ])

        # Write data for each point (convert to physical units using mapper)
        for i in range(nx):
            # Positions [sim → m]
            x_pos = mapper.to_phys_length(positions_sim[i, 0])
            xi_pos = mapper.to_phys_length(positions_sim[i, 3])

            # Velocities [sim → m/s]
            x_vel = mapper.to_phys_velocity(velocities_sim[i, 0])
            xi_vel = mapper.to_phys_velocity(velocities_sim[i, 3])

            # Accelerations [sim → m/s²]
            x_acc = mapper.to_phys_acceleration(accelerations_sim[i, 0])
            xi_acc = mapper.to_phys_acceleration(accelerations_sim[i, 3])

            # Displacements from initial [sim → m]
            delta_x = mapper.to_phys_length(positions_sim[i, 0] - initial_pos_sim[i, 0])
            xi = mapper.to_phys_length(positions_sim[i, 3] - initial_pos_sim[i, 3])

            # Spring properties (to next neighbor if exists)
            L_to_next_phys = 0.0
            delta_L_to_next_phys = 0.0

            if i < nx - 1:  # Not the last point
                # Find next neighbor
                for j_idx in neighbors[i]:
                    j = j_idx.item() if isinstance(j_idx, torch.Tensor) else int(j_idx)
                    if j == i + 1:  # Next neighbor
                        # Current spring vector (sim units)
                        dX_sim = positions_sim[j] - positions_sim[i]
                        L_to_next_sim = np.sqrt(np.sum(dX_sim**2))

                        # Reference spring vector (sim units)
                        dX0_sim = initial_pos_sim[j] - initial_pos_sim[i]
                        L_ref_sim = np.sqrt(np.sum(dX0_sim**2))

                        # Extension beyond reference (sim units)
                        delta_L_to_next_sim = L_to_next_sim - L_ref_sim

                        # Convert to physical [m]
                        L_to_next_phys = mapper.to_phys_length(L_to_next_sim)
                        delta_L_to_next_phys = mapper.to_phys_length(delta_L_to_next_sim)
                        break

            # Forces from neighbors [sim → N]
            f_left_x = mapper.to_phys_force(F_left_sim[i, 0])
            f_left_xi = mapper.to_phys_force(F_left_sim[i, 3])
            f_right_x = mapper.to_phys_force(F_right_sim[i, 0])
            f_right_xi = mapper.to_phys_force(F_right_sim[i, 3])

            # Energy and lateralization (already in J and dimensionless)
            e_amp_kin = E_amp_kin[i]
            e_amp_pot = E_amp_pot[i]
            e_lat_kin = E_lat_kin[i]
            e_lat_pot = E_lat_pot[i]
            r_lat = R_lat[i]

            # Grid spacing in physical units [m]
            h_phys = mapper.to_phys_length(h)

            # Write row (all in SI units)
            writer.writerow([
                i,
                x_pos, xi_pos,
                x_vel, xi_vel,
                x_acc, xi_acc,
                delta_x, xi,
                h_phys,
                L_to_next_phys, delta_L_to_next_phys,
                f_left_x, f_left_xi, f_right_x, f_right_xi,
                e_amp_kin, e_amp_pot, e_lat_kin, e_lat_pot, r_lat
            ])


def run_amplitude_simulation(amplitude_fraction, constants, h_phys, wavelength_phys, nx,
                            num_steps=2000, verbose=False, export_csv=False):
    """
    Run simulation for a given amplitude and return key metrics.

    Parameters
    ----------
    amplitude_fraction : float
        Amplitude as fraction of grid spacing (e.g., 0.1 = 10% of h)
    constants : PhysicalConstants
        Physical constants
    h_phys : float
        Grid spacing in meters
    wavelength_phys : float
        Photon wavelength in meters
    nx : int
        Number of grid points
    num_steps : int
        Number of simulation steps to run (default 2000)
    verbose : bool
        Print detailed output

    Returns
    -------
    results : dict
        Dictionary with simulation results and spatial field data
    """
    if verbose:
        print(f"\n{'='*70}")
        print(f"Running amplitude = {amplitude_fraction:.3f} × h")
        print(f"{'='*70}")

    # Setup (same as main example)
    D = 1
    rho_D = 2.3590e-14  # kg/m
    T_D = rho_D * constants.c**2

    # Rest length from physically calibrated constant
    # L0 = rest_length_frac × a, where rest_length_frac is from continuum calibration
    rest_length_phys = constants.rest_length_frac * h_phys
    c_wave = constants.c
    m_point = rho_D * (h_phys ** D)
    k_spring = T_D * (h_phys ** (D - 2))

    mapper = DimensionalMapper(h_phys=h_phys, c_light=constants.c, mass_reference=m_point)

    h_sim = mapper.to_sim_length(h_phys)
    m_sim = mapper.to_sim_mass(m_point)
    k_sim = mapper.to_sim_spring_constant(k_spring)
    c_wave_sim = mapper.to_sim_velocity(c_wave)
    rest_length_sim = mapper.to_sim_length(rest_length_phys)

    cfl_factor = 0.1
    dt_phys = cfl_factor * h_phys / c_wave
    dt_sim = mapper.to_sim_time(dt_phys)

    domain_length_phys = nx * h_phys

    # Create state
    device = torch.device('cpu')
    dtype = torch.float64
    state = BraneState((nx,), Dimensionality.ONE_D, device, dtype)
    state.initialize_flat_configuration(h_sim)
    initial_positions = state.positions.clone()
    state.set_fixed_boundaries()

    grid = BraneGrid((nx,), Dimensionality.ONE_D, h_sim, device)
    physics = SpringForceComputer(k_sim, rest_length_sim)
    solver = VelocityVerletSolver(dt_sim, m_sim, physics, grid)

    lat_config = LateralizationConfig(amplitude_dim=3, lateral_dims=(0,))
    lateralization = LateralizationMeasurement(
        config=lat_config, grid=grid, m_point=m_point, reference_positions=initial_positions
    )

    # Initialize with specified amplitude
    amplitude_phys = amplitude_fraction * h_phys
    center_position_phys = domain_length_phys / 3.0

    amplitude_sim = mapper.to_sim_length(amplitude_phys)
    wavelength_sim = mapper.to_sim_length(wavelength_phys)
    center_position_sim = mapper.to_sim_length(center_position_phys)

    initialize_wave_shape_1d(state, grid, wavelength_sim, amplitude_sim, center_position_sim)

    initialize_right_moving_velocities_time_reversed(
        state=state, grid=grid, physics=physics,
        m_point=m_sim, wave_speed=c_wave_sim,
        field_component=3, shift_cells=1
    )

    solver.initialize_accelerations(state)
    state.apply_fixed_boundaries()

    initial_energy = solver.compute_energy(state)
    initial_center_sim = track_wave_center(state, grid)

    # Measure initial lateralization
    R_lat_local, R_lat_global_initial, diagnostics = lateralization.measure(state, physics)

    # Run simulation for fixed number of steps
    simulation_time_phys = num_steps * dt_phys

    # Track metrics
    R_lat_global_history = [R_lat_global_initial]
    max_lateral_disp = 0.0

    for step in range(num_steps):
        solver.step(state)

        # Sample every 10% of simulation
        if step % max(1, num_steps // 10) == 0:
            R_lat_local, R_lat_global, diagnostics = lateralization.measure(state, physics)
            R_lat_global_history.append(R_lat_global)

            # Track max lateral displacement
            lateral_disp = (state.positions[:, 0] - initial_positions[:, 0]).abs().max().item()
            lateral_disp_phys = mapper.to_phys_length(lateral_disp)
            max_lateral_disp = max(max_lateral_disp, lateral_disp_phys)

    # Final measurements at step num_steps
    final_energy = solver.compute_energy(state)
    final_center_sim = track_wave_center(state, grid)

    distance_traveled_sim = final_center_sim - initial_center_sim
    distance_traveled_phys = mapper.to_phys_length(distance_traveled_sim)
    measured_speed = distance_traveled_phys / simulation_time_phys
    speed_error = abs(measured_speed - constants.c) / constants.c
    energy_drift = abs(final_energy['total'] - initial_energy['total']) / initial_energy['total']

    R_lat_local_final, R_lat_global_final, diagnostics = lateralization.measure(state, physics)
    R_lat_global_mean = np.mean(R_lat_global_history)
    R_lat_global_max = np.max(R_lat_global_history)

    # Capture spatial field data at final step (for sweep diagrams)
    amplitude_field_sim = state.get_field_component(3).cpu().numpy()  # ξ component
    lateral_disp_sim = (state.positions[:, 0] - initial_positions[:, 0]).cpu().numpy()  # Δx component

    # Convert to physical units
    amplitude_field_phys = mapper.to_phys_length(amplitude_field_sim)
    lateral_disp_phys_array = mapper.to_phys_length(lateral_disp_sim)

    if verbose:
        print(f"  Wavelength: {wavelength_phys:.6e} m | Amplitude: {amplitude_phys:.6e} m ({amplitude_fraction:.3f} × h)")
        print(f"  Speed: {measured_speed:.6e} m/s (error: {speed_error*100:.4f}%)")
        print(f"  Energy drift: {energy_drift:.6e}")
        print(f"  R_lat (final): {R_lat_global_final:.6f}")
        print(f"  R_lat (mean): {R_lat_global_mean:.6f}")
        print(f"  R_lat (max): {R_lat_global_max:.6f}")
        print(f"  Max lateral displacement: {max_lateral_disp:.6e} m")

    # Export CSV snapshot if requested
    if export_csv:
        csv_filename = f'photon_sweep_amp_{amplitude_fraction:.1f}h.csv'
        export_csv_snapshot(csv_filename, state, grid, initial_positions,
                          k_sim, lateralization, physics, h_sim,
                          mapper, rest_length=rest_length_sim)
        if verbose:
            print(f"  ✓ Exported CSV: {csv_filename}")

    return {
        'amplitude_fraction': amplitude_fraction,
        'amplitude_phys': amplitude_phys,
        'measured_speed': measured_speed,
        'speed_error': speed_error,
        'energy_drift': energy_drift,
        'R_lat_initial': R_lat_global_initial,
        'R_lat_final': R_lat_global_final,
        'R_lat_mean': R_lat_global_mean,
        'R_lat_max': R_lat_global_max,
        'R_lat_history': R_lat_global_history,
        'max_lateral_disp': max_lateral_disp,
        'initial_energy': initial_energy['total'],
        'final_energy': final_energy['total'],
        # Spatial field data at final step
        'amplitude_field': amplitude_field_phys,
        'lateral_disp': lateral_disp_phys_array,
        'dt_phys': dt_phys,  # Return for time calculation
    }


def main():
    """Run amplitude sweep experiment."""
    print("="*70)
    print("1D Photon Amplitude Sweep - Lateralization Threshold")
    print("="*70)

    constants = PhysicalConstants()

    print(f"\nPhysical Constants:")
    print(f"  Speed of light c = {constants.c:.6e} m/s")
    print(f"  Compton wavelength λ_C = {constants.lambda_C:.6e} m")

    # Configuration (same as main example)
    wavelength_phys = constants.lambda_C
    points_per_wavelength = 40
    h_phys = wavelength_phys / points_per_wavelength
    nx = 400

    print(f"\nSimulation Configuration:")
    print(f"  Photon wavelength: λ_C = {wavelength_phys:.6e} m")
    print(f"  Grid spacing: h = λ_C / {points_per_wavelength} = {h_phys:.6e} m")
    print(f"  Domain: {nx} points = {nx * h_phys:.6e} m")

    # Theoretical prediction from geometric threshold analysis
    A_th_theory = wavelength_phys / (2 * np.pi)  # λ/(2π)
    A_th_h = A_th_theory / h_phys  # In units of h
    print(f"\nTheoretical Geometric Threshold:")
    print(f"  Critical amplitude A_th ≈ λ/(2π) = {A_th_theory:.6e} m")
    print(f"  In units of h: A_th ≈ {A_th_h:.2f}h")
    print(f"  (50/50 force partition between lateral and amplitude directions)")

    # Amplitude sweep based on geometric threshold analysis:
    # Theory predicts critical amplitude A_th ≈ λ/(2π) ≈ 6.37h
    # (see docs/threshold-analysis-4-12-2025.md)
    # Empirically, threshold appears 2-3 orders of magnitude higher than naive geometric prediction
    # Sweep from 10h to 10000h to capture actual threshold behavior
    amplitude_fractions = np.logspace(1, 4, 30)  # 10.0 to 10000.0 (30 samples)

    # Fixed simulation time
    num_steps = 2000

    print(f"\nAmplitude sweep:")
    print(f"  Testing {len(amplitude_fractions)} amplitudes from {amplitude_fractions[0]:.4f}h to {amplitude_fractions[-1]:.4f}h")
    print(f"  Fixed simulation time: {num_steps} steps")

    results = []

    print(f"\nRunning simulations...")
    # Export CSV for selected amplitudes to analyze forces
    csv_export_indices = [0, 9, 19, 29]  # 10h, ~85h, ~924h, 10000h

    for i, amp_frac in enumerate(amplitude_fractions):
        print(f"\n[{i+1}/{len(amplitude_fractions)}] ", end='')
        # Enable CSV export for selected amplitudes
        export_csv = (i in csv_export_indices)
        result = run_amplitude_simulation(
            amp_frac, constants, h_phys, wavelength_phys, nx,
            num_steps=num_steps,
            verbose=True,
            export_csv=export_csv
        )
        results.append(result)

    # Analysis and plotting
    print(f"\n{'='*70}")
    print("Analysis Complete - Creating Plots")
    print(f"{'='*70}")

    # Extract data for plotting
    amp_fractions = np.array([r['amplitude_fraction'] for r in results])
    amp_phys = np.array([r['amplitude_phys'] for r in results])
    speed_errors = np.array([r['speed_error'] for r in results])
    energy_drifts = np.array([r['energy_drift'] for r in results])
    R_lat_means = np.array([r['R_lat_mean'] for r in results])
    R_lat_maxs = np.array([r['R_lat_max'] for r in results])
    R_lat_finals = np.array([r['R_lat_final'] for r in results])
    max_lateral_disps = np.array([r['max_lateral_disp'] for r in results])

    # Create comprehensive figure
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle('1D Photon Amplitude Sweep - Lateralization Threshold Analysis',
                 fontsize=16, fontweight='bold')

    # 1. Speed error vs amplitude
    ax1 = axes[0, 0]
    ax1.loglog(amp_fractions, speed_errors * 100, 'b-o', linewidth=2, markersize=6)
    ax1.axhline(y=1, color='r', linestyle='--', linewidth=1, alpha=0.5, label='1% error threshold')
    ax1.axhline(y=5, color='orange', linestyle='--', linewidth=1, alpha=0.5, label='5% error threshold')
    ax1.set_xlabel('Amplitude [× h]', fontsize=11)
    ax1.set_ylabel('Wave Speed Error [%]', fontsize=11)
    ax1.set_title('Wave Propagation Speed Deviation', fontsize=12, fontweight='bold')
    ax1.grid(True, alpha=0.3, which='both')
    ax1.legend(fontsize=9)

    # 2. Lateralization ratio vs amplitude
    ax2 = axes[0, 1]
    ax2.semilogx(amp_fractions, R_lat_means, 'g-o', linewidth=2, markersize=6, label='Mean R_lat')
    ax2.semilogx(amp_fractions, R_lat_maxs, 'orange', linestyle='--', linewidth=2, marker='s',
                 markersize=5, label='Max R_lat')
    ax2.axhline(y=0.1, color='r', linestyle='--', linewidth=1, alpha=0.5, label='10% threshold')
    ax2.axhline(y=0.5, color='purple', linestyle='--', linewidth=1, alpha=0.5, label='50% threshold')
    # Naive geometric threshold (A_th ≈ 6.37h) is outside sweep range - empirically too small by 2-3 orders of magnitude
    ax2.set_xlabel('Amplitude [× h]', fontsize=11)
    ax2.set_ylabel('Lateralization Ratio R_lat', fontsize=11)
    ax2.set_title('Energy Coupling to Lateral Modes', fontsize=12, fontweight='bold')
    ax2.set_ylim(0, 1)
    ax2.grid(True, alpha=0.3)
    ax2.legend(fontsize=9)

    # 3. Max lateral displacement vs amplitude
    ax3 = axes[1, 0]
    ax3.loglog(amp_fractions, max_lateral_disps * 1e12, 'purple', linewidth=2, marker='o', markersize=6)
    ax3.set_xlabel('Amplitude [× h]', fontsize=11)
    ax3.set_ylabel('Max Lateral Displacement [pm]', fontsize=11)
    ax3.set_title('Maximum Lateral Distortion', fontsize=12, fontweight='bold')
    ax3.grid(True, alpha=0.3, which='both')

    # 4. Energy conservation vs amplitude
    ax4 = axes[1, 1]
    ax4.loglog(amp_fractions, energy_drifts * 100, 'r-o', linewidth=2, markersize=6)
    ax4.axhline(y=0.01, color='g', linestyle='--', linewidth=1, alpha=0.5, label='0.01% threshold')
    ax4.axhline(y=0.1, color='orange', linestyle='--', linewidth=1, alpha=0.5, label='0.1% threshold')
    ax4.set_xlabel('Amplitude [× h]', fontsize=11)
    ax4.set_ylabel('Energy Drift [%]', fontsize=11)
    ax4.set_title('Energy Conservation', fontsize=12, fontweight='bold')
    ax4.grid(True, alpha=0.3, which='both')
    ax4.legend(fontsize=9)

    plt.tight_layout()
    plt.savefig('photon_1d_amplitude_sweep.png', dpi=150, bbox_inches='tight')
    print(f"  ✓ Saved: photon_1d_amplitude_sweep.png")

    # Find thresholds
    print(f"\n{'='*70}")
    print("Threshold Analysis")
    print(f"{'='*70}")

    # Find amplitude where speed error exceeds 1%
    idx_1pct = np.where(speed_errors > 0.01)[0]
    if len(idx_1pct) > 0:
        threshold_1pct = amp_fractions[idx_1pct[0]]
        print(f"\n1% Speed Error Threshold:")
        print(f"  Amplitude: {threshold_1pct:.4f} × h = {threshold_1pct * h_phys:.6e} m")
        print(f"  R_lat at threshold: {R_lat_means[idx_1pct[0]]:.4f}")

    # Find amplitude where R_lat exceeds 10%
    idx_10pct_lat = np.where(R_lat_means > 0.1)[0]
    if len(idx_10pct_lat) > 0:
        threshold_10pct_lat = amp_fractions[idx_10pct_lat[0]]
        print(f"\n10% Lateralization Threshold:")
        print(f"  Amplitude: {threshold_10pct_lat:.4f} × h = {threshold_10pct_lat * h_phys:.6e} m")
        print(f"  Speed error at threshold: {speed_errors[idx_10pct_lat[0]]*100:.4f}%")

    # Find amplitude where R_lat exceeds 50%
    idx_50pct_lat = np.where(R_lat_means > 0.5)[0]
    if len(idx_50pct_lat) > 0:
        threshold_50pct_lat = amp_fractions[idx_50pct_lat[0]]
        print(f"\n50% Lateralization Threshold (Energy equipartition):")
        print(f"  Amplitude: {threshold_50pct_lat:.4f} × h = {threshold_50pct_lat * h_phys:.6e} m")
        print(f"  Speed error at threshold: {speed_errors[idx_50pct_lat[0]]*100:.4f}%")
    else:
        print(f"\n50% Lateralization Threshold:")
        print(f"  Not reached in sweep (max R_lat = {R_lat_means.max():.4f})")

    print(f"\nSummary:")
    print(f"  Amplitude range tested: {amp_fractions[0]:.4f}h to {amp_fractions[-1]:.4f}h")
    print(f"  Max lateralization reached: {R_lat_maxs.max():.4f}")
    print(f"  Max speed error: {speed_errors.max()*100:.4f}%")

    # Create spatial field diagrams (similar to photon_1d_example but over amplitude instead of time)
    print(f"\nCreating spatial field diagrams...")

    # Select subset of amplitudes to plot (evenly spaced, focusing on threshold region)
    plot_indices = [0, 4, 8, 12, 16, 20, 24, 27, 29]  # 9 amplitudes from 0.1h to 15h
    num_plots = len(plot_indices)

    # Physical x-coordinates
    x_coords_phys = np.arange(nx) * h_phys

    # Calculate physical time from dt_phys (use first result, all have same dt_phys)
    dt_phys = results[0]['dt_phys']
    simulation_time_phys = num_steps * dt_phys
    simulation_time_fs = simulation_time_phys * 1e15  # Convert to femtoseconds

    # 1. Amplitude field ξ(x) for different amplitudes
    fig_amp, axes_amp = plt.subplots(num_plots, 1, figsize=(14, 12))
    fig_amp.suptitle(f'1D Photon Amplitude Field at t={simulation_time_fs:.3f} fs',
                     fontsize=16, fontweight='bold')

    for idx, result_idx in enumerate(plot_indices):
        r = results[result_idx]
        x_nm = x_coords_phys * 1e9
        field_nm = r['amplitude_field'] * 1e9  # Already in physical units

        axes_amp[idx].plot(x_nm, field_nm, 'b-', linewidth=2)
        axes_amp[idx].axhline(y=0, color='k', linestyle='--', linewidth=0.5, alpha=0.3)
        axes_amp[idx].plot([x_nm[0], x_nm[-1]], [0, 0], 'ro',
                          markersize=6, label='Fixed boundaries' if idx == 0 else '')

        axes_amp[idx].set_ylabel('ξ [nm]', fontsize=10)
        axes_amp[idx].set_xlim(x_nm[0], x_nm[-1])
        axes_amp[idx].grid(True, alpha=0.3)

        # Label with amplitude
        amp_frac = r['amplitude_fraction']
        R_lat = r['R_lat_mean']
        axes_amp[idx].text(0.02, 0.95, f'A = {amp_frac:.3f}h, R_lat={R_lat:.3f}',
                          transform=axes_amp[idx].transAxes,
                          fontsize=11, verticalalignment='top',
                          bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

        if idx == 0:
            axes_amp[idx].legend(loc='upper right', fontsize=9)
        if idx == num_plots - 1:
            axes_amp[idx].set_xlabel('Position [nm]', fontsize=11)

    plt.tight_layout()
    plt.savefig('photon_1d_amplitude_sweep_fields.png', dpi=150, bbox_inches='tight')
    print(f"  ✓ Saved: photon_1d_amplitude_sweep_fields.png")

    # 2. Lateral displacement Δx(x) for different amplitudes
    fig_lat, axes_lat = plt.subplots(num_plots, 1, figsize=(14, 12))
    fig_lat.suptitle(f'1D Photon Lateral Displacement at t={simulation_time_fs:.3f} fs',
                     fontsize=16, fontweight='bold')

    for idx, result_idx in enumerate(plot_indices):
        r = results[result_idx]
        x_nm = x_coords_phys * 1e9
        lateral_disp_pm = r['lateral_disp'] * 1e12  # Already in physical units, convert to pm

        axes_lat[idx].plot(x_nm, lateral_disp_pm, 'r-', linewidth=2)
        axes_lat[idx].axhline(y=0, color='k', linestyle='--', linewidth=0.5, alpha=0.3)
        axes_lat[idx].plot([x_nm[0], x_nm[-1]], [0, 0], 'ro',
                          markersize=6, label='Fixed boundaries' if idx == 0 else '')

        axes_lat[idx].set_ylabel('Δx [pm]', fontsize=10)
        axes_lat[idx].set_xlim(x_nm[0], x_nm[-1])
        axes_lat[idx].grid(True, alpha=0.3)

        # Label with amplitude
        amp_frac = r['amplitude_fraction']
        R_lat = r['R_lat_mean']
        axes_lat[idx].text(0.02, 0.95, f'A = {amp_frac:.3f}h, R_lat={R_lat:.3f}',
                          transform=axes_lat[idx].transAxes,
                          fontsize=11, verticalalignment='top',
                          bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

        if idx == 0:
            axes_lat[idx].legend(loc='upper right', fontsize=9)
        if idx == num_plots - 1:
            axes_lat[idx].set_xlabel('Position [nm]', fontsize=11)

    plt.tight_layout()
    plt.savefig('photon_1d_amplitude_sweep_lateral.png', dpi=150, bbox_inches='tight')
    print(f"  ✓ Saved: photon_1d_amplitude_sweep_lateral.png")

    print(f"\n{'='*70}")
    print("Experiment Complete!")
    print(f"{'='*70}")


if __name__ == '__main__':
    main()