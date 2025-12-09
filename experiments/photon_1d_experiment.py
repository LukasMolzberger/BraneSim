"""
1D Photon with Realistic Physical Scales

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
import csv

from branesim.core.state import BraneState, Dimensionality
from branesim.core.grid import BraneGrid
from branesim.core.solver import VelocityVerletSolver
from branesim.physics.forces import SpringForceComputer
from branesim.config.physical_constants import PhysicalConstants
from branesim.physics.dimensional_mapping import DimensionalMapper
from branesim.core.initial_conditions import (
    initialize_right_moving_velocities_time_reversed,
    verify_wave_propagation,
)
from branesim.diagnostics.lateralization import (
    LateralizationMeasurement,
    LateralizationConfig,
)
from branesim.utils import TestRunManager


def initialize_wave_shape_1d(state, grid, wavelength, amplitude, center):
    """
    Initialize ONLY the shape of a 1D wave packet.

    Velocities are initialized separately using time-reversal initialization.
    This uses the actual brane forces to compute consistent initial velocities.
    """
    x = torch.arange(grid.grid_shape[0], device=state.device, dtype=state.dtype) * grid.spacing

    k = 2 * np.pi / wavelength
    sigma = 3 * wavelength / (2 * np.pi)  # Width of wave packet

    # Gaussian envelope
    envelope = amplitude * torch.exp(-((x - center) ** 2) / (2 * sigma ** 2))

    # Hard-truncate Gaussian at 4σ to ensure compact support
    # This prevents far-field regions from showing lateral motion before photon arrival
    cutoff = 4.0 * sigma
    mask = torch.abs(x - center) <= cutoff
    envelope = envelope * mask

    # Position field only - velocities set separately
    state.positions[:, 3] = envelope * torch.cos(k * (x - center))

    print(f"  Wavelength λ = {wavelength:.6e} m ({wavelength/grid.spacing:.1f} × h)")
    print(f"  Width σ = {sigma:.6e} m")
    print(f"  Center = {center:.6e} m")
    print(f"  Amplitude = {amplitude:.6e} m")
    print(f"  Wave number k = {k:.6e} rad/m")


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


def main():
    """Run 1D photon simulation experiment."""
    print("=" * 70)
    print("1D Photon - Experiment")
    print("=" * 70)

    # Initialize test run manager
    run_manager = TestRunManager(experiment_name="photon_1d_experiment")
    print(run_manager.get_summary())

    # Physical constants
    constants = PhysicalConstants()

    print(f"\nPhysical Constants:")
    print(f"  Speed of light c = {constants.c:.6e} m/s")
    print(f"  Compton wavelength λ_C = {constants.lambda_C:.6e} m")
    print(f"  ℏ = {constants.hbar:.6e} J·s")
    print(f"  m_e = {constants.m_e:.6e} kg")

    # Configuration
    # Set photon wavelength to exactly the Compton wavelength
    wavelength_phys = constants.lambda_C  # Photon wavelength = λ_C
    points_per_wavelength = 20  # Grid resolution
    h_phys = wavelength_phys / points_per_wavelength  # Grid spacing
    cfl_factor = 0.1

    D = 1

    # Universal point mass (same for all dimensions - ensures consistent physics)
    m_point = 2.861821e-27  # kg (fixed for all 1D/2D/3D simulations)

    # 1D brane parameters constrained to give wave speed = c
    # Derive mass density from point mass: ρ_D = m_point / h^D
    rho_D = m_point / (h_phys ** D)  # kg/m (linear mass density)
    T_D = rho_D * constants.c**2  # N (tension - computed from c² = T_D/rho_D)

    # Rest length from physically calibrated constant
    # L0 = rest_length_frac × a, where rest_length_frac is from continuum calibration
    rest_length_phys = constants.rest_length_frac * h_phys

    # Wave speed (exactly equals c by construction)
    c_wave = constants.c

    # Axial spring constant (varies with substrate_scale, since T_D varies)
    k_spring = T_D * (h_phys ** (D - 2))



    # Create dimensional mapper for unit conversions
    mapper = DimensionalMapper(
        h_phys=h_phys,
        c_light=constants.c,
        mass_reference=m_point
    )

    # Simulation uses dimensionless units
    # h_sim = 1.0 ALWAYS (by choice of L0 = h_phys)
    # c_wave = c ALWAYS (constraint: T_D/ρ_D = c²)
    h_sim = mapper.to_sim_length(h_phys)  # = 1.0 always
    m_sim = mapper.to_sim_mass(m_point)  # = substrate_scale
    k_sim = mapper.to_sim_spring_constant(k_spring)  # = substrate_scale
    c_wave_sim = mapper.to_sim_velocity(c_wave)  # = 1.0 always (c_wave = c)
    rest_length_sim = mapper.to_sim_length(rest_length_phys)

    # Time step calculation (CFL condition based on ACTUAL wave speed)
    # CRITICAL: Must use c_wave (actual wave speed in brane), not c_light!
    dt_phys = cfl_factor * h_phys / c_wave
    dt_sim = mapper.to_sim_time(dt_phys)

    # Domain size
    nx = 2000
    domain_length_phys = nx * h_phys
    domain_length_sim = nx * h_sim  # = nx * 1.0 = nx

    # Verify wave speed
    expected_wave_speed = np.sqrt(T_D / rho_D)

    print(f"\nPhysical Parameters:")
    print(f"  1D linear mass density ρ_1 = {rho_D:.6e} kg/m")
    print(f"  1D tension T_1 = {T_D:.6e} N")
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
    print(f"  m_point_sim = {m_sim:.6e}  (= substrate_scale)")
    print(f"  k_spring_sim = {k_sim:.6e}  (= substrate_scale)")
    print(f"  dt_sim = {dt_sim:.6e}  (time step)")
    print(f"")
    print(f"  Wave propagation:")
    print(f"    c_wave_sim = √(k_sim/m_sim) = {(k_sim/m_sim)**0.5:.6e}  (= 1.0 always)")
    print(f"    c_wave / c_light = {(k_sim/m_sim)**0.5:.6e}  (FIXED to 1.0)")
    print(f"")

    print(f"\nSimulation Configuration:")
    print(f"  Domain (physical): {nx} points × {h_phys:.3e} m = {domain_length_phys:.6e} m")
    print(f"  Domain (sim units): {nx} points × {h_sim:.1f} = {domain_length_sim:.1f}")
    print(f"  CFL number = {cfl_factor:.3f}")

    # Create components using SIMULATION UNITS
    # Auto-select best available device (GPU if available, otherwise CPU)
    if torch.cuda.is_available():
        device = torch.device('cuda')
        print(f"\n✓ Using NVIDIA GPU: {torch.cuda.get_device_name(0)}")
        dtype = torch.float64
    elif hasattr(torch.backends, 'mps') and torch.backends.mps.is_available():
        device = torch.device('mps')
        print(f"\n✓ Using Apple Silicon GPU (MPS)")
        dtype = torch.float32
        print(f"  Using float32 (MPS doesn't support float64)")
    else:
        device = torch.device('cpu')
        print(f"\n⚠ Using CPU (no GPU detected)")
        dtype = torch.float64

    state = BraneState((nx,), Dimensionality.ONE_D, device, dtype)
    state.initialize_flat_configuration(h_sim)  # Use sim spacing = 1.0

    # Store initial positions for lateral distortion tracking
    initial_positions = state.positions.clone()

    # Set fixed boundaries
    state.set_fixed_boundaries()
    print(f"\nBoundary Conditions:")
    print(f"  Fixed boundaries at x=0 and x={domain_length_phys:.3e} m")
    print(f"  Fixed points: {state.fixed_mask.sum().item()} / {nx}")

    grid = BraneGrid((nx,), Dimensionality.ONE_D, h_sim, device)  # Sim spacing = 1.0


    print(f"\nPretension Implementation (Sim Units):")
    print(f"  Rest length L_0 (sim) = {rest_length_sim:.6e}")
    print(f"  Rest length L_0 (phys) = {rest_length_phys:.6e} m")
    print(f"  Actual spacing (sim) = {h_sim:.6e}")
    print(f"  Spring constant (sim) = {k_sim:.6e}")
    print(f"  Background tension F_0 (sim) = k×(h-L_0) = {k_sim * (h_sim - rest_length_sim):.6e}")

    physics = SpringForceComputer(k_sim, rest_length_sim)
    solver = VelocityVerletSolver(dt_sim, m_sim, physics, grid)

    # Set up lateralization measurement (uses PHYSICAL mass for correct energy scale)
    lat_config = LateralizationConfig(
        amplitude_dim=3,   # ξ = X⁴
        lateral_dims=(0,), # x is the single lateral DOF in 1D brane
    )

    lateralization = LateralizationMeasurement(
        config=lat_config,
        grid=grid,
        m_point=m_point,  # Use physical mass for energy calculations
        reference_positions=initial_positions,
    )

    # Initialize wave packet IN SIMULATION UNITS
    print(f"\nInitializing photon wave packet...")

    # Physical values (wavelength already set to λ_C at configuration)
    amplitude_phys = 10 * h_phys
    center_position_phys = domain_length_phys / 3.0

    # Convert to sim units using mapper
    wavelength_sim = mapper.to_sim_length(wavelength_phys)  # = 40.0 (points per wavelength)
    amplitude_sim = mapper.to_sim_length(amplitude_phys)    # = 0.1
    center_position_sim = mapper.to_sim_length(center_position_phys)  # = nx/3

    print(f"  Physical wavelength: {wavelength_phys:.6e} m (= λ_C)")
    print(f"  Physical wavelength: {wavelength_phys/constants.lambda_C:.2f} × λ_C")
    print(f"  Sim wavelength: {wavelength_sim:.1f} grid units")
    print(f"  Physical amplitude: {amplitude_phys:.6e} m")
    print(f"  Sim amplitude: {amplitude_sim:.3f} grid units")

    # Step 1: Initialize shape only (in sim units)
    print(f"\n[1] Initializing wave shape (sim units)...")
    initialize_wave_shape_1d(state, grid, wavelength_sim, amplitude_sim, center_position_sim)

    # Step 2: Initialize velocities using time-reversal method (in sim units)
    print(f"\n[2] Initializing velocities using time-reversal method (sim units)...")
    # Note: We use c_wave_sim for initialization, which is the actual wave speed in sim units
    initialize_right_moving_velocities_time_reversed(
        state=state,
        grid=grid,
        physics=physics,
        m_point=m_sim,
        wave_speed=c_wave_sim,  # Actual wave speed in sim units (= √(k_sim/m_sim))
        field_component=3,
        shift_cells=1,
    )

    # Step 3: Compute initial accelerations
    solver.initialize_accelerations(state)
    state.apply_fixed_boundaries()

    # Initial measurements (in sim units, but energy uses physical mass)
    initial_energy = solver.compute_energy(state)
    initial_center_sim = track_wave_center(state, grid)  # Sim units
    initial_center_phys = mapper.to_phys_length(initial_center_sim)

    print(f"\nInitial State:")
    print(f"  Energy = {initial_energy['total']:.6e} J")
    print(f"  Wave center (sim) = {initial_center_sim:.3f} grid units")
    print(f"  Wave center (phys) = {initial_center_phys:.6e} m")

    # Run simulation - fixed number of steps
    num_steps = 20000
    simulation_time_sim = num_steps * dt_sim
    simulation_time_phys = mapper.to_phys_time(simulation_time_sim)
    crossing_time_phys = domain_length_phys / constants.c

    print(f"\nRunning simulation...")
    print(f"  Light crossing time (phys) = {crossing_time_phys:.6e} s")
    print(f"  Simulation time (phys) = {simulation_time_phys:.6e} s")
    print(f"  Simulation time (sim) = {simulation_time_sim:.3f} time units")
    print(f"  Number of steps = {num_steps:,}")

    # Tracking
    times_phys = []  # Physical times for plotting
    centers_sim = []  # Centers in sim units
    energies = []

    # Take snapshots at regular intervals (in physical time)
    num_snapshots = 7
    snapshot_times_phys = np.linspace(0, simulation_time_phys, num_snapshots)
    snapshots = {}
    snapshots_lateral = {}  # Store lateral displacements
    snapshots_vel_amplitude = {}  # Store amplitude velocities
    snapshots_vel_lateral = {}  # Store lateral velocities
    snapshots_lateralization = {}  # Store lateralization ratio
    lateralization_global_history = []  # Store global R_lat over time
    snapshot_steps = {int(t / dt_phys): t for t in snapshot_times_phys}

    # Add extra snapshot at step 1 (right after first time step)
    snapshot_steps[1] = dt_phys
    num_snapshots += 1  # Now 8 total snapshots

    # Physical coordinates for plotting (convert from sim to physical)
    x_coords_phys = np.arange(nx) * h_phys

    print_interval = max(1, num_steps // 20)

    for step in range(num_steps + 1):
        if step in snapshot_steps:
            field = state.get_field_component(3).cpu().numpy()
            snapshots[snapshot_steps[step]] = field.copy()

            # Store lateral displacement (x-component)
            lateral_disp = (state.positions[:, 0] - initial_positions[:, 0]).cpu().numpy()
            snapshots_lateral[snapshot_steps[step]] = lateral_disp.copy()

            # Store amplitude velocity (component 3 = X⁴)
            vel_amplitude = state.velocities[:, 3].cpu().numpy()
            snapshots_vel_amplitude[snapshot_steps[step]] = vel_amplitude.copy()

            # Store lateral velocity (component 0 = X)
            vel_lateral = state.velocities[:, 0].cpu().numpy()
            snapshots_vel_lateral[snapshot_steps[step]] = vel_lateral.copy()

            # Measure lateralization ratio
            R_lat_local, R_lat_global, diagnostics = lateralization.measure(state, physics)
            snapshots_lateralization[snapshot_steps[step]] = R_lat_local.cpu().numpy().copy()

            # Export CSV snapshot (converted to SI units)
            csv_filename = f'photon_1d_snapshot_t{step:06d}.csv'
            export_csv_snapshot(csv_filename, state, grid, initial_positions,
                              k_sim, lateralization, physics, h_sim,
                              mapper, rest_length=rest_length_sim)
            if step == 0:
                print(f"  ✓ Exporting CSV snapshots (8 total: t=0, t=1, and 6 regular intervals) in SI units")

            # Debug: print first measurement
            if step == 0:
                print(f"\n[DEBUG] First lateralization measurement:")
                print(f"  R_lat_global: {R_lat_global:.6f}")
                print(f"  R_lat_local range: [{R_lat_local.min().item():.6f}, {R_lat_local.max().item():.6f}]")
                print(f"  E_amp_kin total: {diagnostics['E_amp_kin'].sum().item():.6e}")
                print(f"  E_lat_kin total: {diagnostics['E_lat_kin'].sum().item():.6e}")
                print(f"  E_amp_pot total: {diagnostics['E_amp_pot'].sum().item():.6e}")
                print(f"  E_lat_pot total: {diagnostics['E_lat_pot'].sum().item():.6e}")

        if step % max(1, num_steps // 100) == 0:  # Track 100 points
            center_sim = track_wave_center(state, grid)  # In sim units
            energy = solver.compute_energy(state)
            R_lat_local, R_lat_global, diagnostics = lateralization.measure(state, physics)
            # Convert solver.time (sim units) to physical time
            time_phys = mapper.to_phys_time(solver.time)
            times_phys.append(time_phys)
            centers_sim.append(center_sim)
            energies.append(energy['total'])
            lateralization_global_history.append(R_lat_global)

        if step % print_interval == 0:
            # Convert time and center to physical for printing
            time_phys = mapper.to_phys_time(solver.time)
            center_phys = mapper.to_phys_length(center_sim) if 'center_sim' in locals() else initial_center_phys
            energy_val = energy['total'] if 'energy' in locals() else initial_energy['total']
            print(f"  Step {step:8d}/{num_steps}: t={time_phys:.6e}s, "
                  f"center={center_phys:.3e}m, E={energy_val:.6e}J")

        if step < num_steps:
            solver.step(state)

    # Final analysis
    final_energy = solver.compute_energy(state)
    final_center_sim = track_wave_center(state, grid)  # Sim units
    final_center_phys = mapper.to_phys_length(final_center_sim)

    distance_traveled_sim = final_center_sim - initial_center_sim  # Sim units
    distance_traveled_phys = mapper.to_phys_length(distance_traveled_sim)
    measured_speed = distance_traveled_phys / simulation_time_phys  # Physical units
    speed_error = abs(measured_speed - constants.c) / constants.c

    energy_drift = abs(final_energy['total'] - initial_energy['total']) / initial_energy['total']

    print(f"\n{'=' * 70}")
    print("Results:")
    print(f"{'=' * 70}")

    # Manual wave propagation verification (already converted to physical units)
    expected_displacement = constants.c * simulation_time_phys
    relative_error = abs(measured_speed - constants.c) / constants.c

    print(f"\n--- Wave Propagation Verification ---")
    print(f"Time elapsed: {simulation_time_phys:.6e} s")
    print(f"Initial center: {initial_center_phys:.6e} m")
    print(f"Current center: {final_center_phys:.6e} m")
    print(f"Displacement: {distance_traveled_phys:.6e} m")
    print(f"Expected displacement: {expected_displacement:.6e} m")
    print(f"Expected speed: {constants.c:.6e} m/s")
    print(f"Measured speed: {measured_speed:.6e} m/s")
    print(f"Relative error: {relative_error * 100:.4f}%")

    if relative_error < 0.01:
        print("✓ Wave speed matches c within 1%")
    elif relative_error < 0.05:
        print("⚠ Wave speed matches c within 5%")
    else:
        print(f"✗ Wave speed error {relative_error*100:.1f}% exceeds tolerance")

    print(f"\nEnergy Conservation:")
    print(f"  Initial: {initial_energy['total']:.6e} J")
    print(f"  Final:   {final_energy['total']:.6e} J")
    print(f"  Drift:   {energy_drift:.6e} ({energy_drift*100:.6f}%)")

    if energy_drift < 1e-4:
        print(f"  ✓ Good energy conservation")

    # Visualization
    print(f"\nCreating plots...")

    fig, axes = plt.subplots(num_snapshots, 1, figsize=(14, 12))
    fig.suptitle(f'1D Photon at Realistic Scales (c = {constants.c:.3e} m/s)',
                 fontsize=16, fontweight='bold')

    for idx, (step, t) in enumerate(snapshot_steps.items()):
        if t in snapshots:
            field = snapshots[t]

            # Convert to nanometers for better readability (field in sim units, x_coords_phys in m)
            x_nm = x_coords_phys * 1e9
            field_nm = mapper.to_phys_length(field) * 1e9  # Convert sim → phys → nm

            axes[idx].plot(x_nm, field_nm, 'b-', linewidth=2)
            axes[idx].plot([x_nm[0], x_nm[-1]], [0, 0], 'ro',
                          markersize=8, label='Fixed boundaries' if idx == 0 else '')

            axes[idx].set_ylabel('ξ [nm]', fontsize=11)
            axes[idx].set_xlim(x_nm[0], x_nm[-1])
            axes[idx].set_ylim(-1.5*amplitude_phys*1e9, 1.5*amplitude_phys*1e9)
            axes[idx].grid(True, alpha=0.3)

            # Time in femtoseconds
            t_fs = t * 1e15
            axes[idx].text(0.02, 0.95, f't = {t_fs:.3f} fs',
                          transform=axes[idx].transAxes,
                          fontsize=12, verticalalignment='top',
                          bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

            if idx == 0:
                axes[idx].legend(loc='upper right', fontsize=10)

            if idx == num_snapshots - 1:
                axes[idx].set_xlabel('Position [nm]', fontsize=12)

    plt.tight_layout()
    plt.savefig(run_manager.get_plot_path('photon_1d_example_propagation.png'), dpi=150, bbox_inches='tight')
    print(f"  ✓ Saved: photon_1d_example_propagation.png")

    # Analysis plots
    fig2, axes2 = plt.subplots(2, 1, figsize=(12, 8))

    # Position vs time (in nm and fs)
    times_fs = np.array(times_phys) * 1e15
    centers_phys = mapper.to_phys_length(np.array(centers_sim))  # Convert sim → physical
    centers_nm = centers_phys * 1e9

    axes2[0].plot(times_fs, centers_nm, 'b-', linewidth=2, label='Wave center')
    axes2[0].set_xlabel('Time [fs]', fontsize=12)
    axes2[0].set_ylabel('Wave Center [nm]', fontsize=12)
    axes2[0].set_title('Wave Propagation at Speed of Light', fontsize=14, fontweight='bold')
    axes2[0].grid(True, alpha=0.3)
    axes2[0].legend(fontsize=10)

    # Energy conservation
    energy_array = np.array(energies)
    axes2[1].plot(times_fs, energy_array / initial_energy['total'], 'g-', linewidth=2)
    axes2[1].axhline(y=1.0, color='r', linestyle='--', linewidth=1, alpha=0.5)
    axes2[1].set_xlabel('Time [fs]', fontsize=12)
    axes2[1].set_ylabel('E(t) / E(0)', fontsize=12)
    axes2[1].set_title('Energy Conservation', fontsize=14, fontweight='bold')
    axes2[1].grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(run_manager.get_plot_path('photon_1d_example_analysis.png'), dpi=150, bbox_inches='tight')
    print(f"  ✓ Saved: photon_1d_example_analysis.png")

    # Lateral distortion visualization
    print(f"\nCreating lateral distortion plots...")

    fig3, axes3 = plt.subplots(num_snapshots, 1, figsize=(14, 12))
    fig3.suptitle(f'1D Photon - Lateral Distortion (Left/Right Movement)',
                 fontsize=16, fontweight='bold')

    # Find max lateral displacement for consistent scaling (convert to physical)
    max_lateral_disp = max([mapper.to_phys_length(np.abs(snapshots_lateral[t]).max()) for t in snapshots_lateral.keys()])

    for idx, (step, t) in enumerate(snapshot_steps.items()):
        if t in snapshots_lateral:
            lateral_disp_sim = snapshots_lateral[t]
            lateral_disp_phys = mapper.to_phys_length(lateral_disp_sim)  # Convert to physical

            # Convert to picometers for better readability (lateral displacement is tiny)
            x_nm = x_coords_phys * 1e9
            lateral_disp_pm = lateral_disp_phys * 1e12

            axes3[idx].plot(x_nm, lateral_disp_pm, 'r-', linewidth=2, label='x-displacement')
            axes3[idx].axhline(y=0, color='k', linestyle='--', linewidth=0.5, alpha=0.5)
            axes3[idx].plot([x_nm[0], x_nm[-1]], [0, 0], 'ro',
                          markersize=8, label='Fixed boundaries' if idx == 0 else '')

            axes3[idx].set_ylabel('Δx [pm]', fontsize=11)
            axes3[idx].set_xlim(x_nm[0], x_nm[-1])
            axes3[idx].set_ylim(-1.5*max_lateral_disp*1e12, 1.5*max_lateral_disp*1e12)
            axes3[idx].grid(True, alpha=0.3)

            # Time in femtoseconds
            t_fs = t * 1e15
            axes3[idx].text(0.02, 0.95, f't = {t_fs:.3f} fs',
                          transform=axes3[idx].transAxes,
                          fontsize=12, verticalalignment='top',
                          bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

            if idx == 0:
                axes3[idx].legend(loc='upper right', fontsize=10)

            if idx == num_snapshots - 1:
                axes3[idx].set_xlabel('Position [nm]', fontsize=12)

    plt.tight_layout()
    plt.savefig(run_manager.get_plot_path('photon_1d_example_lateral_distortion.png'), dpi=150, bbox_inches='tight')
    print(f"  ✓ Saved: photon_1d_example_lateral_distortion.png")

    # Amplitude velocity visualization
    print(f"\nCreating amplitude velocity plots...")

    fig4, axes4 = plt.subplots(num_snapshots, 1, figsize=(14, 12))
    fig4.suptitle(f'1D Photon - Amplitude Velocity (v_ξ = ∂ξ/∂t)',
                 fontsize=16, fontweight='bold')

    # Find max amplitude velocity for consistent scaling (convert to physical)
    max_vel_amplitude = max([mapper.to_phys_velocity(np.abs(snapshots_vel_amplitude[t]).max())
                             for t in snapshots_vel_amplitude.keys()])

    for idx, (step, t) in enumerate(snapshot_steps.items()):
        if t in snapshots_vel_amplitude:
            vel_amplitude_sim = snapshots_vel_amplitude[t]
            vel_amplitude_phys = mapper.to_phys_velocity(vel_amplitude_sim)  # Convert to physical

            x_nm = x_coords_phys * 1e9

            axes4[idx].plot(x_nm, vel_amplitude_phys, 'purple', linewidth=2, label='v_ξ')
            axes4[idx].axhline(y=0, color='k', linestyle='--', linewidth=0.5, alpha=0.5)
            axes4[idx].plot([x_nm[0], x_nm[-1]], [0, 0], 'ro',
                          markersize=8, label='Fixed boundaries' if idx == 0 else '')

            axes4[idx].set_ylabel('v_ξ [m/s]', fontsize=11)
            axes4[idx].set_xlim(x_nm[0], x_nm[-1])
            axes4[idx].set_ylim(-1.5*max_vel_amplitude, 1.5*max_vel_amplitude)
            axes4[idx].grid(True, alpha=0.3)

            # Time in femtoseconds
            t_fs = t * 1e15
            axes4[idx].text(0.02, 0.95, f't = {t_fs:.3f} fs',
                          transform=axes4[idx].transAxes,
                          fontsize=12, verticalalignment='top',
                          bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

            if idx == 0:
                axes4[idx].legend(loc='upper right', fontsize=10)

            if idx == num_snapshots - 1:
                axes4[idx].set_xlabel('Position [nm]', fontsize=12)

    plt.tight_layout()
    plt.savefig(run_manager.get_plot_path('photon_1d_example_amplitude_velocity.png'), dpi=150, bbox_inches='tight')
    print(f"  ✓ Saved: photon_1d_example_amplitude_velocity.png")

    # Lateral velocity visualization
    print(f"\nCreating lateral velocity plots...")

    fig5, axes5 = plt.subplots(num_snapshots, 1, figsize=(14, 12))
    fig5.suptitle(f'1D Photon - Lateral Velocity (v_x = ∂x/∂t)',
                 fontsize=16, fontweight='bold')

    # Find max lateral velocity for consistent scaling (convert to physical)
    max_vel_lateral = max([mapper.to_phys_velocity(np.abs(snapshots_vel_lateral[t]).max())
                           for t in snapshots_vel_lateral.keys()])

    for idx, (step, t) in enumerate(snapshot_steps.items()):
        if t in snapshots_vel_lateral:
            vel_lateral_sim = snapshots_vel_lateral[t]
            vel_lateral_phys = mapper.to_phys_velocity(vel_lateral_sim)  # Convert to physical

            x_nm = x_coords_phys * 1e9

            axes5[idx].plot(x_nm, vel_lateral_phys, 'green', linewidth=2, label='v_x')
            axes5[idx].axhline(y=0, color='k', linestyle='--', linewidth=0.5, alpha=0.5)
            axes5[idx].plot([x_nm[0], x_nm[-1]], [0, 0], 'ro',
                          markersize=8, label='Fixed boundaries' if idx == 0 else '')

            axes5[idx].set_ylabel('v_x [m/s]', fontsize=11)
            axes5[idx].set_xlim(x_nm[0], x_nm[-1])
            axes5[idx].set_ylim(-1.5*max_vel_lateral, 1.5*max_vel_lateral)
            axes5[idx].grid(True, alpha=0.3)

            # Time in femtoseconds
            t_fs = t * 1e15
            axes5[idx].text(0.02, 0.95, f't = {t_fs:.3f} fs',
                          transform=axes5[idx].transAxes,
                          fontsize=12, verticalalignment='top',
                          bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

            if idx == 0:
                axes5[idx].legend(loc='upper right', fontsize=10)

            if idx == num_snapshots - 1:
                axes5[idx].set_xlabel('Position [nm]', fontsize=12)

    plt.tight_layout()
    plt.savefig(run_manager.get_plot_path('photon_1d_example_lateral_velocity.png'), dpi=150, bbox_inches='tight')
    print(f"  ✓ Saved: photon_1d_example_lateral_velocity.png")

    # Lateralization ratio visualization
    print(f"\nCreating lateralization ratio plots...")

    fig6, axes6 = plt.subplots(num_snapshots, 1, figsize=(14, 12))
    fig6.suptitle(f'1D Photon - Lateralization Ratio (R_lat = E_lat / E_total)',
                 fontsize=16, fontweight='bold')

    for idx, (step, t) in enumerate(snapshot_steps.items()):
        if t in snapshots_lateralization:
            R_lat = snapshots_lateralization[t]
            x_nm = x_coords_phys * 1e9

            axes6[idx].plot(x_nm, R_lat, 'orange', linewidth=2, label='R_lat')
            axes6[idx].axhline(y=0.5, color='k', linestyle='--', linewidth=0.5, alpha=0.5, label='R_lat=0.5' if idx == 0 else '')
            axes6[idx].set_ylabel('R_lat', fontsize=11)
            axes6[idx].set_xlim(x_nm[0], x_nm[-1])
            axes6[idx].set_ylim(0, 1)
            axes6[idx].grid(True, alpha=0.3)

            t_fs = t * 1e15
            axes6[idx].text(0.02, 0.95, f't = {t_fs:.3f} fs',
                          transform=axes6[idx].transAxes,
                          fontsize=12, verticalalignment='top',
                          bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

            if idx == 0:
                axes6[idx].legend(loc='upper right', fontsize=10)
            if idx == num_snapshots - 1:
                axes6[idx].set_xlabel('Position [nm]', fontsize=12)

    plt.tight_layout()
    plt.savefig(run_manager.get_plot_path('photon_1d_example_lateralization.png'), dpi=150, bbox_inches='tight')
    print(f"  ✓ Saved: photon_1d_example_lateralization.png")

    # Global lateralization ratio vs time
    print(f"\nCreating global lateralization ratio plot...")

    fig7, ax7 = plt.subplots(1, 1, figsize=(12, 6))
    R_lat_global_array = np.array(lateralization_global_history)
    ax7.plot(times_fs, R_lat_global_array, 'orange', linewidth=2, label='Global R_lat')
    ax7.axhline(y=0.5, color='k', linestyle='--', linewidth=1, alpha=0.5, label='R_lat=0.5')
    ax7.set_xlabel('Time [fs]', fontsize=12)
    ax7.set_ylabel('Global R_lat', fontsize=12)
    ax7.set_title('Global Lateralization Ratio vs Time', fontsize=14, fontweight='bold')
    ax7.set_ylim(0, 1)
    ax7.grid(True, alpha=0.3)
    ax7.legend(fontsize=10)

    plt.tight_layout()
    plt.savefig(run_manager.get_plot_path('photon_1d_example_lateralization_global.png'), dpi=150, bbox_inches='tight')
    print(f"  ✓ Saved: photon_1d_example_lateralization_global.png")

    print(f"\n{'=' * 70}")
    print("Simulation complete!")
    print(f"{'=' * 70}")
    print(f"\nPhysical Interpretation:")
    print(f"  Domain size: {domain_length_phys*1e9:.3f} nm = {domain_length_phys/constants.lambda_C:.0f} × λ_C")
    print(f"  Wavelength: {wavelength_phys*1e9:.3f} nm")
    print(f"  Simulation time: {simulation_time_phys*1e15:.3f} femtoseconds")
    print(f"  Distance traveled: {distance_traveled_phys*1e9:.3f} nm")
    print(f"  Number of reflections: ~{int(distance_traveled_phys / domain_length_phys)}")

    # Save configuration
    run_manager.save_config({
        "experiment": "Photon 1D Experiment",
        "device": str(device),
    })

    print(f"\n{'=' * 70}")
    print(f"All outputs saved to: {run_manager.run_dir}")
    print(f"{'=' * 70}")


if __name__ == '__main__':
    main()
