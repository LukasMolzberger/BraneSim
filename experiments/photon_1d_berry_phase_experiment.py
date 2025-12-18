"""
1D Photon with Berry Phase Analysis

Extends the standard photon_1d_experiment.py with Berry phase computation.
Computes the discrete Berry phase profile γ(x) along the 1D chain at snapshot times.

Mathematical Background:
    For normalized states |u_i⟩ at lattice points i, the discrete Berry phase
    increment between neighbors is:
        Δφ_i = arg⟨u_i|u_{i+1}⟩

    The cumulative Berry phase profile is:
        γ_0 = 0,  γ_{i+1} = γ_i + Δφ_i

    The Berry connection (gauge field) is:
        A_x[i] = Δφ_i / h

    where h is the grid spacing.

Complex State Construction:
    The simulation fields are real (ξ = X^4 displacement and ξ̇ velocity).
    To extract phase information, we construct the analytic signal:
        ψ = ξ + i·ξ̇/ω

    where ω = 2π·c/λ is the carrier frequency of the photon.
"""

import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import math
import torch
import numpy as np
import matplotlib.pyplot as plt

from branesim.core.state import BraneState, Dimensionality
from branesim.core.grid import BraneGrid
from branesim.core.solver import VelocityVerletSolver
from branesim.physics.forces import SpringForceComputer
from branesim.config.physical_constants import PhysicalConstants
from branesim.physics.dimensional_mapping import DimensionalMapper
from branesim.core.initial_conditions import (
    initialize_right_moving_velocities_time_reversed,
)
from branesim.analysis import (
    complex_band_state_from_quadrature,
    pointwise_normalize,
    BerryPhase1DConfig,
    berry_phase_profile_along_x,
)
from branesim.visualization import (
    plot_berry_phase_profiles,
    plot_berry_connection_profiles,
)
from branesim.io import export_berry_phase_csv
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


def main():
    """Run 1D photon simulation with Berry phase analysis."""
    print("=" * 70)
    print("1D Photon - Berry Phase Experiment")
    print("=" * 70)

    # Initialize test run manager
    run_manager = TestRunManager(experiment_name="photon_1d_berry_phase_experiment")
    print(run_manager.get_summary())

    # Physical constants
    constants = PhysicalConstants()

    print(f"\nPhysical Constants:")
    print(f"  Speed of light c = {constants.c:.6e} m/s")
    print(f"  Compton wavelength λ_C = {constants.lambda_C:.6e} m")
    print(f"  ℏ = {constants.hbar:.6e} J·s")
    print(f"  m_e = {constants.m_e:.6e} kg")

    # Configuration (same as photon_1d_experiment.py)
    wavelength_phys = constants.lambda_C  # Photon wavelength = λ_C
    points_per_wavelength = 20  # Grid resolution
    h_phys = wavelength_phys / points_per_wavelength  # Grid spacing
    cfl_factor = 0.1

    D = 1

    # Universal point mass
    m_point = 2.861821e-27  # kg

    # 1D brane parameters constrained to give wave speed = c
    rho_D = m_point / (h_phys ** D)  # kg/m (linear mass density)
    T_D = rho_D * constants.c**2  # N (tension)

    # Rest length
    rest_length_phys = constants.rest_length_frac * h_phys

    # Wave speed (exactly equals c by construction)
    c_wave = constants.c

    # Spring constant
    k_spring = T_D * (h_phys ** (D - 2))

    # Create dimensional mapper
    mapper = DimensionalMapper(
        h_phys=h_phys,
        c_light=constants.c,
        mass_reference=m_point
    )

    # Simulation units
    h_sim = mapper.to_sim_length(h_phys)  # = 1.0
    m_sim = mapper.to_sim_mass(m_point)
    k_sim = mapper.to_sim_spring_constant(k_spring)
    c_wave_sim = mapper.to_sim_velocity(c_wave)  # = 1.0
    rest_length_sim = mapper.to_sim_length(rest_length_phys)

    # Time step (CFL condition)
    dt_phys = cfl_factor * h_phys / c_wave
    dt_sim = mapper.to_sim_time(dt_phys)

    # Domain size
    nx = 2000
    domain_length_phys = nx * h_phys
    domain_length_sim = nx * h_sim

    print(f"\nPhysical Parameters:")
    print(f"  1D linear mass density ρ_1 = {rho_D:.6e} kg/m")
    print(f"  1D tension T_1 = {T_D:.6e} N")
    print(f"  Spring constant k = {k_spring:.6e} N/m")
    print(f"  Point mass m = {m_point:.6e} kg")
    print(f"  Time step dt = {dt_phys:.6e} s")

    print(f"\nDimensionless Simulation Parameters:")
    print(f"  h_sim = {h_sim:.6e}")
    print(f"  c_wave_sim = {c_wave_sim:.6e}")
    print(f"  m_sim = {m_sim:.6e}")
    print(f"  k_sim = {k_sim:.6e}")
    print(f"  dt_sim = {dt_sim:.6e}")

    print(f"\nSimulation Configuration:")
    print(f"  Domain (physical): {nx} points × {h_phys:.3e} m = {domain_length_phys:.6e} m")
    print(f"  Domain (sim units): {nx} points × {h_sim:.1f} = {domain_length_sim:.1f}")
    print(f"  CFL number = {cfl_factor:.3f}")

    # Auto-select device
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

    # Create simulation components
    state = BraneState((nx,), Dimensionality.ONE_D, device, dtype)
    state.initialize_flat_configuration(h_sim)
    initial_positions = state.positions.clone()

    # Set fixed boundaries
    state.set_fixed_boundaries()
    print(f"\nBoundary Conditions:")
    print(f"  Fixed boundaries at x=0 and x={domain_length_phys:.3e} m")

    grid = BraneGrid((nx,), Dimensionality.ONE_D, h_sim, device)
    physics = SpringForceComputer(k_sim, rest_length_sim)
    solver = VelocityVerletSolver(dt_sim, m_sim, physics, grid)

    # Initialize wave packet
    print(f"\nInitializing photon wave packet...")

    amplitude_phys = 10 * h_phys
    center_position_phys = domain_length_phys / 3.0

    wavelength_sim = mapper.to_sim_length(wavelength_phys)
    amplitude_sim = mapper.to_sim_length(amplitude_phys)
    center_position_sim = mapper.to_sim_length(center_position_phys)

    print(f"  Physical wavelength: {wavelength_phys:.6e} m (= λ_C)")
    print(f"  Sim wavelength: {wavelength_sim:.1f} grid units")

    # Step 1: Initialize shape
    print(f"\n[1] Initializing wave shape...")
    initialize_wave_shape_1d(state, grid, wavelength_sim, amplitude_sim, center_position_sim)

    # Step 2: Initialize velocities using time-reversal
    print(f"\n[2] Initializing velocities using time-reversal method...")
    initialize_right_moving_velocities_time_reversed(
        state=state,
        grid=grid,
        physics=physics,
        m_point=m_sim,
        wave_speed=c_wave_sim,
        field_component=3,
        shift_cells=1,
    )

    # Step 3: Initialize accelerations
    solver.initialize_accelerations(state)
    state.apply_fixed_boundaries()

    # ========================================================================
    # Berry Phase Configuration
    # ========================================================================

    # Compute carrier frequency for complex state construction
    # ω = 2π·c/λ (in sim units)
    omega_sim = 2.0 * np.pi * float(c_wave_sim) / float(wavelength_sim)

    print(f"\nBerry Phase Configuration:")
    print(f"  Carrier frequency ω_sim = {omega_sim:.6e} rad / sim-time")
    print(f"  Grid spacing h_sim = {h_sim:.6e}")

    berry_cfg = BerryPhase1DConfig(
        spacing=float(h_sim),
        amplitude_threshold=1e-6,  # Mask out low-amplitude noise
        eps=1e-12,
        unwrap=True,
        force_cpu_on_mps=True,  # MPS complex support can be spotty
    )

    # Storage for Berry phase data at snapshot times
    gamma_by_tfs = {}  # Berry phase profiles
    Ax_by_tfs = {}  # Berry connection profiles

    # Precompute position coordinates
    x_coords_sim = grid.get_spatial_coordinates().squeeze()  # [N], sim units
    x_coords_phys = mapper.to_phys_length(x_coords_sim)  # meters
    x_nm = x_coords_phys.cpu().numpy() * 1e9  # nanometers for plotting

    # Edge coordinates for Berry connection (midpoints)
    x_edges_nm = 0.5 * (x_nm[:-1] + x_nm[1:])

    # ========================================================================
    # Run simulation
    # ========================================================================

    num_steps = 20000
    simulation_time_sim = num_steps * dt_sim
    simulation_time_phys = mapper.to_phys_time(simulation_time_sim)

    print(f"\nRunning simulation...")
    print(f"  Simulation time (phys) = {simulation_time_phys:.6e} s")
    print(f"  Number of steps = {num_steps:,}")

    # Snapshots at regular intervals
    num_snapshots = 7
    snapshot_times_phys = np.linspace(0, simulation_time_phys, num_snapshots)
    snapshot_steps = {int(t / dt_phys): t for t in snapshot_times_phys}

    print_interval = max(1, num_steps // 20)

    for step in range(num_steps + 1):
        if step in snapshot_steps:
            t_phys = mapper.to_phys_time(step * dt_sim)
            t_fs = t_phys * 1e15

            # ================================================================
            # Berry Phase Computation
            # ================================================================

            # Extract real fields (sim units)
            xi = state.positions[:, 3]  # Amplitude displacement ξ
            xidot = state.velocities[:, 3]  # Amplitude velocity ξ̇

            # Build complex band state: ψ = ξ + i·ξ̇/ω
            psi = complex_band_state_from_quadrature(
                xi, xidot, omega_sim, eps=berry_cfg.eps
            )

            # Normalize pointwise: |ψ̂⟩ = |ψ⟩ / |ψ|
            psi_hat, amp = pointwise_normalize(psi, eps=berry_cfg.eps)

            # Compute Berry phase profile along x
            result = berry_phase_profile_along_x(psi_hat, amp, berry_cfg)

            # Store results (move to CPU for numpy conversion)
            gamma_by_tfs[t_fs] = result["gamma_wrapped"].detach().cpu().numpy()
            Ax_by_tfs[t_fs] = result["A_x"].detach().cpu().numpy()

            # Export to CSV
            csv_path = run_manager.get_data_path(f"berry_phase_t_{t_fs:.3f}fs.csv")
            export_berry_phase_csv(
                csv_path,
                x_coords_phys.detach().cpu().numpy(),  # Position [m]
                result["gamma_wrapped"].detach().cpu().numpy(),  # Berry phase [rad]
                result["A_x"].detach().cpu().numpy(),  # Berry connection [rad/sim-length]
                amp_m=mapper.to_phys_length(amp).detach().cpu().numpy(),  # Amplitude [m]
            )

            if step == 0:
                print(f"  ✓ Computing Berry phase at {num_snapshots} snapshot times")
                print(f"  ✓ Exporting CSV files (2 files per snapshot: points + edges)")

        if step % print_interval == 0:
            time_phys = mapper.to_phys_time(solver.time)
            print(f"  Step {step:8d}/{num_steps}: t={time_phys:.6e}s")

        if step < num_steps:
            solver.step(state)

    # ========================================================================
    # Visualization
    # ========================================================================

    print(f"\nCreating Berry phase plots...")

    times_fs = sorted(gamma_by_tfs.keys())

    # Plot 1: Berry phase profiles γ(x)
    plot_berry_phase_profiles(
        run_manager,
        x_nm,
        times_fs,
        gamma_by_tfs,
        title="1D Photon - Berry Phase Profile γ(x) (wrapped to [-π, π])",
        filename="photon_1d_berry_phase_profiles.png",
    )
    print(f"  ✓ Saved: photon_1d_berry_phase_profiles.png")

    # Plot 2: Berry connection profiles A_x(x)
    plot_berry_connection_profiles(
        run_manager,
        x_edges_nm,
        times_fs,
        Ax_by_tfs,
        title="1D Photon - Berry Connection A_x(x)",
        filename="photon_1d_berry_connection_profiles.png",
    )
    print(f"  ✓ Saved: photon_1d_berry_connection_profiles.png")

    # ========================================================================
    # Summary
    # ========================================================================

    print(f"\n{'=' * 70}")
    print("Simulation complete!")
    print(f"{'=' * 70}")
    print(f"\nPhysical Interpretation:")
    print(f"  Domain size: {domain_length_phys*1e9:.3f} nm")
    print(f"  Wavelength: {wavelength_phys*1e9:.3f} nm")
    print(f"  Simulation time: {simulation_time_phys*1e15:.3f} femtoseconds")
    print(f"  Berry phase snapshots: {len(times_fs)} times")

    print(f"\nBerry Phase Results:")
    print(f"  Carrier frequency ω = {omega_sim:.6e} rad/sim-time")
    print(f"  Amplitude threshold = {berry_cfg.amplitude_threshold:.6e}")
    print(f"  Grid spacing h = {h_sim:.6e} sim-length")

    # Save configuration
    run_manager.save_config({
        "experiment": "Photon 1D Berry Phase Experiment",
        "device": str(device),
        "wavelength_phys": wavelength_phys,
        "omega_sim": omega_sim,
        "berry_amplitude_threshold": berry_cfg.amplitude_threshold,
    })

    print(f"\n{'=' * 70}")
    print(f"All outputs saved to: {run_manager.run_dir}")
    print(f"{'=' * 70}")


if __name__ == '__main__':
    main()