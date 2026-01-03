"""
1D Photon Simulation with Complete Diagnostics

"""

import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
import torch

from branesim.core.state import BraneState, Dimensionality
from branesim.core.grid import BraneGrid
from branesim.core.solver import VelocityVerletSolver
from branesim.core.dimensions import MassModel
from branesim.initialization.initial_conditions import (
    initialize_wave_shape_1d,
    initialize_right_moving_velocities_time_reversed,
)
from branesim.physics.forces import SpringForceComputer
from branesim.config.physical_constants import PhysicalConstants
from branesim.physics.dimensional_mapping import DimensionalMapper
from branesim.utils import TestRunManager
from branesim.report import FigureSpec, LatexReportGenerator, ReportData

# Visualization
from branesim.visualization import plot_all_brane_1d_standard
from branesim.visualization.displacement_field_viz import (
    displacement_frames_from_positions_frames,
    create_displacement_components_video_1d_in_2d,
    create_displacement_magnitude_angle_video_1d_in_2d,
)
from branesim.berry import compton_omega0, compute_berry_time_series, create_berry_videos


def track_wave_center(state: BraneState, grid: BraneGrid, field_component: int = 3) -> float:
    """Track center of wave energy in simulation units."""
    energy_density = state.velocities[:, field_component] ** 2 + state.positions[:, field_component] ** 2
    total = energy_density.sum()


    if total > 1e-10:
        x_coords = torch.arange(len(energy_density), device=energy_density.device,
                               dtype=energy_density.dtype)
        center = (x_coords * energy_density).sum() / total
        return center.item() * grid.spacing
    return 0.0


def main():
    """Run 1D photon simulation with complete diagnostics."""

    # ========================================================================
    # Configuration
    # ========================================================================

    # Simulation parameters
    num_steps = 20000
    num_snapshots = 7
    add_snapshot_step1 = True  # Add extra snapshot at step 1
    cfl_factor = 0.1

    # Domain and resolution
    points_per_wavelength = 20
    num_wavelengths = 100

    # Wave packet initialization
    amplitude_factor_h = 10.0  # Amplitude = amplitude_factor_h * h_phys
    center_fraction = 1.0 / 3.0  # Wave packet center as fraction of domain

    # Output options
    export_csv_snapshots = True
    experiment_name = "photon_1d_experiment"

    print("=" * 70)
    print(f"1D Photon Simulation with Complete Diagnostics")
    print("=" * 70)

    # ========================================================================
    # Setup
    # ========================================================================

    # Initialize test run manager
    run_manager = TestRunManager(experiment_name=experiment_name)
    print(run_manager.get_summary())

    # Physical constants
    constants = PhysicalConstants()
    print(f"\nPhysical Constants:")
    print(f"  Speed of light c = {constants.c:.6e} m/s")
    print(f"  Compton wavelength λ_C = {constants.lambda_C:.6e} m")

    # Physical parameters
    wavelength_phys = constants.lambda_C
    h_phys = wavelength_phys / points_per_wavelength
    D = 1
    m_point = 2.861821e-27  # kg (universal point mass)

    # 1D brane parameters
    rho_D = m_point / (h_phys ** D)
    T_D = rho_D * constants.c**2
    rest_length_phys = constants.pre_stretch_alpha * h_phys
    c_wave = constants.c
    k_spring = T_D * (h_phys ** (D - 2))

    # Create dimensional mapper
    mapper = DimensionalMapper(
        h_phys=h_phys,
        c_light=constants.c,
        mass_reference=m_point
    )

    # Simulation units
    h_sim = mapper.to_sim_length(h_phys)
    m_sim = mapper.to_sim_mass(m_point)
    k_sim = mapper.to_sim_spring_constant(k_spring)
    c_wave_sim = mapper.to_sim_velocity(c_wave)
    rest_length_sim = mapper.to_sim_length(rest_length_phys)

    # Time step
    dt_phys = cfl_factor * h_phys / c_wave
    dt_sim = mapper.to_sim_time(dt_phys)

    # Domain size
    nx = num_wavelengths * points_per_wavelength
    domain_length_phys = nx * h_phys
    domain_length_sim = nx * h_sim

    print(f"\nPhysical Parameters:")
    print(f"  1D linear mass density ρ_1 = {rho_D:.6e} kg/m")
    print(f"  1D tension T_1 = {T_D:.6e} N")
    print(f"  Point mass m = {m_point:.6e} kg")
    print(f"  Time step dt = {dt_phys:.6e} s")

    print(f"\nSimulation Configuration:")
    print(f"  Domain: {nx} points × {h_phys:.3e} m = {domain_length_phys:.6e} m")
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

    # ========================================================================
    # Create Simulation Components
    # ========================================================================

    # Create state
    state = BraneState((nx,), Dimensionality.ONE_D, device, dtype)
    state.initialize_flat_configuration(h_sim)
    initial_positions = state.positions.clone()

    # Fixed boundaries
    state.set_fixed_boundaries()
    print(f"\nBoundary Conditions:")
    print(f"  Fixed boundaries at x=0 and x={domain_length_phys:.3e} m")

    # Create grid and physics
    grid = BraneGrid((nx,), Dimensionality.ONE_D, h_sim, device)
    physics = SpringForceComputer(k_sim, rest_length_sim)

    # Create mass model
    rho_sim = m_sim / h_sim
    mass_model = MassModel.from_density(
        density=rho_sim,
        intrinsic_dim=1,
        spacing=h_sim,
    )
    solver = VelocityVerletSolver(dt_sim, mass_model, physics, grid)

    expected_c, computed_c, relative_error = solver.verify_wave_speed()

    report = ReportData(
        title="1D Photon Simulation Report",
        experiment_name=experiment_name,
        run_name=run_manager.run_name,
        summary=(
            "1D photon wave packet initialized with a Gaussian envelope and "
            "propagated under a linear spring-mass brane model."
        ),
        metadata={
            "device": str(device),
            "dtype": str(dtype),
            "dimensionality": "1D",
        },
        parameters={
            "num_steps": num_steps,
            "num_snapshots": num_snapshots,
            "cfl_factor": cfl_factor,
            "points_per_wavelength": points_per_wavelength,
            "num_wavelengths": num_wavelengths,
            "wavelength_phys_m": wavelength_phys,
            "h_phys_m": h_phys,
            "h_sim": h_sim,
            "m_point_kg": m_point,
            "rho_1_kg_per_m": rho_D,
            "T_1_newton": T_D,
            "k_spring_N_per_m": k_spring,
            "rest_length_phys_m": rest_length_phys,
            "rest_length_sim": rest_length_sim,
            "pre_stretch_alpha": constants.pre_stretch_alpha,
            "c_wave_sim": c_wave_sim,
            "dt_phys_s": dt_phys,
            "dt_sim": dt_sim,
            "domain_length_phys_m": domain_length_phys,
            "domain_length_sim": domain_length_sim,
            "amplitude_factor_h": amplitude_factor_h,
            "center_fraction": center_fraction,
        },
        choices=[
            "Fixed boundary conditions at both ends.",
            "Velocity Verlet integrator with linear springs.",
            "Right-moving, time-reversed initialization for velocities.",
        ],
        assumptions=[
            "Continuum-to-lattice mapping uses h_sim = 1.0 as unit length.",
            "Wave speed is enforced by T = rho * c^2.",
        ],
        derived={
            "expected_wave_speed_m_per_s": expected_c,
            "computed_wave_speed_m_per_s": computed_c,
            "wave_speed_relative_error": relative_error,
        },
        symbol_map=[
            ("wavelength_phys_m", r"\(\lambda_C\)", wavelength_phys,
             "Compton wavelength used as characteristic scale."),
            ("h_phys_m", r"\(h_\star\)", h_phys,
             "Ground-state geometric spacing (paper v2: coupling/pre-stretch)."),
            ("rest_length_phys_m", r"\(\ell_0\)", rest_length_phys,
             "Spring rest length."),
            ("pre_stretch_alpha", r"\(\alpha=\ell_0/h_\star\)", constants.pre_stretch_alpha,
             "Pre-stretch parameter; alpha<1 yields uniform pre-tension."),
            ("rho_1_kg_per_m", r"\(\rho_m\)", rho_D,
             "Linear mass density in the continuum wave equation."),
            ("T_1_newton", r"\(T\)", T_D,
             "Effective tension in the linearized wave equation."),
            ("k_spring_N_per_m", r"\(k\)", k_spring,
             "Discrete spring constant in the lattice model."),
            ("c_wave_sim", r"\(c_{\mathrm{wave}}\)", c_wave_sim,
             "Wave speed used for calibration (approx. c in physical units)."),
            ("cfl_factor", r"\(\mathrm{CFL}\)", cfl_factor,
             "CFL-like factor used to set the time step."),
            ("dt_phys_s", r"\(\Delta t\)", dt_phys,
             "Physical time step used by the integrator."),
        ],
        notes=[
            "Berry diagnostics use a hardcoded Compton carrier omega0 and the complex amplitude a = sqrt(omega0) u + i v / sqrt(omega0).",
        ],
        figures=[
            FigureSpec(
                path=run_manager.get_plot_path("displacement_1d_components.mp4"),
                caption="Displacement component evolution (1D).",
            ),
            FigureSpec(
                path=run_manager.get_plot_path("displacement_1d_mag_angle.mp4"),
                caption="Displacement magnitude/angle evolution (1D).",
            ),
        ],
    )
    LatexReportGenerator().generate(report, run_manager.get_report_path("report.tex"))

    # ========================================================================
    # Initialize Wave Packet
    # ========================================================================

    print(f"\nInitializing photon wave packet...")

    amplitude_phys = amplitude_factor_h * h_phys
    center_position_phys = domain_length_phys * center_fraction

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

    # Compute carrier frequency for complex state construction
    omega_sim = 2.0 * np.pi * float(c_wave_sim) / float(wavelength_sim)

    print(f"\nWave Parameters:")
    print(f"  Carrier frequency ω_sim = {omega_sim:.6e} rad/sim-time")
    print(f"  Carrier frequency ω_phys = {mapper.to_phys_frequency(omega_sim):.6e} rad/s")

    # ========================================================================
    # Prepare Diagnostics Infrastructure
    # ========================================================================

    # Coordinate arrays
    x_coords_sim = grid.get_spatial_coordinates().squeeze()
    x_coords_phys_m = mapper.to_phys_length(x_coords_sim).cpu().numpy()
    x_nm = x_coords_phys_m * 1e9  # For plotting


    # Prepare snapshot times
    simulation_time_phys = mapper.to_phys_time(num_steps * dt_sim)
    snapshot_times_phys = np.linspace(0, simulation_time_phys, num_snapshots)
    snapshot_steps = {int(t / dt_phys): t for t in snapshot_times_phys}

    if add_snapshot_step1:
        snapshot_steps[1] = dt_phys

    print(f"\nSimulation Setup Complete:")
    print(f"  Total steps: {num_steps:,}")
    print(f"  Snapshots: {len(snapshot_steps)} times")
    print(f"  Simulation time: {simulation_time_phys*1e18:.3f} as")

    # ========================================================================
    # Storage for Results
    # ========================================================================

    # Snapshots (values in simulation units)
    snapshots_xi = {}
    snapshots_v_xi = {}
    snapshots_delta_x = {}
    snapshots_v_x = {}

    # Displacement field video frames (collect full position snapshots)
    frames_X_full = []  # Full 4D positions for each snapshot
    frames_V_full = []  # Full 4D velocities for each snapshot (physical units)
    frames_times_s = []  # Physical times in seconds
    X0_full = initial_positions.detach().cpu().numpy()  # Reference positions (N,4)

    # Tracking
    times_phys_track_s = []
    centers_sim_track = []
    energies_track_J = []

    # ========================================================================
    # Run Simulation
    # ========================================================================

    print(f"\nRunning simulation...")

    print_interval = max(1, num_steps // 20)

    for step in range(num_steps + 1):
        # Capture snapshots
        if step in snapshot_steps:
            t_phys_s = snapshot_steps[step]

            # Store fields (sim units)
            snapshots_xi[t_phys_s] = state.positions[:, 3].cpu().numpy().copy()
            snapshots_v_xi[t_phys_s] = state.velocities[:, 3].cpu().numpy().copy()
            snapshots_delta_x[t_phys_s] = (
                state.positions[:, 0] - initial_positions[:, 0]
            ).cpu().numpy().copy()
            snapshots_v_x[t_phys_s] = state.velocities[:, 0].cpu().numpy().copy()

            # Collect full positions for displacement field videos
            X_phys = mapper.to_phys_length(state.positions).detach().cpu().numpy()
            frames_X_full.append(X_phys)
            frames_times_s.append(t_phys_s)

            V_phys = mapper.to_phys_velocity(state.velocities).detach().cpu().numpy()
            frames_V_full.append(V_phys)

        # Tracking (every 1% of simulation)
        if step % max(1, num_steps // 100) == 0:
            center_sim = track_wave_center(state, grid)
            energy = solver.compute_energy(state)

            time_phys = mapper.to_phys_time(solver.time)
            times_phys_track_s.append(time_phys)
            centers_sim_track.append(center_sim)
            energies_track_J.append(energy['total'])

        if step % print_interval == 0:
            time_phys = mapper.to_phys_time(solver.time)
            print(f"  Step {step:8d}/{num_steps}: t={time_phys:.6e}s")

        if step < num_steps:
            solver.step(state)

    # ========================================================================
    # Visualization
    # ========================================================================

    print(f"\n{'=' * 70}")
    print("Creating Visualizations")
    print(f"{'=' * 70}")

    # Standard brane plots (existing infrastructure)
    # Reconstruct minimal run data for compatibility
    class RunData:
        pass

    run_data = RunData()
    run_data.run_manager = run_manager
    run_data.snapshots_xi = snapshots_xi
    run_data.snapshots_delta_x = snapshots_delta_x
    run_data.snapshots_v_xi = snapshots_v_xi
    run_data.snapshots_v_x = snapshots_v_x
    run_data.x_coords_phys_m = x_coords_phys_m
    run_data.snapshot_times_phys_s = sorted(snapshots_xi.keys())
    run_data.times_phys_track_s = times_phys_track_s
    run_data.centers_sim_track = centers_sim_track
    run_data.energies_track_J = energies_track_J
    run_data.h_sim = h_sim
    run_data.mapper = mapper
    run_data.constants = constants

    plot_all_brane_1d_standard(run_data)
    print(f"  ✓ Standard brane plots")

    # Displacement field videos
    if frames_X_full:
        print(f"\nGenerating displacement field videos...")

        # Convert reference positions to physical units
        X0_phys = mapper.to_phys_length(torch.from_numpy(X0_full)).numpy()

        # Compute displacement frames: u = X - X0
        frames_u_full = displacement_frames_from_positions_frames(frames_X_full, X0_phys)

        # Extract 2D embedding (components 0 and 3: longitudinal and transverse)
        frames_u_2d = [u[:, [0, 3]] for u in frames_u_full]

        # Determine displacement scale (convert to pm for display)
        all_u = np.concatenate(frames_u_2d, axis=0)
        u_max_m = float(np.max(np.abs(all_u)))
        u_max_pm = u_max_m * 1e12

        print(f"  Max displacement: {u_max_pm:.3f} pm")

        # Video 1: Raw components u^0(x,t) and u^3(x,t)
        print(f"  Creating displacement components video...")
        create_displacement_components_video_1d_in_2d(
            frames_u=frames_u_2d,
            times=frames_times_s,
            x_coords=x_nm,  # x in nm
            output_path=run_manager.get_plot_path("displacement_1d_components.mp4"),
            fps=20,
            dpi=140,
            unit_label_x="nm",
            unit_label_u="m",  # Will show in meters (scientific notation)
            title="1D Brane Displacement: Components u⁰ (longitudinal) and u³ (transverse)",
        )
        print(f"    ✓ displacement_1d_components.mp4")

        # Video 2: Magnitude and angle split
        print(f"  Creating magnitude/angle video...")
        create_displacement_magnitude_angle_video_1d_in_2d(
            frames_u=frames_u_2d,
            times=frames_times_s,
            x_coords=x_nm,
            output_path=run_manager.get_plot_path("displacement_1d_mag_angle.mp4"),
            fps=20,
            dpi=140,
            unit_label_x="nm",
            unit_label_a="m",
            title="1D Brane Displacement: Magnitude |u| and Angle θ",
            unwrap_along_x=True,
        )
        print(f"    ✓ displacement_1d_mag_angle.mp4")

        print(f"  ✓ Displacement field videos generated")

        # ====================================================================
        # Berry diagnostics (U(1) phase + connection)
        # ====================================================================

        print(f"\nGenerating Berry diagnostics and videos...")

        # Displacement and velocity frames (physical units)
        frames_u_full = displacement_frames_from_positions_frames(frames_X_full, X0_phys)
        frames_v_full = [v.copy() for v in frames_V_full]  # velocities already in physical units

        # Restrict to the 2D embedding used by the 1D experiment (components 0 and 3)
        frames_u_2d = [u[:, [0, 3]] for u in frames_u_full]
        frames_v_2d = [v[:, [0, 3]] for v in frames_v_full]

        # Hardcoded Compton carrier (as requested)
        omega0 = compton_omega0()
        berry = compute_berry_time_series(
            frames_u=frames_u_2d,
            frames_v=frames_v_2d,
            times_s=frames_times_s,
            omega0=omega0,
            amp_eps_rel=1e-4,
            overlap_eps_rel=1e-3,
            alpha_gamma=0.6,
            alpha_scale=1.0,
            return_psi=False,
        )

        create_berry_videos(
            series=berry,
            intrinsic_dim=1,
            output_dir=run_manager.plots_dir,
            x_coords_1d=x_nm,  # already computed for displacement videos
            filename_prefix="berry_1d",
            fps=20,
            dpi=140,
            unit_label="nm",
        )

        print(f"  ✓ Berry videos generated")

    # ========================================================================
    # Summary
    # ========================================================================

    print(f"\n{'=' * 70}")
    print("Simulation Complete!")
    print(f"{'=' * 70}")

    print(f"\nPhysical Interpretation:")
    print(f"  Domain size: {domain_length_phys*1e9:.3f} nm ({num_wavelengths} wavelengths)")
    print(f"  Wavelength: {wavelength_phys*1e9:.3f} nm (Compton)")
    print(f"  Simulation time: {simulation_time_phys*1e18:.3f} as")
    print(f"  Snapshots: {len(snapshots_xi)} times")


    print(f"\nEnergy Conservation:")
    if energies_track_J:
        E_initial = energies_track_J[0]
        E_final = energies_track_J[-1]
        drift = abs(E_final - E_initial) / E_initial * 100
        print(f"  Initial: {E_initial:.6e} J")
        print(f"  Final:   {E_final:.6e} J")
        print(f"  Drift:   {drift:.3f}%")

    # Save configuration
    run_manager.save_config({
        "experiment": experiment_name,
        "num_steps": num_steps,
        "num_snapshots": len(snapshots_xi),
        "wavelength_phys_m": wavelength_phys,
        "omega_sim": omega_sim,
    })

    print(f"\n{'=' * 70}")
    print(f"All outputs saved to: {run_manager.run_dir}")
    print(f"{'=' * 70}")


if __name__ == '__main__':
    main()
