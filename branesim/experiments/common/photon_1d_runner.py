"""
Common runner for 1D photon simulations.

Eliminates code duplication between photon_1d_experiment.py and
photon_1d_berry_phase_experiment.py by providing a shared simulation runner.
"""

from __future__ import annotations
from dataclasses import dataclass, field
import numpy as np
import torch

from branesim.core.state import BraneState, Dimensionality
from branesim.core.grid import BraneGrid
from branesim.core.solver import VelocityVerletSolver
from branesim.core.dimensions import MassModel
from branesim.core.initial_conditions import (
    initialize_wave_shape_1d,
    initialize_right_moving_velocities_time_reversed,
)
from branesim.physics.forces import SpringForceComputer
from branesim.config.physical_constants import PhysicalConstants
from branesim.physics.dimensional_mapping import DimensionalMapper
from branesim.diagnostics.lateralization import (
    LateralizationMeasurement,
    LateralizationConfig,
)
from branesim.io import export_photon_1d_snapshot_csv
from branesim.utils import TestRunManager


@dataclass
class Photon1DConfig:
    """
    Configuration for 1D photon simulation experiments.

    Attributes
    ----------
    num_steps : int
        Number of simulation time steps
    num_snapshots : int
        Number of snapshot times (evenly spaced)
    add_snapshot_step1 : bool
        If True, add an extra snapshot at step=1 (right after first time step)
    points_per_wavelength : int
        Grid resolution (points per wavelength)
    num_wavelengths : int
        Domain size in wavelengths
    amplitude_factor_h : float
        Amplitude as multiple of grid spacing: A = amplitude_factor_h * h_phys
    center_fraction : float
        Wave packet center as fraction of domain length
    export_csv_snapshots : bool
        If True, export CSV snapshots at snapshot times
    experiment_name : str
        Name for the test run directory
    cfl_factor : float
        CFL number for time step calculation
    """
    num_steps: int = 20000
    num_snapshots: int = 7
    add_snapshot_step1: bool = True
    points_per_wavelength: int = 20
    num_wavelengths: int = 100
    amplitude_factor_h: float = 10.0
    center_fraction: float = 1.0 / 3.0
    export_csv_snapshots: bool = True
    experiment_name: str = "photon_1d_experiment"
    cfl_factor: float = 0.1


@dataclass
class Photon1DRunData:
    """
    Complete data from a 1D photon simulation run.

    Contains all simulation objects, parameters, and collected data.
    """
    # Management
    run_manager: TestRunManager
    constants: PhysicalConstants
    mapper: DimensionalMapper

    # Simulation objects
    grid: BraneGrid
    state: BraneState
    solver: VelocityVerletSolver
    physics: SpringForceComputer
    lateralization: LateralizationMeasurement
    initial_positions: torch.Tensor

    # Core scales (both sim and phys)
    h_sim: float
    dt_sim: float
    h_phys: float
    dt_phys: float
    wavelength_phys: float
    wavelength_sim: float
    c_wave_sim: float
    omega_sim: float  # Carrier frequency 2π·c/λ in sim units

    # Domain info
    nx: int
    domain_length_phys: float

    # Snapshot axes
    x_coords_phys_m: np.ndarray  # Position coordinates [m], shape [N]
    snapshot_times_phys_s: list[float]  # Snapshot times [s]
    snapshot_steps: dict[int, float]  # step -> time_phys mapping

    # Snapshots (values in simulation units, as numpy arrays)
    snapshots_xi: dict[float, np.ndarray]       # ξ(t) [N] (scalar, backward compat)
    snapshots_delta_x: dict[float, np.ndarray]  # (x - x0)(t) [N] (scalar, backward compat)
    snapshots_v_xi: dict[float, np.ndarray]     # vξ(t) [N] (scalar, backward compat)
    snapshots_v_x: dict[float, np.ndarray]      # vx(t) [N] (scalar, backward compat)
    snapshots_R_lat: dict[float, np.ndarray]    # R_lat(t) [N] (dimensionless)
    snapshots_xi_vec: dict[float, np.ndarray]   # Full displacement vector [N, 4]
    snapshots_v_vec: dict[float, np.ndarray]    # Full velocity vector [N, 4]

    # Tracking (values in appropriate units)
    times_phys_track_s: list[float] = field(default_factory=list)
    centers_sim_track: list[float] = field(default_factory=list)
    energies_track_J: list[float] = field(default_factory=list)
    R_lat_global_track: list[float] = field(default_factory=list)


def _track_wave_center(state: BraneState, grid: BraneGrid, field_component: int = 3) -> float:
    """Track center of wave energy in simulation units."""
    energy_density = state.velocities[:, field_component] ** 2 + state.positions[:, field_component] ** 2
    total = energy_density.sum()

    if total > 1e-10:
        x_coords = torch.arange(len(energy_density), device=energy_density.device,
                               dtype=energy_density.dtype)
        center = (x_coords * energy_density).sum() / total
        return center.item() * grid.spacing
    return 0.0


def build_photon_1d_simulation(cfg: Photon1DConfig) -> Photon1DRunData:
    """
    Build and initialize a 1D photon simulation.

    Sets up all simulation objects, initializes wave packet, and prepares
    coordinate arrays and snapshot times.

    Parameters
    ----------
    cfg : Photon1DConfig
        Configuration parameters

    Returns
    -------
    Photon1DRunData
        Initialized simulation data container
    """
    print("=" * 70)
    print(f"1D Photon Simulation - {cfg.experiment_name}")
    print("=" * 70)

    # Initialize test run manager
    run_manager = TestRunManager(experiment_name=cfg.experiment_name)
    print(run_manager.get_summary())

    # Physical constants
    constants = PhysicalConstants()

    print(f"\nPhysical Constants:")
    print(f"  Speed of light c = {constants.c:.6e} m/s")
    print(f"  Compton wavelength λ_C = {constants.lambda_C:.6e} m")

    # Configuration
    wavelength_phys = constants.lambda_C
    h_phys = wavelength_phys / cfg.points_per_wavelength
    cfl_factor = cfg.cfl_factor

    D = 1
    m_point = 2.861821e-27  # kg (universal point mass)

    # 1D brane parameters
    rho_D = m_point / (h_phys ** D)
    T_D = rho_D * constants.c**2
    rest_length_phys = constants.rest_length_frac * h_phys
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
    nx = cfg.num_wavelengths * cfg.points_per_wavelength
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

    # Create simulation components
    state = BraneState((nx,), Dimensionality.ONE_D, device, dtype)
    state.initialize_flat_configuration(h_sim)
    initial_positions = state.positions.clone()

    # Fixed boundaries
    state.set_fixed_boundaries()
    print(f"\nBoundary Conditions:")
    print(f"  Fixed boundaries at x=0 and x={domain_length_phys:.3e} m")

    grid = BraneGrid((nx,), Dimensionality.ONE_D, h_sim, device)
    physics = SpringForceComputer(k_sim, rest_length_sim)

    # Create mass model
    # In sim units: density = m_sim / h_sim (linear density for 1D)
    rho_sim = m_sim / h_sim
    mass_model = MassModel.from_density(
        density=rho_sim,
        intrinsic_dim=1,
        spacing=h_sim,
    )
    solver = VelocityVerletSolver(dt_sim, mass_model, physics, grid)

    # Lateralization measurement
    lat_config = LateralizationConfig(
        amplitude_dim=3,
        lateral_dims=(0,),
    )
    lateralization = LateralizationMeasurement(
        config=lat_config,
        grid=grid,
        m_point=m_point,
        reference_positions=initial_positions,
    )

    # Initialize wave packet
    print(f"\nInitializing photon wave packet...")

    amplitude_phys = cfg.amplitude_factor_h * h_phys
    center_position_phys = domain_length_phys * cfg.center_fraction

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

    # Prepare coordinate arrays
    x_coords_sim = grid.get_spatial_coordinates().squeeze()  # [N]
    x_coords_phys_m = mapper.to_phys_length(x_coords_sim).cpu().numpy()  # meters

    # Prepare snapshot times
    simulation_time_phys = mapper.to_phys_time(cfg.num_steps * dt_sim)
    snapshot_times_phys = np.linspace(0, simulation_time_phys, cfg.num_snapshots)
    snapshot_steps = {int(t / dt_phys): t for t in snapshot_times_phys}

    if cfg.add_snapshot_step1:
        snapshot_steps[1] = dt_phys

    print(f"\nSimulation Setup Complete:")
    print(f"  Total steps: {cfg.num_steps:,}")
    print(f"  Snapshots: {len(snapshot_steps)} times")
    print(f"  Simulation time: {simulation_time_phys*1e15:.3f} fs")

    # Create and return run data container
    return Photon1DRunData(
        run_manager=run_manager,
        constants=constants,
        mapper=mapper,
        grid=grid,
        state=state,
        solver=solver,
        physics=physics,
        lateralization=lateralization,
        initial_positions=initial_positions,
        h_sim=h_sim,
        dt_sim=dt_sim,
        h_phys=h_phys,
        dt_phys=dt_phys,
        wavelength_phys=wavelength_phys,
        wavelength_sim=wavelength_sim,
        c_wave_sim=c_wave_sim,
        omega_sim=omega_sim,
        nx=nx,
        domain_length_phys=domain_length_phys,
        x_coords_phys_m=x_coords_phys_m,
        snapshot_times_phys_s=snapshot_times_phys.tolist(),
        snapshot_steps=snapshot_steps,
        snapshots_xi={},
        snapshots_delta_x={},
        snapshots_v_xi={},
        snapshots_v_x={},
        snapshots_R_lat={},
        snapshots_xi_vec={},
        snapshots_v_vec={},
    )


def run_photon_1d(cfg: Photon1DConfig) -> Photon1DRunData:
    """
    Run a complete 1D photon simulation.

    Builds the simulation, runs the time evolution loop, collects snapshots
    and tracking data, optionally exports CSV files.

    Parameters
    ----------
    cfg : Photon1DConfig
        Configuration parameters

    Returns
    -------
    Photon1DRunData
        Complete simulation results with all snapshots and tracking data

    Examples
    --------
    >>> from branesim.experiments.common import Photon1DConfig, run_photon_1d
    >>> cfg = Photon1DConfig(num_steps=10000, experiment_name="my_photon_test")
    >>> run = run_photon_1d(cfg)
    >>> # Access results
    >>> print(f"Final energy: {run.energies_track_J[-1]:.6e} J")
    """
    # Build simulation
    run = build_photon_1d_simulation(cfg)

    # Run simulation loop
    print(f"\nRunning simulation...")

    num_steps = cfg.num_steps
    print_interval = max(1, num_steps // 20)

    for step in range(num_steps + 1):
        if step in run.snapshot_steps:
            t_phys_s = run.snapshot_steps[step]

            # Capture snapshots (sim units)
            run.snapshots_xi[t_phys_s] = run.state.positions[:, 3].cpu().numpy().copy()
            run.snapshots_delta_x[t_phys_s] = (
                run.state.positions[:, 0] - run.initial_positions[:, 0]
            ).cpu().numpy().copy()
            run.snapshots_v_xi[t_phys_s] = run.state.velocities[:, 3].cpu().numpy().copy()
            run.snapshots_v_x[t_phys_s] = run.state.velocities[:, 0].cpu().numpy().copy()

            # Capture full 4-component vectors (displacements, not absolute positions)
            run.snapshots_xi_vec[t_phys_s] = (
                run.state.positions - run.initial_positions
            ).cpu().numpy().copy()
            run.snapshots_v_vec[t_phys_s] = run.state.velocities.cpu().numpy().copy()

            # Lateralization
            R_lat_local, R_lat_global, diagnostics = run.lateralization.measure(
                run.state, run.physics
            )
            run.snapshots_R_lat[t_phys_s] = R_lat_local.cpu().numpy().copy()

            # Optional CSV export
            if cfg.export_csv_snapshots:
                csv_filename = run.run_manager.get_data_path(
                    f'photon_1d_snapshot_t{step:06d}.csv'
                )
                export_photon_1d_snapshot_csv(
                    csv_filename,
                    run.state,
                    run.grid,
                    run.initial_positions,
                    run.physics.spring_constant,
                    run.lateralization,
                    run.physics,
                    run.h_sim,
                    run.mapper,
                    rest_length_sim=run.physics.rest_length,
                )

                if step == 0:
                    print(f"  ✓ Exporting CSV snapshots ({len(run.snapshot_steps)} total)")

        # Tracking (every 1% of simulation)
        if step % max(1, num_steps // 100) == 0:
            center_sim = _track_wave_center(run.state, run.grid)
            energy = run.solver.compute_energy(run.state)
            R_lat_local, R_lat_global, diagnostics = run.lateralization.measure(
                run.state, run.physics
            )

            time_phys = run.mapper.to_phys_time(run.solver.time)
            run.times_phys_track_s.append(time_phys)
            run.centers_sim_track.append(center_sim)
            run.energies_track_J.append(energy['total'])
            run.R_lat_global_track.append(R_lat_global)

        if step % print_interval == 0:
            time_phys = run.mapper.to_phys_time(run.solver.time)
            print(f"  Step {step:8d}/{num_steps}: t={time_phys:.6e}s")

        if step < num_steps:
            run.solver.step(run.state)

    # Final analysis
    print(f"\n{'=' * 70}")
    print("Simulation complete!")
    print(f"{'=' * 70}")

    # Save configuration
    run.run_manager.save_config({
        "experiment": cfg.experiment_name,
        "num_steps": cfg.num_steps,
        "num_snapshots": len(run.snapshot_steps),
        "wavelength_phys": run.wavelength_phys,
        "omega_sim": run.omega_sim,
    })

    print(f"\nOutputs saved to: {run.run_manager.run_dir}")

    return run