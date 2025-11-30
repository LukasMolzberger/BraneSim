"""
1D Photon: Time-Reversal vs Scalar Gradient Initialization

This example compares two methods for initializing velocities:
1. Scalar gradient method: v = -c * ∇ξ (assumes scalar wave model)
2. Time-reversal method: Uses actual brane forces to compute velocities

The time-reversal method should reduce artificial splitting because it accounts
for the full 4D geometry and lateral distortions via geometric coupling.
"""

import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import torch
import numpy as np
import matplotlib.pyplot as plt

from branesim.core.state import BraneState, Dimensionality
from branesim.core.grid import BraneGrid
from branesim.core.solver import VelocityVerletSolver
from branesim.physics.forces import SpringForceComputer
from branesim.config.simulation_config import PhysicalConstants
from branesim.physics.parameters import compton_calibrated_brane_lattice_params
from branesim.core.initial_conditions import (
    initialize_right_moving_velocities,
    initialize_right_moving_velocities_time_reversed,
    measure_wave_speed,
)


def initialize_wave_shape_1d(state, grid, wavelength, amplitude, center):
    """Initialize ONLY the shape of a 1D wave packet."""
    x = torch.arange(grid.grid_shape[0], device=state.device, dtype=state.dtype) * grid.spacing

    k = 2 * np.pi / wavelength
    sigma = 3 * wavelength / (2 * np.pi)  # Width of wave packet

    # Gaussian envelope
    envelope = amplitude * torch.exp(-((x - center) ** 2) / (2 * sigma ** 2))

    # Position field only - velocities set separately
    state.positions[:, 3] = envelope * torch.cos(k * (x - center))

    print(f"  Wavelength λ = {wavelength:.6e} m ({wavelength/grid.spacing:.1f} × h)")
    print(f"  Amplitude = {amplitude:.6e} m")


def run_simulation(method_name, use_time_reversal, params, constants):
    """Run a single simulation with specified initialization method."""
    print(f"\n{'='*70}")
    print(f"Running simulation: {method_name}")
    print(f"{'='*70}")

    # Configuration
    h = params["h"]
    nx = 400
    domain_length = nx * h
    c = constants.c

    mu = params["rho_D"]
    spring_constant = params["k_spring"]

    # CFL condition
    cfl_factor = 0.1
    dt = cfl_factor * h / c

    # Create components
    device = torch.device('cpu')
    dtype = torch.float64

    state = BraneState((nx,), Dimensionality.ONE_D, device, dtype)
    state.initialize_flat_configuration(h)

    # Store initial positions for lateral distortion tracking
    initial_positions = state.positions.clone()

    # Set fixed boundaries
    state.set_fixed_boundaries()

    grid = BraneGrid((nx,), Dimensionality.ONE_D, h, device)

    # Pretension: rest_length = 0 for κ = ρc²
    rest_length = 0.0
    physics = SpringForceComputer(spring_constant, rest_length)
    solver = VelocityVerletSolver(dt, mu, physics, grid)

    # Initialize wave shape
    wavelength = 40 * h
    amplitude = 0.1 * h
    center_position = domain_length / 3.0

    print(f"\nInitializing wave shape...")
    initialize_wave_shape_1d(state, grid, wavelength, amplitude, center_position)

    # Initialize velocities (METHOD SELECTION)
    print(f"\nInitializing velocities using {method_name}...")
    if use_time_reversal:
        # Time-reversal method: uses actual brane forces
        initialize_right_moving_velocities_time_reversed(
            state=state,
            grid=grid,
            physics=physics,
            m_point=params["m_point"],
            wave_speed=c,
            field_component=3,
            shift_cells=1,
        )
    else:
        # Scalar gradient method: assumes scalar wave model
        initialize_right_moving_velocities(
            state=state,
            grid=grid,
            wave_speed=c,
            direction=None,
            field_component=3,
        )

    # Initialize accelerations
    solver.initialize_accelerations(state)
    state.apply_fixed_boundaries()

    # Run simulation
    crossing_time = domain_length / c
    simulation_time = 1.0 * crossing_time  # One crossing
    num_steps = int(simulation_time / dt)

    print(f"\nRunning {num_steps:,} steps...")

    # Take snapshots
    num_snapshots = 5
    snapshot_times = np.linspace(0, simulation_time, num_snapshots)
    snapshots = {}
    snapshots_lateral = {}
    snapshot_steps = {int(t / dt): t for t in snapshot_times}

    for step in range(num_steps + 1):
        if step in snapshot_steps:
            field = state.get_field_component(3).cpu().numpy()
            snapshots[snapshot_steps[step]] = field.copy()

            # Store lateral displacement (x-component)
            lateral_disp = (state.positions[:, 0] - initial_positions[:, 0]).cpu().numpy()
            snapshots_lateral[snapshot_steps[step]] = lateral_disp.copy()

        if step < num_steps:
            solver.step(state)

    # Compute final metrics
    final_energy = solver.compute_energy(state)

    print(f"✓ Simulation complete")
    print(f"  Final energy: {final_energy['total']:.6e} J")

    return {
        'snapshots': snapshots,
        'snapshots_lateral': snapshots_lateral,
        'snapshot_times': snapshot_times,
        'grid_spacing': h,
        'num_points': nx,
        'amplitude': amplitude,
    }


def main():
    """Compare scalar gradient vs time-reversal initialization."""
    print("="*70)
    print("Time-Reversal vs Scalar Gradient Initialization Comparison")
    print("="*70)

    # Physical constants
    constants = PhysicalConstants()

    # Configuration
    lambda_C_multiplier = 10.0
    h = constants.lambda_C * lambda_C_multiplier

    params = compton_calibrated_brane_lattice_params(
        grid_spacing_m=h,
        dimensionality=1,
        c=constants.c
    )
    params["h"] = h

    print(f"\nPhysical Setup:")
    print(f"  Speed of light c = {constants.c:.6e} m/s")
    print(f"  Grid spacing h = {lambda_C_multiplier:.0f} × λ_C = {h:.6e} m")

    # Run both simulations
    results_scalar = run_simulation(
        "Scalar Gradient Method",
        use_time_reversal=False,
        params=params,
        constants=constants,
    )

    results_time_reversal = run_simulation(
        "Time-Reversal Method",
        use_time_reversal=True,
        params=params,
        constants=constants,
    )

    # Create comparison plots
    print(f"\n{'='*70}")
    print("Creating comparison plots...")
    print(f"{'='*70}")

    num_snapshots = len(results_scalar['snapshot_times'])
    h = results_scalar['grid_spacing']
    nx = results_scalar['num_points']
    x_coords = np.arange(nx) * h
    x_nm = x_coords * 1e9
    amplitude = results_scalar['amplitude']

    # Amplitude comparison
    fig1, axes1 = plt.subplots(num_snapshots, 2, figsize=(16, 12))
    fig1.suptitle('Amplitude Field Comparison: Scalar Gradient vs Time-Reversal',
                  fontsize=16, fontweight='bold')

    for idx, t in enumerate(results_scalar['snapshot_times']):
        # Scalar gradient method (left column)
        field_scalar = results_scalar['snapshots'][t]
        field_nm_scalar = field_scalar * 1e9

        axes1[idx, 0].plot(x_nm, field_nm_scalar, 'b-', linewidth=2)
        axes1[idx, 0].plot([x_nm[0], x_nm[-1]], [0, 0], 'ro', markersize=6)
        axes1[idx, 0].set_ylabel('ξ [nm]', fontsize=11)
        axes1[idx, 0].set_xlim(x_nm[0], x_nm[-1])
        axes1[idx, 0].set_ylim(-1.5*amplitude*1e9, 1.5*amplitude*1e9)
        axes1[idx, 0].grid(True, alpha=0.3)

        t_fs = t * 1e15
        axes1[idx, 0].text(0.02, 0.95, f't = {t_fs:.3f} fs',
                          transform=axes1[idx, 0].transAxes,
                          fontsize=11, verticalalignment='top',
                          bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.8))

        if idx == 0:
            axes1[idx, 0].set_title('Scalar Gradient Method\n(v = -c ∇ξ)',
                                   fontsize=12, fontweight='bold')
        if idx == num_snapshots - 1:
            axes1[idx, 0].set_xlabel('Position [nm]', fontsize=12)

        # Time-reversal method (right column)
        field_tr = results_time_reversal['snapshots'][t]
        field_nm_tr = field_tr * 1e9

        axes1[idx, 1].plot(x_nm, field_nm_tr, 'r-', linewidth=2)
        axes1[idx, 1].plot([x_nm[0], x_nm[-1]], [0, 0], 'ro', markersize=6)
        axes1[idx, 1].set_ylabel('ξ [nm]', fontsize=11)
        axes1[idx, 1].set_xlim(x_nm[0], x_nm[-1])
        axes1[idx, 1].set_ylim(-1.5*amplitude*1e9, 1.5*amplitude*1e9)
        axes1[idx, 1].grid(True, alpha=0.3)

        axes1[idx, 1].text(0.02, 0.95, f't = {t_fs:.3f} fs',
                          transform=axes1[idx, 1].transAxes,
                          fontsize=11, verticalalignment='top',
                          bbox=dict(boxstyle='round', facecolor='lightcoral', alpha=0.8))

        if idx == 0:
            axes1[idx, 1].set_title('Time-Reversal Method\n(Full 4D Forces)',
                                   fontsize=12, fontweight='bold')
        if idx == num_snapshots - 1:
            axes1[idx, 1].set_xlabel('Position [nm]', fontsize=12)

    plt.tight_layout()
    plt.savefig('time_reversal_comparison_amplitude.png', dpi=150, bbox_inches='tight')
    print(f"  ✓ Saved: time_reversal_comparison_amplitude.png")

    # Lateral distortion comparison
    fig2, axes2 = plt.subplots(num_snapshots, 2, figsize=(16, 12))
    fig2.suptitle('Lateral Distortion Comparison: Scalar Gradient vs Time-Reversal',
                  fontsize=16, fontweight='bold')

    # Find max lateral displacement for consistent scaling
    max_lat_scalar = max([np.abs(results_scalar['snapshots_lateral'][t]).max()
                          for t in results_scalar['snapshot_times']])
    max_lat_tr = max([np.abs(results_time_reversal['snapshots_lateral'][t]).max()
                      for t in results_time_reversal['snapshot_times']])
    max_lateral = max(max_lat_scalar, max_lat_tr)

    for idx, t in enumerate(results_scalar['snapshot_times']):
        # Scalar gradient method (left column)
        lateral_scalar = results_scalar['snapshots_lateral'][t]
        lateral_pm_scalar = lateral_scalar * 1e12

        axes2[idx, 0].plot(x_nm, lateral_pm_scalar, 'b-', linewidth=2)
        axes2[idx, 0].axhline(y=0, color='k', linestyle='--', linewidth=0.5, alpha=0.5)
        axes2[idx, 0].plot([x_nm[0], x_nm[-1]], [0, 0], 'ro', markersize=6)
        axes2[idx, 0].set_ylabel('Δx [pm]', fontsize=11)
        axes2[idx, 0].set_xlim(x_nm[0], x_nm[-1])
        axes2[idx, 0].set_ylim(-1.5*max_lateral*1e12, 1.5*max_lateral*1e12)
        axes2[idx, 0].grid(True, alpha=0.3)

        t_fs = t * 1e15
        axes2[idx, 0].text(0.02, 0.95, f't = {t_fs:.3f} fs',
                          transform=axes2[idx, 0].transAxes,
                          fontsize=11, verticalalignment='top',
                          bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.8))

        if idx == 0:
            axes2[idx, 0].set_title('Scalar Gradient Method',
                                   fontsize=12, fontweight='bold')
        if idx == num_snapshots - 1:
            axes2[idx, 0].set_xlabel('Position [nm]', fontsize=12)

        # Time-reversal method (right column)
        lateral_tr = results_time_reversal['snapshots_lateral'][t]
        lateral_pm_tr = lateral_tr * 1e12

        axes2[idx, 1].plot(x_nm, lateral_pm_tr, 'r-', linewidth=2)
        axes2[idx, 1].axhline(y=0, color='k', linestyle='--', linewidth=0.5, alpha=0.5)
        axes2[idx, 1].plot([x_nm[0], x_nm[-1]], [0, 0], 'ro', markersize=6)
        axes2[idx, 1].set_ylabel('Δx [pm]', fontsize=11)
        axes2[idx, 1].set_xlim(x_nm[0], x_nm[-1])
        axes2[idx, 1].set_ylim(-1.5*max_lateral*1e12, 1.5*max_lateral*1e12)
        axes2[idx, 1].grid(True, alpha=0.3)

        axes2[idx, 1].text(0.02, 0.95, f't = {t_fs:.3f} fs',
                          transform=axes2[idx, 1].transAxes,
                          fontsize=11, verticalalignment='top',
                          bbox=dict(boxstyle='round', facecolor='lightcoral', alpha=0.8))

        if idx == 0:
            axes2[idx, 1].set_title('Time-Reversal Method',
                                   fontsize=12, fontweight='bold')
        if idx == num_snapshots - 1:
            axes2[idx, 1].set_xlabel('Position [nm]', fontsize=12)

    plt.tight_layout()
    plt.savefig('time_reversal_comparison_lateral.png', dpi=150, bbox_inches='tight')
    print(f"  ✓ Saved: time_reversal_comparison_lateral.png")

    # Difference plot (final snapshot)
    fig3, axes3 = plt.subplots(2, 1, figsize=(14, 8))
    fig3.suptitle('Difference: Time-Reversal - Scalar Gradient (Final Snapshot)',
                  fontsize=16, fontweight='bold')

    t_final = results_scalar['snapshot_times'][-1]

    # Amplitude difference
    field_diff = (results_time_reversal['snapshots'][t_final] -
                  results_scalar['snapshots'][t_final])
    field_diff_nm = field_diff * 1e9

    axes3[0].plot(x_nm, field_diff_nm, 'g-', linewidth=2)
    axes3[0].axhline(y=0, color='k', linestyle='--', linewidth=1, alpha=0.5)
    axes3[0].set_ylabel('Δξ [nm]', fontsize=12)
    axes3[0].set_xlabel('Position [nm]', fontsize=12)
    axes3[0].set_title('Amplitude Difference', fontsize=13, fontweight='bold')
    axes3[0].grid(True, alpha=0.3)

    # Lateral displacement difference
    lateral_diff = (results_time_reversal['snapshots_lateral'][t_final] -
                    results_scalar['snapshots_lateral'][t_final])
    lateral_diff_pm = lateral_diff * 1e12

    axes3[1].plot(x_nm, lateral_diff_pm, 'purple', linewidth=2)
    axes3[1].axhline(y=0, color='k', linestyle='--', linewidth=1, alpha=0.5)
    axes3[1].set_ylabel('ΔΔx [pm]', fontsize=12)
    axes3[1].set_xlabel('Position [nm]', fontsize=12)
    axes3[1].set_title('Lateral Displacement Difference', fontsize=13, fontweight='bold')
    axes3[1].grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig('time_reversal_comparison_difference.png', dpi=150, bbox_inches='tight')
    print(f"  ✓ Saved: time_reversal_comparison_difference.png")

    print(f"\n{'='*70}")
    print("Analysis complete!")
    print(f"{'='*70}")
    print("\nKey Observations:")
    print("  - Scalar gradient method: Simple v = -c ∇ξ (scalar wave assumption)")
    print("  - Time-reversal method: Uses full 4D brane forces with geometric coupling")
    print("  - Differences highlight the impact of accounting for lateral distortion")
    print("  - Time-reversal should reduce artificial splitting artifacts")


if __name__ == '__main__':
    main()