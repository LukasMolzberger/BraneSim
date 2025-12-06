"""
Electron Stability Test Experiment

Tests the initialization and short-term stability of the W&vdM toroidal
double-loop electron model.

This experiment:
1. Initializes a 3D brane grid at Compton scale
2. Creates an electron using calibrated parameters
3. Runs a short simulation (a few Compton periods)
4. Measures stability metrics (energy leakage, shape drift, mode purity)

The goal is to verify that the analytical ansatz produces a reasonable
initial configuration that doesn't immediately blow up or disperse.
"""

import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import math
import torch
import matplotlib.pyplot as plt
import argparse

from branesim.core.state import BraneState, Dimensionality
from branesim.core.grid import BraneGrid
from branesim.core.solver import VelocityVerletSolver
from branesim.physics.forces import SpringForceComputer
from branesim.config.simulation_config import PhysicalConstants
from branesim.physics.electron_initialization import (
    calibrate_electron_init_params,
    init_electron_state,
)
from branesim.diagnostics.electron_stability import (
    compute_electron_stability_loss,
    build_electron_masks,
    StabilityMetrics,
    compute_total_energy,
    compute_total_momentum,
    compute_total_spin,
)

# Import visualization functions
from electron_visualization import (
    visualize_initial_state,
    create_all_animations,
)


def setup_experiment(
    grid_shape=(20, 20, 20),
    N_periods=3.0,
    eta_cfl=0.2,
    amplitude_scale=1e-13,
):
    """
    Set up experiment parameters.

    Returns:
        Dictionary with all configuration parameters
    """
    # Physical constants
    constants = PhysicalConstants()

    # Grid parameters
    lambda_C = constants.lambda_C
    R_est = lambda_C / (2.0 * math.pi)  # Electron major radius

    # FIX: Keep physical box size constant, vary resolution with grid_shape
    # Box should contain electron + buffer
    # Electron fits in box of ~2R = λ_C/π ≈ 0.32λ_C
    # Add buffer for radiation: use 1.5λ_C box (tighter but still reasonable)
    box_size_phys = 1.5 * lambda_C

    # Now compute grid spacing from box size and grid dimensions
    nx, ny, nz = grid_shape
    h = box_size_phys / nx  # Grid spacing determined by resolution

    print(f"\n=== Grid Configuration ===")
    print(f"Grid size: {nx} × {ny} × {nz} = {nx*ny*nz:,} points")
    print(f"Physical box: {box_size_phys*1e12:.2f} pm = {box_size_phys/lambda_C:.1f} λ_C")
    print(f"Grid spacing h = {h:.6e} m = {h*1e12:.3f} pm")
    print(f"  h/λ_C = {h/lambda_C:.4f}")
    print(f"  λ_C/h = {lambda_C/h:.1f} (points per λ_C)")

    # Key metric: points around electron circumference
    circumference = 2 * math.pi * R_est
    points_around = circumference / h
    print(f"\n=== Electron Resolution ===")
    print(f"Electron major radius R = {R_est*1e12:.3f} pm")
    print(f"R/h = {R_est/h:.2f} grid spacings")
    print(f"Circumference = {circumference*1e12:.3f} pm")
    print(f"Points around circumference: ~{points_around:.1f}")
    if points_around < 10:
        print(f"  ⚠️  WARNING: Need at least 10-20 points for minimal resolution!")
    elif points_around < 20:
        print(f"  ⚠️  Marginal resolution")
    else:
        print(f"  ✓ Good resolution")

    # Sanity check
    if box_size_phys < 4 * R_est:
        print(f"WARNING: Grid may be too small for electron torus!")
        print(f"  Box size: {box_size_phys:.6e} m")
        print(f"  Electron diameter estimate: {2*R_est:.6e} m")

    # Brane parameters - Compton calibration
    # Point mass from density: ρ_m h³ where ρ_m is chosen to give c = √(T/ρ_m)
    # For our target c = 3×10⁸ m/s and tension T = k/h:
    # ρ_m = T/c² = k/(h·c²)
    # m_point = ρ_m h³ = k h²/c²
    k = 1e3  # Spring constant [N/m]
    rest_length_frac = 0.9  # rest_length = rest_length_frac * h

    # Compton-calibrated point mass
    m_point = k * h * h / (constants.c ** 2)

    print(f"\nBrane parameters:")
    print(f"  Spring constant k = {k:.6e} N/m")
    print(f"  Rest length = {rest_length_frac:.2f} h")
    print(f"  Point mass m_p = {m_point:.6e} kg")
    print(f"  Expected wave speed c = sqrt(k/(h·ρ_m)) = {constants.c:.6e} m/s")

    # Time stepping
    # CFL condition: dt < η h / c
    # For FCC lattice, η ≈ 0.33
    dt = eta_cfl * h / constants.c

    # Simulation duration: Run for N_periods Compton periods
    T_compton = 2 * math.pi / (constants.m_e * constants.c ** 2 / constants.hbar)
    T_total = N_periods * T_compton
    n_steps = int(T_total / dt)
    snapshot_interval = max(1, n_steps // 20)  # 20 snapshots

    print(f"\nTime stepping:")
    print(f"  dt = {dt:.6e} s ({dt/T_compton:.6f} T_C)")
    print(f"  CFL η = {eta_cfl:.2f}")
    print(f"  Compton period T_C = {T_compton:.6e} s")
    print(f"  Total time = {N_periods:.1f} T_C = {T_total:.6e} s")
    print(f"  Total steps = {n_steps}")
    print(f"  Snapshot interval = {snapshot_interval}")

    config = {
        'constants': constants,
        'grid_shape': grid_shape,
        'h': h,
        'k': k,
        'rest_length_frac': rest_length_frac,
        'm_point': m_point,
        'dt': dt,
        'n_steps': n_steps,
        'snapshot_interval': snapshot_interval,
        'N_periods': N_periods,
        'amplitude_scale': amplitude_scale,
    }

    return config


def initialize_electron(state, grid, config):
    """
    Initialize electron in the center of the grid.

    Args:
        state: BraneState
        grid: BraneGrid
        config: Configuration dictionary

    Returns:
        ElectronInitParams used for initialization
    """
    constants = config['constants']
    h = config['h']

    # Electron center: center of grid
    nx, ny, nz = config['grid_shape']
    center_grid = (nx // 2, ny // 2, nz // 2)
    center_phys = (
        center_grid[0] * h,
        center_grid[1] * h,
        center_grid[2] * h,
    )

    print(f"\n{'='*60}")
    print(f"Initializing Electron at center = {center_phys}")
    print(f"{'='*60}")

    # Get calibrated parameters
    # Use amplitude scale from config
    amplitude_scale = config['amplitude_scale']

    params = calibrate_electron_init_params(
        constants=constants,
        grid_spacing=h,
        center=center_phys,
        amplitude_scale=amplitude_scale,
    )

    # Initialize electron state
    init_electron_state(state, params)

    # Measure initial physical properties
    m_point = config['m_point']
    E_initial = compute_total_energy(state, m_point)
    P_initial = compute_total_momentum(state, m_point)
    S_initial = compute_total_spin(state, m_point, params.center)

    # Compare with expected values
    E_target = constants.m_e * constants.c ** 2
    P_target = 0.0
    S_target = 0.5 * constants.hbar

    print(f"\n=== Initial Physical Properties ===")
    print(f"  Total energy E:")
    print(f"    Measured: {E_initial:.6e} J")
    print(f"    Target (m_e c²): {E_target:.6e} J")
    print(f"    Ratio E/E_target: {E_initial/E_target:.6f}")
    print(f"\n  Total momentum |P|:")
    print(f"    Measured: {torch.norm(P_initial[:3]).item():.6e} kg·m/s")
    print(f"    Target: {P_target:.6e} kg·m/s")
    print(f"\n  Spin magnitude |S|:")
    print(f"    Measured: {torch.norm(S_initial).item():.6e} J·s")
    print(f"    Target (ℏ/2): {S_target:.6e} J·s")
    print(f"    Ratio |S|/(ℏ/2): {torch.norm(S_initial).item()/S_target:.6f}")

    return params


def run_simulation(state, grid, physics, solver, config):
    """
    Run simulation and collect snapshots.

    Args:
        state: Initial BraneState
        grid: BraneGrid
        physics: SpringForceComputer
        solver: VelocityVerletSolver
        config: Configuration dictionary

    Returns:
        List of BraneState snapshots
    """
    n_steps = config['n_steps']
    snapshot_interval = config['snapshot_interval']

    print(f"\n{'='*60}")
    print(f"Running Simulation")
    print(f"{'='*60}")

    states = [state.clone()]  # Store initial state

    for step in range(n_steps):
        solver.step(state)

        # Collect snapshot
        if (step + 1) % snapshot_interval == 0:
            states.append(state.clone())
            print(f"  Step {step+1}/{n_steps} ({100*(step+1)/n_steps:.1f}%)")

    print(f"\nSimulation complete. Collected {len(states)} snapshots.")
    return states


def analyze_stability(states, params, config):
    """
    Analyze electron stability and compute loss metrics.

    Args:
        states: List of BraneState snapshots
        params: ElectronInitParams
        config: Configuration dictionary

    Returns:
        (loss, metrics): Stability loss and detailed metrics
    """
    print(f"\n{'='*60}")
    print(f"Analyzing Stability")
    print(f"{'='*60}")

    constants = config['constants']
    dt = config['dt']
    m_point = config['m_point']

    # Build masks
    tube_mask, core_mask = build_electron_masks(
        positions=states[0].positions,
        center=params.center,
        R=params.R,
        tube_max_radius=params.tube_max_radius,
        core_radius_frac=0.5,
    )

    total_points = tube_mask.numel()
    tube_count = tube_mask.sum().item()
    core_count = core_mask.sum().item()

    print(f"\n=== Mask Statistics ===")
    print(f"  Tube points: {tube_count}/{total_points} ({100*tube_count/total_points:.2f}%)")
    print(f"  Core points: {core_count}/{total_points} ({100*core_count/total_points:.2f}%)")

    # Compute stability loss
    loss, metrics = compute_electron_stability_loss(
        states=states,
        tube_mask=tube_mask,
        core_mask=core_mask,
        constants=constants,
        m_point=m_point,
        dt=dt * config['snapshot_interval'],
        target_omega=params.compton_omega,
        center=params.center,
        spin_axis=params.spin_axis,
    )

    print(f"\n=== Stability Results ===")
    print(f"Total Loss: {loss:.6e}")
    print(metrics)

    return loss, metrics


def visualize_results(states, params, config):
    """
    Create visualizations of electron evolution.

    Args:
        states: List of BraneState snapshots
        params: ElectronInitParams
        config: Configuration dictionary
    """
    print(f"\n{'='*60}")
    print(f"Creating Visualizations")
    print(f"{'='*60}")

    # Extract amplitude field at different times
    snapshot_indices = [0, len(states)//2, len(states)-1]  # Initial, middle, final

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    for idx, s_idx in enumerate(snapshot_indices):
        state = states[s_idx]
        X4 = state.positions[:, 3].cpu().numpy()

        # For 3D data, take a slice through the center
        nx, ny, nz = config['grid_shape']
        z_center = nz // 2
        X4_grid = X4.reshape((nx, ny, nz))
        X4_slice = X4_grid[:, :, z_center]

        ax = axes[idx]
        im = ax.imshow(
            X4_slice.T,
            origin='lower',
            cmap='RdBu_r',
            extent=[0, nx*config['h']*1e12, 0, ny*config['h']*1e12],
        )
        ax.set_xlabel('X [pm]')
        ax.set_ylabel('Y [pm]')

        # Compute time correctly
        t = s_idx * config['snapshot_interval'] * config['dt']
        ax.set_title(f"t = {t:.6e} s")

        plt.colorbar(im, ax=ax, label='X⁴ [m]')

    plt.tight_layout()
    plt.savefig('electron_evolution.png', dpi=150)
    print(f"  Saved: electron_evolution.png")

    plt.close()


def main():
    """Main experiment runner."""
    # Parse command-line arguments
    parser = argparse.ArgumentParser(
        description='Electron Stability Test Experiment',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        '--periods',
        type=float,
        default=3.0,
        help='Number of Compton periods to simulate'
    )
    parser.add_argument(
        '--amp',
        type=float,
        default=1e-13,
        help='Amplitude scale in meters (e.g., 1e-13 for 0.1 pm)'
    )
    parser.add_argument(
        '--grid',
        type=int,
        nargs=3,
        default=[20, 20, 20],
        metavar=('NX', 'NY', 'NZ'),
        help='Grid dimensions'
    )
    parser.add_argument(
        '--cfl',
        type=float,
        default=0.2,
        help='CFL parameter (eta)'
    )
    args = parser.parse_args()

    print(f"\n{'='*60}")
    print(f"Electron Stability Test Experiment")
    print(f"{'='*60}")

    # Setup with command-line parameters
    config = setup_experiment(
        grid_shape=tuple(args.grid),
        N_periods=args.periods,
        eta_cfl=args.cfl,
        amplitude_scale=args.amp,
    )

    # Create grid and state
    device = torch.device('cpu')  # Use CPU for now
    grid = BraneGrid(
        grid_shape=config['grid_shape'],
        spacing=config['h'],
        dimension=Dimensionality.THREE_D,
        device=device,
    )

    state = BraneState(
        grid_shape=config['grid_shape'],
        dimension=Dimensionality.THREE_D,
        dtype=torch.float64,  # Use double precision for better accuracy
    )

    # Initialize flat configuration
    state.initialize_flat_configuration(config['h'])

    # Create physics and solver
    m_point = config['m_point']
    rest_length = config['rest_length_frac'] * config['h']

    physics = SpringForceComputer(
        spring_constant=config['k'],
        rest_length=rest_length,
    )

    # Convert point mass to mass density
    h = config['h']
    mass_density = m_point / (h ** 3)  # ρ_m = m_point / h³

    solver = VelocityVerletSolver(
        dt=config['dt'],
        mass_density=mass_density,
        physics=physics,
        grid=grid,
    )

    # Initialize electron
    params = initialize_electron(state, grid, config)

    # Visualize initial state (3 orthogonal slices)
    visualize_initial_state(state, params, config)

    # Run simulation
    states = run_simulation(state, grid, physics, solver, config)

    # Analyze stability
    loss, metrics = analyze_stability(states, params, config)

    # Visualize
    visualize_results(states, params, config)

    # Create animations (6 videos: 3 amplitude + 3 distortion)
    create_all_animations(states, config)

    print(f"\n{'='*60}")
    print(f"Experiment Complete")
    print(f"{'='*60}")
    print(f"\nFinal Stability Loss: {loss:.6e}")
    print(f"\nNext steps:")
    print(f"  1. Refine amplitude scale to match charge = -e")
    print(f"  2. Adjust geometry to match energy = m_e c²")
    print(f"  3. Run optimization to minimize stability loss")
    print(f"  4. Increase grid resolution and simulation time")


if __name__ == '__main__':
    main()