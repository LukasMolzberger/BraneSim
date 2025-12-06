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
import numpy as np
import matplotlib.pyplot as plt

from branesim.core.state import BraneState, Dimensionality
from branesim.core.grid import BraneGrid
from branesim.core.solver import VelocityVerletSolver
from branesim.physics.forces import SpringForceComputer
from branesim.config.simulation_config import PhysicalConstants
from branesim.physics.dimensional_mapping import DimensionalMapper
from branesim.physics.electron_initialization import (
    calibrate_electron_init_params,
    init_electron_state,
)
from branesim.diagnostics.electron_stability import (
    compute_electron_stability_loss,
    build_electron_masks,
    StabilityMetrics,
)


def setup_experiment():
    """
    Set up experiment parameters.

    Returns:
        Dictionary with all configuration parameters
    """
    # Physical constants
    constants = PhysicalConstants()

    # Grid parameters
    # Use a smaller grid for faster testing
    # Target: ~5 grid points per Compton wavelength in each dimension
    lambda_C = constants.lambda_C
    h = lambda_C / 5.0  # Grid spacing

    # Grid should be large enough to contain electron + some buffer
    # Electron diameter ~ 2R ~ λ_C/π ~ 0.32 λ_C
    # Use 3 λ_C in each dimension for buffer
    grid_extent = 3.0 * lambda_C
    n_per_side = int(grid_extent / h) + 1

    # For initial testing, use a smaller grid to make it faster
    # You can increase this for production runs
    nx, ny, nz = 20, 20, 20  # Smaller for testing
    print(f"\nGrid size: {nx} × {ny} × {nz} = {nx*ny*nz} points")
    print(f"Grid spacing h = {h:.6e} m ({h/lambda_C:.4f} λ_C)")

    # Brane parameters
    # These should be calibrated from your Compton-scale requirements
    # For now, use reasonable guesses
    k = 1e3  # Spring constant (will be mapped to sim units)
    rest_length_frac = 0.9  # rest_length = rest_length_frac * h
    m_point_frac = 1.0  # Point mass relative to reference

    # Time stepping
    # CFL condition: dt < η h / c
    # For FCC lattice, η ≈ 0.33
    # Use more conservative value for stability
    eta_cfl = 0.2
    dt = eta_cfl * h / constants.c

    # Simulation duration: Run for N_periods Compton periods
    T_compton = 2 * math.pi / (constants.m_e * constants.c ** 2 / constants.hbar)
    N_periods = 3  # Run for 3 Compton periods
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
        'grid_shape': (nx, ny, nz),
        'h': h,
        'k': k,
        'rest_length_frac': rest_length_frac,
        'm_point_frac': m_point_frac,
        'dt': dt,
        'n_steps': n_steps,
        'snapshot_interval': snapshot_interval,
        'N_periods': N_periods,
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
    # Start with a small amplitude scale - this should be refined
    amplitude_scale = 1e-13  # meters (~ 0.1 pm)

    params = calibrate_electron_init_params(
        constants=constants,
        grid_spacing=h,
        center=center_phys,
        amplitude_scale=amplitude_scale,
    )

    # Initialize electron state
    init_electron_state(state, params)

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
        solver.step(state, grid, physics)

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
    h = config['h']
    m_point = config['m_point_frac'] * 1e-30  # Placeholder mass

    # Build masks
    tube_mask, core_mask = build_electron_masks(
        positions=states[0].positions,
        center=params.center,
        R=params.R,
        tube_max_radius=params.tube_max_radius,
        core_radius_frac=0.5,
    )

    print(f"\n=== Mask Statistics ===")
    print(f"  Tube points: {tube_mask.sum().item()}/{len(tube_mask)} ({100*tube_mask.sum().item()/len(tube_mask):.2f}%)")
    print(f"  Core points: {core_mask.sum().item()}/{len(core_mask)} ({100*core_mask.sum().item()/len(core_mask):.2f}%)")

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
    times_to_plot = [0, len(states)//2, -1]  # Initial, middle, final

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    for idx, t_idx in enumerate(times_to_plot):
        state = states[t_idx]
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
        ax.set_title(f"t = {t_idx * config['snapshot_interval'] * config['dt']:.6e} s")
        plt.colorbar(im, ax=ax, label='X⁴ [m]')

    plt.tight_layout()
    plt.savefig('electron_evolution.png', dpi=150)
    print(f"  Saved: electron_evolution.png")

    plt.close()


def main():
    """Main experiment runner."""
    print(f"\n{'='*60}")
    print(f"Electron Stability Test Experiment")
    print(f"{'='*60}")

    # Setup
    config = setup_experiment()

    # Create grid and state
    grid = BraneGrid(
        grid_shape=config['grid_shape'],
        spacing=config['h'],
        dimension=Dimensionality.THREE_D,
    )

    state = BraneState(
        grid_shape=config['grid_shape'],
        dimension=Dimensionality.THREE_D,
        dtype=torch.float64,  # Use double precision for better accuracy
    )

    # Initialize flat configuration
    state.initialize_flat_configuration(config['h'])

    # Create physics and solver
    m_point = config['m_point_frac'] * 1e-30  # Placeholder
    rest_length = config['rest_length_frac'] * config['h']

    physics = SpringForceComputer(
        spring_constant=config['k'],
        rest_length=rest_length,
        mass_point=m_point,
    )

    solver = VelocityVerletSolver()

    # Initialize electron
    params = initialize_electron(state, grid, config)

    # Run simulation
    states = run_simulation(state, grid, physics, solver, config)

    # Analyze stability
    loss, metrics = analyze_stability(states, params, config)

    # Visualize
    visualize_results(states, params, config)

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