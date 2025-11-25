"""
Example 1D Wave Simulation

Demonstrates 1D Gaussian pulse propagation on a brane string.
Validates wave speed and energy conservation.
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
from branesim.config.simulation_config import SimulationConfig
from branesim.visualization.renderer import MatplotlibRenderer


def initialize_gaussian_pulse(
    state: BraneState,
    grid: BraneGrid,
    center: float,
    width: float,
    amplitude: float,
    wave_speed: float
):
    """
    Initialize Gaussian pulse in the 4th dimension (amplitude).

    The pulse is initialized with both position and velocity to create
    a right-moving wave.

    Args:
        state: BraneState to initialize
        grid: BraneGrid with spatial coordinates
        center: Center position [m]
        width: Gaussian width σ [m]
        amplitude: Peak amplitude [m]
        wave_speed: Wave speed for velocity initialization [m/s]
    """
    # Get spatial coordinates
    x = (torch.arange(grid.grid_shape[0], device=state.device, dtype=state.dtype) *
         grid.spacing)

    # Gaussian envelope
    envelope = amplitude * torch.exp(-((x - center) ** 2) / (2 * width ** 2))

    # Set 4th component (amplitude)
    state.positions[:, 3] = envelope

    # Set velocity for right-moving pulse
    # For traveling Gaussian: ∂ξ/∂t = c * ∂ξ/∂x = c * (x-x0)/σ² * ξ
    state.velocities[:, 3] = wave_speed * (x - center) / width ** 2 * envelope


def track_pulse_center(state: BraneState) -> float:
    """
    Track center of mass of the pulse.

    Args:
        state: BraneState

    Returns:
        Center position index (float)
    """
    field = state.positions[:, 3].abs()
    total = field.sum()

    if total > 1e-10:
        indices = torch.arange(field.shape[0], device=field.device, dtype=field.dtype)
        center = (indices * field).sum() / total
        return center.item()
    else:
        return 0.0


def main():
    """Run 1D wave simulation."""
    print("=" * 60)
    print("1D Wave Simulation: Gaussian Pulse")
    print("=" * 60)

    # Configuration
    config = SimulationConfig.for_1d_test(
        nx=200,
        wave_speed=1.0,  # m/s
        spacing=0.01,    # m
        cfl_factor=0.4,
        device='cpu',
        dtype='float32'
    )

    print("\nConfiguration:")
    print(config)

    # Create state, grid, physics, and solver
    device = torch.device(config.device)
    dtype = torch.float32 if config.dtype == 'float32' else torch.float64

    state = BraneState(config.grid_shape, config.dimension, device, dtype)
    state.initialize_flat_configuration(config.grid_spacing)

    grid = BraneGrid(config.grid_shape, config.dimension, config.grid_spacing, device)

    physics = SpringForceComputer(
        config.spring_constant,
        config.rest_length,
        config.critical_strain
    )

    solver = VelocityVerletSolver(
        config.time_step,
        config.mass_density,
        physics,
        grid
    )

    # Initialize Gaussian pulse
    center = 0.3  # m (off-center to allow propagation)
    width = 0.05   # m (narrower pulse)
    amplitude = 0.0005  # m (small displacement)

    initialize_gaussian_pulse(state, grid, center, width, amplitude, config.compute_wave_speed())

    # Initialize accelerations
    solver.initialize_accelerations(state)

    # Track initial energy
    initial_energy = solver.compute_energy(state)
    print(f"\nInitial Energy:")
    print(f"  Kinetic:   {initial_energy['kinetic']:.6e} J")
    print(f"  Potential: {initial_energy['potential']:.6e} J")
    print(f"  Total:     {initial_energy['total']:.6e} J")

    # Track pulse center for wave speed measurement
    initial_center = track_pulse_center(state)
    print(f"\nInitial pulse center: {initial_center * config.grid_spacing:.4f} m")

    # Create renderer
    renderer = MatplotlibRenderer(config.dimension, component_idx=3)

    # Simulation loop
    num_steps = 1000
    render_interval = 50
    energy_history = []
    center_history = []

    print(f"\nRunning {num_steps} steps...")

    plt.ion()  # Interactive mode

    for step in range(num_steps):
        # Time step
        solver.step(state)

        # Track energy
        energy = solver.compute_energy(state)
        energy_history.append(energy['total'])

        # Track pulse center
        center = track_pulse_center(state)
        center_history.append(center)

        # Render
        if step % render_interval == 0:
            renderer.render_field(
                state,
                grid,
                title=f"1D Wave at t = {solver.time:.4f} s (step {step})",
                show=True
            )

            print(f"Step {step:4d}: t={solver.time:.4f}s, "
                  f"E={energy['total']:.6e}J, center={center*config.grid_spacing:.4f}m")

    plt.ioff()

    # Final statistics
    final_energy = solver.compute_energy(state)
    final_center = track_pulse_center(state)

    print("\n" + "=" * 60)
    print("Final Results:")
    print("=" * 60)

    # Energy conservation
    energy_drift = abs(final_energy['total'] - initial_energy['total']) / initial_energy['total']
    print(f"\nEnergy Conservation:")
    print(f"  Initial: {initial_energy['total']:.6e} J")
    print(f"  Final:   {final_energy['total']:.6e} J")
    print(f"  Drift:   {energy_drift:.6e} ({energy_drift*100:.4f}%)")

    if energy_drift < 1e-6:
        print("  ✓ Energy conserved to < 1e-6 (PASS)")
    else:
        print("  ✗ Energy drift > 1e-6 (FAIL)")

    # Wave speed measurement
    distance_traveled = (final_center - initial_center) * config.grid_spacing
    time_elapsed = solver.time
    measured_speed = distance_traveled / time_elapsed
    expected_speed = config.compute_wave_speed()
    speed_error = abs(measured_speed - expected_speed) / expected_speed

    print(f"\nWave Speed:")
    print(f"  Expected:  {expected_speed:.6f} m/s")
    print(f"  Measured:  {measured_speed:.6f} m/s")
    print(f"  Error:     {speed_error:.6e} ({speed_error*100:.4f}%)")

    if speed_error < 0.01:
        print("  ✓ Wave speed within 1% (PASS)")
    else:
        print("  ✗ Wave speed error > 1% (FAIL)")

    # Plot energy history
    plt.figure(figsize=(12, 4))

    plt.subplot(1, 2, 1)
    plt.plot(np.array(energy_history) / initial_energy['total'])
    plt.xlabel('Time Step')
    plt.ylabel('E(t) / E(0)')
    plt.title('Energy Conservation')
    plt.grid(True, alpha=0.3)

    plt.subplot(1, 2, 2)
    plt.plot(np.array(center_history) * config.grid_spacing)
    plt.xlabel('Time Step')
    plt.ylabel('Pulse Center [m]')
    plt.title('Pulse Propagation')
    plt.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig('1d_wave_results.png', dpi=150)
    print("\nResults saved to '1d_wave_results.png'")

    plt.show()


if __name__ == '__main__':
    main()
