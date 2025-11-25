"""
2D Photon with Fixed Boundary Points

Demonstrates circular wave propagation with fixed boundary conditions.
All edge points are held stationary, and waves reflect from them.
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
from branesim.physics.linear_tension_forces import LinearTensionForceComputer
from branesim.config.simulation_config import SimulationConfig


def initialize_circular_wave(state, grid, wavelength, amplitude, wave_speed, center_x, center_y):
    """
    Initialize a circular wave packet centered at (center_x, center_y).

    Args:
        state: BraneState
        grid: BraneGrid
        wavelength: Wavelength of the carrier wave [m]
        amplitude: Peak amplitude [m]
        wave_speed: Wave speed [m/s]
        center_x: Center x coordinate [m]
        center_y: Center y coordinate [m]
    """
    # Get spatial coordinates
    coords = grid.get_spatial_coordinates()
    x = coords[:, 0]
    y = coords[:, 1]

    # Distance from center
    r = torch.sqrt((x - center_x)**2 + (y - center_y)**2)

    k = 2 * np.pi / wavelength
    omega = wave_speed * k
    sigma = 3 * wavelength / (2 * np.pi)  # Width of wave packet

    # Gaussian envelope
    envelope = amplitude * torch.exp(-(r**2) / (2 * sigma**2))

    # Envelope derivative: dA/dr = -(r/σ²) * A(r)
    envelope_derivative = -(r / (sigma**2)) * envelope

    # Position: circular wave pattern
    state.positions[:, 3] = envelope * torch.cos(k * r)

    # Velocity: expanding circular wave
    # v = ∂/∂t[A(r) cos(kr)] for expanding wave
    # For an outgoing wave, we need: v = -ω A(r) sin(kr) - c A'(r) cos(kr)
    state.velocities[:, 3] = (
        -omega * envelope * torch.sin(k * r) +
        (-wave_speed) * envelope_derivative * torch.cos(k * r)
    )

    print(f"  Wavelength: {wavelength:.4f} m")
    print(f"  Center: ({center_x:.4f}, {center_y:.4f}) m")
    print(f"  Width (σ): {sigma:.4f} m")
    print(f"  Amplitude: {amplitude:.6f} m")


def main():
    """Run 2D photon simulation with fixed boundaries."""
    print("=" * 70)
    print("2D Photon with Fixed Boundary Points")
    print("=" * 70)

    # Configuration
    nx, ny = 60, 60
    wave_speed = 1.0  # m/s
    spacing = 0.02  # m
    cfl_factor = 0.05  # Very conservative for 2D with 8 neighbors

    config = SimulationConfig.for_2d_test(
        nx=nx, ny=ny,
        wave_speed=wave_speed,
        spacing=spacing,
        cfl_factor=cfl_factor
    )

    print(f"\nConfiguration:")
    print(f"  Domain: {nx * spacing:.2f} × {ny * spacing:.2f} m ({nx}×{ny} points)")
    print(f"  Grid spacing h = {spacing:.5f} m")
    print(f"  Wave speed c = {wave_speed:.6f} m/s")
    print(f"  Time step dt = {config.time_step:.6f} s")
    print(f"  CFL number = {cfl_factor:.3f}")

    # Create components
    device = torch.device('cpu')
    dtype = torch.float64

    state = BraneState((nx, ny), Dimensionality.TWO_D, device, dtype)
    state.initialize_flat_configuration(spacing)

    # Set fixed boundaries BEFORE initializing the wave
    state.set_fixed_boundaries()
    print(f"\nBoundary Conditions:")
    print(f"  Fixed points: {state.fixed_mask.sum().item()} / {nx * ny}")
    print(f"  All edge points: FIXED")

    grid = BraneGrid((nx, ny), Dimensionality.TWO_D, spacing, device)

    # Use LinearTensionForceComputer for stability
    physics = LinearTensionForceComputer(config.tension, spacing)
    solver = VelocityVerletSolver(config.time_step, config.mass_density, physics, grid)

    # Initialize circular wave packet in the center
    print(f"\nInitializing circular wave packet...")
    wavelength = 0.15  # m
    amplitude = 0.002  # m
    center_x = (nx - 1) * spacing / 2.0
    center_y = (ny - 1) * spacing / 2.0

    initialize_circular_wave(state, grid, wavelength, amplitude, wave_speed, center_x, center_y)

    solver.initialize_accelerations(state)

    # Apply fixed boundaries to initial state
    state.apply_fixed_boundaries()

    # Initial measurements
    initial_energy = solver.compute_energy(state)

    print(f"\nInitial State:")
    print(f"  Energy: {initial_energy['total']:.6e} J")

    # Run simulation
    print(f"\nRunning simulation...")
    simulation_time = 1.5  # seconds
    num_steps = int(simulation_time / config.time_step)

    # Tracking
    times = []
    energies = []

    snapshot_times = [0.0, 0.3, 0.6, 0.9, 1.2, 1.5]
    snapshots = {}
    snapshot_steps = {int(t / config.time_step): t for t in snapshot_times}

    print_interval = num_steps // 10

    for step in range(num_steps + 1):
        if step in snapshot_steps:
            field = state.get_field_component(3).cpu().numpy()
            snapshots[snapshot_steps[step]] = field.copy()

        if step % 10 == 0:
            energy = solver.compute_energy(state)
            times.append(solver.time)
            energies.append(energy['total'])

        if step % print_interval == 0:
            print(f"  Step {step:6d}/{num_steps}: t={solver.time:.4f}s, "
                  f"E={energy['total']:.6e}J")

        if step < num_steps:
            solver.step(state)

    # Final analysis
    final_energy = solver.compute_energy(state)

    energy_drift = abs(final_energy['total'] - initial_energy['total']) / initial_energy['total']

    print(f"\n{'=' * 70}")
    print("Results:")
    print(f"{'=' * 70}")
    print(f"\nWave Propagation:")
    print(f"  Circular wave expands from center and reflects from edges")

    print(f"\nEnergy Conservation:")
    print(f"  Initial: {initial_energy['total']:.6e} J")
    print(f"  Final:   {final_energy['total']:.6e} J")
    print(f"  Drift:   {energy_drift:.6e} ({energy_drift*100:.6f}%)")

    if energy_drift < 1e-3:
        print(f"  ✓ Good energy conservation")
    else:
        print(f"  ⚠ Energy drift significant")

    # Visualization
    print(f"\nCreating plots...")

    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    fig.suptitle('2D Photon with Fixed Boundaries (c = 1.0 m/s)',
                 fontsize=16, fontweight='bold')

    axes = axes.flatten()

    # Create coordinate arrays
    x_coords = np.arange(nx) * spacing
    y_coords = np.arange(ny) * spacing
    X, Y = np.meshgrid(x_coords, y_coords, indexing='ij')

    for idx, t in enumerate(snapshot_times):
        if t in snapshots:
            field = snapshots[t].reshape(nx, ny)

            # Plot heatmap
            im = axes[idx].imshow(field.T, origin='lower',
                                 extent=[0, x_coords[-1], 0, y_coords[-1]],
                                 cmap='RdBu_r', vmin=-amplitude*1.5, vmax=amplitude*1.5)

            axes[idx].set_xlabel('x [m]', fontsize=10)
            axes[idx].set_ylabel('y [m]', fontsize=10)
            axes[idx].set_title(f't = {t:.2f} s', fontsize=12, fontweight='bold')
            axes[idx].set_aspect('equal')

            # Add colorbar
            plt.colorbar(im, ax=axes[idx], label='ξ [m]', fraction=0.046, pad=0.04)

    plt.tight_layout()
    plt.savefig('photon_2d_fixed_boundaries_propagation.png', dpi=150, bbox_inches='tight')
    print(f"  ✓ Saved: photon_2d_fixed_boundaries_propagation.png")

    # Energy conservation plot
    fig2, ax = plt.subplots(figsize=(10, 6))

    energy_array = np.array(energies)
    ax.plot(times, energy_array / initial_energy['total'], 'g-', linewidth=2)
    ax.axhline(y=1.0, color='r', linestyle='--', linewidth=1, alpha=0.5)
    ax.set_xlabel('Time [s]', fontsize=12)
    ax.set_ylabel('E(t) / E(0)', fontsize=12)
    ax.set_title('Energy Conservation', fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig('photon_2d_fixed_boundaries_energy.png', dpi=150, bbox_inches='tight')
    print(f"  ✓ Saved: photon_2d_fixed_boundaries_energy.png")

    print(f"\n{'=' * 70}")
    print("Simulation complete!")
    print(f"{'=' * 70}")


if __name__ == '__main__':
    main()
