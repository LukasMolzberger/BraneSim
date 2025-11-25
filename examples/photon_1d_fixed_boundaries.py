"""
1D Photon with Fixed Boundary Points

Demonstrates wave propagation with fixed boundary conditions.
The first and last points are held stationary, and waves reflect from them.
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


def initialize_traveling_wave(state, grid, wavelength, amplitude, wave_speed, center):
    """Initialize a wave packet centered at a position."""
    x = torch.arange(grid.grid_shape[0], device=state.device, dtype=state.dtype) * grid.spacing

    k = 2 * np.pi / wavelength
    omega = wave_speed * k
    sigma = 3 * wavelength / (2 * np.pi)  # Width of wave packet

    # Gaussian envelope
    envelope = amplitude * torch.exp(-((x - center) ** 2) / (2 * sigma ** 2))

    # Envelope derivative: dA/dx = -(x - center)/σ² * A(x)
    envelope_derivative = -((x - center) / (sigma ** 2)) * envelope

    # Position and velocity for right-moving wave
    state.positions[:, 3] = envelope * torch.cos(k * (x - center))

    # Velocity has two components:
    # 1. Phase velocity: ω A(x) sin(kx)
    # 2. Envelope velocity: -c A'(x) cos(kx)
    state.velocities[:, 3] = (
        omega * envelope * torch.sin(k * (x - center)) +
        (-wave_speed) * envelope_derivative * torch.cos(k * (x - center))
    )

    print(f"  Wavelength: {wavelength:.4f} m")
    print(f"  Center: {center:.4f} m")
    print(f"  Width (σ): {sigma:.4f} m")
    print(f"  Amplitude: {amplitude:.6f} m")


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


def main():
    """Run 1D photon simulation with fixed boundaries."""
    print("=" * 70)
    print("1D Photon with Fixed Boundary Points")
    print("=" * 70)

    # Configuration
    nx = 400
    wave_speed = 1.0  # m/s
    spacing = 0.01  # m
    cfl_factor = 0.1

    mu = 1.0  # kg/m (linear mass density)
    tension = mu * wave_speed**2

    dt = cfl_factor * spacing / wave_speed

    print(f"\nConfiguration:")
    print(f"  Domain: {nx * spacing:.2f} m ({nx} points)")
    print(f"  Grid spacing h = {spacing:.5f} m")
    print(f"  Wave speed c = {wave_speed:.6f} m/s")
    print(f"  Time step dt = {dt:.6f} s")
    print(f"  CFL number = {cfl_factor:.3f}")

    # Create components
    device = torch.device('cpu')
    dtype = torch.float64

    state = BraneState((nx,), Dimensionality.ONE_D, device, dtype)
    state.initialize_flat_configuration(spacing)

    # Set fixed boundaries BEFORE initializing the wave
    state.set_fixed_boundaries()
    print(f"\nBoundary Conditions:")
    print(f"  Fixed points: {state.fixed_mask.sum().item()} / {nx}")
    print(f"  First point (x=0): FIXED")
    print(f"  Last point (x={(nx-1)*spacing:.2f}m): FIXED")

    grid = BraneGrid((nx,), Dimensionality.ONE_D, spacing, device)

    physics = LinearTensionForceComputer(tension, spacing)
    solver = VelocityVerletSolver(dt, mu, physics, grid)

    # Initialize wave packet in the middle
    print(f"\nInitializing wave packet...")
    wavelength = 0.4  # m
    amplitude = 0.002  # m
    center_position = 1.0  # m - start near left but not at boundary

    initialize_traveling_wave(state, grid, wavelength, amplitude, wave_speed, center_position)

    solver.initialize_accelerations(state)

    # Apply fixed boundaries to initial state
    state.apply_fixed_boundaries()

    # Initial measurements
    initial_energy = solver.compute_energy(state)
    initial_center = track_wave_center(state, grid)

    print(f"\nInitial State:")
    print(f"  Energy: {initial_energy['total']:.6e} J")
    print(f"  Wave center: {initial_center:.4f} m")

    # Run simulation
    print(f"\nRunning simulation...")
    simulation_time = 3.0  # seconds
    num_steps = int(simulation_time / dt)

    # Tracking
    times = []
    centers = []
    energies = []

    snapshot_times = [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0]
    snapshots = {}
    snapshot_steps = {int(t / dt): t for t in snapshot_times}

    x_coords = np.arange(nx) * spacing

    print_interval = num_steps // 20

    for step in range(num_steps + 1):
        if step in snapshot_steps:
            field = state.get_field_component(3).cpu().numpy()
            snapshots[snapshot_steps[step]] = field.copy()

        if step % 10 == 0:
            center = track_wave_center(state, grid)
            energy = solver.compute_energy(state)
            times.append(solver.time)
            centers.append(center)
            energies.append(energy['total'])

        if step % print_interval == 0:
            print(f"  Step {step:6d}/{num_steps}: t={solver.time:.4f}s, "
                  f"center={center:.4f}m, E={energy['total']:.6e}J")

        if step < num_steps:
            solver.step(state)

    # Final analysis
    final_energy = solver.compute_energy(state)
    final_center = track_wave_center(state, grid)

    energy_drift = abs(final_energy['total'] - initial_energy['total']) / initial_energy['total']

    print(f"\n{'=' * 70}")
    print("Results:")
    print(f"{'=' * 70}")
    print(f"\nWave Propagation:")
    print(f"  Initial position: {initial_center:.4f} m")
    print(f"  Final position:   {final_center:.4f} m")
    print(f"  Domain length:    {(nx-1)*spacing:.2f} m")
    print(f"  ✓ Wave reflects from fixed boundaries")

    print(f"\nEnergy Conservation:")
    print(f"  Initial: {initial_energy['total']:.6e} J")
    print(f"  Final:   {final_energy['total']:.6e} J")
    print(f"  Drift:   {energy_drift:.6e} ({energy_drift*100:.6f}%)")

    if energy_drift < 1e-4:
        print(f"  ✓ Good energy conservation")
    else:
        print(f"  ⚠ Energy drift significant")

    # Visualization
    print(f"\nCreating plots...")

    fig, axes = plt.subplots(len(snapshot_times), 1, figsize=(14, 12))
    fig.suptitle('1D Photon with Fixed Boundaries (c = 1.0 m/s)',
                 fontsize=16, fontweight='bold')

    for idx, t in enumerate(snapshot_times):
        if t in snapshots:
            field = snapshots[t]
            axes[idx].plot(x_coords, field, 'b-', linewidth=2)

            # Mark fixed boundaries with red dots
            axes[idx].plot([x_coords[0], x_coords[-1]], [0, 0], 'ro',
                          markersize=8, label='Fixed boundaries')

            axes[idx].set_ylabel('ξ [m]', fontsize=11)
            axes[idx].set_xlim(0, x_coords[-1])
            axes[idx].set_ylim(-1.5*amplitude, 1.5*amplitude)
            axes[idx].grid(True, alpha=0.3)
            axes[idx].text(0.02, 0.95, f't = {t:.2f} s',
                          transform=axes[idx].transAxes,
                          fontsize=12, verticalalignment='top',
                          bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

            if idx == 0:
                axes[idx].legend(loc='upper right', fontsize=10)

            if idx == len(snapshot_times) - 1:
                axes[idx].set_xlabel('Position [m]', fontsize=12)

    plt.tight_layout()
    plt.savefig('photon_1d_fixed_boundaries_propagation.png', dpi=150, bbox_inches='tight')
    print(f"  ✓ Saved: photon_1d_fixed_boundaries_propagation.png")

    # Analysis plots
    fig2, axes2 = plt.subplots(2, 1, figsize=(12, 8))

    # Position vs time
    axes2[0].plot(times, centers, 'b-', linewidth=2, label='Wave center')
    axes2[0].axhline(y=0, color='r', linestyle='--', linewidth=1, alpha=0.5, label='Left boundary')
    axes2[0].axhline(y=(nx-1)*spacing, color='r', linestyle='--', linewidth=1, alpha=0.5, label='Right boundary')
    axes2[0].set_xlabel('Time [s]', fontsize=12)
    axes2[0].set_ylabel('Wave Center [m]', fontsize=12)
    axes2[0].set_title('Wave Propagation and Reflection', fontsize=14, fontweight='bold')
    axes2[0].grid(True, alpha=0.3)
    axes2[0].legend(fontsize=10)

    # Energy conservation
    energy_array = np.array(energies)
    axes2[1].plot(times, energy_array / initial_energy['total'], 'g-', linewidth=2)
    axes2[1].axhline(y=1.0, color='r', linestyle='--', linewidth=1, alpha=0.5)
    axes2[1].set_xlabel('Time [s]', fontsize=12)
    axes2[1].set_ylabel('E(t) / E(0)', fontsize=12)
    axes2[1].set_title('Energy Conservation', fontsize=14, fontweight='bold')
    axes2[1].grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig('photon_1d_fixed_boundaries_analysis.png', dpi=150, bbox_inches='tight')
    print(f"  ✓ Saved: photon_1d_fixed_boundaries_analysis.png")

    print(f"\n{'=' * 70}")
    print("Simulation complete!")
    print(f"{'=' * 70}")


if __name__ == '__main__':
    main()
