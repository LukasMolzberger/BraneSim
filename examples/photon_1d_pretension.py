"""
1D Photon with Pre-Tensioned Springs

This example demonstrates proper wave propagation using pre-tensioned springs
to achieve linear restoring force for transverse displacements.
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


def initialize_traveling_sine_wave(
    state: BraneState,
    grid: BraneGrid,
    wavelength: float,
    amplitude: float,
    wave_speed: float,
    phase_offset: float = 0.0
):
    """Initialize a pure traveling sine wave (right-moving)."""
    x = (torch.arange(grid.grid_shape[0], device=state.device, dtype=state.dtype) *
         grid.spacing)

    k = 2 * np.pi / wavelength
    omega = wave_speed * k

    # Position: ξ = A sin(kx + φ)
    state.positions[:, 3] = amplitude * torch.sin(k * x + phase_offset)

    # Velocity: ∂ξ/∂t = -Aω cos(kx + φ) [RIGHT-moving]
    state.velocities[:, 3] = -amplitude * omega * torch.cos(k * x + phase_offset)

    print(f"  Wavelength: {wavelength:.4f} m")
    print(f"  Amplitude: {amplitude:.6f} m")
    print(f"  k = {k:.4f} rad/m")
    print(f"  ω = {omega:.4f} rad/s")
    print(f"  Speed: c = ω/k = {omega/k:.6f} m/s")


def track_wave_center(state: BraneState, grid: BraneGrid) -> float:
    """Track center of energy."""
    energy_density = state.velocities[:, 3] ** 2 + state.positions[:, 3] ** 2
    total_energy = energy_density.sum()

    if total_energy > 1e-10:
        x_coords = torch.arange(energy_density.shape[0],
                               device=energy_density.device,
                               dtype=energy_density.dtype)
        center = (x_coords * energy_density).sum() / total_energy
        return center.item() * grid.spacing
    return 0.0


def main():
    """Run 1D photon simulation with pre-tensioned springs."""
    print("=" * 70)
    print("1D Photon with Pre-Tensioned Springs")
    print("=" * 70)

    # Configuration with small CFL for stability with stiff springs
    config = SimulationConfig.for_1d_test(
        nx=400,
        wave_speed=1.0,  # m/s
        spacing=0.01,    # m
        cfl_factor=0.15,  # Reduced for stability with stiff springs!
        device='cpu',
        dtype='float64'
    )

    print(f"\nConfiguration:")
    print(f"  Domain: {config.grid_shape[0] * config.grid_spacing:.2f} m")
    print(f"  Grid spacing h: {config.grid_spacing:.5f} m")
    print(f"  Time step dt: {config.time_step:.6f} s")
    print(f"  CFL number: {config.compute_cfl_number():.3f}")
    print(f"  Wave speed c: {config.compute_wave_speed():.6f} m/s")

    # Create simulation
    device = torch.device(config.device)
    dtype = torch.float64

    state = BraneState(config.grid_shape, config.dimension, device, dtype)
    state.initialize_flat_configuration(config.grid_spacing)

    grid = BraneGrid(config.grid_shape, config.dimension, config.grid_spacing, device)

    physics = SpringForceComputer(config.spring_constant, config.rest_length)
    solver = VelocityVerletSolver(config.time_step, config.mass_density, physics, grid)

    # Stability check
    omega_local = np.sqrt(config.spring_constant / (config.mass_density * config.grid_spacing))
    stability_param = omega_local * config.time_step
    print(f"\nStability Analysis:")
    print(f"  Local oscillation ω = √(k/m) = {omega_local:.2f} rad/s")
    print(f"  Stability parameter ω·dt = {stability_param:.3f} (should be < 2)")
    if stability_param < 2:
        print(f"  ✓ Stable!")
    else:
        print(f"  ✗ WARNING: May be unstable!")

    # Initialize wave
    print(f"\nInitializing wave...")
    wavelength = 0.4  # m
    amplitude = 0.001  # m (smaller amplitude for linear regime)

    initialize_traveling_sine_wave(
        state, grid, wavelength, amplitude,
        config.compute_wave_speed(), 0.0
    )

    solver.initialize_accelerations(state)

    # Initial measurements
    initial_energy = solver.compute_energy(state)
    initial_center = track_wave_center(state, grid)

    print(f"\nInitial State:")
    print(f"  Energy: {initial_energy['total']:.6e} J")
    print(f"  Wave center: {initial_center:.4f} m")

    # Run simulation
    print(f"\nRunning simulation...")
    simulation_time = 2.0  # seconds
    num_steps = int(simulation_time / config.time_step)

    # Tracking
    times = []
    centers = []
    energies = []

    snapshot_times = [0.0, 0.5, 1.0, 1.5, 2.0]
    snapshots = {}
    snapshot_steps = {int(t / config.time_step): t for t in snapshot_times}

    print_interval = num_steps // 20

    for step in range(num_steps + 1):
        if step in snapshot_steps:
            field = state.get_field_component(3).cpu().numpy()
            snapshots[snapshot_steps[step]] = field.copy()

        if step % 20 == 0:
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

    distance_traveled = final_center - initial_center
    measured_speed = distance_traveled / simulation_time
    expected_speed = config.compute_wave_speed()
    speed_error = abs(measured_speed - expected_speed) / expected_speed

    energy_drift = abs(final_energy['total'] - initial_energy['total']) / initial_energy['total']

    print(f"\n{'=' * 70}")
    print("Results:")
    print(f"{'=' * 70}")
    print(f"\nWave Propagation:")
    print(f"  Initial position: {initial_center:.4f} m")
    print(f"  Final position:   {final_center:.4f} m")
    print(f"  Distance traveled: {distance_traveled:.4f} m")
    print(f"  Measured speed:   {measured_speed:.6f} m/s")
    print(f"  Expected speed:   {expected_speed:.6f} m/s")
    print(f"  Speed error:      {speed_error:.6e} ({speed_error*100:.4f}%)")

    if speed_error < 0.05:
        print(f"  ✓ Speed within 5% (PASS)")
    else:
        print(f"  ✗ Speed error > 5% (FAIL)")

    print(f"\nEnergy Conservation:")
    print(f"  Initial: {initial_energy['total']:.6e} J")
    print(f"  Final:   {final_energy['total']:.6e} J")
    print(f"  Drift:   {energy_drift:.6e} ({energy_drift*100:.6f}%)")

    if energy_drift < 1e-4:
        print(f"  ✓ Energy conserved (PASS)")
    else:
        print(f"  ✗ Energy drift > 1e-4 (FAIL)")

    # Visualization
    print(f"\nCreating plots...")

    fig, axes = plt.subplots(len(snapshot_times), 1, figsize=(14, 10))
    fig.suptitle('1D Photon Propagation (Pre-Tensioned Springs)',
                 fontsize=16, fontweight='bold')

    x_coords = np.arange(grid.grid_shape[0]) * grid.spacing

    for idx, t in enumerate(snapshot_times):
        if t in snapshots:
            field = snapshots[t]
            axes[idx].plot(x_coords, field, 'b-', linewidth=2)
            axes[idx].set_ylabel('Amplitude [m]')
            axes[idx].set_xlim(0, x_coords[-1])
            axes[idx].set_ylim(-1.5*amplitude, 1.5*amplitude)
            axes[idx].grid(True, alpha=0.3)
            axes[idx].text(0.02, 0.95, f't = {t:.2f} s',
                          transform=axes[idx].transAxes,
                          fontsize=12, verticalalignment='top',
                          bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

            if idx == len(snapshot_times) - 1:
                axes[idx].set_xlabel('Position [m]')

    plt.tight_layout()
    plt.savefig('photon_1d_pretension_propagation.png', dpi=150, bbox_inches='tight')
    print(f"  ✓ Saved: photon_1d_pretension_propagation.png")

    # Position vs time plot
    fig2, axes2 = plt.subplots(2, 1, figsize=(12, 8))

    # Plot 1: Position
    axes2[0].plot(times, centers, 'b-', linewidth=2, label='Measured')
    axes2[0].plot(times, [initial_center + expected_speed * t for t in times],
                  'r--', linewidth=2, label=f'Expected (c={expected_speed:.3f} m/s)')
    axes2[0].set_xlabel('Time [s]', fontsize=12)
    axes2[0].set_ylabel('Wave Center [m]', fontsize=12)
    axes2[0].set_title('Wave Propagation', fontsize=14, fontweight='bold')
    axes2[0].grid(True, alpha=0.3)
    axes2[0].legend()

    # Plot 2: Energy
    energy_array = np.array(energies)
    axes2[1].plot(times, energy_array / initial_energy['total'], 'g-', linewidth=2)
    axes2[1].axhline(y=1.0, color='r', linestyle='--', linewidth=1, alpha=0.5)
    axes2[1].set_xlabel('Time [s]', fontsize=12)
    axes2[1].set_ylabel('E(t) / E(0)', fontsize=12)
    axes2[1].set_title('Energy Conservation', fontsize=14, fontweight='bold')
    axes2[1].grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig('photon_1d_pretension_analysis.png', dpi=150, bbox_inches='tight')
    print(f"  ✓ Saved: photon_1d_pretension_analysis.png")

    print(f"\n{'=' * 70}")
    print("Simulation complete!")
    print(f"{'=' * 70}")


if __name__ == '__main__':
    main()
