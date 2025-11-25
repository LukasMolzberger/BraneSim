"""
1D Photon Propagation - Clean Right-Moving Wave

Creates a photon wave packet that propagates cleanly to the right at the speed of light.
Uses proper traveling wave initialization to minimize dispersion.
"""

import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import torch
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, FFMpegWriter

from branesim.core.state import BraneState, Dimensionality
from branesim.core.grid import BraneGrid
from branesim.core.solver import VelocityVerletSolver
from branesim.physics.forces import SpringForceComputer
from branesim.config.simulation_config import SimulationConfig


def initialize_traveling_wave_packet(
    state: BraneState,
    grid: BraneGrid,
    center: float,
    wavelength: float,
    amplitude: float,
    num_cycles: float,
    wave_speed: float
):
    """
    Initialize a clean traveling wave packet (right-moving photon).

    Uses a sinusoidal wave with Gaussian envelope, initialized with both
    position and velocity to create a pure right-moving wave.

    The key is to set both ξ(x,0) and ∂ξ/∂t(x,0) consistently so that
    the wave travels to the right without reflection or dispersion artifacts.

    For a right-moving wave: ξ(x,t) = A(x-ct) exp(i(kx - ωt))
    At t=0: ξ(x,0) = A(x) cos(kx)
            ∂ξ/∂t(x,0) = A'(x)c cos(kx) - A(x)ω sin(kx)

    For slowly varying envelope (A' << kA), we can approximate:
            ∂ξ/∂t(x,0) ≈ -A(x)ω sin(kx) where ω = ck

    Args:
        state: BraneState to initialize
        grid: BraneGrid with spatial coordinates
        center: Center position of wave packet [m]
        wavelength: Wavelength λ [m]
        amplitude: Peak amplitude [m]
        num_cycles: Number of oscillation cycles in the packet
        wave_speed: Wave speed c [m/s]
    """
    # Spatial coordinates
    x = (torch.arange(grid.grid_shape[0], device=state.device, dtype=state.dtype) *
         grid.spacing)

    # Wave parameters
    k = 2 * np.pi / wavelength  # Wave number
    omega = wave_speed * k       # Angular frequency
    sigma = num_cycles * wavelength / (2 * np.pi)  # Gaussian width

    # Gaussian envelope (spatial modulation)
    envelope = amplitude * torch.exp(-((x - center) ** 2) / (2 * sigma ** 2))

    # For a right-moving wave traveling as ξ(x,t) = A(x-ct) cos(k(x-ct))
    # At t=0: ξ(x,0) = A(x) cos(kx)
    #         ∂ξ/∂t(x,0) = +ω A(x) sin(kx)  [RIGHT-moving needs positive!]

    # Position: ξ = A(x) * cos(k(x - x0))
    state.positions[:, 3] = envelope * torch.cos(k * (x - center))

    # Velocity: ∂ξ/∂t = +ω * A(x) * sin(k(x - x0))
    # Positive sign creates wave moving to the RIGHT
    state.velocities[:, 3] = omega * envelope * torch.sin(k * (x - center))

    print(f"  Wave packet initialized:")
    print(f"    Center: {center:.2f} m")
    print(f"    Wavelength: {wavelength:.4f} m")
    print(f"    Width (σ): {sigma:.4f} m")
    print(f"    Cycles: {num_cycles:.1f}")
    print(f"    k = {k:.2f} rad/m")
    print(f"    ω = {omega:.2f} rad/s")


def track_energy_center(state: BraneState, grid: BraneGrid) -> float:
    """
    Track center of energy (more robust than amplitude center).

    Args:
        state: BraneState
        grid: BraneGrid

    Returns:
        Center of energy position [m]
    """
    # Energy density: E = ½ρ(∂ξ/∂t)² + ½k(∂ξ/∂x)²
    # We'll use velocity squared as proxy for energy
    energy_density = state.velocities[:, 3] ** 2 + state.positions[:, 3] ** 2

    total_energy = energy_density.sum()

    if total_energy > 1e-10:
        x_coords = torch.arange(energy_density.shape[0], device=energy_density.device, dtype=energy_density.dtype)
        center = (x_coords * energy_density).sum() / total_energy
        return center.item() * grid.spacing
    else:
        return 0.0


def main():
    """Run clean 1D photon simulation."""
    print("=" * 70)
    print("1D Photon Propagation - Clean Right-Moving Wave")
    print("=" * 70)

    # Configuration - optimized for clean propagation
    config = SimulationConfig.for_1d_test(
        nx=600,           # Longer domain
        wave_speed=1.0,   # m/s (normalized to c=1)
        spacing=0.005,    # Finer grid (better resolution)
        cfl_factor=0.3,   # More conservative CFL for stability
        device='cpu',
        dtype='float64'   # Higher precision
    )

    print(f"\nSimulation Configuration:")
    print(f"  Domain: {config.grid_shape[0]} points × {config.grid_spacing:.5f} m")
    print(f"  Total length: {config.grid_shape[0] * config.grid_spacing:.2f} m")
    print(f"  Wave speed (c): {config.compute_wave_speed():.6f} m/s")
    print(f"  Time step (dt): {config.time_step:.6f} s")
    print(f"  CFL number: {config.compute_cfl_number():.3f}")

    # Create simulation components
    device = torch.device(config.device)
    dtype = torch.float64 if config.dtype == 'float64' else torch.float32

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

    # Initialize photon wave packet
    print("\nInitializing photon wave packet...")
    center_position = 0.5  # m - start position
    wavelength = 0.2       # m - longer wavelength (40x grid spacing for low dispersion)
    amplitude = 0.001      # m - small amplitude
    num_cycles = 3         # Fewer cycles for narrower spectrum

    initialize_traveling_wave_packet(
        state, grid, center_position, wavelength, amplitude, num_cycles, config.compute_wave_speed()
    )

    # Initialize accelerations
    solver.initialize_accelerations(state)

    # Measure initial state
    initial_energy = solver.compute_energy(state)
    initial_center = track_energy_center(state, grid)

    print(f"\nInitial State:")
    print(f"  Total energy: {initial_energy['total']:.6e} J")
    print(f"  Energy center: {initial_center:.4f} m")

    # Run simulation
    print(f"\nRunning simulation...")

    simulation_time = 2.0  # seconds
    num_steps = int(simulation_time / config.time_step)

    # Storage for tracking
    times = []
    centers = []
    energies = []

    print_interval = num_steps // 20

    for step in range(num_steps):
        solver.step(state)

        if step % 10 == 0:
            center = track_energy_center(state, grid)
            energy = solver.compute_energy(state)
            times.append(solver.time)
            centers.append(center)
            energies.append(energy['total'])

        if step % print_interval == 0 and step > 0:
            print(f"  Step {step:5d}/{num_steps}: t={solver.time:.4f}s, "
                  f"center={center:.4f}m, E={energy['total']:.6e}J")

    # Final measurements
    final_energy = solver.compute_energy(state)
    final_center = track_energy_center(state, grid)

    # Analysis
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
    print(f"  Time elapsed:     {simulation_time:.4f} s")
    print(f"  Measured speed:   {measured_speed:.6f} m/s")
    print(f"  Expected speed:   {expected_speed:.6f} m/s")
    print(f"  Speed error:      {speed_error:.6e} ({speed_error*100:.4f}%)")

    if speed_error < 0.01:
        print(f"  ✓ Speed matches c within 1% (PASS)")
    else:
        print(f"  ✗ Speed error > 1% (FAIL)")

    print(f"\nEnergy Conservation:")
    print(f"  Initial energy: {initial_energy['total']:.6e} J")
    print(f"  Final energy:   {final_energy['total']:.6e} J")
    print(f"  Relative drift: {energy_drift:.6e} ({energy_drift*100:.6f}%)")

    if energy_drift < 1e-6:
        print(f"  ✓ Energy conserved to < 1e-6 (PASS)")
    else:
        print(f"  ✗ Energy drift > 1e-6 (FAIL)")

    # Create plots
    print(f"\nGenerating plots...")

    fig, axes = plt.subplots(3, 1, figsize=(14, 10))

    # Plot 1: Position vs time
    axes[0].plot(times, centers, 'b-', linewidth=2)
    axes[0].plot(times, [initial_center + expected_speed * t for t in times], 'r--',
                 linewidth=2, label=f'Expected (c = {expected_speed:.6f} m/s)')
    axes[0].set_xlabel('Time [s]', fontsize=12)
    axes[0].set_ylabel('Wave Packet Center [m]', fontsize=12)
    axes[0].set_title('Photon Propagation: Position vs Time', fontsize=14, fontweight='bold')
    axes[0].grid(True, alpha=0.3)
    axes[0].legend(fontsize=10)

    # Plot 2: Energy vs time
    energy_array = np.array(energies)
    axes[1].plot(times, energy_array / initial_energy['total'], 'g-', linewidth=2)
    axes[1].axhline(y=1.0, color='r', linestyle='--', linewidth=1, alpha=0.5)
    axes[1].set_xlabel('Time [s]', fontsize=12)
    axes[1].set_ylabel('E(t) / E(0)', fontsize=12)
    axes[1].set_title('Energy Conservation', fontsize=14, fontweight='bold')
    axes[1].grid(True, alpha=0.3)

    # Plot 3: Current wave profile
    x_coords = np.arange(grid.grid_shape[0]) * grid.spacing
    field = state.get_field_component(3).cpu().numpy()
    axes[2].plot(x_coords, field, 'b-', linewidth=2)
    axes[2].axvline(x=final_center, color='r', linestyle='--', linewidth=1,
                    label=f'Center: {final_center:.4f} m')
    axes[2].set_xlabel('Position [m]', fontsize=12)
    axes[2].set_ylabel('Amplitude $\\xi^3$ [m]', fontsize=12)
    axes[2].set_title(f'Wave Profile at t = {simulation_time:.2f} s', fontsize=14, fontweight='bold')
    axes[2].grid(True, alpha=0.3)
    axes[2].legend(fontsize=10)

    plt.tight_layout()
    plt.savefig('photon_1d_clean_results.png', dpi=150, bbox_inches='tight')
    print(f"  ✓ Saved: photon_1d_clean_results.png")

    plt.show()

    print(f"\n{'=' * 70}")
    print("Simulation complete!")
    print(f"{'=' * 70}")


if __name__ == '__main__':
    main()
