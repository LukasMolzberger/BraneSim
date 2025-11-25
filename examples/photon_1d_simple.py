"""
1D Photon as Pure Traveling Wave

Uses a simple traveling sine wave to demonstrate clean propagation at speed of light.
No Gaussian envelope - just a pure sinusoidal wave moving to the right.
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


def initialize_traveling_sine_wave(
    state: BraneState,
    grid: BraneGrid,
    wavelength: float,
    amplitude: float,
    wave_speed: float,
    phase_offset: float = 0.0
):
    """
    Initialize a pure traveling sine wave (right-moving).

    For a right-moving wave: ξ(x,t) = A sin(kx - ωt + φ)
    At t=0:
        ξ(x,0) = A sin(kx + φ)
        ∂ξ/∂t(x,0) = -Aω cos(kx + φ) = Aω sin(kx + φ + π/2)

    Args:
        state: BraneState to initialize
        grid: BraneGrid
        wavelength: Wave wavelength λ [m]
        amplitude: Wave amplitude A [m]
        wave_speed: Wave speed c [m/s]
        phase_offset: Initial phase φ [rad]
    """
    x = (torch.arange(grid.grid_shape[0], device=state.device, dtype=state.dtype) *
         grid.spacing)

    k = 2 * np.pi / wavelength
    omega = wave_speed * k

    # Position: ξ = A sin(kx + φ)
    state.positions[:, 3] = amplitude * torch.sin(k * x + phase_offset)

    # Velocity: ∂ξ/∂t = -Aω cos(kx + φ)
    # This creates a RIGHT-moving wave
    state.velocities[:, 3] = -amplitude * omega * torch.cos(k * x + phase_offset)

    print(f"  Wave initialized:")
    print(f"    Wavelength: {wavelength:.4f} m")
    print(f"    Amplitude: {amplitude:.6f} m")
    print(f"    k = {k:.4f} rad/m")
    print(f"    ω = {omega:.4f} rad/s")
    print(f"    Expected speed: c = ω/k = {omega/k:.6f} m/s")


def measure_wave_speed_fft(state: BraneState, grid: BraneGrid, old_state: BraneState, dt: float) -> float:
    """
    Measure wave speed using phase shift between timesteps.

    Args:
        state: Current state
        grid: Grid
        old_state: Previous state
        dt: Time step

    Returns:
        Measured wave speed [m/s]
    """
    field_now = state.positions[:, 3].cpu().numpy()
    field_old = old_state.positions[:, 3].cpu().numpy()

    # Cross-correlation to find phase shift
    correlation = np.correlate(field_now, field_old, mode='same')
    shift_idx = np.argmax(correlation) - len(field_now) // 2
    shift_distance = shift_idx * grid.spacing

    if abs(dt) > 1e-10:
        speed = shift_distance / dt
        return speed
    else:
        return 0.0


def main():
    """Run simple traveling wave simulation."""
    print("=" * 70)
    print("1D Photon as Pure Traveling Sine Wave")
    print("=" * 70)

    # Configuration
    config = SimulationConfig.for_1d_test(
        nx=400,
        wave_speed=1.0,
        spacing=0.01,      # 100 points per meter
        cfl_factor=0.4,
        device='cpu',
        dtype='float64'
    )

    print(f"\nConfiguration:")
    print(f"  Domain: {config.grid_shape[0] * config.grid_spacing:.2f} m")
    print(f"  Grid spacing: {config.grid_spacing:.4f} m")
    print(f"  Wave speed (c): {config.compute_wave_speed():.6f} m/s")
    print(f"  Time step: {config.time_step:.6f} s")
    print(f"  CFL: {config.compute_cfl_number():.3f}")

    # Create simulation
    device = torch.device(config.device)
    dtype = torch.float64

    state = BraneState(config.grid_shape, config.dimension, device, dtype)
    state.initialize_flat_configuration(config.grid_spacing)

    grid = BraneGrid(config.grid_shape, config.dimension, config.grid_spacing, device)

    physics = SpringForceComputer(config.spring_constant, config.rest_length)
    solver = VelocityVerletSolver(config.time_step, config.mass_density, physics, grid)

    # Debug: Print actual parameters
    print(f"\n  Spring constant k = {config.spring_constant:.2f} N/m")
    print(f"  Rest length L_0 = {config.rest_length:.6f} m")
    print(f"  Mass density ρ_m = {config.mass_density:.6f} kg/m³")
    print(f"  Mass per point m = ρ_m * h = {config.mass_density * config.grid_spacing:.6f} kg")
    print(f"  Theoretical wave speed = √(k/m) = {np.sqrt(config.spring_constant / (config.mass_density * config.grid_spacing)):.6f} m/s")

    # Initialize wave
    print("\nInitializing traveling wave...")
    wavelength = 0.4  # m - well resolved (40 points per wavelength)
    amplitude = 0.002  # m
    phase_offset = 0.0

    initialize_traveling_sine_wave(state, grid, wavelength, amplitude, config.compute_wave_speed(), phase_offset)

    solver.initialize_accelerations(state)

    # Store initial state
    initial_state = state.clone()
    initial_energy = solver.compute_energy(state)

    print(f"\nInitial energy: {initial_energy['total']:.6e} J")

    # Run simulation with snapshots
    print(f"\nRunning simulation...")
    simulation_time = 1.0  # 1 second
    num_steps = int(simulation_time / config.time_step)

    snapshot_times = [0.0, 0.25, 0.5, 0.75, 1.0]
    snapshots = {}
    snapshot_step = {}

    for t in snapshot_times:
        step = int(t / config.time_step)
        snapshot_step[step] = t

    x_coords = np.arange(grid.grid_shape[0]) * grid.spacing

    # Store old state for velocity measurement
    old_state = state.clone()
    measurement_interval = 50

    measured_speeds = []
    times_for_speed = []

    for step in range(num_steps + 1):
        if step in snapshot_step:
            field = state.get_field_component(3).cpu().numpy()
            snapshots[snapshot_step[step]] = field.copy()

        if step % measurement_interval == 0 and step > 0:
            speed = measure_wave_speed_fft(state, grid, old_state, measurement_interval * config.time_step)
            measured_speeds.append(speed)
            times_for_speed.append(solver.time)
            old_state = state.clone()

        if step % 100 == 0:
            energy = solver.compute_energy(state)
            print(f"  Step {step:5d}/{num_steps}: t={solver.time:.4f}s, E={energy['total']:.6e}J")

        if step < num_steps:
            solver.step(state)

    # Final analysis
    final_energy = solver.compute_energy(state)
    energy_drift = abs(final_energy['total'] - initial_energy['total']) / initial_energy['total']

    avg_speed = np.mean(measured_speeds) if measured_speeds else 0.0
    speed_error = abs(avg_speed - config.compute_wave_speed()) / config.compute_wave_speed()

    print(f"\n{'=' * 70}")
    print("Results:")
    print(f"{'=' * 70}")
    print(f"\nWave Propagation:")
    print(f"  Expected speed: {config.compute_wave_speed():.6f} m/s")
    print(f"  Measured speed: {avg_speed:.6f} m/s")
    print(f"  Speed error: {speed_error:.6e} ({speed_error*100:.4f}%)")

    if speed_error < 0.01:
        print(f"  ✓ Speed within 1% (PASS)")
    else:
        print(f"  ✗ Speed error > 1% (FAIL)")

    print(f"\nEnergy Conservation:")
    print(f"  Initial: {initial_energy['total']:.6e} J")
    print(f"  Final:   {final_energy['total']:.6e} J")
    print(f"  Drift:   {energy_drift:.6e} ({energy_drift*100:.6f}%)")

    if energy_drift < 1e-6:
        print(f"  ✓ Energy conserved (PASS)")
    else:
        print(f"  ✗ Energy drift > 1e-6 (FAIL)")

    # Create visualization
    print(f"\nCreating plots...")

    fig, axes = plt.subplots(len(snapshot_times), 1, figsize=(14, 10))
    fig.suptitle('1D Traveling Wave Propagation (Photon)', fontsize=16, fontweight='bold')

    for idx, t in enumerate(snapshot_times):
        if t in snapshots:
            field = snapshots[t]
            axes[idx].plot(x_coords, field, 'b-', linewidth=2)
            axes[idx].set_ylabel('Amplitude [m]')
            axes[idx].set_xlim(0, x_coords[-1])
            axes[idx].set_ylim(-1.2*amplitude, 1.2*amplitude)
            axes[idx].grid(True, alpha=0.3)
            axes[idx].text(0.02, 0.95, f't = {t:.2f} s', transform=axes[idx].transAxes,
                          fontsize=12, verticalalignment='top',
                          bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

            if idx == len(snapshot_times) - 1:
                axes[idx].set_xlabel('Position [m]')

    plt.tight_layout()
    plt.savefig('photon_1d_simple_propagation.png', dpi=150, bbox_inches='tight')
    print(f"  ✓ Saved: photon_1d_simple_propagation.png")

    # Speed vs time plot
    if measured_speeds:
        fig2, ax2 = plt.subplots(figsize=(10, 6))
        ax2.plot(times_for_speed, measured_speeds, 'bo-', linewidth=2, markersize=6, label='Measured')
        ax2.axhline(y=config.compute_wave_speed(), color='r', linestyle='--', linewidth=2, label='Expected')
        ax2.set_xlabel('Time [s]', fontsize=12)
        ax2.set_ylabel('Wave Speed [m/s]', fontsize=12)
        ax2.set_title('Wave Speed Over Time', fontsize=14, fontweight='bold')
        ax2.grid(True, alpha=0.3)
        ax2.legend(fontsize=11)
        plt.tight_layout()
        plt.savefig('photon_1d_speed_measurement.png', dpi=150, bbox_inches='tight')
        print(f"  ✓ Saved: photon_1d_speed_measurement.png")

    print(f"\n{'=' * 70}")
    print("Complete!")
    print(f"{'=' * 70}")


if __name__ == '__main__':
    main()
