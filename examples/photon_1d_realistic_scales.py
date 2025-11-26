"""
1D Photon with Realistic Physical Scales

Uses actual speed of light c = 299,792,458 m/s and physical length scales
based on the Compton wavelength.
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
from branesim.config.simulation_config import PhysicalConstants
from branesim.core.initial_conditions import (
    initialize_right_moving_velocities,
    measure_wave_speed,
    verify_wave_propagation,
)


def initialize_wave_shape_1d(state, grid, wavelength, amplitude, center):
    """
    Initialize ONLY the shape of a 1D wave packet.

    Velocities are initialized separately using initialize_right_moving_velocities().
    This allows dimension-independent velocity initialization.
    """
    x = torch.arange(grid.grid_shape[0], device=state.device, dtype=state.dtype) * grid.spacing

    k = 2 * np.pi / wavelength
    sigma = 3 * wavelength / (2 * np.pi)  # Width of wave packet

    # Gaussian envelope
    envelope = amplitude * torch.exp(-((x - center) ** 2) / (2 * sigma ** 2))

    # Position field only - velocities set separately
    state.positions[:, 3] = envelope * torch.cos(k * (x - center))

    print(f"  Wavelength λ = {wavelength:.6e} m ({wavelength/grid.spacing:.1f} × h)")
    print(f"  Width σ = {sigma:.6e} m")
    print(f"  Center = {center:.6e} m")
    print(f"  Amplitude = {amplitude:.6e} m")
    print(f"  Wave number k = {k:.6e} rad/m")


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
    """Run 1D photon simulation with realistic physical scales."""
    print("=" * 70)
    print("1D Photon - Realistic Physical Scales")
    print("=" * 70)

    # Physical constants
    constants = PhysicalConstants()

    print(f"\nPhysical Constants:")
    print(f"  Speed of light c = {constants.c:.6e} m/s")
    print(f"  Compton wavelength λ_C = {constants.lambda_C:.6e} m")
    print(f"  ℏ = {constants.hbar:.6e} J·s")
    print(f"  m_e = {constants.m_e:.6e} kg")

    # Configuration with realistic scales
    # Grid spacing as multiple of Compton wavelength
    lambda_C_multiplier = 10.0  # Grid spacing = 10 × λ_C
    h = constants.lambda_C * lambda_C_multiplier

    # Domain size
    nx = 400
    domain_length = nx * h

    # Wave speed (speed of light)
    c = constants.c

    # CFL condition for stability
    cfl_factor = 0.1
    dt = cfl_factor * h / c

    # Mass density: choose to give desired wave speed
    # For c = √(T/μ), with T = 1.0 N, we need μ = T/c²
    tension = 1.0  # N
    mu = tension / c**2  # kg/m (linear mass density)

    print(f"\nSimulation Configuration:")
    print(f"  Grid spacing h = {h:.6e} m ({lambda_C_multiplier:.0f} × λ_C)")
    print(f"  Domain: {nx} points × {h:.3e} m = {domain_length:.6e} m")
    print(f"  Time step dt = {dt:.6e} s")
    print(f"  CFL number = {cfl_factor:.3f}")
    print(f"  Tension T = {tension:.3f} N")
    print(f"  Linear mass density μ = {mu:.6e} kg/m")
    print(f"  Expected wave speed = √(T/μ) = {np.sqrt(tension/mu):.6e} m/s")

    # Create components
    device = torch.device('cpu')
    dtype = torch.float64

    state = BraneState((nx,), Dimensionality.ONE_D, device, dtype)
    state.initialize_flat_configuration(h)

    # Store initial positions for lateral distortion tracking
    initial_positions = state.positions.clone()

    # Set fixed boundaries
    state.set_fixed_boundaries()
    print(f"\nBoundary Conditions:")
    print(f"  Fixed boundaries at x=0 and x={domain_length:.3e} m")
    print(f"  Fixed points: {state.fixed_mask.sum().item()} / {nx}")

    grid = BraneGrid((nx,), Dimensionality.ONE_D, h, device)

    physics = LinearTensionForceComputer(tension, h)
    solver = VelocityVerletSolver(dt, mu, physics, grid)

    # Initialize wave packet (two-step process)
    print(f"\nInitializing photon wave packet...")

    # Wavelength: multiple of grid spacing for good resolution
    wavelength = 40 * h  # 40 points per wavelength

    # Amplitude: small compared to domain
    amplitude = 0.1 * h

    # Center: in the left third of domain
    center_position = domain_length / 3.0

    # Step 1: Initialize shape only
    print(f"\n[1] Initializing wave shape...")
    initialize_wave_shape_1d(state, grid, wavelength, amplitude, center_position)

    # Step 2: Initialize velocities for right-moving wave at speed c
    print(f"\n[2] Initializing velocities for right-moving wave...")
    initialize_right_moving_velocities(
        state=state,
        grid=grid,
        wave_speed=c,
        direction=None,  # Default: +x
        field_component=3
    )

    # Step 3: Compute initial accelerations
    solver.initialize_accelerations(state)
    state.apply_fixed_boundaries()

    # Initial measurements
    initial_energy = solver.compute_energy(state)
    initial_center = track_wave_center(state, grid)

    print(f"\nInitial State:")
    print(f"  Energy = {initial_energy['total']:.6e} J")
    print(f"  Wave center = {initial_center:.6e} m")

    # Run simulation
    # For realistic scales, simulate a very short time!
    # Time for light to cross domain: t = L/c
    crossing_time = domain_length / c
    simulation_time = 3.0 * crossing_time  # 3 crossings

    num_steps = int(simulation_time / dt)

    print(f"\nRunning simulation...")
    print(f"  Light crossing time = {crossing_time:.6e} s")
    print(f"  Simulation time = {simulation_time:.6e} s")
    print(f"  Number of steps = {num_steps:,}")

    # Tracking
    times = []
    centers = []
    energies = []

    # Take snapshots at regular intervals
    num_snapshots = 7
    snapshot_times = np.linspace(0, simulation_time, num_snapshots)
    snapshots = {}
    snapshots_lateral = {}  # Store lateral displacements
    snapshot_steps = {int(t / dt): t for t in snapshot_times}

    x_coords = np.arange(nx) * h

    print_interval = max(1, num_steps // 20)

    for step in range(num_steps + 1):
        if step in snapshot_steps:
            field = state.get_field_component(3).cpu().numpy()
            snapshots[snapshot_steps[step]] = field.copy()

            # Store lateral displacement (x-component)
            lateral_disp = (state.positions[:, 0] - initial_positions[:, 0]).cpu().numpy()
            snapshots_lateral[snapshot_steps[step]] = lateral_disp.copy()

        if step % max(1, num_steps // 100) == 0:  # Track 100 points
            center = track_wave_center(state, grid)
            energy = solver.compute_energy(state)
            times.append(solver.time)
            centers.append(center)
            energies.append(energy['total'])

        if step % print_interval == 0:
            print(f"  Step {step:8d}/{num_steps}: t={solver.time:.6e}s, "
                  f"center={center:.3e}m, E={energy['total']:.6e}J")

        if step < num_steps:
            solver.step(state)

    # Final analysis
    final_energy = solver.compute_energy(state)
    final_center = track_wave_center(state, grid)

    distance_traveled = final_center - initial_center
    measured_speed = distance_traveled / simulation_time
    speed_error = abs(measured_speed - c) / c

    energy_drift = abs(final_energy['total'] - initial_energy['total']) / initial_energy['total']

    print(f"\n{'=' * 70}")
    print("Results:")
    print(f"{'=' * 70}")

    # Use the verification function
    verify_wave_propagation(
        state=state,
        grid=grid,
        wave_speed_expected=c,
        time_elapsed=simulation_time,
        initial_center_x=initial_center,
        field_component=3
    )

    print(f"\nEnergy Conservation:")
    print(f"  Initial: {initial_energy['total']:.6e} J")
    print(f"  Final:   {final_energy['total']:.6e} J")
    print(f"  Drift:   {energy_drift:.6e} ({energy_drift*100:.6f}%)")

    if energy_drift < 1e-4:
        print(f"  ✓ Good energy conservation")

    # Visualization
    print(f"\nCreating plots...")

    fig, axes = plt.subplots(num_snapshots, 1, figsize=(14, 12))
    fig.suptitle(f'1D Photon at Realistic Scales (c = {c:.3e} m/s)',
                 fontsize=16, fontweight='bold')

    for idx, (step, t) in enumerate(snapshot_steps.items()):
        if t in snapshots:
            field = snapshots[t]

            # Convert to nanometers for better readability
            x_nm = x_coords * 1e9
            field_nm = field * 1e9

            axes[idx].plot(x_nm, field_nm, 'b-', linewidth=2)
            axes[idx].plot([x_nm[0], x_nm[-1]], [0, 0], 'ro',
                          markersize=8, label='Fixed boundaries' if idx == 0 else '')

            axes[idx].set_ylabel('ξ [nm]', fontsize=11)
            axes[idx].set_xlim(x_nm[0], x_nm[-1])
            axes[idx].set_ylim(-1.5*amplitude*1e9, 1.5*amplitude*1e9)
            axes[idx].grid(True, alpha=0.3)

            # Time in femtoseconds
            t_fs = t * 1e15
            axes[idx].text(0.02, 0.95, f't = {t_fs:.3f} fs',
                          transform=axes[idx].transAxes,
                          fontsize=12, verticalalignment='top',
                          bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

            if idx == 0:
                axes[idx].legend(loc='upper right', fontsize=10)

            if idx == num_snapshots - 1:
                axes[idx].set_xlabel('Position [nm]', fontsize=12)

    plt.tight_layout()
    plt.savefig('photon_1d_realistic_scales_propagation.png', dpi=150, bbox_inches='tight')
    print(f"  ✓ Saved: photon_1d_realistic_scales_propagation.png")

    # Analysis plots
    fig2, axes2 = plt.subplots(2, 1, figsize=(12, 8))

    # Position vs time (in nm and fs)
    times_fs = np.array(times) * 1e15
    centers_nm = np.array(centers) * 1e9

    axes2[0].plot(times_fs, centers_nm, 'b-', linewidth=2, label='Wave center')
    axes2[0].set_xlabel('Time [fs]', fontsize=12)
    axes2[0].set_ylabel('Wave Center [nm]', fontsize=12)
    axes2[0].set_title('Wave Propagation at Speed of Light', fontsize=14, fontweight='bold')
    axes2[0].grid(True, alpha=0.3)
    axes2[0].legend(fontsize=10)

    # Energy conservation
    energy_array = np.array(energies)
    axes2[1].plot(times_fs, energy_array / initial_energy['total'], 'g-', linewidth=2)
    axes2[1].axhline(y=1.0, color='r', linestyle='--', linewidth=1, alpha=0.5)
    axes2[1].set_xlabel('Time [fs]', fontsize=12)
    axes2[1].set_ylabel('E(t) / E(0)', fontsize=12)
    axes2[1].set_title('Energy Conservation', fontsize=14, fontweight='bold')
    axes2[1].grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig('photon_1d_realistic_scales_analysis.png', dpi=150, bbox_inches='tight')
    print(f"  ✓ Saved: photon_1d_realistic_scales_analysis.png")

    # Lateral distortion visualization
    print(f"\nCreating lateral distortion plots...")

    fig3, axes3 = plt.subplots(num_snapshots, 1, figsize=(14, 12))
    fig3.suptitle(f'1D Photon - Lateral Distortion (Left/Right Movement)',
                 fontsize=16, fontweight='bold')

    # Find max lateral displacement for consistent scaling
    max_lateral_disp = max([np.abs(snapshots_lateral[t]).max() for t in snapshots_lateral.keys()])

    for idx, (step, t) in enumerate(snapshot_steps.items()):
        if t in snapshots_lateral:
            lateral_disp = snapshots_lateral[t]

            # Convert to picometers for better readability (lateral displacement is tiny)
            x_nm = x_coords * 1e9
            lateral_disp_pm = lateral_disp * 1e12

            axes3[idx].plot(x_nm, lateral_disp_pm, 'r-', linewidth=2, label='x-displacement')
            axes3[idx].axhline(y=0, color='k', linestyle='--', linewidth=0.5, alpha=0.5)
            axes3[idx].plot([x_nm[0], x_nm[-1]], [0, 0], 'ro',
                          markersize=8, label='Fixed boundaries' if idx == 0 else '')

            axes3[idx].set_ylabel('Δx [pm]', fontsize=11)
            axes3[idx].set_xlim(x_nm[0], x_nm[-1])
            axes3[idx].set_ylim(-1.5*max_lateral_disp*1e12, 1.5*max_lateral_disp*1e12)
            axes3[idx].grid(True, alpha=0.3)

            # Time in femtoseconds
            t_fs = t * 1e15
            axes3[idx].text(0.02, 0.95, f't = {t_fs:.3f} fs',
                          transform=axes3[idx].transAxes,
                          fontsize=12, verticalalignment='top',
                          bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

            if idx == 0:
                axes3[idx].legend(loc='upper right', fontsize=10)

            if idx == num_snapshots - 1:
                axes3[idx].set_xlabel('Position [nm]', fontsize=12)

    plt.tight_layout()
    plt.savefig('photon_1d_realistic_scales_lateral_distortion.png', dpi=150, bbox_inches='tight')
    print(f"  ✓ Saved: photon_1d_realistic_scales_lateral_distortion.png")

    print(f"\n{'=' * 70}")
    print("Simulation complete!")
    print(f"{'=' * 70}")
    print(f"\nPhysical Interpretation:")
    print(f"  Domain size: {domain_length*1e9:.3f} nm = {domain_length/constants.lambda_C:.0f} × λ_C")
    print(f"  Wavelength: {wavelength*1e9:.3f} nm")
    print(f"  Simulation time: {simulation_time*1e15:.3f} femtoseconds")
    print(f"  Distance traveled: {distance_traveled*1e9:.3f} nm")
    print(f"  Number of reflections: ~{int(distance_traveled / domain_length)}")


if __name__ == '__main__':
    main()
