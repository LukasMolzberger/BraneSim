"""
2D Photon in Tunnel Geometry with Realistic Physical Scales

Uses actual speed of light c = 299,792,458 m/s and physical length scales
based on the Compton wavelength.
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
from branesim.physics.linear_tension_forces import LinearTensionForceComputer
from branesim.config.simulation_config import PhysicalConstants


def initialize_plane_wave_packet_x(state, grid, amplitude, wave_speed, center_x, width_x, wavelength):
    """
    Initialize a wave packet moving in the +x direction with standing wave mode in y.

    Uses fundamental standing wave mode sin(πy/L_y) in y-direction to naturally
    satisfy fixed boundary conditions at y=0 and y=L_y. Localized Gaussian in x.

    Args:
        state: BraneState
        grid: BraneGrid
        amplitude: Peak amplitude [m]
        wave_speed: Wave speed [m/s]
        center_x: Center x coordinate [m]
        width_x: Gaussian width σ_x [m]
        wavelength: Carrier wavelength [m]
    """
    # Get spatial coordinates
    coords = grid.get_spatial_coordinates()
    x = coords[:, 0]
    y = coords[:, 1]

    # Domain dimensions
    nx, ny = grid.grid_shape
    domain_length_y = (ny - 1) * grid.spacing

    # Gaussian envelope in x
    envelope_x = amplitude * torch.exp(-((x - center_x) ** 2) / (2 * width_x ** 2))

    # Standing wave mode in y: sin(πy/L_y) satisfies ξ(y=0)=0 and ξ(y=L_y)=0
    y_mode = torch.sin(np.pi * y / domain_length_y)

    # Full envelope: localized in x, standing wave in y
    envelope = envelope_x * y_mode

    # Wave numbers and frequency (same as 1D case)
    k = 2 * np.pi / wavelength
    omega = wave_speed * k

    # Envelope derivative in x: ∂A/∂x = -(x - center_x)/σ_x² * A(x,y)
    envelope_derivative_x = -((x - center_x) / (width_x ** 2)) * envelope

    # Position: ξ(x,y,0) = A(x) * sin(πy/L_y) * cos(k(x-x₀))
    state.positions[:, 3] = envelope * torch.cos(k * (x - center_x))

    # Velocity: For right-moving wave packet (same formula as 1D case)
    # ∂ξ/∂t = sin(πy/L_y) * [ω*A(x)*sin(k(x-x₀)) - c*A'(x)*cos(k(x-x₀))]
    state.velocities[:, 3] = (
        omega * envelope * torch.sin(k * (x - center_x)) +
        (-wave_speed) * envelope_derivative_x * torch.cos(k * (x - center_x))
    )

    print(f"  Wavelength λ = {wavelength:.6e} m ({wavelength/grid.spacing:.1f} × h)")
    print(f"  Width σ_x = {width_x:.6e} m")
    print(f"  Center x = {center_x:.6e} m")
    print(f"  Amplitude = {amplitude:.6e} m")
    print(f"  Wave number k = {k:.6e} rad/m")
    print(f"  Angular frequency ω = {omega:.6e} rad/s")
    print(f"  Y-mode: sin(πy/{domain_length_y:.3e}) - fundamental standing wave")
    print(f"  Satisfies fixed boundary conditions at y=0 and y=L_y")


def main():
    """Run 2D photon simulation in tunnel geometry with realistic scales."""
    print("=" * 70)
    print("2D Photon in Tunnel - Realistic Physical Scales")
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

    # Domain size - tunnel geometry (long in x, narrow in y)
    nx = 400  # Long tunnel
    ny = 50   # Narrow tunnel
    domain_length_x = nx * h
    domain_length_y = ny * h

    # Wave speed (speed of light)
    c = constants.c

    # CFL condition for stability
    cfl_factor = 0.1
    dt = cfl_factor * h / c

    # Mass density: For 2D, we need surface mass density σ = T/c²
    # to get wave speed c = √(T/σ)
    tension = 1.0  # N
    sigma = tension / c**2  # kg/m² (surface mass density)

    print(f"\nSimulation Configuration:")
    print(f"  Grid spacing h = {h:.6e} m ({lambda_C_multiplier:.0f} × λ_C)")
    print(f"  Domain: {nx} × {ny} points")
    print(f"  Domain size: {domain_length_x:.6e} × {domain_length_y:.6e} m")
    print(f"  Aspect ratio: {nx/ny:.1f}:1 (tunnel geometry)")
    print(f"  Time step dt = {dt:.6e} s")
    print(f"  CFL number = {cfl_factor:.3f}")
    print(f"  Tension T = {tension:.3f} N")
    print(f"  Surface mass density σ = {sigma:.6e} kg/m²")
    print(f"  Expected wave speed = √(T/σ) = {np.sqrt(tension/sigma):.6e} m/s")

    # Create components
    device = torch.device('cpu')
    dtype = torch.float64

    state = BraneState((nx, ny), Dimensionality.TWO_D, device, dtype)
    state.initialize_flat_configuration(h)

    # Set fixed boundaries (all edges)
    state.set_fixed_boundaries()
    print(f"\nBoundary Conditions:")
    print(f"  Fixed points: {state.fixed_mask.sum().item()} / {nx * ny}")
    print(f"  Top/bottom walls: FIXED (tunnel)")
    print(f"  Left/right walls: FIXED (end caps)")

    grid = BraneGrid((nx, ny), Dimensionality.TWO_D, h, device)

    physics = LinearTensionForceComputer(tension, h)
    solver = VelocityVerletSolver(dt, sigma, physics, grid)

    # Initialize wave packet
    print(f"\nInitializing photon wave packet...")

    # Wavelength: multiple of grid spacing for good resolution
    wavelength = 40 * h  # 40 points per wavelength

    # Amplitude: small compared to domain
    amplitude = 0.1 * h

    # Width: several wavelengths (only in x, plane wave in y)
    width_x = 3 * wavelength / (2 * np.pi)  # Similar to 1D example

    # Center: in the left third of domain
    center_x = domain_length_x / 3.0

    initialize_plane_wave_packet_x(state, grid, amplitude, c, center_x, width_x, wavelength)

    solver.initialize_accelerations(state)
    state.apply_fixed_boundaries()

    # Initial measurements
    initial_energy = solver.compute_energy(state)

    print(f"\nInitial State:")
    print(f"  Energy = {initial_energy['total']:.6e} J")

    # Run simulation
    # Time for light to cross domain: t = L/c
    crossing_time = domain_length_x / c
    simulation_time = 3.0 * crossing_time  # 3 crossings

    num_steps = int(simulation_time / dt)

    print(f"\nRunning simulation...")
    print(f"  Light crossing time = {crossing_time:.6e} s")
    print(f"  Simulation time = {simulation_time:.6e} s")
    print(f"  Number of steps = {num_steps:,}")

    # Tracking
    times = []
    energies = []

    # Take snapshots at regular intervals
    num_snapshots = 7
    snapshot_times = np.linspace(0, simulation_time, num_snapshots)
    snapshots = {}
    snapshot_steps = {int(t / dt): t for t in snapshot_times}

    # For animation - save every N steps
    animation_frames = []
    animation_times = []
    frame_interval = max(1, num_steps // 300)  # ~300 frames total

    x_coords = np.arange(nx) * h
    y_coords = np.arange(ny) * h

    print_interval = max(1, num_steps // 20)

    for step in range(num_steps + 1):
        if step in snapshot_steps:
            field = state.get_field_component(3).cpu().numpy()
            snapshots[snapshot_steps[step]] = field.copy()

        # Save frames for animation
        if step % frame_interval == 0:
            field = state.get_field_component(3).cpu().numpy()
            animation_frames.append(field.copy())
            animation_times.append(solver.time)

        if step % max(1, num_steps // 100) == 0:  # Track 100 points
            energy = solver.compute_energy(state)
            times.append(solver.time)
            energies.append(energy['total'])

        if step % print_interval == 0:
            print(f"  Step {step:8d}/{num_steps}: t={solver.time:.6e}s, "
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
    print(f"  Photon propagates at c = {c:.6e} m/s")
    print(f"  Tunnel length: {domain_length_x:.6e} m")
    print(f"  Tunnel height: {domain_length_y:.6e} m")

    print(f"\nEnergy Conservation:")
    print(f"  Initial: {initial_energy['total']:.6e} J")
    print(f"  Final:   {final_energy['total']:.6e} J")
    print(f"  Drift:   {energy_drift:.6e} ({energy_drift*100:.6f}%)")

    if energy_drift < 1e-2:
        print(f"  ✓ Good energy conservation")
    else:
        print(f"  ⚠ Energy drift: {energy_drift*100:.2f}%")

    # Visualization
    print(f"\nCreating plots...")

    fig, axes = plt.subplots(num_snapshots, 1, figsize=(16, 12))
    fig.suptitle(f'2D Photon in Tunnel (c = {c:.3e} m/s)',
                 fontsize=16, fontweight='bold')

    # Convert to nanometers for better readability
    x_nm = x_coords * 1e9
    y_nm = y_coords * 1e9
    amplitude_nm = amplitude * 1e9

    for idx, (step, t) in enumerate(snapshot_steps.items()):
        if t in snapshots:
            field = snapshots[t].reshape(nx, ny)
            field_nm = field * 1e9

            # Plot heatmap
            im = axes[idx].imshow(field_nm.T, origin='lower',
                                 extent=[x_nm[0], x_nm[-1], y_nm[0], y_nm[-1]],
                                 cmap='RdBu_r', vmin=-amplitude_nm*1.2, vmax=amplitude_nm*1.2,
                                 aspect='auto')

            axes[idx].set_ylabel('y [nm]', fontsize=11)
            axes[idx].set_xlim(x_nm[0], x_nm[-1])
            axes[idx].set_ylim(y_nm[0], y_nm[-1])

            # Time in femtoseconds
            t_fs = t * 1e15
            axes[idx].text(0.02, 0.95, f't = {t_fs:.3f} fs',
                          transform=axes[idx].transAxes,
                          fontsize=12, verticalalignment='top',
                          bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

            # Add colorbar
            plt.colorbar(im, ax=axes[idx], label='ξ [nm]', fraction=0.046, pad=0.04)

            if idx == num_snapshots - 1:
                axes[idx].set_xlabel('x [nm]', fontsize=12)

    plt.tight_layout()
    plt.savefig('photon_2d_tunnel_realistic_scales_propagation.png', dpi=150, bbox_inches='tight')
    print(f"  ✓ Saved: photon_2d_tunnel_realistic_scales_propagation.png")

    # Energy conservation plot
    fig2, ax = plt.subplots(figsize=(10, 6))

    times_fs = np.array(times) * 1e15
    energy_array = np.array(energies)
    ax.plot(times_fs, energy_array / initial_energy['total'], 'g-', linewidth=2)
    ax.axhline(y=1.0, color='r', linestyle='--', linewidth=1, alpha=0.5)
    ax.set_xlabel('Time [fs]', fontsize=12)
    ax.set_ylabel('E(t) / E(0)', fontsize=12)
    ax.set_title('Energy Conservation', fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig('photon_2d_tunnel_realistic_scales_energy.png', dpi=150, bbox_inches='tight')
    print(f"  ✓ Saved: photon_2d_tunnel_realistic_scales_energy.png")

    # Create animation
    print(f"\nCreating animation...")
    print(f"  Total frames: {len(animation_frames)}")

    fig_anim, ax_anim = plt.subplots(figsize=(12, 3))

    # Initial frame
    field_init = animation_frames[0].reshape(nx, ny) * 1e9
    im_anim = ax_anim.imshow(field_init.T, origin='lower',
                             extent=[x_nm[0], x_nm[-1], y_nm[0], y_nm[-1]],
                             cmap='RdBu_r', vmin=-amplitude_nm*1.2, vmax=amplitude_nm*1.2,
                             aspect='auto', animated=True)

    ax_anim.set_xlabel('x [nm]', fontsize=12)
    ax_anim.set_ylabel('y [nm]', fontsize=12)
    time_text = ax_anim.text(0.02, 0.95, '', transform=ax_anim.transAxes,
                            fontsize=12, verticalalignment='top',
                            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

    plt.colorbar(im_anim, ax=ax_anim, label='ξ [nm]', fraction=0.046, pad=0.04)
    ax_anim.set_title('2D Photon in Tunnel (c = 299,792,458 m/s)', fontsize=14, fontweight='bold')

    def animate(frame_idx):
        """Update function for animation."""
        field = animation_frames[frame_idx].reshape(nx, ny) * 1e9
        im_anim.set_array(field.T)
        t_fs = animation_times[frame_idx] * 1e15
        time_text.set_text(f't = {t_fs:.3f} fs')
        return [im_anim, time_text]

    anim = FuncAnimation(fig_anim, animate, frames=len(animation_frames),
                        interval=50, blit=True, repeat=True)

    # Save animation
    writer = FFMpegWriter(fps=20, bitrate=2000)
    anim.save('photon_2d_tunnel_realistic_scales.mp4', writer=writer, dpi=100)
    print(f"  ✓ Saved: photon_2d_tunnel_realistic_scales.mp4")

    plt.close(fig_anim)

    print(f"\n{'=' * 70}")
    print("Simulation complete!")
    print(f"{'=' * 70}")
    print(f"\nPhysical Interpretation:")
    print(f"  Domain size: {domain_length_x*1e9:.3f} × {domain_length_y*1e9:.3f} nm")
    print(f"  Domain size: {domain_length_x/constants.lambda_C:.0f} × {domain_length_y/constants.lambda_C:.0f} λ_C")
    print(f"  Wavelength: {wavelength*1e9:.3f} nm")
    print(f"  Simulation time: {simulation_time*1e15:.3f} femtoseconds")


if __name__ == '__main__':
    main()