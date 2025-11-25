"""
2D Photon in Tunnel Geometry

Demonstrates photon propagating horizontally through a tunnel.
Top and bottom boundaries are fixed (tunnel walls), left and right are fixed (end walls).
The photon moves at the speed of light in the +x direction.
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
from branesim.config.simulation_config import SimulationConfig


def initialize_wave_packet_2d(state, grid, amplitude, wave_speed, center_x, center_y, width_x, width_y, wavelength):
    """
    Initialize a localized 2D wave packet moving in the +x direction.

    The packet is localized in both x and y (2D Gaussian envelope) with a sinusoidal carrier.
    This matches the 1D animate_1d_wave.py implementation.

    Args:
        state: BraneState
        grid: BraneGrid
        amplitude: Peak amplitude [m]
        wave_speed: Wave speed [m/s]
        center_x: Center x coordinate [m]
        center_y: Center y coordinate [m]
        width_x: Gaussian width σ_x [m]
        width_y: Gaussian width σ_y [m]
        wavelength: Carrier wavelength [m]
    """
    # Get spatial coordinates
    coords = grid.get_spatial_coordinates()
    x = coords[:, 0]
    y = coords[:, 1]

    # 2D Gaussian envelope localized in both x and y
    envelope = amplitude * torch.exp(
        -((x - center_x) ** 2) / (2 * width_x ** 2)
        -((y - center_y) ** 2) / (2 * width_y ** 2)
    )

    # Wave number and frequency
    k = 2 * np.pi / wavelength
    omega = wave_speed * k

    # Envelope derivative in x: dA/dx = -(x - center_x)/σ_x² * A(x,y)
    envelope_derivative = -((x - center_x) / (width_x ** 2)) * envelope

    # Position: ξ(x,y,0) = A(x,y) * cos(kx)
    state.positions[:, 3] = envelope * torch.cos(k * x)

    # Velocity: For right-moving wave packet
    # ∂ξ/∂t = -c * dA/dx * cos(kx) + ω * A(x,y) * sin(kx)
    state.velocities[:, 3] = (
        -wave_speed * envelope_derivative * torch.cos(k * x) +
        omega * envelope * torch.sin(k * x)
    )

    print(f"  Center: ({center_x:.4f}, {center_y:.4f}) m")
    print(f"  Width (σ_x, σ_y): ({width_x:.4f}, {width_y:.4f}) m")
    print(f"  Wavelength: {wavelength:.4f} m")
    print(f"  Amplitude: {amplitude:.6f} m")


def main():
    """Run 2D photon simulation in tunnel geometry."""
    print("=" * 70)
    print("2D Photon in Tunnel Geometry")
    print("=" * 70)

    # Tunnel configuration: wide but narrow
    nx = 200  # Long tunnel
    ny = 50   # Increased height to see particle better
    wave_speed = 1.0  # m/s
    spacing = 0.01  # m
    cfl_factor = 0.05  # Conservative for 2D

    config = SimulationConfig.for_2d_test(
        nx=nx, ny=ny,
        wave_speed=wave_speed,
        spacing=spacing,
        cfl_factor=cfl_factor
    )

    print(f"\nTunnel Configuration:")
    print(f"  Domain: {nx * spacing:.2f} × {ny * spacing:.2f} m ({nx}×{ny} points)")
    print(f"  Aspect ratio: {nx/ny:.1f}:1 (tunnel geometry)")
    print(f"  Grid spacing h = {spacing:.5f} m")
    print(f"  Wave speed c = {wave_speed:.6f} m/s")
    print(f"  Time step dt = {config.time_step:.6f} s")
    print(f"  CFL number = {cfl_factor:.3f}")

    # Create components
    device = torch.device('cpu')
    dtype = torch.float64

    state = BraneState((nx, ny), Dimensionality.TWO_D, device, dtype)
    state.initialize_flat_configuration(spacing)

    # Set fixed boundaries (all edges)
    state.set_fixed_boundaries()
    print(f"\nBoundary Conditions:")
    print(f"  Fixed points: {state.fixed_mask.sum().item()} / {nx * ny}")
    print(f"  Top/bottom walls: FIXED (tunnel)")
    print(f"  Left/right walls: FIXED (end caps)")

    grid = BraneGrid((nx, ny), Dimensionality.TWO_D, spacing, device)

    # Use LinearTensionForceComputer for stability
    physics = LinearTensionForceComputer(config.tension, spacing)
    solver = VelocityVerletSolver(config.time_step, config.mass_density, physics, grid)

    # Initialize localized wave packet (particle-like) on the left
    print(f"\nInitializing localized 2D wave packet (particle)...")
    amplitude = 0.002  # m
    center_x = 0.4  # Start at 0.4m from left edge
    center_y = (ny - 1) * spacing / 2.0  # Center in tunnel height
    width_x = 0.3  # m (pulse width in x - wider like 1D example)
    width_y = 0.05  # m (pulse width in y - narrow to see localization)
    wavelength = 0.2  # m (carrier wavelength)

    initialize_wave_packet_2d(state, grid, amplitude, wave_speed, center_x, center_y, width_x, width_y, wavelength)

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

    # For static plots
    snapshot_times = [0.0, 0.3, 0.6, 0.9, 1.2, 1.5]
    snapshots = {}
    snapshot_steps = {int(t / config.time_step): t for t in snapshot_times}

    # For animation - save every 10 steps
    animation_frames = []
    animation_times = []
    frame_interval = 10

    print_interval = num_steps // 10

    for step in range(num_steps + 1):
        if step in snapshot_steps:
            field = state.get_field_component(3).cpu().numpy()
            snapshots[snapshot_steps[step]] = field.copy()

        # Save frames for animation
        if step % frame_interval == 0:
            field = state.get_field_component(3).cpu().numpy()
            animation_frames.append(field.copy())
            animation_times.append(solver.time)

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
    print(f"  Photon propagates horizontally through tunnel at c = {wave_speed} m/s")
    print(f"  Tunnel length: {(nx-1)*spacing:.2f} m")
    print(f"  Tunnel height: {(ny-1)*spacing:.2f} m")

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

    fig, axes = plt.subplots(len(snapshot_times), 1, figsize=(16, 12))
    fig.suptitle('2D Photon in Tunnel (c = 1.0 m/s)',
                 fontsize=16, fontweight='bold')

    # Create coordinate arrays
    x_coords = np.arange(nx) * spacing
    y_coords = np.arange(ny) * spacing
    X, Y = np.meshgrid(x_coords, y_coords, indexing='ij')

    for idx, t in enumerate(snapshot_times):
        if t in snapshots:
            field = snapshots[t].reshape(nx, ny)

            # Plot heatmap (using 'RdBu_r' colormap for oscillating wave packet)
            im = axes[idx].imshow(field.T, origin='lower',
                                 extent=[0, x_coords[-1], 0, y_coords[-1]],
                                 cmap='RdBu_r', vmin=-amplitude*1.2, vmax=amplitude*1.2,
                                 aspect='auto')

            axes[idx].set_ylabel('y [m]', fontsize=11)
            axes[idx].set_xlim(0, x_coords[-1])
            axes[idx].set_ylim(0, y_coords[-1])

            # Add text label
            axes[idx].text(0.02, 0.95, f't = {t:.2f} s',
                          transform=axes[idx].transAxes,
                          fontsize=12, verticalalignment='top',
                          bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

            # Add colorbar
            plt.colorbar(im, ax=axes[idx], label='ξ [m]', fraction=0.046, pad=0.04)

            if idx == len(snapshot_times) - 1:
                axes[idx].set_xlabel('x [m]', fontsize=12)

    plt.tight_layout()
    plt.savefig('photon_2d_tunnel_propagation.png', dpi=150, bbox_inches='tight')
    print(f"  ✓ Saved: photon_2d_tunnel_propagation.png")

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
    plt.savefig('photon_2d_tunnel_energy.png', dpi=150, bbox_inches='tight')
    print(f"  ✓ Saved: photon_2d_tunnel_energy.png")

    # Create animation
    print(f"\nCreating animation...")
    print(f"  Total frames: {len(animation_frames)}")

    fig_anim, ax_anim = plt.subplots(figsize=(12, 3))

    # Initial frame
    field_init = animation_frames[0].reshape(nx, ny)
    im_anim = ax_anim.imshow(field_init.T, origin='lower',
                             extent=[0, x_coords[-1], 0, y_coords[-1]],
                             cmap='RdBu_r', vmin=-amplitude*1.2, vmax=amplitude*1.2,
                             aspect='auto', animated=True)

    ax_anim.set_xlabel('x [m]', fontsize=12)
    ax_anim.set_ylabel('y [m]', fontsize=12)
    time_text = ax_anim.text(0.02, 0.95, '', transform=ax_anim.transAxes,
                            fontsize=12, verticalalignment='top',
                            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

    plt.colorbar(im_anim, ax=ax_anim, label='ξ [m]', fraction=0.046, pad=0.04)
    ax_anim.set_title('2D Photon in Tunnel (c = 1.0 m/s)', fontsize=14, fontweight='bold')

    def animate(frame_idx):
        """Update function for animation."""
        field = animation_frames[frame_idx].reshape(nx, ny)
        im_anim.set_array(field.T)
        time_text.set_text(f't = {animation_times[frame_idx]:.3f} s')
        return [im_anim, time_text]

    anim = FuncAnimation(fig_anim, animate, frames=len(animation_frames),
                        interval=50, blit=True, repeat=True)

    # Save animation
    writer = FFMpegWriter(fps=20, bitrate=2000)
    anim.save('photon_2d_tunnel.mp4', writer=writer, dpi=100)
    print(f"  ✓ Saved: photon_2d_tunnel.mp4")

    plt.close(fig_anim)

    print(f"\n{'=' * 70}")
    print("Simulation complete!")
    print(f"{'=' * 70}")


if __name__ == '__main__':
    main()
