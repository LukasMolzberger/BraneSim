"""
Create 10-second video of 1D photon propagation.
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


def initialize_traveling_wave(state, grid, wavelength, amplitude, wave_speed):
    """Initialize a right-moving sine wave."""
    x = torch.arange(grid.grid_shape[0], device=state.device, dtype=state.dtype) * grid.spacing
    k = 2 * np.pi / wavelength
    omega = wave_speed * k

    state.positions[:, 3] = amplitude * torch.sin(k * x)
    state.velocities[:, 3] = -amplitude * omega * torch.cos(k * x)


def main():
    """Create 10-second video of photon propagation."""
    print("=" * 70)
    print("Creating 10-Second Video: 1D Photon Propagation")
    print("=" * 70)

    # Configuration
    nx = 600  # Longer domain for better visualization
    wave_speed = 1.0  # m/s
    spacing = 0.01  # m
    cfl_factor = 0.1

    mu = 1.0  # kg/m
    tension = mu * wave_speed**2
    dt = cfl_factor * spacing / wave_speed

    print(f"\nConfiguration:")
    print(f"  Domain: {nx * spacing:.2f} m ({nx} points)")
    print(f"  Wave speed: {wave_speed} m/s")
    print(f"  Time step: {dt:.6f} s")
    print(f"  CFL: {cfl_factor}")

    # Create components
    device = torch.device('cpu')
    dtype = torch.float64

    state = BraneState((nx,), Dimensionality.ONE_D, device, dtype)
    state.initialize_flat_configuration(spacing)

    grid = BraneGrid((nx,), Dimensionality.ONE_D, spacing, device)
    physics = LinearTensionForceComputer(tension, spacing)
    solver = VelocityVerletSolver(dt, mu, physics, grid)

    # Initialize wave
    print(f"\nInitializing photon wave packet...")
    wavelength = 0.5  # m
    amplitude = 0.003  # m - larger for visibility

    initialize_traveling_wave(state, grid, wavelength, amplitude, wave_speed)
    solver.initialize_accelerations(state)

    print(f"  Wavelength: {wavelength} m")
    print(f"  Amplitude: {amplitude} m")

    # Video parameters
    fps = 30
    duration = 10.0  # seconds
    total_frames = int(fps * duration)
    simulation_time = 4.0  # Run for 4 seconds of simulation time

    # Calculate steps per frame
    total_steps = int(simulation_time / dt)
    steps_per_frame = max(1, total_steps // total_frames)

    print(f"\nVideo settings:")
    print(f"  Duration: {duration} s")
    print(f"  Frame rate: {fps} fps")
    print(f"  Total frames: {total_frames}")
    print(f"  Simulation time: {simulation_time} s")
    print(f"  Steps per frame: {steps_per_frame}")

    # Prepare for animation
    x_coords = np.arange(nx) * spacing

    # Store frames
    print(f"\nRunning simulation and capturing frames...")
    frames = []
    times = []

    for frame_idx in range(total_frames):
        # Store current state
        field = state.get_field_component(3).cpu().numpy()
        frames.append(field.copy())
        times.append(solver.time)

        # Progress indicator
        if frame_idx % 30 == 0:
            energy = solver.compute_energy(state)
            print(f"  Frame {frame_idx}/{total_frames}: t={solver.time:.3f}s, E={energy['total']:.6e}J")

        # Advance simulation
        for _ in range(steps_per_frame):
            solver.step(state)

    print(f"\nCreating video animation...")

    # Create figure
    fig, ax = plt.subplots(figsize=(16, 6))

    # Initial plot
    line, = ax.plot(x_coords, frames[0], 'b-', linewidth=2)
    ax.set_xlim(0, x_coords[-1])
    ax.set_ylim(-1.5*amplitude, 1.5*amplitude)
    ax.set_xlabel('Position [m]', fontsize=14)
    ax.set_ylabel('Displacement ξ [m]', fontsize=14)
    ax.set_title('1D Photon Propagation (c = 1.0 m/s)', fontsize=16, fontweight='bold')
    ax.grid(True, alpha=0.3)

    # Time text
    time_text = ax.text(0.02, 0.95, '', transform=ax.transAxes,
                       fontsize=14, verticalalignment='top',
                       bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

    # Update function
    def update(frame):
        line.set_ydata(frames[frame])
        time_text.set_text(f't = {times[frame]:.3f} s')
        return line, time_text

    # Create animation
    anim = FuncAnimation(fig, update, frames=total_frames,
                        interval=1000/fps, blit=True)

    # Save video
    output_file = 'photon_1d_propagation_10s.mp4'
    print(f"\nSaving video to {output_file}...")

    writer = FFMpegWriter(fps=fps, bitrate=2000,
                         codec='libx264',
                         extra_args=['-pix_fmt', 'yuv420p'])

    anim.save(output_file, writer=writer, dpi=120)

    file_size = os.path.getsize(output_file) / (1024 * 1024)
    print(f"  ✓ Video saved: {output_file}")
    print(f"  ✓ File size: {file_size:.2f} MB")
    print(f"  ✓ Resolution: 1920x720")
    print(f"  ✓ Duration: {duration} seconds")
    print(f"  ✓ Frame rate: {fps} fps")

    plt.close()

    print(f"\n{'=' * 70}")
    print("Video creation complete!")
    print(f"{'=' * 70}")


if __name__ == '__main__':
    main()
