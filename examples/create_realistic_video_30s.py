"""
Create 30-second video of 1D photon at realistic physical scales.
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


def initialize_traveling_wave(state, grid, wavelength, amplitude, wave_speed, center):
    """Initialize a wave packet centered at a position."""
    x = torch.arange(grid.grid_shape[0], device=state.device, dtype=state.dtype) * grid.spacing

    k = 2 * np.pi / wavelength
    omega = wave_speed * k
    sigma = 3 * wavelength / (2 * np.pi)

    envelope = amplitude * torch.exp(-((x - center) ** 2) / (2 * sigma ** 2))

    # Envelope derivative: dA/dx = -(x - center)/σ² * A(x)
    envelope_derivative = -((x - center) / (sigma ** 2)) * envelope

    state.positions[:, 3] = envelope * torch.cos(k * (x - center))

    # Velocity with both phase and envelope terms
    state.velocities[:, 3] = (
        omega * envelope * torch.sin(k * (x - center)) +
        (-wave_speed) * envelope_derivative * torch.cos(k * (x - center))
    )


def main():
    """Create 30-second video at realistic scales."""
    print("=" * 70)
    print("Creating 30-Second Video: 1D Photon at Realistic Scales")
    print("=" * 70)

    # Physical constants
    constants = PhysicalConstants()

    print(f"\nPhysical Constants:")
    print(f"  Speed of light c = {constants.c:.6e} m/s")
    print(f"  Compton wavelength λ_C = {constants.lambda_C:.6e} m")

    # Configuration
    lambda_C_multiplier = 10.0
    h = constants.lambda_C * lambda_C_multiplier

    nx = 600  # More points for smoother visualization
    domain_length = nx * h

    c = constants.c
    cfl_factor = 0.1
    dt = cfl_factor * h / c

    tension = 1.0
    mu = tension / c**2

    print(f"\nConfiguration:")
    print(f"  Domain: {nx} points × {h:.3e} m = {domain_length:.6e} m ({domain_length*1e9:.3f} nm)")
    print(f"  Time step: {dt:.6e} s")
    print(f"  Speed of light: {c:.6e} m/s")

    # Create components
    device = torch.device('cpu')
    dtype = torch.float64

    state = BraneState((nx,), Dimensionality.ONE_D, device, dtype)
    state.initialize_flat_configuration(h)

    # Set fixed boundaries
    state.set_fixed_boundaries()
    print(f"  Fixed boundaries at x=0 and x={domain_length*1e9:.3f} nm")

    grid = BraneGrid((nx,), Dimensionality.ONE_D, h, device)
    physics = LinearTensionForceComputer(tension, h)
    solver = VelocityVerletSolver(dt, mu, physics, grid)

    # Initialize wave
    print(f"\nInitializing photon wave packet...")
    wavelength = 40 * h
    amplitude = 0.1 * h
    center_position = domain_length / 3.0

    initialize_traveling_wave(state, grid, wavelength, amplitude, c, center_position)
    solver.initialize_accelerations(state)
    state.apply_fixed_boundaries()

    print(f"  Wavelength: {wavelength*1e9:.3f} nm")
    print(f"  Amplitude: {amplitude*1e12:.3f} pm")
    print(f"  Center: {center_position*1e9:.3f} nm")

    # Video parameters
    fps = 30
    duration = 30.0  # 30 seconds
    total_frames = int(fps * duration)

    # Simulation time: show multiple bounces
    crossing_time = domain_length / c
    simulation_time = 6.0 * crossing_time  # 6 crossings for good visualization

    total_steps = int(simulation_time / dt)
    steps_per_frame = max(1, total_steps // total_frames)

    print(f"\nVideo settings:")
    print(f"  Duration: {duration} s")
    print(f"  Frame rate: {fps} fps")
    print(f"  Total frames: {total_frames}")
    print(f"  Simulation time: {simulation_time:.6e} s ({simulation_time*1e15:.3f} fs)")
    print(f"  Light crossing time: {crossing_time:.6e} s ({crossing_time*1e15:.3f} fs)")
    print(f"  Total steps: {total_steps:,}")
    print(f"  Steps per frame: {steps_per_frame}")

    # Prepare
    x_coords = np.arange(nx) * h

    # Store frames
    print(f"\nRunning simulation and capturing frames...")
    frames = []
    times = []

    for frame_idx in range(total_frames):
        field = state.get_field_component(3).cpu().numpy()
        frames.append(field.copy())
        times.append(solver.time)

        if frame_idx % 100 == 0:
            energy = solver.compute_energy(state)
            print(f"  Frame {frame_idx:4d}/{total_frames}: t={solver.time:.3e}s ({solver.time*1e15:.3f} fs), E={energy['total']:.6e}J")

        for _ in range(steps_per_frame):
            solver.step(state)

    print(f"\nCreating video animation...")

    # Create figure
    fig, ax = plt.subplots(figsize=(16, 6))

    # Convert to nanometers for display
    x_nm = x_coords * 1e9

    # Initial plot
    field_nm = frames[0] * 1e9
    line, = ax.plot(x_nm, field_nm, 'b-', linewidth=2)

    # Mark fixed boundaries
    boundary_markers = ax.plot([x_nm[0], x_nm[-1]], [0, 0], 'ro',
                               markersize=10, label='Fixed boundaries', zorder=5)[0]

    ax.set_xlim(x_nm[0], x_nm[-1])
    ax.set_ylim(-1.5*amplitude*1e9, 1.5*amplitude*1e9)
    ax.set_xlabel('Position [nm]', fontsize=14)
    ax.set_ylabel('Displacement ξ [nm]', fontsize=14)
    ax.set_title(f'1D Photon at Realistic Scales (c = {c:.3e} m/s, λ_C = {constants.lambda_C*1e15:.0f} fm)',
                 fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3)
    ax.legend(loc='upper right', fontsize=12)

    # Time text (in femtoseconds)
    time_text = ax.text(0.02, 0.95, '', transform=ax.transAxes,
                       fontsize=14, verticalalignment='top',
                       bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

    # Update function
    def update(frame):
        field_nm = frames[frame] * 1e9
        line.set_ydata(field_nm)
        t_fs = times[frame] * 1e15
        time_text.set_text(f't = {t_fs:.3f} fs')
        return line, boundary_markers, time_text

    # Create animation
    anim = FuncAnimation(fig, update, frames=total_frames,
                        interval=1000/fps, blit=True)

    # Save video
    output_file = 'photon_1d_realistic_scales_30s.mp4'
    print(f"\nSaving video to {output_file}...")

    writer = FFMpegWriter(fps=fps, bitrate=3000,
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
    print(f"\nPhysical Summary:")
    print(f"  The video shows a photon wave packet bouncing between")
    print(f"  fixed boundaries separated by {domain_length*1e9:.2f} nanometers.")
    print(f"  Over {simulation_time*1e15:.1f} femtoseconds, the photon completes")
    print(f"  ~6 round trips at the speed of light.")
    print(f"  Grid spacing: {h*1e12:.2f} pm = 10 × Compton wavelength")
    print(f"  Wavelength: {wavelength*1e9:.2f} nm")


if __name__ == '__main__':
    main()
