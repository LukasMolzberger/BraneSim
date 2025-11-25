"""
Create 10-second video of 1D wave propagation.

Exports as high-quality MP4 video at 30 fps.
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


def initialize_wave_packet(
    state: BraneState,
    grid: BraneGrid,
    center: float,
    width: float,
    amplitude: float,
    wavelength: float,
    wave_speed: float
):
    """
    Initialize wave packet with sinusoidal carrier.

    Args:
        state: BraneState to initialize
        grid: BraneGrid with spatial coordinates
        center: Center position [m]
        width: Gaussian width σ [m]
        amplitude: Peak amplitude [m]
        wavelength: Carrier wavelength [m]
        wave_speed: Wave speed for velocity initialization [m/s]
    """
    x = (torch.arange(grid.grid_shape[0], device=state.device, dtype=state.dtype) *
         grid.spacing)

    envelope = amplitude * torch.exp(-((x - center) ** 2) / (2 * width ** 2))

    k = 2 * np.pi / wavelength
    omega = wave_speed * k

    state.positions[:, 3] = envelope * torch.cos(k * x)
    state.velocities[:, 3] = -envelope * omega * torch.sin(k * x)


def main():
    """Create 10-second video of 1D wave propagation."""
    print("=" * 70)
    print("Creating 10-Second 1D Wave Propagation Video")
    print("=" * 70)

    # Video parameters
    video_duration = 10.0  # seconds
    fps = 30  # frames per second
    total_frames = int(video_duration * fps)  # 300 frames

    print(f"\nVideo settings:")
    print(f"  Duration: {video_duration} seconds")
    print(f"  Frame rate: {fps} fps")
    print(f"  Total frames: {total_frames}")

    # Configuration - longer domain and more steps
    config = SimulationConfig.for_1d_test(
        nx=500,
        wave_speed=1.0,  # m/s
        spacing=0.01,    # m
        cfl_factor=0.4,
        device='cpu',
        dtype='float32'
    )

    print(f"\nSimulation setup:")
    print(f"  Domain length: {config.grid_shape[0] * config.grid_spacing:.2f} m")
    print(f"  Wave speed: {config.compute_wave_speed():.2f} m/s")
    print(f"  CFL number: {config.compute_cfl_number():.3f}")
    print(f"  Time step: {config.time_step:.6f} s")

    # Create simulation components
    device = torch.device(config.device)
    dtype = torch.float32 if config.dtype == 'float32' else torch.float64

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

    # Initialize wave packet
    center = 1.0  # m
    width = 0.4   # m (wider for better visualization)
    amplitude = 0.003  # m (larger amplitude)
    wavelength = 0.15  # m

    print(f"\nWave packet:")
    print(f"  Initial position: {center:.2f} m")
    print(f"  Width (σ): {width:.2f} m")
    print(f"  Wavelength (λ): {wavelength:.2f} m")
    print(f"  Amplitude: {amplitude:.4f} m")

    initialize_wave_packet(state, grid, center, width, amplitude, wavelength, config.compute_wave_speed())
    solver.initialize_accelerations(state)

    # Calculate simulation steps needed
    simulation_time = video_duration  # Simulate 10 seconds
    num_steps = int(simulation_time / config.time_step)
    capture_interval = max(1, num_steps // total_frames)  # Steps between frames

    print(f"\nSimulation:")
    print(f"  Total simulation time: {simulation_time} s")
    print(f"  Total steps: {num_steps}")
    print(f"  Capture interval: {capture_interval} steps/frame")

    # Storage for frames
    x_coords = np.arange(grid.grid_shape[0]) * grid.spacing
    frames_data = []
    times = []

    print(f"\nRunning simulation and capturing {total_frames} frames...")
    print("Progress: ", end='', flush=True)

    # Run simulation and capture frames
    frame_count = 0
    progress_marks = 20
    progress_interval = num_steps // progress_marks

    for step in range(num_steps):
        solver.step(state)

        if step % capture_interval == 0 and frame_count < total_frames:
            field = state.get_field_component(3).cpu().numpy()
            frames_data.append(field.copy())
            times.append(solver.time)
            frame_count += 1

        if step % progress_interval == 0:
            print("▓", end='', flush=True)

    # Ensure we have exactly the right number of frames
    if frame_count < total_frames:
        # Capture final frame if needed
        field = state.get_field_component(3).cpu().numpy()
        frames_data.append(field.copy())
        times.append(solver.time)

    print(" Done!")
    print(f"Captured {len(frames_data)} frames")

    # Create animation
    print("\nCreating video animation...")

    fig, ax = plt.subplots(figsize=(16, 6))

    # Find global min/max for consistent axis limits
    all_data = np.array(frames_data)
    y_min, y_max = all_data.min(), all_data.max()
    y_margin = 0.15 * (y_max - y_min) if y_max > y_min else 0.001

    line, = ax.plot([], [], 'b-', linewidth=2.5, label='Wave amplitude $\\xi^3$')
    time_text = ax.text(0.02, 0.96, '', transform=ax.transAxes,
                        fontsize=14, verticalalignment='top',
                        bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.9))

    ax.set_xlim(0, x_coords[-1])
    ax.set_ylim(y_min - y_margin, y_max + y_margin)
    ax.set_xlabel('Position [m]', fontsize=14, fontweight='bold')
    ax.set_ylabel('Amplitude $\\xi^3$ [m]', fontsize=14, fontweight='bold')
    ax.set_title('1D Wave Packet Propagation on Brane String', fontsize=16, fontweight='bold')
    ax.grid(True, alpha=0.3, linewidth=0.5)
    ax.legend(loc='upper right', fontsize=12)

    # Add zero line for reference
    ax.axhline(y=0, color='k', linestyle='--', linewidth=0.8, alpha=0.3)

    def init():
        """Initialize animation."""
        line.set_data([], [])
        time_text.set_text('')
        return line, time_text

    def animate(frame_idx):
        """Update animation frame."""
        line.set_data(x_coords, frames_data[frame_idx])
        time_text.set_text(
            f'Time: {times[frame_idx]:.2f} s\n'
            f'Frame: {frame_idx + 1}/{len(frames_data)}'
        )
        return line, time_text

    anim = FuncAnimation(
        fig,
        animate,
        init_func=init,
        frames=len(frames_data),
        interval=1000/fps,  # milliseconds between frames
        blit=True
    )

    # Save as MP4 video
    output_file = '1d_wave_animation_10s.mp4'
    print(f"\nSaving video as: {output_file}")
    print("This may take a minute...")

    writer = FFMpegWriter(
        fps=fps,
        bitrate=3000,  # Higher bitrate for better quality
        codec='libx264',
        extra_args=['-pix_fmt', 'yuv420p']  # For compatibility
    )

    anim.save(output_file, writer=writer, dpi=120)

    print(f"\n{'=' * 70}")
    print("Video creation complete!")
    print(f"{'=' * 70}")

    # Get file info
    file_size = os.path.getsize(output_file) / (1024 * 1024)  # Convert to MB
    print(f"\nOutput file: {output_file}")
    print(f"File size: {file_size:.2f} MB")
    print(f"Duration: {video_duration} seconds")
    print(f"Resolution: 1920x720 pixels")
    print(f"Frame rate: {fps} fps")
    print(f"Format: MP4 (H.264)")

    plt.close()


if __name__ == '__main__':
    main()
