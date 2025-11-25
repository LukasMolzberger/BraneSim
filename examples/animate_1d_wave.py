"""
Animated 1D Wave Simulation

Creates an animation showing wave propagation on a 1D brane string.
Exports as both MP4 video and GIF.
"""

import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import torch
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter, FFMpegWriter

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

    This creates a Gaussian envelope modulated by a sine wave,
    which propagates more cleanly than a pure Gaussian.

    Args:
        state: BraneState to initialize
        grid: BraneGrid with spatial coordinates
        center: Center position [m]
        width: Gaussian width σ [m]
        amplitude: Peak amplitude [m]
        wavelength: Carrier wavelength [m]
        wave_speed: Wave speed for velocity initialization [m/s]
    """
    # Get spatial coordinates
    x = (torch.arange(grid.grid_shape[0], device=state.device, dtype=state.dtype) *
         grid.spacing)

    # Gaussian envelope
    envelope = amplitude * torch.exp(-((x - center) ** 2) / (2 * width ** 2))

    # Wave number
    k = 2 * np.pi / wavelength
    omega = wave_speed * k

    # Position: ξ(x,0) = A(x) * cos(kx)
    state.positions[:, 3] = envelope * torch.cos(k * x)

    # Velocity: ∂ξ/∂t = -A(x) * ω * sin(kx)
    # For right-moving wave, we also need the envelope derivative term
    state.velocities[:, 3] = -envelope * omega * torch.sin(k * x)


def main():
    """Run animated 1D wave simulation."""
    print("=" * 60)
    print("1D Wave Animation: Wave Packet Propagation")
    print("=" * 60)

    # Configuration - longer domain to see propagation
    config = SimulationConfig.for_1d_test(
        nx=400,
        wave_speed=1.0,  # m/s
        spacing=0.01,    # m
        cfl_factor=0.4,
        device='cpu',
        dtype='float32'
    )

    print("\nConfiguration:")
    print(f"  Domain length: {config.grid_shape[0] * config.grid_spacing:.2f} m")
    print(f"  Wave speed: {config.compute_wave_speed():.2f} m/s")
    print(f"  CFL number: {config.compute_cfl_number():.3f}")

    # Create state, grid, physics, and solver
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
    width = 0.3   # m (wider packet)
    amplitude = 0.002  # m
    wavelength = 0.2  # m (carrier wavelength)

    print(f"\nWave packet:")
    print(f"  Center: {center:.2f} m")
    print(f"  Width: {width:.2f} m")
    print(f"  Wavelength: {wavelength:.2f} m")
    print(f"  Amplitude: {amplitude:.4f} m")

    initialize_wave_packet(state, grid, center, width, amplitude, wavelength, config.compute_wave_speed())

    # Initialize accelerations
    solver.initialize_accelerations(state)

    # Storage for animation
    num_steps = 2000
    capture_interval = 5  # Capture every 5 steps
    num_frames = num_steps // capture_interval

    x_coords = np.arange(grid.grid_shape[0]) * grid.spacing
    frames_data = []
    times = []

    print(f"\nRunning {num_steps} steps, capturing {num_frames} frames...")

    # Run simulation and capture frames
    for step in range(num_steps):
        solver.step(state)

        if step % capture_interval == 0:
            # Extract field and store
            field = state.get_field_component(3).cpu().numpy()
            frames_data.append(field.copy())
            times.append(solver.time)

            if step % 200 == 0:
                print(f"  Step {step:4d}/{num_steps}, t={solver.time:.3f}s")

    print("\nCreating animation...")

    # Create animation
    fig, ax = plt.subplots(figsize=(14, 6))

    # Find global min/max for consistent axis limits
    all_data = np.array(frames_data)
    y_min, y_max = all_data.min(), all_data.max()
    y_margin = 0.1 * (y_max - y_min)

    line, = ax.plot([], [], 'b-', linewidth=2, label='Wave amplitude')
    time_text = ax.text(0.02, 0.95, '', transform=ax.transAxes,
                        fontsize=12, verticalalignment='top',
                        bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    ax.set_xlim(0, x_coords[-1])
    ax.set_ylim(y_min - y_margin, y_max + y_margin)
    ax.set_xlabel('Position [m]', fontsize=12)
    ax.set_ylabel('Amplitude $\\xi^3$ [m]', fontsize=12)
    ax.set_title('1D Wave Packet Propagation on Brane String', fontsize=14)
    ax.grid(True, alpha=0.3)
    ax.legend(loc='upper right')

    def init():
        """Initialize animation."""
        line.set_data([], [])
        time_text.set_text('')
        return line, time_text

    def animate(frame_idx):
        """Update animation frame."""
        line.set_data(x_coords, frames_data[frame_idx])
        time_text.set_text(f'Time: {times[frame_idx]:.3f} s\nFrame: {frame_idx}/{len(frames_data)}')
        return line, time_text

    anim = FuncAnimation(
        fig,
        animate,
        init_func=init,
        frames=len(frames_data),
        interval=50,  # 50ms between frames = 20 fps
        blit=True
    )

    # Save as GIF
    print("\nSaving animation as GIF...")
    gif_path = '1d_wave_animation.gif'
    writer_gif = PillowWriter(fps=20)
    anim.save(gif_path, writer=writer_gif)
    print(f"  ✓ Saved: {gif_path}")

    # Try to save as MP4 (requires ffmpeg)
    try:
        print("\nSaving animation as MP4...")
        mp4_path = '1d_wave_animation.mp4'
        writer_mp4 = FFMpegWriter(fps=20, bitrate=1800)
        anim.save(mp4_path, writer=writer_mp4)
        print(f"  ✓ Saved: {mp4_path}")
    except Exception as e:
        print(f"  ⚠ Could not save MP4 (ffmpeg not installed?): {e}")
        print("    GIF version is available!")

    print("\n" + "=" * 60)
    print("Animation complete!")
    print("=" * 60)
    print(f"\nOutput files:")
    print(f"  - {gif_path}")
    if os.path.exists('1d_wave_animation.mp4'):
        print(f"  - 1d_wave_animation.mp4")

    plt.close()


if __name__ == '__main__':
    main()
