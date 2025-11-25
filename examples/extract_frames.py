"""
Extract key frames from the simulation for preview.
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


def initialize_wave_packet(state, grid, center, width, amplitude, wavelength, wave_speed):
    """Initialize wave packet."""
    x = (torch.arange(grid.grid_shape[0], device=state.device, dtype=state.dtype) *
         grid.spacing)
    envelope = amplitude * torch.exp(-((x - center) ** 2) / (2 * width ** 2))
    k = 2 * np.pi / wavelength
    omega = wave_speed * k
    state.positions[:, 3] = envelope * torch.cos(k * x)
    state.velocities[:, 3] = -envelope * omega * torch.sin(k * x)


def main():
    """Extract key frames from simulation."""
    # Configuration
    config = SimulationConfig.for_1d_test(
        nx=400, wave_speed=1.0, spacing=0.01, cfl_factor=0.4,
        device='cpu', dtype='float32'
    )

    device = torch.device(config.device)
    dtype = torch.float32

    state = BraneState(config.grid_shape, config.dimension, device, dtype)
    state.initialize_flat_configuration(config.grid_spacing)
    grid = BraneGrid(config.grid_shape, config.dimension, config.grid_spacing, device)

    physics = SpringForceComputer(config.spring_constant, config.rest_length, config.critical_strain)
    solver = VelocityVerletSolver(config.time_step, config.mass_density, physics, grid)

    # Initialize
    initialize_wave_packet(state, grid, 1.0, 0.3, 0.002, 0.2, config.compute_wave_speed())
    solver.initialize_accelerations(state)

    # Extract frames at different times
    frame_times = [0, 500, 1000, 1500, 2000]  # Steps
    frame_labels = ['t=0.0s', 't=2.0s', 't=4.0s', 't=6.0s', 't=8.0s']

    x_coords = np.arange(grid.grid_shape[0]) * grid.spacing

    fig, axes = plt.subplots(5, 1, figsize=(14, 10))
    fig.suptitle('1D Wave Packet Propagation - Key Frames', fontsize=16, fontweight='bold')

    for idx, (target_step, label) in enumerate(zip(frame_times, frame_labels)):
        # Run simulation to target step
        while solver.step_count < target_step:
            solver.step(state)

        # Extract and plot
        field = state.get_field_component(3).cpu().numpy()
        axes[idx].plot(x_coords, field, 'b-', linewidth=2)
        axes[idx].set_ylabel('Amplitude [m]')
        axes[idx].set_xlim(0, x_coords[-1])
        axes[idx].set_ylim(-0.0025, 0.0025)
        axes[idx].grid(True, alpha=0.3)
        axes[idx].text(0.02, 0.95, label, transform=axes[idx].transAxes,
                      fontsize=11, verticalalignment='top',
                      bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

        if idx == 4:
            axes[idx].set_xlabel('Position [m]')

    plt.tight_layout()
    plt.savefig('1d_wave_keyframes.png', dpi=150, bbox_inches='tight')
    print("✓ Saved key frames to: 1d_wave_keyframes.png")
    plt.close()


if __name__ == '__main__':
    main()
