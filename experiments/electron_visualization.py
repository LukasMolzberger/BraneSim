"""
Enhanced Electron Visualization

Creates:
1. Initial state visualization (3 orthogonal slices)
2. 6 videos: 3 amplitude evolution + 3 lateral distortion evolution
"""

import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import torch
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, FFMpegWriter

def visualize_initial_state(state, params, config):
    """
    Visualize initial electron state in 3 orthogonal slices.

    Args:
        state: BraneState after initialization
        params: ElectronInitParams
        config: Configuration dictionary
    """
    print(f"\n{'='*60}")
    print(f"Creating Initial State Visualization")
    print(f"{'='*60}")

    nx, ny, nz = config['grid_shape']
    h = config['h']

    # Get amplitude field
    X4 = state.positions[:, 3].cpu().numpy()
    X4_grid = X4.reshape((nx, ny, nz))

    # Get lateral positions for distortion visualization
    X_lat = state.positions[:, :3].cpu().numpy()

    # Create flat grid for comparison
    x = np.arange(nx) * h
    y = np.arange(ny) * h
    z = np.arange(nz) * h
    X_flat, Y_flat, Z_flat = np.meshgrid(x, y, z, indexing='ij')

    # Compute lateral distortion
    X_lat_grid = X_lat.reshape((nx, ny, nz, 3))
    distortion_x = X_lat_grid[:, :, :, 0] - X_flat
    distortion_y = X_lat_grid[:, :, :, 1] - Y_flat
    distortion_z = X_lat_grid[:, :, :, 2] - Z_flat
    distortion_mag = np.sqrt(distortion_x**2 + distortion_y**2 + distortion_z**2)

    # Center indices
    cx, cy, cz = nx//2, ny//2, nz//2

    # Create figure with 3x2 layout
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    fig.suptitle('Initial Electron State - 3 Orthogonal Slices', fontsize=16, fontweight='bold')

    # XY slice (through center in z)
    ax = axes[0, 0]
    im = ax.imshow(X4_grid[:, :, cz].T, origin='lower', cmap='RdBu_r',
                   extent=[0, nx*h*1e12, 0, ny*h*1e12])
    ax.set_xlabel('X [pm]')
    ax.set_ylabel('Y [pm]')
    ax.set_title('XY Slice - Amplitude X⁴')
    plt.colorbar(im, ax=ax, label='X⁴ [m]')

    ax = axes[1, 0]
    im = ax.imshow(distortion_mag[:, :, cz].T, origin='lower', cmap='hot',
                   extent=[0, nx*h*1e12, 0, ny*h*1e12])
    ax.set_xlabel('X [pm]')
    ax.set_ylabel('Y [pm]')
    ax.set_title('XY Slice - Lateral Distortion')
    plt.colorbar(im, ax=ax, label='|Δr| [m]')

    # XZ slice (through center in y)
    ax = axes[0, 1]
    im = ax.imshow(X4_grid[:, cy, :].T, origin='lower', cmap='RdBu_r',
                   extent=[0, nx*h*1e12, 0, nz*h*1e12])
    ax.set_xlabel('X [pm]')
    ax.set_ylabel('Z [pm]')
    ax.set_title('XZ Slice - Amplitude X⁴')
    plt.colorbar(im, ax=ax, label='X⁴ [m]')

    ax = axes[1, 1]
    im = ax.imshow(distortion_mag[:, cy, :].T, origin='lower', cmap='hot',
                   extent=[0, nx*h*1e12, 0, nz*h*1e12])
    ax.set_xlabel('X [pm]')
    ax.set_ylabel('Z [pm]')
    ax.set_title('XZ Slice - Lateral Distortion')
    plt.colorbar(im, ax=ax, label='|Δr| [m]')

    # YZ slice (through center in x)
    ax = axes[0, 2]
    im = ax.imshow(X4_grid[cx, :, :].T, origin='lower', cmap='RdBu_r',
                   extent=[0, ny*h*1e12, 0, nz*h*1e12])
    ax.set_xlabel('Y [pm]')
    ax.set_ylabel('Z [pm]')
    ax.set_title('YZ Slice - Amplitude X⁴')
    plt.colorbar(im, ax=ax, label='X⁴ [m]')

    ax = axes[1, 2]
    im = ax.imshow(distortion_mag[cx, :, :].T, origin='lower', cmap='hot',
                   extent=[0, ny*h*1e12, 0, nz*h*1e12])
    ax.set_xlabel('Y [pm]')
    ax.set_ylabel('Z [pm]')
    ax.set_title('YZ Slice - Lateral Distortion')
    plt.colorbar(im, ax=ax, label='|Δr| [m]')

    plt.tight_layout()
    plt.savefig('electron_initial_state.png', dpi=150, bbox_inches='tight')
    print(f"  ✓ Saved: electron_initial_state.png")
    plt.close()


def collect_animation_frames(states, config):
    """
    Collect animation frames from simulation states.

    Args:
        states: List of BraneState snapshots
        config: Configuration dictionary

    Returns:
        Dictionary with frames for each slice and type
    """
    print(f"\n{'='*60}")
    print(f"Collecting Animation Frames")
    print(f"{'='*60}")

    nx, ny, nz = config['grid_shape']
    h = config['h']

    # Center indices
    cx, cy, cz = nx//2, ny//2, nz//2

    # Create flat grid for distortion calculation
    x = np.arange(nx) * h
    y = np.arange(ny) * h
    z = np.arange(nz) * h
    X_flat, Y_flat, Z_flat = np.meshgrid(x, y, z, indexing='ij')

    frames = {
        'amplitude_xy': [],
        'amplitude_xz': [],
        'amplitude_yz': [],
        'distortion_xy': [],
        'distortion_xz': [],
        'distortion_yz': [],
        'times': []
    }

    for idx, state in enumerate(states):
        # Amplitude field
        X4 = state.positions[:, 3].cpu().numpy()
        X4_grid = X4.reshape((nx, ny, nz))

        frames['amplitude_xy'].append(X4_grid[:, :, cz])
        frames['amplitude_xz'].append(X4_grid[:, cy, :])
        frames['amplitude_yz'].append(X4_grid[cx, :, :])

        # Lateral distortion
        X_lat = state.positions[:, :3].cpu().numpy()
        X_lat_grid = X_lat.reshape((nx, ny, nz, 3))
        distortion_x = X_lat_grid[:, :, :, 0] - X_flat
        distortion_y = X_lat_grid[:, :, :, 1] - Y_flat
        distortion_z = X_lat_grid[:, :, :, 2] - Z_flat
        distortion_mag = np.sqrt(distortion_x**2 + distortion_y**2 + distortion_z**2)

        frames['distortion_xy'].append(distortion_mag[:, :, cz])
        frames['distortion_xz'].append(distortion_mag[:, cy, :])
        frames['distortion_yz'].append(distortion_mag[cx, :, :])

        # Time
        t = idx * config['snapshot_interval'] * config['dt']
        frames['times'].append(t)

    print(f"  Collected {len(frames['times'])} frames")
    return frames


def create_animation(frames, slice_name, field_type, config, filename):
    """
    Create and save animation for a specific slice and field type.

    Args:
        frames: List of 2D arrays (frames)
        slice_name: 'xy', 'xz', or 'yz'
        field_type: 'amplitude' or 'distortion'
        config: Configuration dictionary
        filename: Output filename
    """
    nx, ny, nz = config['grid_shape']
    h = config['h']
    T_compton = 2 * np.pi / (config['constants'].m_e * config['constants'].c ** 2 / config['constants'].hbar)

    # Determine extents and labels based on slice
    if slice_name == 'xy':
        extent = [0, nx*h*1e12, 0, ny*h*1e12]
        xlabel, ylabel = 'X [pm]', 'Y [pm]'
    elif slice_name == 'xz':
        extent = [0, nx*h*1e12, 0, nz*h*1e12]
        xlabel, ylabel = 'X [pm]', 'Z [pm]'
    else:  # yz
        extent = [0, ny*h*1e12, 0, nz*h*1e12]
        xlabel, ylabel = 'Y [pm]', 'Z [pm]'

    # Determine colormap and scaling
    if field_type == 'amplitude':
        cmap = 'RdBu_r'
        max_val = max(np.abs(f).max() for f in frames)
        vmin, vmax = -max_val, max_val
        label = 'X⁴ [m]'
        title = f'{slice_name.upper()} Slice - Amplitude Evolution'
    else:
        cmap = 'hot'
        vmin = 0
        vmax = max(f.max() for f in frames)
        label = '|Δr| [m]'
        title = f'{slice_name.upper()} Slice - Lateral Distortion'

    # Create figure
    fig, ax = plt.subplots(figsize=(10, 8))

    # Initial frame
    im = ax.imshow(frames[0].T, origin='lower', extent=extent,
                   cmap=cmap, vmin=vmin, vmax=vmax, animated=True)
    ax.set_xlabel(xlabel, fontsize=12)
    ax.set_ylabel(ylabel, fontsize=12)
    time_text = ax.text(0.02, 0.95, '', transform=ax.transAxes,
                       fontsize=12, verticalalignment='top',
                       bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    plt.colorbar(im, ax=ax, label=label, fraction=0.046, pad=0.04)
    ax.set_title(title, fontsize=14, fontweight='bold')

    def animate(frame_idx):
        """Update function for animation."""
        im.set_array(frames[frame_idx].T)
        t = frames_dict['times'][frame_idx]
        time_text.set_text(f't = {t/T_compton:.3f} T_C')
        return [im, time_text]

    # Create animation
    anim = FuncAnimation(fig, animate, frames=len(frames),
                        interval=100, blit=True, repeat=True)

    # Save
    writer = FFMpegWriter(fps=10, bitrate=2000)
    anim.save(filename, writer=writer, dpi=100)
    print(f"  ✓ Saved: {filename}")

    plt.close(fig)


# Store frames globally for animate function
frames_dict = {}

def create_all_animations(states, config):
    """
    Create all 6 animations (3 amplitude + 3 distortion).

    Args:
        states: List of BraneState snapshots
        config: Configuration dictionary
    """
    global frames_dict

    print(f"\n{'='*60}")
    print(f"Creating Animations")
    print(f"{'='*60}")

    # Collect frames
    frames_dict = collect_animation_frames(states, config)

    # Create amplitude animations
    create_animation(frames_dict['amplitude_xy'], 'xy', 'amplitude', config, 'electron_amplitude_xy.mp4')
    create_animation(frames_dict['amplitude_xz'], 'xz', 'amplitude', config, 'electron_amplitude_xz.mp4')
    create_animation(frames_dict['amplitude_yz'], 'yz', 'amplitude', config, 'electron_amplitude_yz.mp4')

    # Create distortion animations
    create_animation(frames_dict['distortion_xy'], 'xy', 'distortion', config, 'electron_distortion_xy.mp4')
    create_animation(frames_dict['distortion_xz'], 'xz', 'distortion', config, 'electron_distortion_xz.mp4')
    create_animation(frames_dict['distortion_yz'], 'yz', 'distortion', config, 'electron_distortion_yz.mp4')

    print(f"\n✓ All animations created successfully!")