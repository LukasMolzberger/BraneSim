"""
Volumetric 3D rendering for brane fields.

Uses matplotlib voxels with RGBA mapping and camera orbit for animations.
"""

from typing import Callable, List, Optional, Tuple

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, FFMpegWriter

from branesim.visualization.brane_3d_viz import camera_orbit


def downsample_volume(
    field_3d: np.ndarray,
    grid_shape: Tuple[int, int, int],
    subsample_factor: int = 2,
) -> np.ndarray:
    """
    Downsample a flattened 3D field into a regular 3D volume.

    Args:
        field_3d: Flattened field values (N,)
        grid_shape: (nx, ny, nz)
        subsample_factor: Take every Nth point along each axis

    Returns:
        volume: (nx_sub, ny_sub, nz_sub) array
    """
    nx, ny, nz = grid_shape
    volume = field_3d.reshape(nx, ny, nz)
    return volume[::subsample_factor, ::subsample_factor, ::subsample_factor]


def _volume_to_rgba(
    values: np.ndarray,
    cmap_name: str,
    vmin: float,
    vmax: float,
    gamma: float,
    alpha_scale: float,
    density_threshold: float,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Map volume values to RGBA colors and a mask.
    """
    denom = max(vmax - vmin, 1e-12)
    normalized = np.clip((values - vmin) / denom, 0.0, 1.0)
    cmap = plt.get_cmap(cmap_name)
    colors = cmap(normalized)

    abs_max = max(np.abs(vmin), np.abs(vmax), 1e-12)
    magnitude = np.abs(values) / abs_max
    colors[..., 3] = alpha_scale * (magnitude ** gamma)

    mask = magnitude >= density_threshold
    return mask, colors


def _setup_volume_axes(
    ax,
    grid_shape: Tuple[int, int, int],
    spacing: float,
    xlabel: str,
    ylabel: str,
    zlabel: str,
) -> None:
    """
    Configure 3D axes for volume rendering.
    """
    nx, ny, nz = grid_shape
    ax.set_xlim(0.0, nx * spacing)
    ax.set_ylim(0.0, ny * spacing)
    ax.set_zlim(0.0, nz * spacing)
    ax.set_box_aspect([nx, ny, nz])
    ax.set_xlabel(xlabel, fontsize=11)
    ax.set_ylabel(ylabel, fontsize=11)
    ax.set_zlabel(zlabel, fontsize=11)
    ax.grid(False)


def create_3d_volume_animation(
    frames: List[np.ndarray],
    times: List[float],
    spacing: float,
    output_path: str,
    cmap_name: str = "RdBu_r",
    fps: int = 20,
    dpi: int = 120,
    camera_motion: Optional[Callable[[int, int], Tuple[float, float]]] = None,
    density_threshold: float = 0.08,
    alpha_scale: float = 0.9,
    gamma: float = 1.1,
    figsize: Tuple[int, int] = (10, 8),
    xlabel: str = "x [nm]",
    ylabel: str = "y [nm]",
    zlabel: str = "z [nm]",
    title_template: str = "3D Electron Volume (t = {:.2f} fs)",
) -> None:
    """
    Create an animated volumetric rendering using voxels.

    Args:
        frames: List of 3D volume arrays
        times: Physical times for each frame [s]
        spacing: Voxel edge length in the same units as labels
        output_path: Output video file path
        cmap_name: Matplotlib colormap name
        fps: Frames per second
        dpi: Output resolution
        camera_motion: Optional function(frame_idx, num_frames) -> (elev, azim)
        density_threshold: Normalized magnitude threshold for rendering voxels
        alpha_scale: Overall opacity multiplier
        gamma: Opacity gamma correction
        figsize: Figure size in inches
        xlabel, ylabel, zlabel: Axis labels
        title_template: Title template with {:.2f} placeholder for time (fs)
    """
    if not frames:
        return

    grid_shape = frames[0].shape
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111, projection="3d")
    _setup_volume_axes(ax, grid_shape, spacing, xlabel, ylabel, zlabel)

    # Precompute voxel corner coordinates
    x_edges = np.arange(grid_shape[0] + 1) * spacing
    y_edges = np.arange(grid_shape[1] + 1) * spacing
    z_edges = np.arange(grid_shape[2] + 1) * spacing
    X, Y, Z = np.meshgrid(x_edges, y_edges, z_edges, indexing="ij")

    # Global color scale
    all_values = np.concatenate([frame.ravel() for frame in frames])
    abs_max = np.abs(all_values).max()
    if abs_max < 1e-12:
        abs_max = 1.0
    vmin, vmax = -abs_max, abs_max

    if camera_motion is None:
        camera_motion = camera_orbit

    # Title text
    time_0_fs = times[0] * 1e15
    time_text = ax.text2D(
        0.5,
        0.95,
        title_template.format(time_0_fs),
        transform=ax.transAxes,
        ha="center",
        fontsize=13,
        fontweight="bold",
    )

    def update(frame_idx: int):
        values = frames[frame_idx]
        time_fs = times[frame_idx] * 1e15

        ax.collections.clear()
        ax.patches.clear()

        mask, colors = _volume_to_rgba(
            values,
            cmap_name=cmap_name,
            vmin=vmin,
            vmax=vmax,
            gamma=gamma,
            alpha_scale=alpha_scale,
            density_threshold=density_threshold,
        )
        ax.voxels(X, Y, Z, mask, facecolors=colors, edgecolor="none")

        elev, azim = camera_motion(frame_idx, len(frames))
        ax.view_init(elev=elev, azim=azim)
        time_text.set_text(title_template.format(time_fs))
        return ax.collections + [time_text]

    anim = FuncAnimation(
        fig,
        update,
        frames=len(frames),
        interval=1000 / fps,
        blit=False,
        repeat=True,
    )

    writer = FFMpegWriter(fps=fps, bitrate=2400)
    anim.save(output_path, writer=writer, dpi=dpi)
    plt.close(fig)
