"""
3D Point Cloud Visualization for Brane Fields

Provides reusable 3D visualization using scatter plots with color/opacity mapping.
Supports camera animation and video generation.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, FFMpegWriter
from mpl_toolkits.mplot3d import Axes3D
from typing import Optional, Tuple, List, Callable
import torch


def subsample_3d_field(
    field_3d: np.ndarray,
    grid_shape: Tuple[int, int, int],
    h_phys: float,
    subsample_factor: int = 2,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Subsample a 3D field for visualization.

    Args:
        field_3d: Flattened field values (N,)
        grid_shape: (nx, ny, nz)
        h_phys: Physical lattice spacing [m]
        subsample_factor: Take every Nth point along each axis

    Returns:
        coords: (M, 3) physical coordinates [m]
        values: (M,) field values at those coordinates
    """
    nx, ny, nz = grid_shape

    # Reshape to 3D
    field = field_3d.reshape(nx, ny, nz)

    # Subsample
    field_sub = field[::subsample_factor, ::subsample_factor, ::subsample_factor]
    nx_sub, ny_sub, nz_sub = field_sub.shape

    # Generate coordinate grid
    x = np.arange(nx_sub) * subsample_factor * h_phys
    y = np.arange(ny_sub) * subsample_factor * h_phys
    z = np.arange(nz_sub) * subsample_factor * h_phys

    X, Y, Z = np.meshgrid(x, y, z, indexing='ij')

    # Flatten
    coords = np.stack([X.ravel(), Y.ravel(), Z.ravel()], axis=-1)
    values = field_sub.ravel()

    return coords, values


def field_to_colors(
    values: np.ndarray,
    cmap_name: str = 'RdBu_r',
    vmin: Optional[float] = None,
    vmax: Optional[float] = None,
    gamma: float = 1.5,
    alpha_scale: float = 0.8,
) -> np.ndarray:
    """
    Map field values to RGBA colors with opacity based on magnitude.

    Args:
        values: Field values
        cmap_name: Matplotlib colormap name
        vmin, vmax: Value range for color mapping (None = auto from data)
        gamma: Gamma correction for opacity (higher = more transparent low values)
        alpha_scale: Overall opacity multiplier

    Returns:
        colors: (N, 4) RGBA array
    """
    # Auto-determine symmetric range if not provided
    if vmin is None or vmax is None:
        abs_max = np.abs(values).max()
        if abs_max < 1e-12:
            abs_max = 1.0
        vmin = -abs_max if vmin is None else vmin
        vmax = abs_max if vmax is None else vmax

    # Normalize to [0, 1]
    denom = max(vmax - vmin, 1e-12)
    normalized = np.clip((values - vmin) / denom, 0, 1)

    # Map to colors
    cmap = plt.get_cmap(cmap_name)
    colors = cmap(normalized)

    # Set opacity based on absolute magnitude (with gamma correction)
    magnitude_norm = np.abs(values) / (np.abs(values).max() + 1e-12)
    colors[:, 3] = alpha_scale * (magnitude_norm ** gamma)

    return colors


def render_3d_field(
    coords: np.ndarray,
    values: np.ndarray,
    ax: Axes3D,
    cmap_name: str = 'RdBu_r',
    point_size: float = 4.0,
    gamma: float = 1.5,
    alpha_scale: float = 0.8,
    vmin: Optional[float] = None,
    vmax: Optional[float] = None,
) -> None:
    """
    Render 3D field as point cloud.

    Args:
        coords: (N, 3) coordinates
        values: (N,) field values
        ax: Matplotlib 3D axis
        cmap_name: Colormap name
        point_size: Scatter point size
        gamma: Opacity gamma correction
        alpha_scale: Overall opacity
        vmin, vmax: Color scale limits
    """
    colors = field_to_colors(values, cmap_name, vmin, vmax, gamma, alpha_scale)

    ax.scatter(
        coords[:, 0],
        coords[:, 1],
        coords[:, 2],
        s=point_size,
        c=colors,
        marker='o',
        depthshade=False,  # Disable depth shading to preserve our opacity
    )


def setup_3d_axes(
    ax: Axes3D,
    coords: np.ndarray,
    xlabel: str = 'x [nm]',
    ylabel: str = 'y [nm]',
    zlabel: str = 'z [nm]',
    title: str = '',
    equal_aspect: bool = True,
    margin_factor: float = 1.1,
) -> None:
    """
    Configure 3D axes with proper limits and labels.

    Args:
        ax: Matplotlib 3D axis
        coords: (N, 3) coordinates to determine extent
        xlabel, ylabel, zlabel: Axis labels
        title: Plot title
        equal_aspect: Whether to use equal aspect ratio
        margin_factor: Expand limits by this factor
    """
    # Compute extent
    x_min, y_min, z_min = coords.min(axis=0)
    x_max, y_max, z_max = coords.max(axis=0)

    x_range = x_max - x_min
    y_range = y_max - y_min
    z_range = z_max - z_min

    if equal_aspect:
        max_range = max(x_range, y_range, z_range)
        x_mid = 0.5 * (x_max + x_min)
        y_mid = 0.5 * (y_max + y_min)
        z_mid = 0.5 * (z_max + z_min)

        half_range = 0.5 * max_range * margin_factor
        ax.set_xlim(x_mid - half_range, x_mid + half_range)
        ax.set_ylim(y_mid - half_range, y_mid + half_range)
        ax.set_zlim(z_mid - half_range, z_mid + half_range)
        ax.set_box_aspect([1, 1, 1])
    else:
        ax.set_xlim(x_min * margin_factor, x_max * margin_factor)
        ax.set_ylim(y_min * margin_factor, y_max * margin_factor)
        ax.set_zlim(z_min * margin_factor, z_max * margin_factor)

    ax.set_xlabel(xlabel, fontsize=11)
    ax.set_ylabel(ylabel, fontsize=11)
    ax.set_zlabel(zlabel, fontsize=11)

    if title:
        ax.set_title(title, fontsize=13, fontweight='bold')


def camera_orbit(
    frame: int,
    num_frames: int,
    elev_start: float = 20,
    elev_end: float = 20,
    azim_start: float = -60,
    azim_revolutions: float = 1.0,
) -> Tuple[float, float]:
    """
    Calculate camera position for orbital animation.

    Args:
        frame: Current frame number
        num_frames: Total number of frames
        elev_start: Starting elevation angle [degrees]
        elev_end: Ending elevation angle [degrees]
        azim_start: Starting azimuth angle [degrees]
        azim_revolutions: Number of full rotations

    Returns:
        elev, azim: Camera angles [degrees]
    """
    t = frame / max(num_frames - 1, 1)
    elev = elev_start + t * (elev_end - elev_start)
    azim = azim_start + t * azim_revolutions * 360.0
    return elev, azim


def create_3d_animation(
    frames_data: List[Tuple[np.ndarray, np.ndarray]],
    times: List[float],
    output_path: str,
    cmap_name: str = 'RdBu_r',
    point_size: float = 4.0,
    gamma: float = 1.5,
    alpha_scale: float = 0.8,
    xlabel: str = 'x [nm]',
    ylabel: str = 'y [nm]',
    zlabel: str = 'z [nm]',
    title_template: str = '3D Field - t = {:.3f} fs',
    fps: int = 20,
    dpi: int = 100,
    camera_motion: Optional[Callable[[int, int], Tuple[float, float]]] = None,
    figsize: Tuple[int, int] = (10, 8),
) -> None:
    """
    Create animated 3D visualization.

    Args:
        frames_data: List of (coords, values) tuples for each frame
        times: Time values for each frame [s]
        output_path: Output video file path
        cmap_name: Colormap name
        point_size: Scatter point size
        gamma: Opacity gamma correction
        alpha_scale: Overall opacity
        xlabel, ylabel, zlabel: Axis labels
        title_template: Title template with {:.3f} placeholder for time
        fps: Frames per second
        dpi: Output resolution
        camera_motion: Optional function(frame_idx, num_frames) -> (elev, azim)
        figsize: Figure size in inches
    """
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111, projection='3d')

    # Use first frame to setup axes
    coords_0, values_0 = frames_data[0]
    setup_3d_axes(ax, coords_0, xlabel, ylabel, zlabel)

    # Determine global color scale from all frames
    all_values = np.concatenate([vals for _, vals in frames_data])
    abs_max = np.abs(all_values).max()
    if abs_max < 1e-12:
        abs_max = 1.0
    vmin, vmax = -abs_max, abs_max

    # Initial view
    if camera_motion is not None:
        elev, azim = camera_motion(0, len(frames_data))
        ax.view_init(elev=elev, azim=azim)
    else:
        ax.view_init(elev=20, azim=-60)

    # Create initial scatter
    colors_0 = field_to_colors(values_0, cmap_name, vmin, vmax, gamma, alpha_scale)
    scatter = ax.scatter(
        coords_0[:, 0],
        coords_0[:, 1],
        coords_0[:, 2],
        s=point_size,
        c=colors_0,
        marker='o',
        depthshade=False,
    )

    # Title text - initialize with first frame time
    time_0_fs = times[0] * 1e15
    time_text = ax.text2D(0.5, 0.95, title_template.format(time_0_fs),
                         transform=ax.transAxes,
                         ha='center', fontsize=13, fontweight='bold')

    def update(frame_idx):
        """Update function for animation."""
        coords, values = frames_data[frame_idx]
        time_fs = times[frame_idx] * 1e15

        # Update scatter data
        scatter._offsets3d = (coords[:, 0], coords[:, 1], coords[:, 2])
        colors = field_to_colors(values, cmap_name, vmin, vmax, gamma, alpha_scale)
        scatter.set_color(colors)

        # Update camera if motion specified
        if camera_motion is not None:
            elev, azim = camera_motion(frame_idx, len(frames_data))
            ax.view_init(elev=elev, azim=azim)

        # Update title
        time_text.set_text(title_template.format(time_fs))

        return scatter, time_text

    anim = FuncAnimation(
        fig,
        update,
        frames=len(frames_data),
        interval=1000/fps,
        blit=False,  # 3D doesn't support blitting
        repeat=True
    )

    # Save animation
    writer = FFMpegWriter(fps=fps, bitrate=2000)
    anim.save(output_path, writer=writer, dpi=dpi)
    plt.close(fig)