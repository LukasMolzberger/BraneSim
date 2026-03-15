"""Volumetric rendering utilities owned by visualization component."""

from __future__ import annotations

from typing import Callable, List, Optional, Tuple

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, FFMpegWriter


def downsample_volume(field_3d: np.ndarray, grid_shape: Tuple[int, int, int], subsample_factor: int = 2) -> np.ndarray:
    nx, ny, nz = grid_shape
    volume = field_3d.reshape(nx, ny, nz)
    return volume[::subsample_factor, ::subsample_factor, ::subsample_factor]


def camera_orbit(frame_idx: int, num_frames: int) -> Tuple[float, float]:
    t = frame_idx / max(num_frames - 1, 1)
    elev = 22.0 + 8.0 * np.sin(2 * np.pi * t)
    azim = -60.0 + 360.0 * 0.5 * t
    return elev, azim


def _volume_to_rgba(
    values: np.ndarray,
    cmap_name: str,
    vmin: float,
    vmax: float,
    gamma: float,
    color_gamma: float,
    alpha_scale: float,
    density_threshold: float,
) -> tuple[np.ndarray, np.ndarray]:
    abs_max = max(np.abs(vmin), np.abs(vmax), 1e-12)
    signed = np.clip(values / abs_max, -1.0, 1.0)
    if color_gamma != 1.0:
        signed = np.sign(signed) * (np.abs(signed) ** color_gamma)
    normalized = 0.5 + 0.5 * signed

    cmap = plt.get_cmap(cmap_name)
    colors = cmap(normalized)

    magnitude = np.abs(values) / abs_max
    colors[..., 3] = alpha_scale * (magnitude ** gamma)
    mask = magnitude >= density_threshold
    return mask, colors


def create_3d_volume_animation(
    frames: List[np.ndarray],
    times: List[float],
    spacing: float,
    output_path: str,
    cmap_name: str = "RdBu_r",
    fps: int = 20,
    dpi: int = 120,
    camera_motion: Optional[Callable[[int, int], Tuple[float, float]]] = None,
    density_threshold: float = 0.005,
    alpha_scale: float = 1.0,
    gamma: float = 0.8,
    color_gamma: float = 1.0,
) -> None:
    if not frames:
        return

    grid_shape = frames[0].shape
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection="3d")

    x_edges = np.arange(grid_shape[0] + 1) * spacing
    y_edges = np.arange(grid_shape[1] + 1) * spacing
    z_edges = np.arange(grid_shape[2] + 1) * spacing
    X, Y, Z = np.meshgrid(x_edges, y_edges, z_edges, indexing="ij")

    all_values = np.concatenate([frame.ravel() for frame in frames])
    abs_max = np.abs(all_values).max()
    if abs_max < 1e-12:
        abs_max = 1.0
    vmin, vmax = -abs_max, abs_max

    cam = camera_motion or camera_orbit

    def setup_axes() -> None:
        nx, ny, nz = grid_shape
        ax.set_xlim(0.0, nx * spacing)
        ax.set_ylim(0.0, ny * spacing)
        ax.set_zlim(0.0, nz * spacing)
        ax.set_box_aspect([nx, ny, nz])
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_zticks([])
        ax.set_axis_off()

    def update(frame_idx: int):
        ax.cla()
        setup_axes()
        values = frames[frame_idx]

        mask, colors = _volume_to_rgba(
            values,
            cmap_name=cmap_name,
            vmin=vmin,
            vmax=vmax,
            gamma=gamma,
            color_gamma=color_gamma,
            alpha_scale=alpha_scale,
            density_threshold=density_threshold,
        )
        ax.voxels(X, Y, Z, mask, facecolors=colors, edgecolor="none")

        elev, azim = cam(frame_idx, len(frames))
        ax.view_init(elev=elev, azim=azim)
        return ax.collections

    anim = FuncAnimation(fig, update, frames=len(frames), interval=1000 / fps, blit=False, repeat=True)
    writer = FFMpegWriter(fps=fps, bitrate=2400)
    anim.save(output_path, writer=writer, dpi=dpi)
    plt.close(fig)
