"""branesim.berry.viz

Dimension-specific visualization for Berry diagnostics.

We visualize scalar fields defined on the brane (Berry phase, Berry connection)
using an HSV→RGB hue map for phase-like values and **alpha masking** to show
only regions where the polarization ray is defined.

Outputs
-------
- 1D brane: scatter/line video (x vs value) with per-point RGBA colors.
- 2D brane: RGBA imshow video.
- 3D brane: XY/XZ/YZ slice RGBA videos.

The analytical part (berry.analysis) is independent of intrinsic dimension; the
caller chooses how to reshape N into grid_shape for visualization.
"""

from __future__ import annotations

import os
from typing import List, Optional, Sequence, Tuple

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, FFMpegWriter
from matplotlib.colors import hsv_to_rgb

from .analysis import BerryTimeSeries


# -----------------------------------------------------------------------------
# Color mapping helpers
# -----------------------------------------------------------------------------


def phase_to_rgb(phase: np.ndarray) -> np.ndarray:
    """Map phase in radians to RGB via hue."""
    phase = np.asarray(phase, dtype=np.float64)
    hue = (phase + np.pi) / (2.0 * np.pi)
    hue = np.mod(hue, 1.0)
    sat = np.ones_like(hue)
    val = np.ones_like(hue)
    hsv = np.stack([hue, sat, val], axis=-1)
    return hsv_to_rgb(hsv)


def scalar_to_rgba(phase: np.ndarray, alpha: np.ndarray) -> np.ndarray:
    """Convert a scalar 'phase-like' field to RGBA image."""
    rgb = phase_to_rgb(phase)
    a = np.clip(np.asarray(alpha, dtype=np.float64), 0.0, 1.0)
    return np.concatenate([rgb, a[..., None]], axis=-1)


# -----------------------------------------------------------------------------
# 1D brane in 2D embedding: scatter video(s)
# -----------------------------------------------------------------------------


def create_berry_video_1d(
    values_frames: np.ndarray,
    alpha_frames: np.ndarray,
    times_s: Sequence[float],
    x_coords: np.ndarray,
    output_path: str,
    title_template: str,
    ylabel: str,
    fps: int = 20,
    dpi: int = 140,
    unit_label_x: str = "nm",
    ylim: Tuple[float, float] = (-np.pi, np.pi),
) -> None:
    """Create a 1D scatter video with per-point RGBA colors.

    values_frames: (T,N)
    alpha_frames:  (T,N)
    x_coords:      (N,) in display units
    """
    os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)

    values = np.asarray(values_frames, dtype=np.float64)
    alpha = np.asarray(alpha_frames, dtype=np.float64)
    times = np.asarray(times_s, dtype=np.float64)
    x = np.asarray(x_coords, dtype=np.float64)

    if values.ndim != 2:
        raise ValueError(f"values_frames must be (T,N). Got {values.shape}.")
    if alpha.shape != values.shape:
        raise ValueError(f"alpha_frames must match values_frames. Got {alpha.shape} vs {values.shape}.")
    if x.shape != (values.shape[1],):
        raise ValueError(f"x_coords must be (N,). Got {x.shape}.")

    T, N = values.shape

    fig, ax = plt.subplots(figsize=(10, 5))
    ax.set_xlabel(f"x [{unit_label_x}]")
    ax.set_ylabel(ylabel)
    ax.set_ylim(*ylim)
    ax.grid(True, alpha=0.25)

    rgba0 = np.concatenate([phase_to_rgb(values[0]), alpha[0, :, None]], axis=-1)
    sc = ax.scatter(x, values[0], s=10, c=rgba0, marker='s', linewidths=0)

    time_text = ax.text(
        0.02,
        0.95,
        "",
        transform=ax.transAxes,
        fontsize=12,
        va="top",
        bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.8),
    )

    def update(k: int):
        rgba = np.concatenate([phase_to_rgb(values[k]), alpha[k, :, None]], axis=-1)
        sc.set_offsets(np.c_[x, values[k]])
        sc.set_color(rgba)
        t_as = times[k] * 1e18
        time_text.set_text(f"t = {t_as:.3f} as")
        ax.set_title(title_template.format(t=t_as))
        return (sc, time_text)

    anim = FuncAnimation(fig, update, frames=T, interval=1000 / fps, blit=False)
    writer = FFMpegWriter(fps=fps, bitrate=2500)
    anim.save(output_path, writer=writer, dpi=dpi)
    plt.close(fig)


def create_berry_videos_1d(
    series: BerryTimeSeries,
    x_coords: np.ndarray,
    output_dir: str,
    filename_prefix: str = "berry_1d",
    fps: int = 20,
    dpi: int = 140,
    unit_label_x: str = "nm",
) -> List[str]:
    os.makedirs(output_dir, exist_ok=True)

    out_phase = os.path.join(output_dir, f"{filename_prefix}_phase.mp4")
    out_conn = os.path.join(output_dir, f"{filename_prefix}_connection.mp4")

    create_berry_video_1d(
        values_frames=series.phase,
        alpha_frames=series.alpha,
        times_s=series.times_s,
        x_coords=x_coords,
        output_path=out_phase,
        title_template="Berry phase γ (t = {t:.3f} as)",
        ylabel="γ [rad]",
        fps=fps,
        dpi=dpi,
        unit_label_x=unit_label_x,
        ylim=(-np.pi, np.pi),
    )

    create_berry_video_1d(
        values_frames=series.connection,
        alpha_frames=series.alpha,
        times_s=series.times_s,
        x_coords=x_coords,
        output_path=out_conn,
        title_template="Berry connection A_t (t = {t:.3f} as)",
        ylabel="A_t [rad/s]",
        fps=fps,
        dpi=dpi,
        unit_label_x=unit_label_x,
        ylim=(
            float(np.nanmin(series.connection)),
            float(np.nanmax(series.connection)),
        ) if np.isfinite(series.connection).all() else (-1.0, 1.0),
    )

    return [out_phase, out_conn]


# -----------------------------------------------------------------------------
# 2D brane in 3D embedding: RGBA image video(s)
# -----------------------------------------------------------------------------


def create_berry_video_2d(
    values_frames: np.ndarray,
    alpha_frames: np.ndarray,
    times_s: Sequence[float],
    grid_shape: Tuple[int, int],
    spacing: float,
    output_path: str,
    title: str,
    fps: int = 20,
    dpi: int = 120,
    display_scale: float = 1e9,
    unit_label: str = "nm",
    aspect: str = "equal",
) -> None:
    """Create a 2D RGBA imshow video.

    values_frames: (T,N)
    alpha_frames:  (T,N)
    grid_shape:    (nx,ny)
    spacing:       physical spacing (m)
    """
    os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)

    values = np.asarray(values_frames, dtype=np.float64)
    alpha = np.asarray(alpha_frames, dtype=np.float64)
    times = np.asarray(times_s, dtype=np.float64)

    nx, ny = grid_shape
    T, N = values.shape
    if N != nx * ny:
        raise ValueError(f"values_frames has N={N}, but grid_shape implies {nx*ny}.")
    if alpha.shape != values.shape:
        raise ValueError("alpha_frames must match values_frames.")

    extent = [0, (nx - 1) * spacing * display_scale, 0, (ny - 1) * spacing * display_scale]

    fig, ax = plt.subplots(figsize=(8, 7))
    ax.set_aspect(aspect)
    ax.set_xlabel(f"x [{unit_label}]")
    ax.set_ylabel(f"y [{unit_label}]")

    v0 = values[0].reshape(nx, ny)
    a0 = alpha[0].reshape(nx, ny)
    rgba0 = scalar_to_rgba(v0, a0)

    im = ax.imshow(
        rgba0.swapaxes(0, 1),
        origin="lower",
        extent=extent,
        interpolation="nearest",
        animated=True,
    )

    time_text = ax.text(
        0.02,
        0.95,
        "",
        transform=ax.transAxes,
        fontsize=12,
        va="top",
        bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.8),
    )
    ax.set_title(title)

    def update(k: int):
        v = values[k].reshape(nx, ny)
        a = alpha[k].reshape(nx, ny)
        rgba = scalar_to_rgba(v, a)
        im.set_data(rgba.swapaxes(0, 1))
        t_as = times[k] * 1e18
        time_text.set_text(f"t = {t_as:.3f} as")
        return (im, time_text)

    anim = FuncAnimation(fig, update, frames=T, interval=1000 / fps, blit=True)
    writer = FFMpegWriter(fps=fps, bitrate=2500)
    anim.save(output_path, writer=writer, dpi=dpi)
    plt.close(fig)


def create_berry_videos_2d(
    series: BerryTimeSeries,
    grid_shape: Tuple[int, int],
    spacing: float,
    output_dir: str,
    filename_prefix: str = "berry_2d",
    fps: int = 20,
    dpi: int = 120,
    display_scale: float = 1e9,
    unit_label: str = "nm",
) -> List[str]:
    os.makedirs(output_dir, exist_ok=True)

    out_phase = os.path.join(output_dir, f"{filename_prefix}_phase.mp4")
    out_conn = os.path.join(output_dir, f"{filename_prefix}_connection.mp4")

    create_berry_video_2d(
        values_frames=series.phase,
        alpha_frames=series.alpha,
        times_s=series.times_s,
        grid_shape=grid_shape,
        spacing=spacing,
        output_path=out_phase,
        title="Berry phase γ (hue) with amplitude mask (alpha)",
        fps=fps,
        dpi=dpi,
        display_scale=display_scale,
        unit_label=unit_label,
    )

    create_berry_video_2d(
        values_frames=series.connection,
        alpha_frames=series.alpha,
        times_s=series.times_s,
        grid_shape=grid_shape,
        spacing=spacing,
        output_path=out_conn,
        title="Berry connection A_t (hue) with amplitude mask (alpha)",
        fps=fps,
        dpi=dpi,
        display_scale=display_scale,
        unit_label=unit_label,
    )

    return [out_phase, out_conn]


# -----------------------------------------------------------------------------
# 3D brane in 4D embedding: slice RGBA videos
# -----------------------------------------------------------------------------


def _extract_slice_scalar(
    flat: np.ndarray,
    grid_shape: Tuple[int, int, int],
    plane: str,
    index: Optional[int] = None,
) -> np.ndarray:
    nx, ny, nz = grid_shape
    F = flat.reshape(nx, ny, nz)
    plane = plane.lower()
    if plane == "xy":
        if index is None:
            index = nz // 2
        return F[:, :, index]
    if plane == "xz":
        if index is None:
            index = ny // 2
        return F[:, index, :]
    if plane == "yz":
        if index is None:
            index = nx // 2
        return F[index, :, :]
    raise ValueError(f"Unknown plane '{plane}', expected 'xy', 'xz', or 'yz'.")


def create_berry_slices_videos_3d(
    values_frames: np.ndarray,
    alpha_frames: np.ndarray,
    times_s: Sequence[float],
    grid_shape: Tuple[int, int, int],
    spacing: float,
    output_dir: str,
    filename_prefix: str,
    title_prefix: str,
    planes: Sequence[str] = ("xy", "xz", "yz"),
    fps: int = 20,
    dpi: int = 120,
    display_scale: float = 1e9,
    unit_label: str = "nm",
) -> List[str]:
    os.makedirs(output_dir, exist_ok=True)

    values = np.asarray(values_frames, dtype=np.float64)
    alpha = np.asarray(alpha_frames, dtype=np.float64)
    times = np.asarray(times_s, dtype=np.float64)

    nx, ny, nz = grid_shape
    T, N = values.shape
    if N != nx * ny * nz:
        raise ValueError(f"values_frames has N={N}, but grid_shape implies {nx*ny*nz}.")
    if alpha.shape != values.shape:
        raise ValueError("alpha_frames must match values_frames.")

    extent_xy = [0, (nx - 1) * spacing * display_scale, 0, (ny - 1) * spacing * display_scale]
    extent_xz = [0, (nx - 1) * spacing * display_scale, 0, (nz - 1) * spacing * display_scale]
    extent_yz = [0, (ny - 1) * spacing * display_scale, 0, (nz - 1) * spacing * display_scale]
    extents = {"xy": extent_xy, "xz": extent_xz, "yz": extent_yz}
    axis_labels = {
        "xy": (f"x [{unit_label}]", f"y [{unit_label}]"),
        "xz": (f"x [{unit_label}]", f"z [{unit_label}]"),
        "yz": (f"y [{unit_label}]", f"z [{unit_label}]"),
    }

    saved: List[str] = []

    for plane in planes:
        plane = plane.lower()
        fig, ax = plt.subplots(figsize=(9, 7))
        ax.set_aspect("equal" if plane == "xy" else "auto")
        ax.set_xlabel(axis_labels[plane][0])
        ax.set_ylabel(axis_labels[plane][1])

        v0 = _extract_slice_scalar(values[0], grid_shape, plane)
        a0 = _extract_slice_scalar(alpha[0], grid_shape, plane)
        rgba0 = scalar_to_rgba(v0, a0)

        im = ax.imshow(
            rgba0.swapaxes(0, 1),
            origin="lower",
            extent=extents[plane],
            interpolation="nearest",
            animated=True,
        )

        time_text = ax.text(
            0.02,
            0.95,
            "",
            transform=ax.transAxes,
            fontsize=12,
            va="top",
            bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.8),
        )
        ax.set_title(f"{title_prefix} — {plane.upper()} slice")

        def update(k: int):
            v = _extract_slice_scalar(values[k], grid_shape, plane)
            a = _extract_slice_scalar(alpha[k], grid_shape, plane)
            rgba = scalar_to_rgba(v, a)
            im.set_data(rgba.swapaxes(0, 1))
            t_as = times[k] * 1e18
            time_text.set_text(f"t = {t_as:.3f} as")
            return (im, time_text)

        anim = FuncAnimation(fig, update, frames=T, interval=1000 / fps, blit=True)
        out_path = os.path.join(output_dir, f"{filename_prefix}_{plane}.mp4")
        writer = FFMpegWriter(fps=fps, bitrate=2500)
        anim.save(out_path, writer=writer, dpi=dpi)
        plt.close(fig)
        saved.append(out_path)

    return saved


def create_berry_videos_3d(
    series: BerryTimeSeries,
    grid_shape: Tuple[int, int, int],
    spacing: float,
    output_dir: str,
    filename_prefix: str = "berry_3d",
    planes: Sequence[str] = ("xy", "xz", "yz"),
    fps: int = 20,
    dpi: int = 120,
    display_scale: float = 1e9,
    unit_label: str = "nm",
) -> List[str]:
    os.makedirs(output_dir, exist_ok=True)

    saved = []
    saved += create_berry_slices_videos_3d(
        values_frames=series.phase,
        alpha_frames=series.alpha,
        times_s=series.times_s,
        grid_shape=grid_shape,
        spacing=spacing,
        output_dir=output_dir,
        filename_prefix=f"{filename_prefix}_phase",
        title_prefix="Berry phase γ (hue) with amplitude mask (alpha)",
        planes=planes,
        fps=fps,
        dpi=dpi,
        display_scale=display_scale,
        unit_label=unit_label,
    )
    saved += create_berry_slices_videos_3d(
        values_frames=series.connection,
        alpha_frames=series.alpha,
        times_s=series.times_s,
        grid_shape=grid_shape,
        spacing=spacing,
        output_dir=output_dir,
        filename_prefix=f"{filename_prefix}_connection",
        title_prefix="Berry connection A_t (hue) with amplitude mask (alpha)",
        planes=planes,
        fps=fps,
        dpi=dpi,
        display_scale=display_scale,
        unit_label=unit_label,
    )
    return saved


# -----------------------------------------------------------------------------
# Convenience dispatcher
# -----------------------------------------------------------------------------


def create_berry_videos(
    series: BerryTimeSeries,
    intrinsic_dim: int,
    output_dir: str,
    *,
    x_coords_1d: Optional[np.ndarray] = None,
    grid_shape_2d: Optional[Tuple[int, int]] = None,
    grid_shape_3d: Optional[Tuple[int, int, int]] = None,
    spacing: Optional[float] = None,
    filename_prefix: str = "berry",
    fps: int = 20,
    dpi: int = 120,
    display_scale: float = 1e9,
    unit_label: str = "nm",
) -> List[str]:
    """Create Berry videos for intrinsic dimension 1/2/3.

    Parameters
    ----------
    intrinsic_dim:
        1, 2, or 3 (brane dimension).
    x_coords_1d:
        Required if intrinsic_dim==1. Display coordinates.
    grid_shape_2d, grid_shape_3d, spacing:
        Required for intrinsic_dim==2/3.
    """
    if intrinsic_dim == 1:
        if x_coords_1d is None:
            raise ValueError("x_coords_1d is required for intrinsic_dim==1")
        return create_berry_videos_1d(
            series=series,
            x_coords=x_coords_1d,
            output_dir=output_dir,
            filename_prefix=filename_prefix,
            fps=fps,
            dpi=max(dpi, 140),
            unit_label_x=unit_label,
        )

    if intrinsic_dim == 2:
        if grid_shape_2d is None or spacing is None:
            raise ValueError("grid_shape_2d and spacing are required for intrinsic_dim==2")
        return create_berry_videos_2d(
            series=series,
            grid_shape=grid_shape_2d,
            spacing=float(spacing),
            output_dir=output_dir,
            filename_prefix=filename_prefix,
            fps=fps,
            dpi=dpi,
            display_scale=display_scale,
            unit_label=unit_label,
        )

    if intrinsic_dim == 3:
        if grid_shape_3d is None or spacing is None:
            raise ValueError("grid_shape_3d and spacing are required for intrinsic_dim==3")
        return create_berry_videos_3d(
            series=series,
            grid_shape=grid_shape_3d,
            spacing=float(spacing),
            output_dir=output_dir,
            filename_prefix=filename_prefix,
            planes=("xy", "xz", "yz"),
            fps=fps,
            dpi=dpi,
            display_scale=display_scale,
            unit_label=unit_label,
        )

    raise ValueError(f"Unsupported intrinsic_dim={intrinsic_dim}, expected 1,2,3")
