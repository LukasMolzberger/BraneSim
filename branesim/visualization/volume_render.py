"""3D volume renderer for spacelike slices of the worldvolume.

Adapted from git c0f1aaa7^ components/visualization/volume_render.py for the
current worldvolume format (branesim-block-v1) and the U(1) vortex seed.

Two render modes:
  - "opacity"  : voxel opacity = excess-energy density or |displacement|
  - "phase_rgb": voxel colour = U(1) carrier phase (HSV hue), opacity = |u|

Both produce 3D voxel animations with an orbiting camera via FFMpeg.

No physics logic.  Layer F (visualization) — see principles §2.1.
"""

from __future__ import annotations

import os
import tempfile
from pathlib import Path
from typing import Callable, List, Optional, Sequence, Tuple

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, FFMpegWriter
from matplotlib.colors import hsv_to_rgb

os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "matplotlib"))


# ---------------------------------------------------------------------------
# Camera helpers (unchanged from legacy)
# ---------------------------------------------------------------------------


def camera_orbit(frame_idx: int, num_frames: int) -> Tuple[float, float]:
    """Slowly orbiting camera: elevation oscillates, azimuth increases."""
    t = frame_idx / max(num_frames - 1, 1)
    elev = 22.0 + 8.0 * np.sin(2 * np.pi * t)
    azim = -60.0 + 360.0 * 0.5 * t
    return elev, azim


# ---------------------------------------------------------------------------
# Phase -> RGB (adapted from legacy berry.py _phase_to_rgb)
# ---------------------------------------------------------------------------


def phase_to_rgb(phase: np.ndarray) -> np.ndarray:
    """Convert phase (radians, any range) to RGB via HSV hue wheel.

    Parameters
    ----------
    phase : ndarray, arbitrary shape
        Phase values in radians.

    Returns
    -------
    rgb : ndarray, same shape + trailing dimension 3, float in [0,1].
    """
    hue = np.mod((phase + np.pi) / (2.0 * np.pi), 1.0)
    sat = np.ones_like(hue)
    val = np.ones_like(hue)
    return hsv_to_rgb(np.stack([hue, sat, val], axis=-1))


# ---------------------------------------------------------------------------
# RGBA assembly for a volume frame
# ---------------------------------------------------------------------------


def _frame_to_rgba(
    amplitude: np.ndarray,
    phase: Optional[np.ndarray],
    gamma: float,
    alpha_scale: float,
    density_threshold: float,
    cmap_name: str,
    amp_max: float,
) -> Tuple[np.ndarray, np.ndarray]:
    """Convert an (nx, ny, nz) amplitude + optional phase into (mask, rgba).

    Parameters
    ----------
    amplitude : ndarray, shape (nx, ny, nz)
        Opacity-driving field (|displacement| or energy density).
    phase : ndarray or None, shape (nx, ny, nz)
        If given, use HSV phase colouring.  Otherwise use ``cmap_name``.
    gamma : float
        Opacity gamma (< 1 boosts faint features).
    alpha_scale : float
        Global opacity multiplier.
    density_threshold : float
        Fractional threshold below which opacity = 0 (hides vacuum noise).
    cmap_name : str
        Matplotlib colormap name used when phase is None.
    amp_max : float
        Reference amplitude for normalisation.

    Returns
    -------
    mask : bool ndarray, shape (nx, ny, nz)
    rgba : float ndarray, shape (nx, ny, nz, 4) in [0, 1]
    """
    norm = np.clip(amplitude / max(amp_max, 1e-30), 0.0, 1.0)
    mask = norm >= density_threshold

    if phase is not None:
        rgb = phase_to_rgb(phase)  # (nx, ny, nz, 3)
    else:
        cmap = plt.get_cmap(cmap_name)
        rgb = cmap(norm)[..., :3]  # (nx, ny, nz, 3)

    alpha = alpha_scale * (norm ** gamma)  # (nx, ny, nz)
    rgba = np.concatenate([rgb, alpha[..., None]], axis=-1)  # (nx, ny, nz, 4)
    return mask, rgba


# ---------------------------------------------------------------------------
# Main 3-D volume animation
# ---------------------------------------------------------------------------


def create_volume_animation(
    frames_amplitude: List[np.ndarray],
    frames_phase: Optional[List[np.ndarray]],
    grid_shape: Tuple[int, int, int],
    spacing: float,
    output_path: str,
    times: Optional[Sequence[float]] = None,
    cmap_name: str = "inferno",
    fps: int = 20,
    dpi: int = 100,
    camera_motion: Optional[Callable[[int, int], Tuple[float, float]]] = None,
    density_threshold: float = 0.005,
    alpha_scale: float = 1.0,
    gamma: float = 0.7,
    title_prefix: str = "",
) -> None:
    """Render a 3-D voxel animation from amplitude (and optional phase) frames.

    Parameters
    ----------
    frames_amplitude : list of ndarray, each shape (nx, ny, nz)
        Per-time-slice opacity-driving field.
    frames_phase : list of ndarray or None, each shape (nx, ny, nz)
        If provided, hue = U(1) phase; otherwise coloured by cmap.
    grid_shape : (nx, ny, nz)
    spacing : float
        Lattice spacing for axis labeling.
    output_path : str
        Output .mp4 file path.
    times : sequence of float or None
        Time values for title.
    cmap_name : str
        Fallback colormap when phase is None.
    fps, dpi : int
    camera_motion : callable or None
    density_threshold, alpha_scale, gamma : float
        Opacity tuning.
    title_prefix : str
    """
    if not frames_amplitude:
        raise ValueError("frames_amplitude must be non-empty")

    nx, ny, nz = grid_shape
    n_frames = len(frames_amplitude)
    times_arr = list(times) if times is not None else list(range(n_frames))

    # Global amplitude max for consistent normalisation across frames
    amp_max = max(float(np.max(np.abs(f))) for f in frames_amplitude)
    if amp_max < 1e-30:
        amp_max = 1.0

    # Voxel edge arrays
    x_edges = np.arange(nx + 1) * spacing
    y_edges = np.arange(ny + 1) * spacing
    z_edges = np.arange(nz + 1) * spacing
    X, Y, Z = np.meshgrid(x_edges, y_edges, z_edges, indexing="ij")

    cam = camera_motion or camera_orbit

    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection="3d")

    def _setup_axes() -> None:
        ax.set_xlim(0.0, nx * spacing)
        ax.set_ylim(0.0, ny * spacing)
        ax.set_zlim(0.0, nz * spacing)
        ax.set_box_aspect([nx, ny, nz])
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_zticks([])
        ax.set_axis_off()

    def _update(frame_idx: int):
        ax.cla()
        _setup_axes()

        amp = frames_amplitude[frame_idx]
        ph = frames_phase[frame_idx] if frames_phase else None

        mask, rgba = _frame_to_rgba(
            amp, ph,
            gamma=gamma,
            alpha_scale=alpha_scale,
            density_threshold=density_threshold,
            cmap_name=cmap_name,
            amp_max=amp_max,
        )
        if np.any(mask):
            ax.voxels(X, Y, Z, mask, facecolors=rgba, edgecolor="none")

        t_label = f"{float(times_arr[frame_idx]):.2f}"
        colour_label = "phase→RGB" if frames_phase else cmap_name
        ax.set_title(
            f"{title_prefix}  t={t_label}  [{colour_label}]",
            fontsize=9,
        )
        elev, azim = cam(frame_idx, n_frames)
        ax.view_init(elev=elev, azim=azim)
        return ax.collections

    anim = FuncAnimation(fig, _update, frames=n_frames, interval=1000 // fps, blit=False)
    out = Path(output_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    writer = FFMpegWriter(fps=fps, bitrate=2400)
    anim.save(str(out), writer=writer, dpi=dpi)
    plt.close(fig)


# ---------------------------------------------------------------------------
# 2-D midplane slice animation
# ---------------------------------------------------------------------------


def create_slice_animation(
    frames_field: List[np.ndarray],
    frames_phase: Optional[List[np.ndarray]],
    grid_shape: Tuple[int, int, int],
    spacing: float,
    plane: str,
    output_path: str,
    times: Optional[Sequence[float]] = None,
    cmap_name: str = "inferno",
    fps: int = 20,
    dpi: int = 120,
    title_prefix: str = "",
) -> None:
    """Render a 2-D midplane slice animation.

    Parameters
    ----------
    frames_field : list of ndarray, each shape (nx, ny, nz)
        Scalar field (e.g. |displacement|).
    frames_phase : list of ndarray or None, each shape (nx, ny, nz)
        If provided, renders phase→RGB; otherwise uses cmap.
    plane : str
        One of ``"xy"``, ``"xz"``, ``"yz"``.
    output_path : str
    """
    if not frames_field:
        raise ValueError("frames_field must be non-empty")
    nx, ny, nz = grid_shape
    n_frames = len(frames_field)
    times_arr = list(times) if times is not None else list(range(n_frames))

    def _extract_slice(vol: np.ndarray) -> np.ndarray:
        if plane == "xy":
            return vol[:, :, nz // 2]
        elif plane == "xz":
            return vol[:, ny // 2, :]
        elif plane == "yz":
            return vol[nx // 2, :, :]
        else:
            raise ValueError(f"plane must be xy/xz/yz; got {plane!r}")

    if plane == "xy":
        extent = [0, spacing * nx, 0, spacing * ny]
        xlabel, ylabel = "x", "y"
    elif plane == "xz":
        extent = [0, spacing * nx, 0, spacing * nz]
        xlabel, ylabel = "x", "z"
    else:
        extent = [0, spacing * ny, 0, spacing * nz]
        xlabel, ylabel = "y", "z"

    use_phase = frames_phase is not None
    amp_max = max(float(np.max(np.abs(_extract_slice(f)))) for f in frames_field)
    if amp_max < 1e-30:
        amp_max = 1.0

    fig, ax = plt.subplots(figsize=(8, 7))
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    slice0_amp = _extract_slice(frames_field[0])
    if use_phase:
        slice0_ph = _extract_slice(frames_phase[0])
        rgb0 = phase_to_rgb(slice0_ph)
        alpha0 = np.clip(np.abs(slice0_amp) / amp_max, 0, 1)[..., None]
        rgba0 = np.concatenate([rgb0, alpha0], axis=-1)
        im = ax.imshow(
            rgba0.swapaxes(0, 1),
            origin="lower",
            extent=extent,
            interpolation="nearest",
            animated=True,
        )
        colour_label = "phase→RGB"
    else:
        im = ax.imshow(
            slice0_amp.T,
            origin="lower",
            cmap=cmap_name,
            extent=extent,
            vmin=0.0,
            vmax=amp_max,
            interpolation="nearest",
            animated=True,
        )
        plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        colour_label = cmap_name

    title = ax.set_title(
        f"{title_prefix}  {plane}-slice  t={float(times_arr[0]):.2f}  [{colour_label}]",
        fontsize=9,
    )

    def _update(frame_idx: int):
        sl_amp = _extract_slice(frames_field[frame_idx])
        if use_phase:
            sl_ph = _extract_slice(frames_phase[frame_idx])
            rgb = phase_to_rgb(sl_ph)
            alpha = np.clip(np.abs(sl_amp) / amp_max, 0, 1)[..., None]
            rgba = np.concatenate([rgb, alpha], axis=-1)
            im.set_data(rgba.swapaxes(0, 1))
        else:
            im.set_data(sl_amp.T)
        title.set_text(
            f"{title_prefix}  {plane}-slice  "
            f"t={float(times_arr[frame_idx]):.2f}  [{colour_label}]"
        )
        return im, title

    anim = FuncAnimation(fig, _update, frames=n_frames, interval=1000 // fps, blit=False)
    out = Path(output_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    writer = FFMpegWriter(fps=fps, bitrate=2600)
    anim.save(str(out), writer=writer, dpi=dpi)
    plt.close(fig)
