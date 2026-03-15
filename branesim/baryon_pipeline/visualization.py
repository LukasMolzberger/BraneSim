"""Visualization component: render trajectory files into independent output products."""

from __future__ import annotations

import os
import tempfile
from pathlib import Path
from typing import Any

os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "matplotlib"))
os.environ.setdefault("XDG_CACHE_HOME", tempfile.gettempdir())

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation, FFMpegWriter

from branesim.visualization.brane_volume_viz import create_3d_volume_animation, downsample_volume

from .io import iter_trajectory_frames, load_trajectory_array, load_trajectory_manifest


def _extract_plane(field_3d: np.ndarray, plane: str, index: int | None = None) -> np.ndarray:
    nx, ny, nz = field_3d.shape
    if plane == "xy":
        idx = nz // 2 if index is None else int(index)
        return field_3d[:, :, idx]
    if plane == "xz":
        idx = ny // 2 if index is None else int(index)
        return field_3d[:, idx, :]
    if plane == "yz":
        idx = nx // 2 if index is None else int(index)
        return field_3d[idx, :, :]
    raise ValueError(f"Unsupported plane {plane!r}; expected one of 'xy', 'xz', 'yz'")


def _frame_field(
    positions: np.ndarray,
    rest_positions: np.ndarray,
    component: int,
    use_displacement: bool,
) -> np.ndarray:
    if component < 0 or component >= positions.shape[1]:
        raise ValueError(f"component must be in [0,{positions.shape[1]-1}], got {component}")
    field = positions[:, component]
    if use_displacement:
        field = field - rest_positions[:, component]
    return field


def render_volume(
    trajectory_path: str | Path,
    output_path: str | Path,
    *,
    component: int = 3,
    use_displacement: bool = True,
    frame_stride: int = 1,
    subsample: int = 2,
    fps: int = 20,
    dpi: int = 120,
    title_template: str = "Baryon Volume (t = {:.3f})",
) -> dict[str, Any]:
    """Render a volumetric movie from compressed simulation checkpoints."""

    manifest = load_trajectory_manifest(trajectory_path)
    grid_shape = tuple(int(v) for v in manifest["lattice"]["grid_shape"])
    spacing = float(manifest["lattice"]["spacing"])

    rest_positions = load_trajectory_array(trajectory_path, "aux/rest_positions.npy")

    frames: list[np.ndarray] = []
    times: list[float] = []

    for frame in iter_trajectory_frames(trajectory_path, frame_stride=frame_stride):
        field_flat = _frame_field(frame.positions, rest_positions, component, use_displacement)
        frames.append(downsample_volume(field_flat, grid_shape, subsample_factor=subsample))
        times.append(frame.time)

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    create_3d_volume_animation(
        frames=frames,
        times=times,
        spacing=spacing * subsample,
        output_path=str(output_path),
        fps=fps,
        dpi=dpi,
        title_template=title_template,
    )

    return {
        "output": str(output_path),
        "num_frames": len(frames),
        "mode": "volume",
        "component": component,
        "use_displacement": use_displacement,
    }


def render_slice(
    trajectory_path: str | Path,
    output_path: str | Path,
    *,
    plane: str = "xy",
    index: int | None = None,
    component: int = 3,
    use_displacement: bool = True,
    frame_stride: int = 1,
    fps: int = 20,
    dpi: int = 140,
    cmap: str = "RdBu_r",
) -> dict[str, Any]:
    """Render a 2D slice movie from compressed simulation checkpoints."""

    manifest = load_trajectory_manifest(trajectory_path)
    grid_shape = tuple(int(v) for v in manifest["lattice"]["grid_shape"])
    spacing = float(manifest["lattice"]["spacing"])

    rest_positions = load_trajectory_array(trajectory_path, "aux/rest_positions.npy")

    slices: list[np.ndarray] = []
    times: list[float] = []

    for frame in iter_trajectory_frames(trajectory_path, frame_stride=frame_stride):
        field_flat = _frame_field(frame.positions, rest_positions, component, use_displacement)
        field_3d = field_flat.reshape(grid_shape)
        slices.append(_extract_plane(field_3d, plane=plane, index=index))
        times.append(frame.time)

    if not slices:
        raise ValueError("No frames found in trajectory")

    vmax = max(float(np.max(np.abs(s))) for s in slices)
    if vmax < 1e-12:
        vmax = 1.0

    if plane == "xy":
        extent = [0.0, spacing * (grid_shape[0] - 1), 0.0, spacing * (grid_shape[1] - 1)]
        xlabel, ylabel = "x", "y"
    elif plane == "xz":
        extent = [0.0, spacing * (grid_shape[0] - 1), 0.0, spacing * (grid_shape[2] - 1)]
        xlabel, ylabel = "x", "z"
    else:
        extent = [0.0, spacing * (grid_shape[1] - 1), 0.0, spacing * (grid_shape[2] - 1)]
        xlabel, ylabel = "y", "z"

    fig, ax = plt.subplots(figsize=(8, 7))
    im = ax.imshow(
        slices[0].T,
        origin="lower",
        cmap=cmap,
        vmin=-vmax,
        vmax=vmax,
        extent=extent,
        interpolation="nearest",
        animated=True,
    )
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    title = ax.set_title(f"Baryon slice ({plane}) t={times[0]:.3f}")
    plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

    def update(frame_idx: int):
        im.set_data(slices[frame_idx].T)
        title.set_text(f"Baryon slice ({plane}) t={times[frame_idx]:.3f}")
        return im, title

    anim = FuncAnimation(fig, update, frames=len(slices), interval=1000 / fps, blit=False)

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    writer = FFMpegWriter(fps=fps, bitrate=2600)
    anim.save(str(output_path), writer=writer, dpi=dpi)
    plt.close(fig)

    return {
        "output": str(output_path),
        "num_frames": len(slices),
        "mode": "slice",
        "plane": plane,
        "component": component,
        "use_displacement": use_displacement,
    }


def run_visualization_component(
    trajectory_path: str | Path,
    output_path: str | Path,
    *,
    mode: str,
    **kwargs: Any,
) -> dict[str, Any]:
    """Dispatch visualization mode independent of simulation execution."""

    if mode == "volume":
        return render_volume(trajectory_path, output_path, **kwargs)
    if mode == "slice":
        return render_slice(trajectory_path, output_path, **kwargs)
    raise ValueError(f"Unsupported visualization mode {mode!r}; expected 'volume' or 'slice'")
