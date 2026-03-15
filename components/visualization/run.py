"""Component 3: visualization.

Renders volume or slice videos from simulation output file.
"""

from __future__ import annotations

import argparse
import os
import tempfile
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "matplotlib"))
os.environ.setdefault("XDG_CACHE_HOME", tempfile.gettempdir())

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation, FFMpegWriter

from components.shared import iter_frames, load_manifest, load_npy
from components.visualization.volume_render import create_3d_volume_animation, downsample_volume


def _field_values(
    positions: np.ndarray,
    rest_positions: np.ndarray,
    field_mode: str,
    component: int,
    use_displacement: bool,
) -> np.ndarray:
    disp = positions - rest_positions

    if field_mode == "component":
        src = disp if use_displacement else positions
        if component < 0 or component >= src.shape[1]:
            raise ValueError("component index out of range")
        return src[:, component]

    if field_mode == "magnitude_xyz":
        src = disp[:, :3] if use_displacement else positions[:, :3]
        return np.linalg.norm(src, axis=1)

    if field_mode == "magnitude_all":
        src = disp if use_displacement else positions
        return np.linalg.norm(src, axis=1)

    raise ValueError(f"Unknown field_mode {field_mode!r}")


def _extract_plane(field_3d: np.ndarray, plane: str, index: int | None) -> np.ndarray:
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
    raise ValueError("plane must be one of: xy,xz,yz")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Component 3: visualization")
    p.add_argument("--input", required=True)
    p.add_argument("--output", required=True)
    p.add_argument("--mode", required=True, choices=("volume", "slice"))

    p.add_argument("--field-mode", default="magnitude_all", choices=("component", "magnitude_xyz", "magnitude_all"))
    p.add_argument("--component", type=int, default=3)
    p.add_argument("--absolute-field", action="store_true")

    p.add_argument("--frame-stride", type=int, default=1)
    p.add_argument("--fps", type=int, default=20)
    p.add_argument("--dpi", type=int, default=120)

    p.add_argument("--subsample", type=int, default=2)
    p.add_argument("--density-threshold", type=float, default=0.003)
    p.add_argument("--alpha-scale", type=float, default=1.0)
    p.add_argument("--gamma", type=float, default=0.8)
    p.add_argument("--color-gamma", type=float, default=1.0)
    p.add_argument("--cmap", type=str, default=None)

    p.add_argument("--plane", default="xy", choices=("xy", "xz", "yz"))
    p.add_argument("--index", type=int, default=None)
    return p.parse_args()


def _run_volume(args: argparse.Namespace, grid_shape: tuple[int, int, int], spacing: float, rest: np.ndarray) -> None:
    frames = []
    times = []
    for fr in iter_frames(args.input, frame_stride=args.frame_stride):
        field_flat = _field_values(
            fr.positions,
            rest,
            field_mode=args.field_mode,
            component=args.component,
            use_displacement=not args.absolute_field,
        )
        frames.append(downsample_volume(field_flat, grid_shape, subsample_factor=args.subsample))
        times.append(fr.time)

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)

    cmap = args.cmap if args.cmap else ("inferno" if args.field_mode.startswith("magnitude") else "RdBu_r")
    create_3d_volume_animation(
        frames=frames,
        times=times,
        spacing=spacing * args.subsample,
        output_path=str(out),
        cmap_name=cmap,
        fps=args.fps,
        dpi=args.dpi,
        density_threshold=args.density_threshold,
        alpha_scale=args.alpha_scale,
        gamma=args.gamma,
        color_gamma=args.color_gamma,
    )

    print("Visualization complete")
    print("  mode: volume")
    print(f"  output: {out}")
    print(f"  frames: {len(frames)}")


def _run_slice(args: argparse.Namespace, grid_shape: tuple[int, int, int], spacing: float, rest: np.ndarray) -> None:
    slices = []
    times = []

    for fr in iter_frames(args.input, frame_stride=args.frame_stride):
        field_flat = _field_values(
            fr.positions,
            rest,
            field_mode=args.field_mode,
            component=args.component,
            use_displacement=not args.absolute_field,
        )
        slices.append(_extract_plane(field_flat.reshape(grid_shape), args.plane, args.index))
        times.append(fr.time)

    vmax = max(float(np.max(np.abs(s))) for s in slices)
    if vmax < 1e-12:
        vmax = 1.0

    if args.plane == "xy":
        extent = [0.0, spacing * (grid_shape[0] - 1), 0.0, spacing * (grid_shape[1] - 1)]
        xlabel, ylabel = "x", "y"
    elif args.plane == "xz":
        extent = [0.0, spacing * (grid_shape[0] - 1), 0.0, spacing * (grid_shape[2] - 1)]
        xlabel, ylabel = "x", "z"
    else:
        extent = [0.0, spacing * (grid_shape[1] - 1), 0.0, spacing * (grid_shape[2] - 1)]
        xlabel, ylabel = "y", "z"

    cmap = args.cmap if args.cmap else ("inferno" if args.field_mode.startswith("magnitude") else "RdBu_r")
    fig, ax = plt.subplots(figsize=(8, 7))
    im = ax.imshow(
        slices[0].T,
        origin="lower",
        cmap=cmap,
        extent=extent,
        vmin=0.0 if args.field_mode.startswith("magnitude") else -vmax,
        vmax=vmax,
        interpolation="nearest",
        animated=True,
    )
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    title = ax.set_title(f"Baryon Slice ({args.plane}) t={times[0]:.3f}")
    plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

    def update(frame_idx: int):
        im.set_data(slices[frame_idx].T)
        title.set_text(f"Baryon Slice ({args.plane}) t={times[frame_idx]:.3f}")
        return im, title

    anim = FuncAnimation(fig, update, frames=len(slices), interval=1000 / args.fps, blit=False)

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    writer = FFMpegWriter(fps=args.fps, bitrate=2600)
    anim.save(str(out), writer=writer, dpi=args.dpi)
    plt.close(fig)

    print("Visualization complete")
    print("  mode: slice")
    print(f"  output: {out}")
    print(f"  frames: {len(slices)}")


def main() -> None:
    args = parse_args()
    manifest = load_manifest(args.input)
    lattice = manifest["lattice"]
    grid_shape = tuple(int(v) for v in lattice["grid_shape"])
    spacing = float(lattice["spacing"])
    rest = load_npy(args.input, "aux/rest_positions.npy")

    if args.mode == "volume":
        _run_volume(args, grid_shape, spacing, rest)
    else:
        _run_slice(args, grid_shape, spacing, rest)


if __name__ == "__main__":
    main()
