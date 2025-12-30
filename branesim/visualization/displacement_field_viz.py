"""displacement_field_viz.py

Displacement-field visualizations for brane simulations.

This module focuses on the *displacement field*

    u(x,t) = X(x,t) - X0(x)

where X is the current embedding position of each brane lattice point and X0 is
the reference (rest) embedding.

It provides reusable video generators for

* 3D brane in 4D embedding:
  1) 3D arrow-field animation (projected to 3D)
  2) Direction/magnitude split as RGBA slices (XY/XZ/YZ) with:
       - color  = 4D unit-direction (S^3) encoded by HSV
       - alpha  = displacement magnitude

* 2D brane in 3D embedding:
  1) 3D arrow-field animation
  2) Direction/magnitude split as RGBA image:
       - color  = 3D unit-direction (S^2) encoded by HSV
       - alpha  = displacement magnitude

* 1D brane in 2D embedding:
  1) Raw components u^1(x,t), u^2(x,t) (line animation)
  2) Split magnitude a(x,t) and angle theta(x,t) (two-panel line animation)

Notes
-----
* The “arrow style” is implemented as line segments (fast, robust) and uses
  Line3DCollection for 3D.
* For the 3D brane-in-4D arrow visualization, we *project* u ∈ R^4 to a 3-vector.
  By default we take (u^1,u^2,u^3). If you prefer a different projection, pass
  `proj=(0,1,3)` etc.
"""

from __future__ import annotations

import os
from dataclasses import dataclass
from typing import Callable, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, FFMpegWriter
from matplotlib.colors import hsv_to_rgb


def displacement_from_positions(X: np.ndarray, X0: np.ndarray) -> np.ndarray:
    """Compute displacement u = X - X0.

    This is a tiny helper to keep experiment scripts clean.

    Parameters
    ----------
    X:
        (N,D) array of current embedding coordinates.
    X0:
        (N,D) array of reference embedding coordinates.

    Returns
    -------
    u:
        (N,D) displacement array.
    """
    X = np.asarray(X)
    X0 = np.asarray(X0)
    if X.shape != X0.shape:
        raise ValueError(f"X and X0 must have the same shape, got {X.shape} vs {X0.shape}.")
    return X - X0


def displacement_frames_from_positions_frames(
    frames_X: Sequence[np.ndarray],
    X0: np.ndarray,
) -> List[np.ndarray]:
    """Compute displacement frames from position frames and a fixed reference X0."""
    X0 = np.asarray(X0)
    return [displacement_from_positions(np.asarray(X), X0) for X in frames_X]


# -----------------------------------------------------------------------------
# Core math helpers
# -----------------------------------------------------------------------------

def displacement_magnitude(u: np.ndarray, axis: int = -1, eps: float = 1e-12) -> np.ndarray:
    """Return |u| with numerical safety."""
    return np.sqrt(np.maximum(np.sum(u * u, axis=axis), 0.0) + eps * 0.0)


def displacement_direction(u: np.ndarray, axis: int = -1, eps: float = 1e-12) -> np.ndarray:
    """Return u/|u|, defined as 0 where |u| ~ 0."""
    a = displacement_magnitude(u, axis=axis, eps=eps)
    a_safe = np.maximum(a, eps)
    return u / np.expand_dims(a_safe, axis=axis)


def alpha_from_magnitude(
    a: np.ndarray,
    a_max: float,
    gamma: float = 1.3,
    alpha_scale: float = 0.9,
    eps: float = 1e-12,
) -> np.ndarray:
    """Map magnitude to opacity (alpha)."""
    denom = max(a_max, eps)
    x = np.clip(a / denom, 0.0, 1.0)
    return np.clip(alpha_scale * (x ** gamma), 0.0, 1.0)


# -----------------------------------------------------------------------------
# Direction → color encodings
# -----------------------------------------------------------------------------

def s3_direction_to_rgb(n: np.ndarray) -> np.ndarray:
    """Encode a 4D unit direction n ∈ S^3 as RGB via an HSV scheme.

    Scheme (uses all 4 components):
      * hue        ← atan2(n2, n1)
      * saturation ← remapped n3
      * value      ← remapped n4

    This is continuous on S^3 except for the standard hue branch cut.

    Parameters
    ----------
    n : (...,4) array
        Unit directions.

    Returns
    -------
    rgb : (...,3) array in [0,1]
    """
    n1, n2, n3, n4 = n[..., 0], n[..., 1], n[..., 2], n[..., 3]
    hue = (np.arctan2(n2, n1) + np.pi) / (2.0 * np.pi)  # [0,1]
    sat = 0.25 + 0.75 * (0.5 * (n3 + 1.0))              # [0.25,1]
    val = 0.25 + 0.75 * (0.5 * (n4 + 1.0))              # [0.25,1]
    hsv = np.stack([hue, np.clip(sat, 0, 1), np.clip(val, 0, 1)], axis=-1)
    return hsv_to_rgb(hsv)


def s2_direction_to_rgb(n: np.ndarray) -> np.ndarray:
    """Encode a 3D unit direction n ∈ S^2 as RGB via HSV.

    Scheme:
      * hue   ← atan2(n_y, n_x)
      * value ← remapped n_z
      * saturation = 1
    """
    nx, ny, nz = n[..., 0], n[..., 1], n[..., 2]
    hue = (np.arctan2(ny, nx) + np.pi) / (2.0 * np.pi)
    sat = np.ones_like(hue)
    val = 0.25 + 0.75 * (0.5 * (nz + 1.0))
    hsv = np.stack([hue, sat, np.clip(val, 0, 1)], axis=-1)
    return hsv_to_rgb(hsv)


def s1_angle_to_rgb(theta: np.ndarray) -> np.ndarray:
    """Encode an angle theta ∈ (-pi,pi] as RGB using hue."""
    hue = (theta + np.pi) / (2.0 * np.pi)
    hsv = np.stack([hue, np.ones_like(hue), np.ones_like(hue)], axis=-1)
    return hsv_to_rgb(hsv)


# -----------------------------------------------------------------------------
# Coordinate helpers
# -----------------------------------------------------------------------------

def coords_from_grid_shape_3d(
    grid_shape: Tuple[int, int, int],
    spacing: float,
    display_scale: float = 1.0,
) -> np.ndarray:
    """Return flattened (N,3) coordinates for an (nx,ny,nz) lattice."""
    nx, ny, nz = grid_shape
    x = np.arange(nx) * spacing * display_scale
    y = np.arange(ny) * spacing * display_scale
    z = np.arange(nz) * spacing * display_scale
    X, Y, Z = np.meshgrid(x, y, z, indexing="ij")
    return np.stack([X.ravel(), Y.ravel(), Z.ravel()], axis=-1)


def coords_from_grid_shape_2d(
    grid_shape: Tuple[int, int],
    spacing: float,
    display_scale: float = 1.0,
    z0: float = 0.0,
) -> np.ndarray:
    """Return flattened (N,3) coordinates for an (nx,ny) lattice embedded in z=z0."""
    nx, ny = grid_shape
    x = np.arange(nx) * spacing * display_scale
    y = np.arange(ny) * spacing * display_scale
    X, Y = np.meshgrid(x, y, indexing="ij")
    Z = np.full_like(X, z0)
    return np.stack([X.ravel(), Y.ravel(), Z.ravel()], axis=-1)


def subsample_indices_3d(grid_shape: Tuple[int, int, int], subsample: int) -> np.ndarray:
    """Return flat indices for subsampling a 3D lattice."""
    nx, ny, nz = grid_shape
    ii = np.arange(0, nx, subsample)
    jj = np.arange(0, ny, subsample)
    kk = np.arange(0, nz, subsample)
    I, J, K = np.meshgrid(ii, jj, kk, indexing="ij")
    return (I * (ny * nz) + J * nz + K).ravel()


def subsample_indices_2d(grid_shape: Tuple[int, int], subsample: int) -> np.ndarray:
    """Return flat indices for subsampling a 2D lattice."""
    nx, ny = grid_shape
    ii = np.arange(0, nx, subsample)
    jj = np.arange(0, ny, subsample)
    I, J = np.meshgrid(ii, jj, indexing="ij")
    return (I * ny + J).ravel()


def set_equal_aspect_3d(ax, coords: np.ndarray, margin: float = 1.08) -> None:
    """Set equal aspect 3D limits from coords."""
    mins = coords.min(axis=0)
    maxs = coords.max(axis=0)
    ranges = maxs - mins
    max_range = float(np.max(ranges)) if np.max(ranges) > 0 else 1.0
    mid = 0.5 * (mins + maxs)
    half = 0.5 * max_range * margin
    ax.set_xlim(mid[0] - half, mid[0] + half)
    ax.set_ylim(mid[1] - half, mid[1] + half)
    ax.set_zlim(mid[2] - half, mid[2] + half)
    try:
        ax.set_box_aspect([1, 1, 1])
    except Exception:
        pass


def default_camera_orbit(frame_idx: int, num_frames: int) -> Tuple[float, float]:
    """A gentle orbital camera motion."""
    t = frame_idx / max(num_frames - 1, 1)
    elev = 18.0 + 10.0 * np.sin(2.0 * np.pi * t)  # small bob
    azim = -60.0 + 360.0 * 0.6 * t
    return elev, azim


# -----------------------------------------------------------------------------
# 3D brane in 4D: arrow-field video
# -----------------------------------------------------------------------------

def create_displacement_arrows_video_3d_in_4d(
    frames_u: Sequence[np.ndarray],
    times: Sequence[float],
    grid_shape: Tuple[int, int, int],
    spacing: float,
    output_path: str,
    proj: Tuple[int, int, int] = (0, 1, 2),
    subsample: int = 3,
    arrow_fraction_of_extent: float = 0.06,
    fps: int = 20,
    dpi: int = 120,
    display_scale: float = 1e9,
    unit_label: str = "nm",
    title_template: str = "Displacement arrows (t = {t:.2f} as)",
    camera_motion: Optional[Callable[[int, int], Tuple[float, float]]] = None,
) -> None:
    """Video: 3D arrows for a 3D brane embedded in 4D (u projected to 3D).

    Parameters
    ----------
    frames_u:
        Iterable of displacement arrays, each shaped (N,4) or (N,D>=max(proj)+1).
        N must equal nx*ny*nz.
    times:
        Physical times in seconds (only used for labeling).
    spacing:
        Lattice spacing in meters.
    display_scale:
        Multiply coordinates by this factor for display (default: 1e9 => nm).
    """
    from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
    from mpl_toolkits.mplot3d.art3d import Line3DCollection

    os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)

    # Precompute coordinates and subsampling indices.
    coords_all = coords_from_grid_shape_3d(grid_shape, spacing, display_scale=display_scale)
    idx = subsample_indices_3d(grid_shape, subsample=subsample)
    coords = coords_all[idx]

    # Global arrow length scale based on extent.
    extent = np.max(coords_all.max(axis=0) - coords_all.min(axis=0))
    arrow_len = float(arrow_fraction_of_extent * (extent if extent > 0 else 1.0))

    # Prepare figure.
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection="3d")
    set_equal_aspect_3d(ax, coords_all)
    ax.set_xlabel(f"x [{unit_label}]")
    ax.set_ylabel(f"y [{unit_label}]")
    ax.set_zlabel(f"z [{unit_label}]")

    cam = camera_motion or default_camera_orbit

    def make_segments(u_proj: np.ndarray) -> np.ndarray:
        """Create line segments (M,2,3) from coords and projected vectors."""
        v = u_proj[idx]
        norms = np.linalg.norm(v, axis=1)
        valid = norms > 1e-15  # Much lower threshold to capture full Gaussian envelope
        v_hat = np.zeros_like(v)
        v_hat[valid] = v[valid] / norms[valid, None]
        tips = coords + arrow_len * v_hat
        seg = np.stack([coords, tips], axis=1)
        return seg

    # Initial segments.
    u0 = np.asarray(frames_u[0])
    u0p = u0[:, list(proj)]
    seg0 = make_segments(u0p)
    lc = Line3DCollection(seg0, colors="C0", linewidths=1.0, alpha=0.85)
    ax.add_collection3d(lc)

    # Title text.
    t0_as = float(times[0]) * 1e18
    title = ax.text2D(
        0.5,
        0.95,
        title_template.format(t=t0_as),
        transform=ax.transAxes,
        ha="center",
        fontsize=13,
        fontweight="bold",
    )

    def update(frame_idx: int):
        u = np.asarray(frames_u[frame_idx])
        up = u[:, list(proj)]
        seg = make_segments(up)
        lc.set_segments(seg)
        elev, azim = cam(frame_idx, len(frames_u))
        ax.view_init(elev=elev, azim=azim)
        t_as = float(times[frame_idx]) * 1e18
        title.set_text(title_template.format(t=t_as))
        return (lc, title)

    anim = FuncAnimation(fig, update, frames=len(frames_u), interval=1000 / fps, blit=False)
    writer = FFMpegWriter(fps=fps, bitrate=2500)
    anim.save(output_path, writer=writer, dpi=dpi)
    plt.close(fig)


# -----------------------------------------------------------------------------
# 3D brane in 4D: direction (S^3) + magnitude (alpha) slice videos
# -----------------------------------------------------------------------------

def _extract_slice_u4(
    u_flat: np.ndarray,
    grid_shape: Tuple[int, int, int],
    plane: str,
    index: Optional[int] = None,
) -> np.ndarray:
    """Extract a 2D slice from u(x) shaped (N,4) on a 3D lattice."""
    nx, ny, nz = grid_shape
    u = u_flat.reshape(nx, ny, nz, -1)
    plane = plane.lower()
    if plane == "xy":
        if index is None:
            index = nz // 2
        return u[:, :, index, :]
    if plane == "xz":
        if index is None:
            index = ny // 2
        return u[:, index, :, :]
    if plane == "yz":
        if index is None:
            index = nx // 2
        return u[index, :, :, :]
    raise ValueError(f"Unknown plane '{plane}', expected 'xy', 'xz', or 'yz'.")


def create_displacement_diralpha_slices_videos_3d_in_4d(
    frames_u: Sequence[np.ndarray],
    times: Sequence[float],
    grid_shape: Tuple[int, int, int],
    spacing: float,
    output_dir: str,
    filename_prefix: str = "displacement_3d_diralpha",
    planes: Sequence[str] = ("xy", "xz", "yz"),
    fps: int = 20,
    dpi: int = 150,
    display_scale: float = 1e9,
    unit_label: str = "nm",
    alpha_gamma: float = 1.3,
    alpha_scale: float = 0.95,
) -> List[str]:
    """Create 3 slice videos (XY/XZ/YZ) with RGBA encoding:

      color = S^3 direction of u (unit vector in R^4)
      alpha = |u|
    """
    os.makedirs(output_dir, exist_ok=True)

    nx, ny, nz = grid_shape
    extent_xy = [0, (nx - 1) * spacing * display_scale, 0, (ny - 1) * spacing * display_scale]
    extent_xz = [0, (nx - 1) * spacing * display_scale, 0, (nz - 1) * spacing * display_scale]
    extent_yz = [0, (ny - 1) * spacing * display_scale, 0, (nz - 1) * spacing * display_scale]
    extents = {"xy": extent_xy, "xz": extent_xz, "yz": extent_yz}
    axis_labels = {
        "xy": (f"x [{unit_label}]", f"y [{unit_label}]"),
        "xz": (f"x [{unit_label}]", f"z [{unit_label}]"),
        "yz": (f"y [{unit_label}]", f"z [{unit_label}]"),
    }

    # Determine global max magnitude across all frames and requested planes.
    a_max = 0.0
    for u in frames_u:
        for plane in planes:
            sl = _extract_slice_u4(np.asarray(u), grid_shape, plane)
            a = np.linalg.norm(sl, axis=-1)
            a_max = max(a_max, float(np.max(a)))
    if a_max <= 1e-30:
        a_max = 1.0

    saved = []

    for plane in planes:
        plane = plane.lower()
        fig, ax = plt.subplots(figsize=(9, 7))
        ax.set_aspect("equal" if plane == "xy" else "auto")
        ax.set_xlabel(axis_labels[plane][0])
        ax.set_ylabel(axis_labels[plane][1])

        # Initialize image.
        sl0 = _extract_slice_u4(np.asarray(frames_u[0]), grid_shape, plane)
        a0 = np.linalg.norm(sl0, axis=-1)
        n0 = sl0 / np.maximum(a0[..., None], 1e-12)
        rgb0 = s3_direction_to_rgb(n0)
        alpha0 = alpha_from_magnitude(a0, a_max=a_max, gamma=alpha_gamma, alpha_scale=alpha_scale)
        rgba0 = np.concatenate([rgb0, alpha0[..., None]], axis=-1)

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
        ax.set_title(f"Displacement direction (color) + magnitude (alpha) — {plane.upper()} slice")

        def update(frame_idx: int):
            sl = _extract_slice_u4(np.asarray(frames_u[frame_idx]), grid_shape, plane)
            a = np.linalg.norm(sl, axis=-1)
            n = sl / np.maximum(a[..., None], 1e-12)
            rgb = s3_direction_to_rgb(n)
            alpha = alpha_from_magnitude(a, a_max=a_max, gamma=alpha_gamma, alpha_scale=alpha_scale)
            rgba = np.concatenate([rgb, alpha[..., None]], axis=-1)
            im.set_data(rgba.swapaxes(0, 1))
            t_as = float(times[frame_idx]) * 1e18
            time_text.set_text(f"t = {t_as:.3f} as")
            return (im, time_text)

        anim = FuncAnimation(fig, update, frames=len(frames_u), interval=1000 / fps, blit=True)
        out_path = os.path.join(output_dir, f"{filename_prefix}_{plane}.mp4")
        writer = FFMpegWriter(fps=fps, bitrate=2500)
        anim.save(out_path, writer=writer, dpi=dpi)
        plt.close(fig)
        saved.append(out_path)

    return saved


# -----------------------------------------------------------------------------
# 2D brane in 3D: arrow-field video
# -----------------------------------------------------------------------------

def create_displacement_arrows_video_2d_in_3d(
    frames_u: Sequence[np.ndarray],
    times: Sequence[float],
    grid_shape: Tuple[int, int],
    spacing: float,
    output_path: str,
    subsample: int = 2,
    arrow_fraction_of_extent: float = 0.08,
    fps: int = 20,
    dpi: int = 140,
    display_scale: float = 1e9,
    unit_label: str = "nm",
    title_template: str = "Displacement arrows (t = {t:.2f} as)",
    camera_motion: Optional[Callable[[int, int], Tuple[float, float]]] = None,
) -> None:
    """Video: 3D arrows for a 2D brane embedded in 3D."""
    from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
    from mpl_toolkits.mplot3d.art3d import Line3DCollection

    os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)

    coords_all = coords_from_grid_shape_2d(grid_shape, spacing, display_scale=display_scale, z0=0.0)
    idx = subsample_indices_2d(grid_shape, subsample=subsample)
    coords = coords_all[idx]

    extent = np.max(coords_all.max(axis=0) - coords_all.min(axis=0))
    arrow_scale = float(arrow_fraction_of_extent * (extent if extent > 0 else 1.0))

    # Find global max displacement magnitude for scaling
    u_max = 0.0
    for u in frames_u:
        u_arr = np.asarray(u)
        u_max = max(u_max, float(np.max(np.linalg.norm(u_arr, axis=1))))
    if u_max <= 1e-30:
        u_max = 1.0

    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection="3d")
    set_equal_aspect_3d(ax, coords_all)
    ax.set_xlabel(f"x [{unit_label}]")
    ax.set_ylabel(f"y [{unit_label}]")
    ax.set_zlabel(f"z [{unit_label}]")

    cam = camera_motion or (lambda i, n: (25.0, -60.0 + 360.0 * 0.35 * (i / max(n - 1, 1))))

    def make_segments(u3: np.ndarray) -> np.ndarray:
        v = u3[idx]
        norms = np.linalg.norm(v, axis=1)
        valid = norms > 1e-15  # Much lower threshold to capture full Gaussian envelope
        # Scale arrows by actual magnitude (not normalized)
        scale_factors = np.zeros(len(v))
        scale_factors[valid] = (norms[valid] / u_max) * arrow_scale
        tips = coords + v * (scale_factors / np.maximum(norms, 1e-30))[:, None]
        return np.stack([coords, tips], axis=1)

    u0 = np.asarray(frames_u[0])
    seg0 = make_segments(u0)
    lc = Line3DCollection(seg0, colors="C0", linewidths=1.1, alpha=0.85)
    ax.add_collection3d(lc)

    t0_as = float(times[0]) * 1e18
    title = ax.text2D(0.5, 0.95, title_template.format(t=t0_as), transform=ax.transAxes,
                      ha="center", fontsize=13, fontweight="bold")

    def update(frame_idx: int):
        u = np.asarray(frames_u[frame_idx])
        lc.set_segments(make_segments(u))
        elev, azim = cam(frame_idx, len(frames_u))
        ax.view_init(elev=elev, azim=azim)
        t_as = float(times[frame_idx]) * 1e18
        title.set_text(title_template.format(t=t_as))
        return (lc, title)

    anim = FuncAnimation(fig, update, frames=len(frames_u), interval=1000 / fps, blit=False)
    writer = FFMpegWriter(fps=fps, bitrate=2500)
    anim.save(output_path, writer=writer, dpi=dpi)
    plt.close(fig)


# -----------------------------------------------------------------------------
# 2D brane in 3D: direction (S^2) + magnitude (alpha) video
# -----------------------------------------------------------------------------

def create_displacement_diralpha_video_2d_in_3d(
    frames_u: Sequence[np.ndarray],
    times: Sequence[float],
    grid_shape: Tuple[int, int],
    spacing: float,
    output_path: str,
    fps: int = 20,
    dpi: int = 160,
    display_scale: float = 1e9,
    unit_label: str = "nm",
    alpha_gamma: float = 1.3,
    alpha_scale: float = 0.95,
) -> None:
    """Video: RGBA encoding on the full 2D brane.

      color = S^2 direction (unit vector in R^3)
      alpha = |u|
    """
    os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)

    nx, ny = grid_shape
    extent = [0, (nx - 1) * spacing * display_scale, 0, (ny - 1) * spacing * display_scale]

    # Global max magnitude.
    a_max = 0.0
    for u in frames_u:
        uu = np.asarray(u).reshape(nx, ny, -1)
        a_max = max(a_max, float(np.max(np.linalg.norm(uu, axis=-1))))
    if a_max <= 1e-30:
        a_max = 1.0

    fig, ax = plt.subplots(figsize=(9, 7))
    ax.set_xlabel(f"x [{unit_label}]")
    ax.set_ylabel(f"y [{unit_label}]")
    ax.set_aspect("equal")
    ax.set_title("Displacement direction (color) + magnitude (alpha)")

    u0 = np.asarray(frames_u[0]).reshape(nx, ny, -1)
    a0 = np.linalg.norm(u0, axis=-1)
    n0 = u0 / np.maximum(a0[..., None], 1e-12)
    rgb0 = s2_direction_to_rgb(n0)
    alpha0 = alpha_from_magnitude(a0, a_max=a_max, gamma=alpha_gamma, alpha_scale=alpha_scale)
    rgba0 = np.concatenate([rgb0, alpha0[..., None]], axis=-1)

    im = ax.imshow(rgba0.swapaxes(0, 1), origin="lower", extent=extent, interpolation="nearest", animated=True)
    time_text = ax.text(
        0.02,
        0.95,
        "",
        transform=ax.transAxes,
        fontsize=12,
        va="top",
        bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.8),
    )

    def update(frame_idx: int):
        u = np.asarray(frames_u[frame_idx]).reshape(nx, ny, -1)
        a = np.linalg.norm(u, axis=-1)
        n = u / np.maximum(a[..., None], 1e-12)
        rgb = s2_direction_to_rgb(n)
        alpha = alpha_from_magnitude(a, a_max=a_max, gamma=alpha_gamma, alpha_scale=alpha_scale)
        rgba = np.concatenate([rgb, alpha[..., None]], axis=-1)
        im.set_data(rgba.swapaxes(0, 1))
        t_as = float(times[frame_idx]) * 1e18
        time_text.set_text(f"t = {t_as:.3f} as")
        return (im, time_text)

    anim = FuncAnimation(fig, update, frames=len(frames_u), interval=1000 / fps, blit=True)
    writer = FFMpegWriter(fps=fps, bitrate=2500)
    anim.save(output_path, writer=writer, dpi=dpi)
    plt.close(fig)


# -----------------------------------------------------------------------------
# 1D brane in 2D: component + polar split videos
# -----------------------------------------------------------------------------

def create_displacement_components_video_1d_in_2d(
    frames_u: Sequence[np.ndarray],
    times: Sequence[float],
    x_coords: np.ndarray,
    output_path: str,
    fps: int = 30,
    dpi: int = 140,
    unit_label_x: str = "nm",
    unit_label_u: str = "pm",
    title: str = "1D displacement components",
) -> None:
    """Video: plot u^1(x,t) and u^2(x,t) as two lines."""
    os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)

    x = np.asarray(x_coords)
    u0 = np.asarray(frames_u[0])
    if u0.ndim != 2 or u0.shape[1] < 2:
        raise ValueError("frames_u must contain arrays shaped (N,2) (or more columns).")

    # Determine global y-limits.
    all_u = np.concatenate([np.asarray(u)[:, :2] for u in frames_u], axis=0)
    ymax = float(np.max(np.abs(all_u)))
    if ymax <= 1e-30:
        ymax = 1.0
    ylim = 1.15 * ymax

    fig, ax = plt.subplots(figsize=(12, 6))
    ax.set_title(title, fontsize=14, fontweight="bold")
    ax.set_xlabel(f"x [{unit_label_x}]")
    ax.set_ylabel(f"u [{unit_label_u}]")
    ax.grid(True, alpha=0.3)
    ax.set_xlim(float(x.min()), float(x.max()))
    ax.set_ylim(-ylim, ylim)

    line1, = ax.plot(x, u0[:, 0], linewidth=2.0, label="u¹")
    line2, = ax.plot(x, u0[:, 1], linewidth=2.0, label="u²")
    ax.legend(loc="upper right")

    time_text = ax.text(
        0.02,
        0.95,
        "",
        transform=ax.transAxes,
        fontsize=12,
        va="top",
        bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.8),
    )

    def update(frame_idx: int):
        u = np.asarray(frames_u[frame_idx])
        line1.set_ydata(u[:, 0])
        line2.set_ydata(u[:, 1])
        t_as = float(times[frame_idx]) * 1e18
        time_text.set_text(f"t = {t_as:.3f} as")
        return (line1, line2, time_text)

    anim = FuncAnimation(fig, update, frames=len(frames_u), interval=1000 / fps, blit=True)
    writer = FFMpegWriter(fps=fps, bitrate=2500)
    anim.save(output_path, writer=writer, dpi=dpi)
    plt.close(fig)


def create_displacement_magnitude_angle_video_1d_in_2d(
    frames_u: Sequence[np.ndarray],
    times: Sequence[float],
    x_coords: np.ndarray,
    output_path: str,
    fps: int = 30,
    dpi: int = 140,
    unit_label_x: str = "nm",
    unit_label_a: str = "pm",
    title: str = "1D displacement: magnitude + angle",
    unwrap_along_x: bool = True,
    magnitude_threshold: float = 0.01,  # Relative threshold for angle masking
) -> None:
    """Video: two-panel plot of magnitude a(x,t) and angle theta(x,t).

    The angle is only plotted where the magnitude exceeds magnitude_threshold * amax,
    since the angle becomes numerically unstable when |u| ≈ 0.
    """
    os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)

    x = np.asarray(x_coords)
    u0 = np.asarray(frames_u[0])
    if u0.ndim != 2 or u0.shape[1] < 2:
        raise ValueError("frames_u must contain arrays shaped (N,2) (or more columns).")

    # Precompute global limits.
    mags = [np.linalg.norm(np.asarray(u)[:, :2], axis=1) for u in frames_u]
    amax = float(np.max(np.concatenate(mags)))
    if amax <= 1e-30:
        amax = 1.0

    # Threshold for valid angle (relative to maximum magnitude)
    a_threshold = magnitude_threshold * amax

    fig, (ax_a, ax_th) = plt.subplots(2, 1, figsize=(12, 8), sharex=True)
    fig.suptitle(title, fontsize=14, fontweight="bold")

    # Initial data
    a0 = mags[0]
    th0 = np.arctan2(u0[:, 1], u0[:, 0])

    # Mask where magnitude is too small
    mask0 = a0 > a_threshold
    x_masked0 = x[mask0]
    th0_masked = th0[mask0]

    if unwrap_along_x and len(th0_masked) > 0:
        th0_masked = np.unwrap(th0_masked)

    ax_a.set_ylabel(f"|u| [{unit_label_a}]")
    ax_a.grid(True, alpha=0.3)
    ax_a.set_ylim(0.0, 1.15 * amax)
    line_a, = ax_a.plot(x, a0, linewidth=2.0, color='C0')

    ax_th.set_ylabel("θ [rad]")
    ax_th.set_xlabel(f"x [{unit_label_x}]")
    ax_th.grid(True, alpha=0.3)
    ax_th.set_ylim(-np.pi, np.pi)  # Standard angle range
    line_th, = ax_th.plot(x_masked0, th0_masked, linewidth=2.0, color='C1', marker='.')

    # Add threshold indicator on magnitude plot
    ax_a.axhline(y=a_threshold, color='red', linestyle='--', alpha=0.5, linewidth=1,
                 label=f'Angle threshold ({magnitude_threshold*100:.0f}% of max)')
    ax_a.legend(loc='upper right', fontsize=9)

    time_text = ax_a.text(
        0.02,
        0.85,
        "",
        transform=ax_a.transAxes,
        fontsize=12,
        va="top",
        bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.8),
    )

    def update(frame_idx: int):
        u = np.asarray(frames_u[frame_idx])
        a = np.linalg.norm(u[:, :2], axis=1)
        th = np.arctan2(u[:, 1], u[:, 0])

        # Mask where magnitude is too small
        mask = a > a_threshold
        x_masked = x[mask]
        th_masked = th[mask]

        if unwrap_along_x and len(th_masked) > 0:
            th_masked = np.unwrap(th_masked)

        # Update magnitude line
        line_a.set_ydata(a)

        # Update angle line with masked data
        line_th.set_data(x_masked, th_masked)

        t_as = float(times[frame_idx]) * 1e18
        time_text.set_text(f"t = {t_as:.3f} as")
        return (line_a, line_th, time_text)

    anim = FuncAnimation(fig, update, frames=len(frames_u), interval=1000 / fps, blit=True)
    writer = FFMpegWriter(fps=fps, bitrate=2500)
    anim.save(output_path, writer=writer, dpi=dpi)
    plt.close(fig)
