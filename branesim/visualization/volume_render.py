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
from matplotlib.collections import LineCollection
from matplotlib.colors import hsv_to_rgb, to_rgba, PowerNorm

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
# Optional 2-D overlays (add-ons; off by default, used by create_slice_animation)
#
# Idiom borrowed from R. Behiel's GL-vortex animations (Vortex.py / PsiThetaPsi.py):
#   - short line "dashes" oriented by the local U(1) phase make the winding and the
#     core defect legible where flat hue alone is ambiguous (and they fade out in the
#     dark |Psi|->0 region, where the hue is meaningless);
#   - dots advected along the supercurrent j ~ |Psi|^2 grad(theta) visualise the
#     swirl around the vortex core.
# Pure visualization, no physics logic (principles §2.1).
# ---------------------------------------------------------------------------


def _inplane_current(phase_slice: np.ndarray, amp_slice: np.ndarray,
                     spacing: float) -> Tuple[np.ndarray, np.ndarray]:
    """Branch-cut-safe in-plane supercurrent j = Im(conj(Psi) grad Psi) = |Psi|^2 grad(theta).

    Reconstructs Psi = amp * exp(i*phase) so the gradient never crosses the phase
    branch cut.  Returns (j0, j1), components along the two slice axes.
    """
    psi = amp_slice * np.exp(1j * phase_slice)
    g0 = np.gradient(psi, spacing, axis=0)
    g1 = np.gradient(psi, spacing, axis=1)
    j0 = np.imag(np.conj(psi) * g0)
    j1 = np.imag(np.conj(psi) * g1)
    return j0, j1


def _phase_dash_segments(
    phase_slice: np.ndarray,
    amp_slice: np.ndarray,
    spacing: float,
    stride: int,
    length_frac: float,
    amp_max: float,
    amp_threshold: float,
    base_rgba: Tuple[float, float, float, float],
) -> Tuple[np.ndarray, np.ndarray]:
    """Build LineCollection segments + per-segment RGBA for the phase-dash overlay.

    Each sampled grid point gets a short segment centred on the cell, oriented along
    the local phase direction (cos theta, sin theta).  Segment opacity tracks the
    normalised amplitude (faded/hidden where |Psi| is small).
    """
    n0, n1 = phase_slice.shape
    seg_len = length_frac * stride * spacing
    ii = np.arange(0, n0, stride)
    jj = np.arange(0, n1, stride)
    I, J = np.meshgrid(ii, jj, indexing="ij")
    xc = (I + 0.5) * spacing
    yc = (J + 0.5) * spacing
    th = phase_slice[I, J]
    dx = 0.5 * seg_len * np.cos(th)
    dy = 0.5 * seg_len * np.sin(th)
    p0 = np.stack([xc - dx, yc - dy], axis=-1)
    p1 = np.stack([xc + dx, yc + dy], axis=-1)
    segs = np.stack([p0, p1], axis=-2).reshape(-1, 2, 2)

    amp_norm = np.clip(amp_slice[I, J].ravel() / max(amp_max, 1e-30), 0.0, 1.0)
    rgba = np.tile(np.asarray(base_rgba, dtype=float), (segs.shape[0], 1))
    rgba[:, 3] = np.where(amp_norm >= amp_threshold, base_rgba[3] * amp_norm, 0.0)
    return segs, rgba


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


def _write_placeholder_animation(
    output_path: str, fps: int, dpi: int, text: str
) -> None:
    """Write a valid 1-frame mp4 carrying ``text`` (E15 skip placeholder).

    Used when a volume render is skipped (degenerate field), so downstream
    consumers and ``.md`` pointers still find an mp4 rather than a dangling path.
    """
    fig = plt.figure(figsize=(10, 8))
    fig.text(0.5, 0.5, text, ha="center", va="center", fontsize=12, wrap=True)
    out = Path(output_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    writer = FFMpegWriter(fps=fps, bitrate=600)
    with writer.saving(fig, str(out), dpi):
        writer.grab_frame()
    plt.close(fig)


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
    min_signal: float = 1e-9,
    max_voxels_per_axis: int = 48,
) -> None:
    """Render a 3-D voxel animation from amplitude (and optional phase) frames.

    Robustness guards (E15; see EXPERIMENT_OPEN_PROBLEMS.md)
    -------------------------------------------------------
    ``min_signal`` : if the global ``max|field|`` is below this absolute floor the
        field carries no signal (e.g. the SU(3) channel of a *pure* U(1) seed,
        ~1e-11).  Relative normalization would then amplify numerical noise until
        the voxel mask fills, and matplotlib ``voxels()`` / ``autoscale_view``
        becomes pathological — this is the bug that hung a 96³ run for ~11 h.
        Such fields are rendered as a 1-frame "skipped" placeholder instead.
    ``max_voxels_per_axis`` : ``voxels()`` is O(voxels); a 96³ frame is intractable
        (~8–9 min EACH at best).  If any grid axis exceeds this, the field is
        coarsened by an integer stride before rendering so the work stays bounded.

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

    # Global amplitude max for consistent normalisation across frames.
    amp_max = max(float(np.max(np.abs(f))) for f in frames_amplitude)

    # --- E15 guard 1: degenerate near-zero field -> skip (placeholder frame) ---
    if amp_max < min_signal:
        _write_placeholder_animation(
            output_path, fps, dpi,
            f"{title_prefix}\n\nvolume render skipped\nmax|field| = {amp_max:.2e}"
            f"  <  min_signal = {min_signal:.0e}\n(field is numerically zero)",
        )
        return

    # --- E15 guard 2: cap voxel count (voxels() is O(voxels); 96^3 intractable) ---
    stride = max(1, int(np.ceil(max(grid_shape) / max_voxels_per_axis)))
    if stride > 1:
        frames_amplitude = [f[::stride, ::stride, ::stride] for f in frames_amplitude]
        if frames_phase is not None:
            frames_phase = [p[::stride, ::stride, ::stride] for p in frames_phase]
        nx, ny, nz = frames_amplitude[0].shape
        spacing = spacing * stride  # keep physical extent correct after coarsening

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
    phase_dashes: bool = False,
    dash_stride: int = 3,
    dash_length_frac: float = 0.85,
    dash_amp_threshold: float = 0.12,
    dash_color: str = "k",
    dash_alpha: float = 0.55,
    dash_linewidth: float = 0.9,
    flow_dots: bool = False,
    n_flow_dots: int = 350,
    flow_speed: float = 1.5,
    dot_lifetime_frames: int = 16,
    dot_color: str = "w",
    dot_size: float = 5.0,
    flow_seed: int = 0,
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

    Optional overlays (add-ons; all default off, base render unchanged)
    -------------------------------------------------------------------
    phase_dashes : bool
        Overlay short line segments oriented by the local U(1) phase (requires
        ``frames_phase``).  Opacity tracks amplitude so dashes fade where |Psi|→0.
    dash_stride, dash_length_frac, dash_amp_threshold, dash_color, dash_alpha,
    dash_linewidth :
        Dash sampling stride, length (fraction of stride·spacing), relative-amplitude
        cut below which a dash is hidden, colour, base opacity, line width.
    flow_dots : bool
        Overlay dots advected along the supercurrent j ~ |Psi|² ∇θ (requires
        ``frames_phase``), spawned with density ∝ amplitude and aged out — the
        swirl around the vortex core.
    n_flow_dots, flow_speed, dot_lifetime_frames, dot_color, dot_size, flow_seed :
        Target live-dot count, advection speed (cells/frame scale), dot lifetime in
        frames, colour, marker size, and the RNG seed (reproducible spawning).
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

    # --- optional overlays (add-ons; require phase, base render untouched) ---
    overlay_dashes = phase_dashes and use_phase
    overlay_dots = flow_dots and use_phase
    ext_x, ext_y = extent[1], extent[3]
    dash_lc = None
    if overlay_dashes:
        base_rgba = to_rgba(dash_color, dash_alpha)
        segs0, seg_rgba0 = _phase_dash_segments(
            slice0_ph, np.abs(slice0_amp), spacing, dash_stride,
            dash_length_frac, amp_max, dash_amp_threshold, base_rgba,
        )
        dash_lc = LineCollection(segs0, colors=seg_rgba0,
                                 linewidths=dash_linewidth, zorder=5)
        ax.add_collection(dash_lc)
    dot_scat = None
    dot_rgb = to_rgba(dot_color)[:3]
    dot_state = {"pos": np.zeros((0, 2)), "age": np.zeros((0,), dtype=int)}
    rng = np.random.default_rng(flow_seed)
    if overlay_dots:
        dot_scat = ax.scatter([], [], s=dot_size, facecolors="none",
                              edgecolors="none", zorder=8)

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
        artists = [im, title]

        if overlay_dashes:
            segs, seg_rgba = _phase_dash_segments(
                sl_ph, np.abs(sl_amp), spacing, dash_stride,
                dash_length_frac, amp_max, dash_amp_threshold,
                to_rgba(dash_color, dash_alpha),
            )
            dash_lc.set_segments(segs)
            dash_lc.set_color(seg_rgba)
            artists.append(dash_lc)

        if overlay_dots:
            pos, age = dot_state["pos"], dot_state["age"]
            if pos.shape[0]:
                j0, j1 = _inplane_current(sl_ph, np.abs(sl_amp), spacing)
                ix = np.clip((pos[:, 0] / spacing).astype(int), 0, sl_ph.shape[0] - 1)
                iy = np.clip((pos[:, 1] / spacing).astype(int), 0, sl_ph.shape[1] - 1)
                vx, vy = j0[ix, iy], j1[ix, iy]
                mag = np.hypot(vx, vy) + 1e-30
                pos = pos.copy()
                pos[:, 0] += flow_speed * spacing * vx / mag
                pos[:, 1] += flow_speed * spacing * vy / mag
                age = age + 1
            alive = (
                (age < dot_lifetime_frames)
                & (pos[:, 0] >= 0) & (pos[:, 0] <= ext_x)
                & (pos[:, 1] >= 0) & (pos[:, 1] <= ext_y)
            )
            pos, age = pos[alive], age[alive]
            n_spawn = max(0, n_flow_dots - pos.shape[0])
            if n_spawn:
                w = np.clip(np.abs(sl_amp) / amp_max, 0, 1).ravel()
                tot = w.sum()
                if tot > 0:
                    idx = rng.choice(w.size, size=n_spawn, p=w / tot)
                    n1 = sl_amp.shape[1]
                    sx = (idx // n1 + rng.random(n_spawn)) * spacing
                    sy = (idx % n1 + rng.random(n_spawn)) * spacing
                    pos = np.concatenate([pos, np.stack([sx, sy], axis=-1)], axis=0)
                    age = np.concatenate([age, np.zeros(n_spawn, dtype=int)])
            dot_state["pos"], dot_state["age"] = pos, age
            if pos.shape[0]:
                fade = np.sin(np.pi * age / max(dot_lifetime_frames, 1)) ** 2
                rgba_d = np.tile(np.append(dot_rgb, 1.0), (pos.shape[0], 1))
                rgba_d[:, 3] = fade
                dot_scat.set_offsets(pos)
                dot_scat.set_facecolors(rgba_d)
            else:
                dot_scat.set_offsets(np.empty((0, 2)))
            artists.append(dot_scat)

        return tuple(artists)

    anim = FuncAnimation(fig, _update, frames=n_frames, interval=1000 // fps, blit=False)
    out = Path(output_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    writer = FFMpegWriter(fps=fps, bitrate=2600)
    anim.save(str(out), writer=writer, dpi=dpi)
    plt.close(fig)


# ---------------------------------------------------------------------------
# U(1) vs SU(3) energy-channel comparison animation
# ---------------------------------------------------------------------------


def create_channel_energy_animation(
    frames_u1: List[np.ndarray],
    frames_su3: List[np.ndarray],
    grid_shape: Tuple[int, int, int],
    spacing: float,
    plane: str,
    output_path: str,
    u1_energy: Optional[Sequence[float]] = None,
    su3_energy: Optional[Sequence[float]] = None,
    times: Optional[Sequence[float]] = None,
    fps: int = 15,
    dpi: int = 120,
    shared_scale: bool = True,
    gamma: float = 1.0,
    title_prefix: str = "",
) -> None:
    """Side-by-side U(1) vs SU(3) energy-density slice movie + integrated curves.

    Three panels per frame:
      - Left   : U(1) (trace / EM) energy density on the midplane slice (Blues).
      - Centre : SU(3) (traceless / colour) energy density on the slice (Oranges).
      - Right  : box-integrated U(1) and SU(3) energy vs time, with a moving
                 time marker and the running U(1):SU(3) energy fraction.

    The two density fields are the per-node squared norms of the lateral
    displacement projected with P_U1 and P_SU3 (see
    ``diagnostics.alpha_separability.projection_operators``).  With
    ``shared_scale=True`` (default) both maps use one common colour normalisation
    so the eye reads the *relative* energy content honestly — a pure-trace seed
    shows a bright U(1) map and an empty SU(3) map.

    Parameters
    ----------
    frames_u1, frames_su3 : list of ndarray, each (nx, ny, nz)
        Per-time-slice U(1) and SU(3) energy-density volumes.
    grid_shape : (nx, ny, nz)
    spacing : float
    plane : str
        One of ``"xy"``, ``"xz"``, ``"yz"``.
    output_path : str
        Output .mp4 path.
    u1_energy, su3_energy : sequence of float or None
        Pre-integrated box energies per slice for the time-series panel.  If
        None, they are computed as the per-frame sum of the density volumes.
    times : sequence of float or None
    fps, dpi : int
    shared_scale : bool
        If True, both spatial maps share one colour normalisation (honest relative
        magnitude — but a 3%-of-U(1) SU(3) channel renders ~black).  If False, each
        map autoscales to its OWN max, so faint SU(3) structure is visible; the
        honest magnitude ratio is then carried by the integrated-energy panel.
    gamma : float
        Power-law display stretch (``PowerNorm``).  ``gamma<1`` brightens low
        densities → more contrast for faint structure.  ``1.0`` = linear.
    title_prefix : str
    """
    if not frames_u1 or not frames_su3:
        raise ValueError("frames_u1 and frames_su3 must be non-empty")
    if len(frames_u1) != len(frames_su3):
        raise ValueError("frames_u1 and frames_su3 must have equal length")

    nx, ny, nz = grid_shape
    n_frames = len(frames_u1)
    times_arr = list(times) if times is not None else list(range(n_frames))

    def _extract_slice(vol: np.ndarray) -> np.ndarray:
        if plane == "xy":
            return vol[:, :, nz // 2]
        elif plane == "xz":
            return vol[:, ny // 2, :]
        elif plane == "yz":
            return vol[nx // 2, :, :]
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

    # Integrated box energies for the time-series panel.
    u1_e = (np.asarray(u1_energy, dtype=float) if u1_energy is not None
            else np.array([float(f.sum()) for f in frames_u1]))
    su3_e = (np.asarray(su3_energy, dtype=float) if su3_energy is not None
             else np.array([float(f.sum()) for f in frames_su3]))

    # Colour normalisation across all frames (shared so magnitudes compare).
    vmax_u1 = max(float(np.max(_extract_slice(f))) for f in frames_u1)
    vmax_su3 = max(float(np.max(_extract_slice(f))) for f in frames_su3)
    if shared_scale:
        vmax_u1 = vmax_su3 = max(vmax_u1, vmax_su3)
    vmax_u1 = vmax_u1 or 1.0
    vmax_su3 = vmax_su3 or 1.0

    def _norm(vmax: float):
        """PowerNorm contrast stretch (gamma<1 brightens faint structure)."""
        return PowerNorm(gamma=gamma, vmin=0.0, vmax=vmax)

    fig, (ax_u1, ax_su3, ax_ts) = plt.subplots(1, 3, figsize=(20, 6.5))
    fig.suptitle(
        f"{title_prefix}  —  U(1) (trace/EM) vs SU(3) (traceless/colour) energy content",
        fontsize=11, fontweight="bold",
    )

    im_u1 = ax_u1.imshow(
        _extract_slice(frames_u1[0]).T, origin="lower", extent=extent,
        cmap="Blues", norm=_norm(vmax_u1), interpolation="nearest", animated=True,
    )
    ax_u1.set_xlabel(xlabel); ax_u1.set_ylabel(ylabel)
    ax_u1.set_title(f"U(1) trace energy density  ({plane} midplane)")
    plt.colorbar(im_u1, ax=ax_u1, fraction=0.046, pad=0.04)

    im_su3 = ax_su3.imshow(
        _extract_slice(frames_su3[0]).T, origin="lower", extent=extent,
        cmap="Oranges", norm=_norm(vmax_su3), interpolation="nearest", animated=True,
    )
    ax_su3.set_xlabel(xlabel); ax_su3.set_ylabel(ylabel)
    ax_su3.set_title(f"SU(3) traceless energy density  ({plane} midplane)")
    plt.colorbar(im_su3, ax=ax_su3, fraction=0.046, pad=0.04)

    # Time-series panel: integrated channel energies + moving marker.
    ax_ts.plot(times_arr, u1_e, color="tab:blue", label="U(1) trace (EM)")
    ax_ts.plot(times_arr, su3_e, color="tab:orange", label="SU(3) traceless (colour)")
    ax_ts.set_xlabel("time"); ax_ts.set_ylabel("box-integrated channel energy")
    ax_ts.set_title("Integrated U(1) vs SU(3) energy over the loop")
    ax_ts.legend(loc="upper right", fontsize=9)
    y_top = float(max(u1_e.max(), su3_e.max())) or 1.0
    ax_ts.set_ylim(0.0, 1.08 * y_top)
    vline = ax_ts.axvline(times_arr[0], color="0.3", lw=1.0, ls="--")
    dot_u1, = ax_ts.plot([times_arr[0]], [u1_e[0]], "o", color="tab:blue", ms=6)
    dot_su3, = ax_ts.plot([times_arr[0]], [su3_e[0]], "o", color="tab:orange", ms=6)
    frac_txt = ax_ts.text(
        0.03, 0.95, "", transform=ax_ts.transAxes, fontsize=9, va="top",
        bbox=dict(boxstyle="round", fc="white", ec="0.7", alpha=0.85),
    )

    def _frac_label(idx: int) -> str:
        tot = u1_e[idx] + su3_e[idx] + 1e-40
        return (f"t = {float(times_arr[idx]):.2f}\n"
                f"U(1): {100 * u1_e[idx] / tot:5.1f}%\n"
                f"SU(3): {100 * su3_e[idx] / tot:5.1f}%")

    frac_txt.set_text(_frac_label(0))

    def _update(frame_idx: int):
        im_u1.set_data(_extract_slice(frames_u1[frame_idx]).T)
        im_su3.set_data(_extract_slice(frames_su3[frame_idx]).T)
        t = times_arr[frame_idx]
        vline.set_xdata([t, t])
        dot_u1.set_data([t], [u1_e[frame_idx]])
        dot_su3.set_data([t], [su3_e[frame_idx]])
        frac_txt.set_text(_frac_label(frame_idx))
        return im_u1, im_su3, vline, dot_u1, dot_su3, frac_txt

    anim = FuncAnimation(fig, _update, frames=n_frames, interval=1000 // fps, blit=False)
    out = Path(output_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    writer = FFMpegWriter(fps=fps, bitrate=2800)
    anim.save(str(out), writer=writer, dpi=dpi)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Berry / emergent-EM animation: B-flux threading the axis + running Berry phase
# ---------------------------------------------------------------------------


def create_berry_em_animation(
    frames_amp: List[np.ndarray],
    frames_phase: List[np.ndarray],
    grid_shape: Tuple[int, int, int],
    spacing: float,
    output_path: str,
    times: Optional[Sequence[float]] = None,
    n_t: Optional[int] = None,
    fps: int = 15,
    dpi: int = 120,
    title_prefix: str = "",
    quiver_stride: int = 4,
) -> None:
    """Animate the emergent-EM / Berry fields of the carrier over the time loop.

    Two panels per frame:
      - Left: the magnetic flux ``B_z`` threading the vortex axis on the xy-midplane
        (``B_z = ∂_x A_y − ∂_y A_x`` of the U(1) Berry connection
        ``A_i = Im(ψ̂* ∂_i ψ̂)``, ``ψ̂ = Ψ/|Ψ|``), with the in-plane connection as a
        quiver — the "B field threads the donut hole" picture.
      - Right: the running **accumulated Berry phase** ∮A_t over the loop (the
        geometric/EM-charge channel), with a moving marker and the ``2π·n_t``
        closure target if ``n_t`` is given.

    Same Berry-connection math the D4/D5 diagnostics use; this is its *visualization*.

    CAVEAT (read with the trust ledger): on a seed or a NON-converged state the
    connection is a faithful picture of the field but a read-back of the
    prescribed/transient carrier (E8) — only on a CONVERGED solution does it depict
    emergent EM with a quantized (``2π·n_t``) loop phase.
    """
    if not frames_amp:
        raise ValueError("frames_amp must be non-empty")
    nx, ny, nz = grid_shape
    n_frames = len(frames_amp)
    times_arr = list(times) if times is not None else list(range(n_frames))
    zc = nz // 2

    def _mid(vol: np.ndarray) -> np.ndarray:
        return vol[:, :, zc]

    def _connection_and_Bz(amp2d: np.ndarray, ph2d: np.ndarray):
        """In-plane Berry connection (A_x,A_y) and out-of-plane B_z on a 2D slice."""
        psi = amp2d * np.exp(1j * ph2d)
        amax = float(np.max(np.abs(psi))) or 1.0
        eps = 1e-3 * amax
        psihat = np.where(np.abs(psi) > eps, psi / (np.abs(psi) + 1e-300), 0.0 + 0.0j)
        gx = np.gradient(psihat, spacing, axis=0)
        gy = np.gradient(psihat, spacing, axis=1)
        A_x = np.imag(np.conj(psihat) * gx)
        A_y = np.imag(np.conj(psihat) * gy)
        B_z = np.gradient(A_y, spacing, axis=0) - np.gradient(A_x, spacing, axis=1)
        return A_x, A_y, B_z

    # Per-slice B_z + connection (precompute), and the accumulated Berry phase.
    Bz_frames, Ax_frames, Ay_frames = [], [], []
    dphase = np.zeros(n_frames)
    for l in range(n_frames):
        amp2d = _mid(frames_amp[l]); ph2d = _mid(frames_phase[l])
        A_x, A_y, B_z = _connection_and_Bz(amp2d, ph2d)
        Bz_frames.append(B_z); Ax_frames.append(A_x); Ay_frames.append(A_y)
        # amplitude-weighted wrapped temporal phase advance -> accumulated Berry phase
        if l > 0:
            d = _mid(frames_phase[l]) - _mid(frames_phase[l - 1])
            d = np.mod(d + np.pi, 2 * np.pi) - np.pi
            w = _mid(frames_amp[l]) ** 2
            dphase[l] = float(np.sum(w * d) / (np.sum(w) + 1e-300))
    berry_accum = np.cumsum(dphase)

    bz_max = max(float(np.max(np.abs(b))) for b in Bz_frames) or 1.0
    extent = [0, nx * spacing, 0, ny * spacing]
    xs = np.arange(0, nx, quiver_stride) * spacing
    ys = np.arange(0, ny, quiver_stride) * spacing
    QX, QY = np.meshgrid(xs, ys, indexing="ij")

    _apply_style() if "_apply_style" in globals() else None
    fig, (ax_b, ax_p) = plt.subplots(1, 2, figsize=(15, 6.5))
    fig.suptitle(f"{title_prefix}  —  emergent EM: B-flux on axis + accumulated Berry phase",
                 fontsize=11, fontweight="bold")

    im = ax_b.imshow(Bz_frames[0].T, origin="lower", extent=extent, cmap="RdBu_r",
                     vmin=-bz_max, vmax=bz_max, interpolation="nearest", animated=True)
    plt.colorbar(im, ax=ax_b, fraction=0.046, pad=0.04, label="B_z (flux ∥ axis)")
    qv = ax_b.quiver(QX, QY, Ax_frames[0][::quiver_stride, ::quiver_stride],
                     Ay_frames[0][::quiver_stride, ::quiver_stride],
                     color="k", alpha=0.5, scale=None)
    ax_b.set_xlabel("x"); ax_b.set_ylabel("y")
    ax_b.set_title("B_z threading the vortex axis (xy-midplane) + Berry connection A")

    ax_p.plot(times_arr, berry_accum, color="tab:purple", label="∮A_t accumulated")
    if n_t is not None:
        for s in (+1, -1):
            ax_p.axhline(s * 2 * np.pi * n_t, color="0.6", ls="--", lw=0.9)
        ax_p.text(0.02, 0.95, f"closure target ±2π·n_t = ±{2*np.pi*n_t:.2f}",
                  transform=ax_p.transAxes, fontsize=8, va="top")
    marker, = ax_p.plot([times_arr[0]], [berry_accum[0]], "o", color="tab:purple", ms=7)
    ax_p.set_xlabel("time"); ax_p.set_ylabel("accumulated Berry phase [rad]")
    ax_p.set_title("Accumulated carrier (Berry) phase over the loop")
    ax_p.legend(loc="lower left", fontsize=9)

    def _update(i: int):
        im.set_data(Bz_frames[i].T)
        qv.set_UVC(Ax_frames[i][::quiver_stride, ::quiver_stride],
                   Ay_frames[i][::quiver_stride, ::quiver_stride])
        marker.set_data([times_arr[i]], [berry_accum[i]])
        ax_b.set_title(f"B_z + Berry connection  (t={float(times_arr[i]):.2f})")
        return im, qv, marker

    anim = FuncAnimation(fig, _update, frames=n_frames, interval=1000 // fps, blit=False)
    out = Path(output_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    writer = FFMpegWriter(fps=fps, bitrate=2800)
    anim.save(str(out), writer=writer, dpi=dpi)
    plt.close(fig)
