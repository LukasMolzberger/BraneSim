"""Diagnostics-owned Berry analysis and optional video output."""

from __future__ import annotations

import os
from dataclasses import dataclass
from typing import Optional, Sequence, Tuple

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, FFMpegWriter
from matplotlib.colors import hsv_to_rgb


@dataclass(frozen=True)
class BerryTimeSeries:
    phase: np.ndarray
    connection: np.ndarray
    delta_phase: np.ndarray
    amplitude: np.ndarray
    alpha: np.ndarray
    overlap_abs: np.ndarray
    omega0: float
    times_s: np.ndarray


def _as_numpy_frames(frames: Sequence[np.ndarray]) -> list[np.ndarray]:
    out = []
    for f in frames:
        arr = np.asarray(f)
        if arr.ndim != 2:
            raise ValueError(f"Each frame must be 2D (N,E), got {arr.shape}")
        out.append(arr)
    return out


def complex_amplitude_from_u_v(u: np.ndarray, v: np.ndarray, omega0: float) -> np.ndarray:
    if u.shape != v.shape:
        raise ValueError("u and v must have same shape")
    sw = np.sqrt(float(omega0))
    return sw * u.astype(np.float64) + 1j * (v.astype(np.float64) / sw)


def _norm_rows(x: np.ndarray, eps: float = 1e-30) -> np.ndarray:
    return np.sqrt(np.maximum(np.sum(np.abs(x) ** 2, axis=-1), 0.0) + eps * 0.0)


def compute_berry_time_series(
    frames_u: Sequence[np.ndarray],
    frames_v: Sequence[np.ndarray],
    times_s: Sequence[float],
    omega0: Optional[float] = None,
    amp_eps_rel: float = 1e-4,
    overlap_eps_rel: float = 1e-6,
    alpha_gamma: float = 1.0,
    alpha_scale: float = 0.95,
) -> BerryTimeSeries:
    U = _as_numpy_frames(frames_u)
    V = _as_numpy_frames(frames_v)
    if len(U) != len(V):
        raise ValueError("frames_u and frames_v must have same length")
    if len(U) < 2:
        raise ValueError("Need at least 2 frames")

    T = len(U)
    times = np.asarray(times_s, dtype=np.float64)
    if times.shape != (T,):
        raise ValueError(f"times_s must have shape ({T},)")

    N, E = U[0].shape
    for k in range(T):
        if U[k].shape != (N, E) or V[k].shape != (N, E):
            raise ValueError("All frames must share shape (N,E)")

    w0 = float(omega0) if omega0 is not None else 1.0
    if w0 <= 0:
        raise ValueError("omega0 must be positive")

    A = np.empty((T, N, E), dtype=np.complex128)
    amp = np.empty((T, N), dtype=np.float64)
    for k in range(T):
        A[k] = complex_amplitude_from_u_v(U[k], V[k], w0)
        amp[k] = _norm_rows(A[k])

    amp_max = float(np.max(amp)) if np.max(amp) > 0 else 1.0
    amp_eps = float(amp_eps_rel) * amp_max

    PSI = np.empty((T, N, E), dtype=np.complex128)
    for k in range(T):
        denom = np.maximum(amp[k], amp_eps)
        PSI[k] = A[k] / denom[:, None]

    delta = np.zeros((T, N), dtype=np.float64)
    conn = np.zeros((T, N), dtype=np.float64)
    overlap_abs = np.zeros((T, N), dtype=np.float64)

    for k in range(T - 1):
        dt = float(times[k + 1] - times[k])
        if dt <= 0:
            raise ValueError("times_s must be strictly increasing")
        ov = np.sum(np.conj(PSI[k]) * PSI[k + 1], axis=-1)
        overlap_abs[k] = np.abs(ov)
        delta_k = np.angle(ov)
        delta[k] = delta_k
        conn[k] = delta_k / dt

    overlap_abs[T - 1] = overlap_abs[T - 2]
    conn[T - 1] = conn[T - 2]
    delta[T - 1] = 0.0
    phase = np.cumsum(delta, axis=0)

    amp_norm = np.clip(amp / max(amp_max, 1e-30), 0.0, 1.0)
    alpha = alpha_scale * (amp_norm ** float(alpha_gamma))
    alpha = np.where(amp >= amp_eps, alpha, 0.0)
    alpha = np.where(overlap_abs >= float(overlap_eps_rel), alpha, 0.0)
    alpha = np.clip(alpha, 0.0, 1.0)

    return BerryTimeSeries(
        phase=phase,
        connection=conn,
        delta_phase=delta,
        amplitude=amp,
        alpha=alpha,
        overlap_abs=overlap_abs,
        omega0=w0,
        times_s=times,
    )


def _phase_to_rgb(phase: np.ndarray) -> np.ndarray:
    hue = np.mod((phase + np.pi) / (2.0 * np.pi), 1.0)
    sat = np.ones_like(hue)
    val = np.ones_like(hue)
    return hsv_to_rgb(np.stack([hue, sat, val], axis=-1))


def create_berry_videos(
    series: BerryTimeSeries,
    output_dir: str,
    grid_shape_3d: Tuple[int, int, int],
    spacing: float,
    filename_prefix: str = "berry",
    fps: int = 20,
    dpi: int = 120,
) -> list[str]:
    os.makedirs(output_dir, exist_ok=True)

    nx, ny, nz = grid_shape_3d
    z_idx = nz // 2

    def render(values: np.ndarray, alpha: np.ndarray, out_path: str, title: str):
        fig, ax = plt.subplots(figsize=(8, 7))

        v0 = values[0].reshape(nx, ny, nz)[:, :, z_idx]
        a0 = alpha[0].reshape(nx, ny, nz)[:, :, z_idx]
        rgba0 = np.concatenate([_phase_to_rgb(v0), np.clip(a0, 0, 1)[..., None]], axis=-1)

        im = ax.imshow(
            rgba0.swapaxes(0, 1),
            origin="lower",
            extent=[0, spacing * (nx - 1), 0, spacing * (ny - 1)],
            interpolation="nearest",
            animated=True,
        )
        ax.set_title(title)
        ax.set_xlabel("x")
        ax.set_ylabel("y")

        def update(k: int):
            v = values[k].reshape(nx, ny, nz)[:, :, z_idx]
            a = alpha[k].reshape(nx, ny, nz)[:, :, z_idx]
            rgba = np.concatenate([_phase_to_rgb(v), np.clip(a, 0, 1)[..., None]], axis=-1)
            im.set_data(rgba.swapaxes(0, 1))
            return (im,)

        anim = FuncAnimation(fig, update, frames=values.shape[0], interval=1000 / fps, blit=False)
        writer = FFMpegWriter(fps=fps, bitrate=2500)
        anim.save(out_path, writer=writer, dpi=dpi)
        plt.close(fig)

    out_phase = os.path.join(output_dir, f"{filename_prefix}_phase.mp4")
    out_conn = os.path.join(output_dir, f"{filename_prefix}_connection.mp4")
    render(series.phase, series.alpha, out_phase, "Berry phase (XY mid-slice)")
    render(series.connection, series.alpha, out_conn, "Berry connection (XY mid-slice)")

    return [out_phase, out_conn]
