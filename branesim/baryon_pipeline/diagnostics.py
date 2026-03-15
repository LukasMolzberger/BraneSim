"""Diagnostics component: Berry analysis and QCD-inspired baryon metrics."""

from __future__ import annotations

import csv
import json
import os
import tempfile
from pathlib import Path
from typing import Any

os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "matplotlib"))
os.environ.setdefault("XDG_CACHE_HOME", tempfile.gettempdir())

import numpy as np

from branesim.berry import compute_berry_time_series, create_berry_videos

from .io import iter_trajectory_frames, load_trajectory_array, load_trajectory_manifest


def _safe_corrcoef(x: np.ndarray) -> np.ndarray:
    if x.shape[0] < 3:
        return np.full((x.shape[1], x.shape[1]), np.nan)
    if np.allclose(x, x[0]):
        return np.full((x.shape[1], x.shape[1]), np.nan)
    return np.corrcoef(x, rowvar=False)


def _frame_qcd_metrics(
    disp_xyz: np.ndarray,
    vel_xyz: np.ndarray,
    coords_xyz: np.ndarray,
    omega_ref: float,
    leakage_radius_factor: float,
) -> dict[str, float]:
    energy_proxy = 0.5 * (disp_xyz ** 2 + (vel_xyz / max(omega_ref, 1e-12)) ** 2)
    channel_energy = np.mean(energy_proxy, axis=0)
    total_energy = float(np.sum(channel_energy))
    if total_energy < 1e-30:
        fractions = np.array([0.0, 0.0, 0.0], dtype=np.float64)
    else:
        fractions = channel_energy / total_energy

    trace_mode = np.sum(disp_xyz, axis=1) / np.sqrt(3.0)
    traceless = disp_xyz - np.mean(disp_xyz, axis=1, keepdims=True)

    amp = np.linalg.norm(disp_xyz, axis=1)
    weights = amp * amp + 1e-20
    weighted_center = np.sum(coords_xyz * weights[:, None], axis=0) / np.sum(weights)
    radial = np.linalg.norm(coords_xyz - weighted_center[None, :], axis=1)
    radius_rms = np.sqrt(np.sum(weights * radial * radial) / np.sum(weights))
    leak_threshold = leakage_radius_factor * radius_rms
    leakage = float(np.sum(weights[radial > leak_threshold]) / np.sum(weights))

    corr = _safe_corrcoef(disp_xyz)
    if np.isnan(corr).all():
        mixing_strength = float("nan")
    else:
        off_diag = np.abs(corr[np.triu_indices(3, k=1)])
        mixing_strength = float(np.nanmean(off_diag))

    return {
        "energy_total_proxy": total_energy,
        "energy_frac_x": float(fractions[0]),
        "energy_frac_y": float(fractions[1]),
        "energy_frac_z": float(fractions[2]),
        "trace_rms": float(np.sqrt(np.mean(trace_mode * trace_mode))),
        "traceless_rms": float(np.sqrt(np.mean(np.sum(traceless * traceless, axis=1)))),
        "radius_rms": float(radius_rms),
        "leakage_fraction": leakage,
        "mixing_strength": mixing_strength,
    }


def run_diagnostics_component(
    trajectory_path: str | Path,
    output_dir: str | Path,
    *,
    frame_stride: int = 1,
    max_frames: int | None = None,
    berry_point_stride: int = 4,
    berry_omega0: float | None = None,
    omega_ref: float = 1.0,
    leakage_radius_factor: float = 2.0,
    render_berry_videos: bool = False,
) -> dict[str, Any]:
    """Compute offline diagnostics from stored simulation trajectory."""

    manifest = load_trajectory_manifest(trajectory_path)
    lattice = manifest["lattice"]
    grid_shape = tuple(int(v) for v in lattice["grid_shape"])
    spacing = float(lattice["spacing"])

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    rest_positions = load_trajectory_array(trajectory_path, "aux/rest_positions.npy")

    berry_idx = np.arange(0, rest_positions.shape[0], max(int(berry_point_stride), 1), dtype=np.int64)
    berry_u_frames: list[np.ndarray] = []
    berry_v_frames: list[np.ndarray] = []
    berry_times: list[float] = []

    metrics_rows: list[dict[str, float]] = []

    frame_counter = 0
    for frame in iter_trajectory_frames(trajectory_path, frame_stride=frame_stride):
        if max_frames is not None and frame_counter >= max_frames:
            break

        disp = frame.positions - rest_positions
        disp_xyz = disp[:, 0:3]
        vel_xyz = frame.velocities[:, 0:3]
        coords_xyz = rest_positions[:, 0:3]

        row = {
            "frame_index": float(frame.index),
            "step": float(frame.step),
            "time": float(frame.time),
        }
        row.update(_frame_qcd_metrics(disp_xyz, vel_xyz, coords_xyz, omega_ref, leakage_radius_factor))
        metrics_rows.append(row)

        berry_u_frames.append(disp[berry_idx, :])
        berry_v_frames.append(frame.velocities[berry_idx, :])
        berry_times.append(frame.time)

        frame_counter += 1

    if len(metrics_rows) == 0:
        raise ValueError("No frames selected for diagnostics")

    berry_series = None
    if len(berry_u_frames) >= 2:
        berry_series = compute_berry_time_series(
            berry_u_frames,
            berry_v_frames,
            berry_times,
            omega0=berry_omega0,
            return_psi=False,
        )

    csv_path = output_dir / "qcd_metrics.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(metrics_rows[0].keys()))
        writer.writeheader()
        writer.writerows(metrics_rows)

    summary: dict[str, Any] = {
        "trajectory": str(trajectory_path),
        "num_frames_used": len(metrics_rows),
        "frame_stride": frame_stride,
        "berry_point_stride": int(berry_point_stride),
        "qcd_metrics_csv": str(csv_path),
        "final_metrics": metrics_rows[-1],
        "mean_metrics": {
            key: float(np.nanmean([row[key] for row in metrics_rows]))
            for key in metrics_rows[0].keys()
            if key not in {"frame_index", "step", "time"}
        },
    }

    if berry_series is not None:
        berry_npz = output_dir / "berry_series.npz"
        np.savez_compressed(
            berry_npz,
            phase=berry_series.phase,
            connection=berry_series.connection,
            delta_phase=berry_series.delta_phase,
            alpha=berry_series.alpha,
            times_s=berry_series.times_s,
            overlap_abs=berry_series.overlap_abs,
            omega0=np.array([berry_series.omega0], dtype=np.float64),
        )
        summary["berry_series"] = str(berry_npz)
        summary["berry_phase_std_final"] = float(np.std(berry_series.phase[-1]))
        summary["berry_connection_mean_abs"] = float(np.mean(np.abs(berry_series.connection)))

        if render_berry_videos:
            berry_video_dir = output_dir / "berry_videos"
            berry_video_dir.mkdir(parents=True, exist_ok=True)
            videos = create_berry_videos(
                berry_series,
                intrinsic_dim=3,
                output_dir=str(berry_video_dir),
                grid_shape_3d=grid_shape,
                spacing=spacing,
                filename_prefix="baryon_berry",
            )
            summary["berry_videos"] = videos

    summary_path = output_dir / "diagnostics_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    summary["summary_json"] = str(summary_path)

    return summary
