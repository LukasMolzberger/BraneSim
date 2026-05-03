"""Dispersion / isotropy diagnostic (sprint 1, subtask 2/3).

Reads a trajectory produced by `components.simulation.run` from a plane-wave
initial condition and extracts the angular frequency `ω(k)` of the standing
wave by projecting each frame onto the seeded wavevector.

Independently callable:
    python -m components.diagnostics.dispersion --input traj.zip --output-dir out/

Output:
    dispersion_summary.json — extracted ω, c=ω/|k|, and the χ²-of-fit
    dispersion_series.npz   — time series of the projection coefficients
"""

from __future__ import annotations

import argparse
import json
import os
import tempfile
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "matplotlib"))
os.environ.setdefault("XDG_CACHE_HOME", tempfile.gettempdir())

import numpy as np

from components.shared import iter_frames, load_manifest, load_npy


def _project_onto_wavevector(
    disp_xyz: np.ndarray,
    coords_xyz: np.ndarray,
    k_vec: np.ndarray,
    polarization: np.ndarray,
) -> tuple[float, float]:
    """Return (cos-projection, sin-projection) of the displacement field onto
    the seeded plane-wave mode.

    For a standing wave `u(x, t) = ε p̂ cos(k·x) cos(ω t)`, the cos-projection
    `(2/N) Σ p̂·u(x_n) cos(k·x_n)` evaluates to `ε cos(ω t)` (modulo
    discrete-Fourier normalization), and the sin-projection vanishes.
    """
    phase = coords_xyz @ k_vec
    cosp = np.cos(phase)
    sinp = np.sin(phase)
    proj = disp_xyz @ polarization
    norm = 2.0 / coords_xyz.shape[0]
    return norm * float(np.sum(proj * cosp)), norm * float(np.sum(proj * sinp))


def _fit_oscillation(times: np.ndarray, signal: np.ndarray) -> dict[str, float]:
    """Fit `signal(t) = A cos(ω t + φ) + B` by minimising least-squares over
    a coarse grid in ω, refined locally. Returns frequency, amplitude, phase,
    offset, and R² of the fit.
    """
    # Estimate the frequency by counting zero crossings of the centered signal.
    s = signal - signal.mean()
    sign_changes = np.where(np.diff(np.sign(s)))[0]
    if len(sign_changes) < 2:
        return {"omega": 0.0, "amplitude": 0.0, "phase": 0.0, "offset": float(signal.mean()), "r2": 0.0}
    half_periods = np.diff(times[sign_changes])
    period_est = 2.0 * float(np.median(half_periods))
    omega_est = 2.0 * np.pi / max(period_est, 1e-12)

    # Refine over a 1D grid centered on the estimate.
    omegas = np.linspace(omega_est * 0.7, omega_est * 1.3, 400)
    best = None
    for omega in omegas:
        c = np.cos(omega * times)
        sn = np.sin(omega * times)
        design = np.stack([c, sn, np.ones_like(times)], axis=1)
        coeff, *_ = np.linalg.lstsq(design, signal, rcond=None)
        residual = signal - design @ coeff
        ss_res = float(np.sum(residual ** 2))
        if best is None or ss_res < best[0]:
            best = (ss_res, omega, coeff)

    ss_res, omega_fit, (a_cos, a_sin, offset) = best
    amplitude = float(np.hypot(a_cos, a_sin))
    phase = float(np.arctan2(-a_sin, a_cos))
    ss_tot = float(np.sum((signal - signal.mean()) ** 2)) + 1e-30
    r2 = 1.0 - ss_res / ss_tot

    return {
        "omega": float(omega_fit),
        "amplitude": amplitude,
        "phase": phase,
        "offset": float(offset),
        "r2": float(r2),
    }


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Dispersion / isotropy diagnostic")
    p.add_argument("--input", required=True)
    p.add_argument("--output-dir", required=True)
    p.add_argument("--frame-stride", type=int, default=1)
    p.add_argument("--max-frames", type=int, default=None)
    return p.parse_args()


def main() -> None:
    args = parse_args()

    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)

    manifest = load_manifest(args.input)
    lattice = manifest["lattice"]
    spacing = float(lattice["spacing"])
    nx, ny, nz = (int(v) for v in lattice["grid_shape"])

    metadata = manifest.get("metadata", {})
    if metadata.get("seed_kind") != "plane_wave":
        raise ValueError(
            f"This diagnostic expects a plane-wave seed; manifest reports seed_kind={metadata.get('seed_kind')!r}."
        )
    debug = metadata.get("debug", {})
    k_vec = np.asarray(debug["k_vec"], dtype=np.float64)
    polarization = np.asarray(debug["polarization"], dtype=np.float64)
    k_index = tuple(int(v) for v in debug["k_index"])
    seed_amplitude = float(metadata.get("seed", {}).get("amplitude", 0.0)) or float(debug.get("amplitude", 0.0))

    rest = load_npy(args.input, "aux/rest_positions.npy")
    coords = rest[:, :3].astype(np.float64, copy=False)

    times: list[float] = []
    cos_proj: list[float] = []
    sin_proj: list[float] = []

    for i, fr in enumerate(iter_frames(args.input, frame_stride=args.frame_stride)):
        if args.max_frames is not None and i >= args.max_frames:
            break
        disp = (fr.positions[:, :3] - coords).astype(np.float64, copy=False)
        c, s = _project_onto_wavevector(disp, coords, k_vec, polarization)
        times.append(float(fr.time))
        cos_proj.append(c)
        sin_proj.append(s)

    if len(times) < 8:
        raise ValueError(f"Need at least 8 frames for a stable fit, got {len(times)}.")

    t_arr = np.asarray(times, dtype=np.float64)
    cos_arr = np.asarray(cos_proj, dtype=np.float64)
    sin_arr = np.asarray(sin_proj, dtype=np.float64)

    fit = _fit_oscillation(t_arr, cos_arr)

    k_mag = float(np.linalg.norm(k_vec))
    k_a = k_mag * spacing
    c_meas = fit["omega"] / k_mag if k_mag > 0 else 0.0
    direction = k_vec / max(k_mag, 1e-30)

    summary = {
        "trajectory": str(Path(args.input).resolve()),
        "k_index": list(k_index),
        "k_vec": [float(v) for v in k_vec],
        "k_magnitude": k_mag,
        "k_a": k_a,
        "k_direction": [float(v) for v in direction],
        "polarization": [float(v) for v in polarization],
        "seed_amplitude": seed_amplitude,
        "num_frames": len(times),
        "t_span": [float(t_arr[0]), float(t_arr[-1])],
        "fit": fit,
        "c_measured": float(c_meas),
        "lattice": {"grid_shape": [nx, ny, nz], "spacing": spacing},
    }

    series_path = out / "dispersion_series.npz"
    np.savez_compressed(
        series_path,
        times=t_arr,
        cos_projection=cos_arr,
        sin_projection=sin_arr,
        k_vec=k_vec,
        polarization=polarization,
    )
    summary["series_npz"] = str(series_path)

    summary_path = out / "dispersion_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    summary["summary_json"] = str(summary_path)

    print("Dispersion diagnostic complete")
    print(f"  k_index = {k_index}, |k|·a = {k_a:.4f}")
    print(f"  fit ω = {fit['omega']:.6f}  (R² = {fit['r2']:.4f})")
    print(f"  c = ω/|k| = {c_meas:.6f}")
    print(f"  → {summary_path}")


if __name__ == "__main__":
    main()
