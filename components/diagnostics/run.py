"""Component 4: diagnostics.

Computes QCD-inspired metrics + Berry diagnostics from trajectory file.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import tempfile
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "matplotlib"))
os.environ.setdefault("XDG_CACHE_HOME", tempfile.gettempdir())

import numpy as np

from components.shared import iter_frames, load_manifest, load_npy
from components.diagnostics.berry import compute_berry_time_series, create_berry_videos


def _safe_corrcoef(x: np.ndarray) -> np.ndarray:
    if x.shape[0] < 3:
        return np.full((x.shape[1], x.shape[1]), np.nan)
    if np.allclose(x, x[0]):
        return np.full((x.shape[1], x.shape[1]), np.nan)
    return np.corrcoef(x, rowvar=False)


def _box_fill_radius(coords_xyz: np.ndarray) -> float:
    """RMS node radius about the geometric centre (uniform weight).

    Frame-invariant (depends only on the rest geometry). This is the value the
    energy-weighted ``radius_rms`` approaches if the field disperses to fill the
    box uniformly, so it is the natural fixed reference for confinement.
    """
    geom_center = np.mean(coords_xyz, axis=0)
    r_geom = np.linalg.norm(coords_xyz - geom_center[None, :], axis=1)
    return float(np.sqrt(np.mean(r_geom * r_geom)))


def _frame_metrics(
    disp_xyz: np.ndarray,
    vel_xyz: np.ndarray,
    coords_xyz: np.ndarray,
    omega_ref: float,
    leakage_radius_factor: float,
    confinement_radius_factor: float,
    box_fill_radius: float,
) -> dict[str, float]:
    # NOTE: this function is intentionally 3-channel. The 3 columns are the
    # spacelike lateral triplet (the SU(3) gauge / colour sector, principles
    # §1.1a), not an arbitrary spatial dimension; the caller slices [:, :3].
    energy_proxy = 0.5 * (disp_xyz ** 2 + (vel_xyz / max(omega_ref, 1e-12)) ** 2)
    channel_energy = np.mean(energy_proxy, axis=0)
    total_energy = float(np.sum(channel_energy))
    fractions = channel_energy / total_energy if total_energy > 1e-30 else np.zeros(3)

    trace_mode = np.sum(disp_xyz, axis=1) / np.sqrt(3.0)
    traceless = disp_xyz - np.mean(disp_xyz, axis=1, keepdims=True)

    amp = np.linalg.norm(disp_xyz, axis=1)
    weights = amp * amp + 1e-20
    total_weight = np.sum(weights)
    center = np.sum(coords_xyz * weights[:, None], axis=0) / total_weight
    radial = np.linalg.norm(coords_xyz - center[None, :], axis=1)
    radius_rms = float(np.sqrt(np.sum(weights * radial * radial) / total_weight))
    leak_threshold = leakage_radius_factor * radius_rms
    leakage = float(np.sum(weights[radial > leak_threshold]) / total_weight)

    # Confinement (non-self-referential). The leakage metric above uses a
    # threshold proportional to the *current* radius_rms, so it stays ~0 even
    # when the packet disperses; it cannot detect spreading. These metrics
    # instead compare against the fixed box geometry (box_fill_radius, passed
    # in and frame-invariant):
    #   spread_ratio      -- radius_rms / box_fill_radius. -> 1 means fully
    #     dispersed (box-fill); << 1 means localized/confined.
    #   confined_fraction -- weighted-energy fraction within a radius that is a
    #     FIXED fraction of the box. The threshold is a box scale, but the
    #     distance `radial` is measured from the energy-weighted centre, so
    #     this reads "fraction of energy within confinement_radius of wherever
    #     the packet currently sits" (robust to an off-centre packet).
    spread_ratio = radius_rms / box_fill_radius if box_fill_radius > 1e-30 else float("nan")
    confinement_radius = confinement_radius_factor * box_fill_radius
    confined_fraction = float(np.sum(weights[radial <= confinement_radius]) / total_weight)

    corr = _safe_corrcoef(disp_xyz)
    mixing = float(np.nanmean(np.abs(corr[np.triu_indices(3, k=1)]))) if not np.isnan(corr).all() else float("nan")

    return {
        "energy_total_proxy": total_energy,
        "energy_frac_x": float(fractions[0]),
        "energy_frac_y": float(fractions[1]),
        "energy_frac_z": float(fractions[2]),
        "trace_rms": float(np.sqrt(np.mean(trace_mode * trace_mode))),
        "traceless_rms": float(np.sqrt(np.mean(np.sum(traceless * traceless, axis=1)))),
        "radius_rms": radius_rms,
        "box_fill_radius": box_fill_radius,
        "spread_ratio": spread_ratio,
        "confined_fraction": confined_fraction,
        "leakage_fraction": leakage,
        "mixing_strength": mixing,
    }


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Component 4: diagnostics")
    p.add_argument("--input", required=True)
    p.add_argument("--output-dir", required=True)
    p.add_argument("--frame-stride", type=int, default=1)
    p.add_argument("--max-frames", type=int, default=None)
    p.add_argument("--berry-point-stride", type=int, default=4)
    p.add_argument("--berry-omega0", type=float, default=None)
    p.add_argument("--omega-ref", type=float, default=1.0)
    p.add_argument("--leakage-radius-factor", type=float, default=2.0)
    p.add_argument("--confinement-radius-factor", type=float, default=0.5,
                   help="confined_fraction counts energy within this fraction of "
                        "box_fill_radius (a fixed box scale, not the packet's own spread).")
    p.add_argument("--render-berry-videos", action="store_true")
    return p.parse_args()


def main() -> None:
    args = parse_args()

    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)

    manifest = load_manifest(args.input)
    lattice = manifest["lattice"]
    grid_shape = tuple(int(v) for v in lattice["grid_shape"])
    spacing = float(lattice["spacing"])

    rest = load_npy(args.input, "aux/rest_positions.npy")
    berry_idx = np.arange(0, rest.shape[0], max(int(args.berry_point_stride), 1), dtype=np.int64)

    # Frame-invariant confinement reference (depends only on rest geometry).
    box_fill_radius = _box_fill_radius(rest[:, :3])

    rows: list[dict[str, float]] = []
    berry_u: list[np.ndarray] = []
    berry_v: list[np.ndarray] = []
    berry_t: list[float] = []

    for i, fr in enumerate(iter_frames(args.input, frame_stride=args.frame_stride)):
        if args.max_frames is not None and i >= args.max_frames:
            break

        disp = fr.positions - rest
        row = {
            "frame_index": float(fr.index),
            "step": float(fr.step),
            "time": float(fr.time),
        }
        row.update(
            _frame_metrics(
                disp[:, :3],
                fr.velocities[:, :3],
                rest[:, :3],
                omega_ref=float(args.omega_ref),
                leakage_radius_factor=float(args.leakage_radius_factor),
                confinement_radius_factor=float(args.confinement_radius_factor),
                box_fill_radius=box_fill_radius,
            )
        )
        rows.append(row)

        berry_u.append(disp[berry_idx, :])
        berry_v.append(fr.velocities[berry_idx, :])
        berry_t.append(fr.time)

    if not rows:
        raise ValueError("No frames selected for diagnostics")

    csv_path = out / "qcd_metrics.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    summary: dict[str, object] = {
        "trajectory": str(Path(args.input).resolve()),
        "num_frames_used": len(rows),
        "frame_stride": int(args.frame_stride),
        "berry_point_stride": int(args.berry_point_stride),
        "qcd_metrics_csv": str(csv_path),
        "final_metrics": rows[-1],
        "mean_metrics": {
            k: float(np.nanmean([row[k] for row in rows]))
            for k in rows[0].keys()
            if k not in {"frame_index", "step", "time"}
        },
        # Spreading over the run: radius_rms(final) / radius_rms(first sampled
        # frame; note this is the first frame at --frame-stride, not necessarily
        # t=0). > 1 means the mode expanded from its seed; combined with a final
        # spread_ratio near 1 this is the signature of dispersion to box-fill
        # (the opposite of confinement).
        "radius_growth": (
            rows[-1]["radius_rms"] / rows[0]["radius_rms"]
            if rows[0]["radius_rms"] > 1e-30 else float("nan")
        ),
    }

    if len(berry_u) >= 2:
        series = compute_berry_time_series(
            berry_u,
            berry_v,
            berry_t,
            omega0=args.berry_omega0,
        )
        berry_npz = out / "berry_series.npz"
        np.savez_compressed(
            berry_npz,
            phase=series.phase,
            connection=series.connection,
            delta_phase=series.delta_phase,
            alpha=series.alpha,
            times_s=series.times_s,
            overlap_abs=series.overlap_abs,
            omega0=np.array([series.omega0], dtype=np.float64),
        )
        summary["berry_series"] = str(berry_npz)
        summary["berry_phase_std_final"] = float(np.std(series.phase[-1]))
        summary["berry_connection_mean_abs"] = float(np.mean(np.abs(series.connection)))

        if args.render_berry_videos:
            berry_dir = out / "berry_videos"
            berry_dir.mkdir(parents=True, exist_ok=True)
            videos = create_berry_videos(
                series,
                output_dir=str(berry_dir),
                grid_shape_3d=grid_shape,
                spacing=spacing,
                filename_prefix="baryon_berry",
            )
            summary["berry_videos"] = videos

    summary_path = out / "diagnostics_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    summary["summary_json"] = str(summary_path)

    print("Diagnostics complete")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
