"""Post-process helper for sprint4b Skyrme-corrected confinement analysis.

Usage (called by run_sweep.sh after each branesim-run invocation):

    python analyze_run.py <run_dir> [--weight-mode strain|lateral]

Reads <run_dir>/worldvolume.zip and <run_dir>/summary.json, computes
confinement metrics via branesim.diagnostics.confinement.confinement_from_worldvolume,
and prints one CSV row to stdout with the fields:

    run_stem,alpha,u0,w,profile_shape,tanh_steepness,n_slices,n_nodes,
    m_ambient,grid_shape,t_phys,walltime_s,max_abs_displacement,
    box_fill_radius,spread_ratio_mean,spread_ratio_final,
    confined_fraction_mean,confined_fraction_final,radius_growth,
    interior_residual_norm,weight_mode,status

The confinement_from_worldvolume call uses weight_mode="strain" by default
(correct for Skyrme/topological solitons; see confinement.py docstring on the
Skyrme-aware note -- lateral weights have algebraic tails for hedgehog seeds).

The worldvolume.zip written by branesim.run_experiment now includes a
manifest["lattice"]["dim"] field (fixed in run_experiment.py r2026-06-04),
so confinement_from_worldvolume resolves dim without external config.

Principles compliance:
  - Read-only: does not modify worldvolume or solver state.
  - No physics: pure IO/diagnostics wiring.
  - Dimension-agnostic: dim is read from manifest.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np


def _parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Analyze a branesim sprint4b run directory.")
    p.add_argument("run_dir", nargs="?", default=None,
                   help="Path to a run directory (contains worldvolume.zip + summary.json). "
                        "Not required when --header is used.")
    p.add_argument("--weight-mode", default="strain", choices=["strain", "lateral"],
                   help="Confinement weight mode. 'strain' is Skyrme-correct (default).")
    p.add_argument("--header", action="store_true", help="Print CSV header line and exit.")
    args = p.parse_args()
    if not args.header and args.run_dir is None:
        p.error("run_dir is required unless --header is given.")
    return args


CSV_HEADER = (
    "run_stem,alpha,u0,w,profile_shape,tanh_steepness,"
    "n_slices,n_nodes,m_ambient,grid_shape,t_phys,"
    "walltime_s,max_abs_displacement,"
    "box_fill_radius,spread_ratio_mean,spread_ratio_final,"
    "confined_fraction_mean,confined_fraction_final,radius_growth,"
    "interior_residual_norm,weight_mode,status"
)


def main() -> None:
    args = _parse_args()

    if args.header:
        print(CSV_HEADER)
        return

    run_dir = Path(args.run_dir)
    wv_path = run_dir / "worldvolume.zip"
    summary_path = run_dir / "summary.json"
    run_stem = run_dir.name

    def _row_error(msg: str) -> None:
        # Emit a minimal row with error status so sweep.csv stays well-formed.
        print(f"{run_stem},,,,,,,,,,,,,,,,,,,,,{msg}", flush=True)
        sys.exit(0)  # non-fatal: sweep continues

    if not wv_path.exists():
        _row_error("error:no-worldvolume-zip")
    if not summary_path.exists():
        _row_error("error:no-summary-json")

    try:
        summary = json.loads(summary_path.read_text(encoding="utf-8"))
    except Exception as exc:
        _row_error(f"error:summary-parse:{exc!r}")
        return

    # Extract config fields for the CSV.
    config = summary.get("config", {})
    seed_cfg = config.get("seed", {})
    action_cfg = config.get("action", {})
    lattice_cfg = config.get("lattice", {})

    alpha = float(action_cfg.get("alpha", float("nan")))
    u0 = float(seed_cfg.get("u0", float("nan")))
    w = float(seed_cfg.get("w", float("nan")))
    profile_shape = str(seed_cfg.get("profile_shape", ""))
    tanh_steepness = float(seed_cfg.get("tanh_steepness", float("nan")))
    n_slices = int(summary.get("n_slices", -1))
    n_nodes = int(summary.get("n_nodes", -1))
    m_ambient = int(summary.get("m_ambient", -1))
    grid_shape = str(lattice_cfg.get("grid_shape", ""))
    dt = float(action_cfg.get("dt", float("nan")))
    t_phys = n_slices * dt
    walltime_s = float(summary.get("walltime_s", float("nan")))
    max_abs_disp = float(summary.get("max_abs_displacement", float("nan")))
    res_norm = float(summary.get("interior_residual_norm", float("nan")))

    # Run confinement diagnostics.
    weight_mode = args.weight_mode
    try:
        from branesim.diagnostics.confinement import confinement_from_worldvolume

        if weight_mode == "strain":
            # confinement_from_worldvolume uses lateral by default; for strain we need
            # the lattice neighbor table, which confinement_from_worldvolume does not
            # currently accept.  Fall back to the summary-based confinement_summary call
            # which accepts a lattice object directly.
            from branesim.core.conventions import LatticeParams
            from branesim.core.lattice import SpacelikeLattice
            from branesim.diagnostics.confinement import confinement_summary
            from branesim.io.contracts import iter_slices, load_npy, load_manifest

            manifest = load_manifest(wv_path)
            dim = int(manifest["lattice"]["dim"])
            ref = load_npy(wv_path, "aux/ref_positions.npy")

            all_positions = []
            for _index, _time, positions in iter_slices(wv_path):
                all_positions.append(positions)
            slices = np.stack(all_positions, axis=0)

            # Rebuild lattice from manifest for neighbor table.
            gs = manifest["lattice"]["grid_shape"]
            spacing = float(manifest["lattice"].get("spacing", 1.0))
            periodic_axes = manifest["lattice"].get("periodic_axes", [True] * dim)
            axial_weight = float(manifest["lattice"].get("axial_weight", 1.0))
            lp = LatticeParams(
                grid_shape=tuple(int(v) for v in gs),
                spacing=spacing,
                periodic_axes=tuple(bool(v) for v in periodic_axes),
                axial_weight=axial_weight,
            )
            lat = SpacelikeLattice(lp)

            result = confinement_summary(
                slices, ref, dim,
                weight_mode="strain",
                lattice=lat,
                k_s=float(action_cfg.get("k_s", 1.0)),
                alpha=alpha,
            )
        else:
            result = confinement_from_worldvolume(wv_path)

        box_fill_radius = float(result["box_fill_radius"])
        spread_ratio_mean = float(np.mean(result["spread_ratio"]))
        spread_ratio_final = float(result["final"]["spread_ratio"])
        confined_fraction_mean = float(np.mean(result["confined_fraction"]))
        confined_fraction_final = float(result["final"]["confined_fraction"])
        radius_growth = float(result["radius_growth"])
        status = "ok"

    except Exception as exc:
        box_fill_radius = float("nan")
        spread_ratio_mean = float("nan")
        spread_ratio_final = float("nan")
        confined_fraction_mean = float("nan")
        confined_fraction_final = float("nan")
        radius_growth = float("nan")
        status = f"diag-error:{exc!r}"[:80]

    # Emit CSV row.
    row = (
        f"{run_stem},{alpha},{u0},{w},{profile_shape},{tanh_steepness},"
        f"{n_slices},{n_nodes},{m_ambient},{grid_shape!r},{t_phys:.2f},"
        f"{walltime_s:.1f},{max_abs_disp:.6g},"
        f"{box_fill_radius:.6g},{spread_ratio_mean:.6f},{spread_ratio_final:.6f},"
        f"{confined_fraction_mean:.6f},{confined_fraction_final:.6f},{radius_growth:.6f},"
        f"{res_norm:.6g},{weight_mode},{status}"
    )
    print(row, flush=True)


if __name__ == "__main__":
    main()